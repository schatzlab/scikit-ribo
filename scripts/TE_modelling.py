#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for fitting the GLM of TE
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 3/31/2016
## ----------------------------------------

from __future__ import print_function
import sys
import argparse
import numpy as np
sys.path.remove('/sonas-hs/lyon/hpc/home/hfang/.local/lib/python3.4/site-packages/statsmodels-0.6.1-py3.4-linux-x86_64.egg') # to be removed
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
from patsy import dmatrices
import pybedtools as pbt
from statsmodels.base.model import GenericLikelihoodModel


class glmTE(object):
    ''' class to perform glm for TE estimation
    '''
    def __init__(self, fn, unmap_fn):
        self.fn = fn
        self.df_all = pd.read_table(self.fn,  header=0)
        self.df_col_names = list(self.df_all.columns.values)
        self.bt_all = pbt.BedTool.from_dataframe(self.df_all)
        self.unmap = pbt.BedTool(unmap_fn)
        self.bt_map = self.bt_all.intersect(self.unmap, v=True)
        self.df = self.bt_map.to_dataframe(names=self.df_col_names)

    def filter(self):
        ## get a summary of the starting df
        print ("[status]\tStarting number of observations:\t" + str(self.df.shape[0]), flush=True)
        print ("[status]\tNumber of variables:\t" + str(self.df.shape[1]-1), flush=True)
        #inds = pd.isnull(self.df).any(1).nonzero()[0]
        #print(inds)
        #print(pd.isnull(self.df).any(1).nonzero()[0])
        #print(self.df.isnull().values.any())

        ## remove missing values
        self.df = self.df.dropna()

        ## group by gene names
        grouped = self.df.groupby("gene_name")
        #ribosome_count_sum = np.array(grouped.aggregate(np.sum).TPM)

        ## codons to exclude , start:'ATG', stop: TAG, TAA, TGA
        self.df = self.df.loc[ (self.df.codon != 'ATG') & (self.df.codon != 'TAG') & (self.df.codon != 'TAA') &
                               (self.df.codon != 'TGA') , :]

        ## define TPM lower bound, fitler df, require TPM > 10, ribosome count > 100
        TPM_lb = 10
        self.df = self.df.loc[self.df.TPM >= TPM_lb, :]
        num_genes = np.unique(self.df.gene_name).shape[0]
        
        ## calculate log(TPM)
        log_TPM = np.log(self.df['TPM'])
        self.df = self.df.assign(log_TPM = log_TPM)

        ## add 1 to ribosome_count, divide self.df.ribosome_count by TPM
        self.df.ribosome_count = self.df.ribosome_count + 1
        #self.df["norm_count"] = (self.df.ribosome_count ) / self.df.TPM

        ## print status
        print ("[status]\tTPM lower bound:\t" + str(TPM_lb), flush=True)
        print ("[status]\tNumber of observations after filtering:\t" + str(self.df.shape[0]), flush=True)
        print ("[status]\tNumber of genes after filtering:\t" + str(num_genes), flush=True)
        
    def nbGlm(self):
        ## define model formula,
        formula = 'ribosome_count ~ C(gene_name) + C(codon) + pair_prob'
        print("[status]\tFormula: " + str(formula), flush=True)
        
        ## define model fitting options
        sovler = "lbfgs" # "IRLS" # "lbfgs"
        tolerence = 1e-8
        num_iter = 100
        print("[status]\tSolver: " + sovler, flush=True)
        print("[status]\tConvergence tolerance: " + str(tolerence), flush=True)
        print("[status]\tMaxiter: " + str(num_iter), flush=True)
        
        ## model fitting NegativeBinomial GLM
        print("[status]\tModel: smf.glm(formula, self.df, family=sm.families.NegativeBinomial(), offset=self.df['log_TPM']", flush=True)
        mod = smf.glm(formula, self.df, family=sm.families.NegativeBinomial(), offset=self.df['log_TPM'])
        res = mod.fit(method=sovler, tol=tolerence, maxiter=num_iter)

        ## print model output
        print ( res.summary() )

        ## alternative fit a gamma dist
        #formula = 'norm_count ~ C(gene_name) + C(codon)  + pair_prob '
        #print("[status]\tFormula: " + str(formula), flush=True)
        #mod = smf.glm(formula, self.df, family=sm.families.Gamma(link=sm.families.links.log))

    def genericGLM(self):
        y, X = dmatrices('ribosome_count ~ C(gene_name) + C(codon) + pair_prob', self.df)
        mod = NBin(y, X)
        res = mod.fit('lbfgs')
        print(res)

class NBin(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(NBin, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        alph = params[-1]
        beta = params[:-1]
        ll = _ll_nb2(self.endog, self.exog, beta, alph)
        return -ll
    def fit(self, start_params=None, maxiter=1000, maxfun=5000, **kwds):
        if start_params == None:
            # Reasonable starting values
            start_params = np.append(np.zeros(self.exog.shape[1]), .5)
            start_params[0] = np.log(self.endog.mean())
        return super(NBin, self).fit(start_params=start_params,
                                     maxiter=maxiter, maxfun=maxfun,
                                     **kwds)



## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input data-frame for modelling")
    parser.add_argument("-u", help="unmappable regions")
    
    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.u!=None):
        ## Load the content of the table
        print("[status]\tstatsmodel version:\t"+ str(statsmodels.__version__), flush=True)
        print("[status]\tReading the ribosome counts: " + str(args.i), flush=True)
        print("[status]\tReading the un-mappable regions: " + str(args.u), flush=True)
        df_fn = args.i
        unmap_fn = args.u

        ## start model fitting
        print("[execute]\tStart the modelling of TE", flush=True)
        model = glmTE(df_fn, unmap_fn)
        print("[execute]\tFilter the df", flush=True)
        model.filter()
        print("[execute]\tFitting the GLM", flush=True)
        model.nbGlm()
        #model.genericGLM()
    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
