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
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np
from patsy import dmatrices
import pybedtools as pbt
from statsmodels.base.model import GenericLikelihoodModel
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as n2r
#from glmnet import PoissonNet
#from memory_profiler import profile


class glmTE(object):
    """ perform glm for TE estimation
    """
    def __init__(self, fn=None, unMappableFn=None, tpmLB=1):
        """
        :param fn:
        :param unMappableFn:
        :param tpmLB:
        """
        self.fn = fn
        self.unMappableFn = unMappableFn
        self.tpmLB = tpmLB

    def loadDat(self):
        """
        loading data: self.df, self.unMappable
        :return: None
        """
        # read the df and un-mappable regions
        self.df = pd.read_table(self.fn,  header=0)
        self.unMappable = pbt.BedTool(self.unMappableFn)
        # compute the codon length
        codons = self.df[["gene", "codon_idx"]]
        codonLen = codons.groupby('gene').max().reset_index()
        codonLen.columns = ['gene', "numCodons"]
        codonLen["numCodons"] = codonLen["numCodons"] - 8
        self.df = pd.merge(self.df, codonLen, how = "left", on = ["gene"])
        # rolling mean of 2' structure pairing probabilities
        paring = self.df[["gene", "codon_idx", "pair_prob"]]
        rollingAvg = paring.groupby("gene").rolling(window=11, center=True).mean().reset_index(0, drop=True)
        rollingAvg.columns = ["gene", "codon_idx", "avg_prob"]
        # merge 2' prob
        self.df = pd.merge(self.df, rollingAvg, how='left', on=["gene", 'codon_idx']).drop(["pair_prob"], axis=1).dropna()

    def filterDf(self):
        """
        filter dataframe, self.df
        :return: None
        """
        # get a summary of the starting df
        sys.stderr.write("[status]\tStarting number of observations:\t" + str(self.df.shape[0]) + "\n")
        sys.stderr.write("[status]\tNumber of variables:\t" + str(self.df.shape[1]-1) + "\n")
        # 1. exclude 5/3 utr codons
        self.df = self.df[(self.df["codon_idx"] <= self.df["numCodons"]) & (self.df["codon_idx"] >= 1)]
        # 2. exclude genes w/ low Ribo TPM
        gene = self.df[["gene", "codon_idx", "ribosome_count"]].groupby("gene").sum().ribosome_count.reset_index(name="riboCnt")
        codonLen = self.df[["gene", "numCodons"]].drop_duplicates(subset=["gene"])
        gene = pd.merge(gene, codonLen, on="gene")
        gene["rpkm"] = gene["riboCnt"] * 1000000000 / (np.sum(gene["riboCnt"]) * gene["numCodons"] * 3)
        gene["tpm"] = gene["rpkm"]  * 1000000 / np.sum(gene["rpkm"])
        gene = gene[(gene["tpm"]>=self.tpmLB)][["gene"]]
        self.df = pd.merge(self.df, gene, on="gene")
        # 3. exclude genes w/ low RNA TPM, ribosome count > 100
        self.df = self.df[(self.df["TPM"] >= self.tpmLB)]
        # 4. filter out unmappable regions
        colNames = list(self.df.columns.values)
        bt = pbt.BedTool.from_dataframe(self.df).intersect(self.unMappable, v=True)
        self.df = bt.to_dataframe(names=colNames)
        # 5. remove genes that have too few codons remained (<10)
        geneFilter = self.df.groupby('gene').size().reset_index(name="remained")
        geneFilter["keep"] = np.where(geneFilter.remained>10, True, False)
        geneFilter = geneFilter[(geneFilter["keep"] == True)][['gene']]
        self.df = pd.merge(self.df, geneFilter, how = "inner", on = ["gene"])
        # 6. exclude stop codons: TAG, TAA, TGA
        self.df = self.df[((self.df["codon"] != 'TAG') & (self.df["codon"] != 'TAA') & (self.df["codon"] != 'TGA'))]
        numGenes = np.unique(self.df.gene).shape[0]
        # summary statistics of the data-frame
        sys.stderr.write("[status]\tTPM lower bound:\t" + str(self.tpmLB) + "\n")
        sys.stderr.write("[status]\tNumber of observations after filtering:\t" + str(self.df.shape[0]) + "\n")
        sys.stderr.write("[status]\tNumber of genes after filtering:\t" + str(numGenes) + "\n")

    def varScaling(self):
        """
        variable scaling
        :return:
        """
        # calculate log(TPM)
        tmp = self.df[["gene", "TPM"]].drop_duplicates()
        tmp['logTPM'] = np.log(tmp['TPM'])
        tmp["logTPM_scaled"] = (tmp['logTPM'] - np.mean(tmp['logTPM'])) / np.std(tmp['logTPM'])
        tmp.drop(["logTPM", "TPM"], axis=1, inplace=True)
        self.df = pd.merge(self.df, tmp, how='left', on=["gene"])
        # variable scaling
        self.df["avgProb_scaled"] = (self.df["avg_prob"] - np.mean(self.df["avg_prob"])) / np.std(self.df["avg_prob"])
        # remove unnecessary cols and NAs
        self.df = self.df.drop(["TPM", "avg_prob", "numCodons"], axis=1).dropna()
        self.df.to_csv(path_or_buf='filtered.txt', sep='\t', header=True, index=False, float_format='%.4f')
        #self.df = self.df[["chrom", "start", "end", "gene", "codon_idx", "gene_strand", "codon", "logTPM_scaled", "avgProb_scaled", "ribosome_count"]]

    def nbGlm(self):
        ## define model formula
        self.df = self.df[['ribosome_count', 'gene', 'codon', 'avgProb_scaled', 'logTPM_scaled']]
        formula = 'ribosome_count ~ C(gene) + C(codon) + avgProb_scaled'
        sys.stderr.write("[status]\tFormula: " + str(formula) + "\n")
        ## define model fitting options
        sovler = "IRLS" # "lbfgs"
        tolerence = 1e-4
        numIter = 100
        sys.stderr.write("[status]\tSolver: " + sovler + "\n")
        sys.stderr.write("[status]\tConvergence tolerance: " + str(tolerence) + "\n")
        sys.stderr.write("[status]\tMaxiter: " + str(numIter) + "\n")
        ## model fitting NegativeBinomial GLM
        sys.stderr.write("[status]\tModel: smf.glm(formula, self.df, family=sm.families.NegativeBinomial(), offset=self.df['logTPM_scaled']" + "\n")
        mod = smf.glm(formula, self.df, family=sm.families.NegativeBinomial(), offset=self.df['logTPM_scaled'])
        res = mod.fit(method=sovler, tol=tolerence, maxiter=numIter)
        ## print model output
        print(res.summary(), flush=True)

    def glmnetFit(self):
        X = np.array(pd.get_dummies(self.df[["gene", "codon"]]))
        y = np.array(self.df["ribosome_count"])
        lambdas = [1600, 400, 200] + list(range(100, 0, -1))
        m = PoissonNet(n_splits=0, alpha=0,lambda_path=lambdas)
        m = m.fit(X, y)
        coefs = m.coef_path_
        print(coefs)

    def checkrpy(self):
        # import R's utility package
        utils = rpackages.importr('utils')
        # select a mirror for R packages
        utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
        # R package names
        packnames = ('glmnet')
        # R vector of strings
        from rpy2.robjects.vectors import StrVector
        # Selectively install what needs to be install.
        # We are fancy, just because we can.
        names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
        if len(names_to_install) > 0:
            utils.install_packages(StrVector(packnames), dependencies=True, repos='http://cran.rstudio.com/')

    def rpyFit(self):
        n2r.activate()
        r = ro.r
        r.library('glmnet')
        #
        #from rpy2.robjects import pandas2ri
        #pandas2ri.activate()
        #
        self.df = pd.read_table("filtered.txt",  header=0)
        offsets = ro.FloatVector(list(self.df.logTPM_scaled.transpose()))
        X = np.array(pd.get_dummies(self.df[["gene", "codon"]]))
        y = np.array(self.df["ribosome_count"])
        y = ro.IntVector(list(y.transpose()))  # use factors
        fit = r['cv.glmnet'](X, y, family="poisson", alpha=0, parallel=True)
        #fit = r['cv.glmnet'](X, y, family="poisson", alpha=0, offset=offsets, parallel=True)
        fit = np.asanyarray(fit.rx2('lambda'))
        print(fit)
        #lambdas = [1600, 400, 200] + list(range(100, 0, -1))


## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input data-frame for modelling")
    parser.add_argument("-u", help="unmappable regions")
    parser.add_argument("-l", help="lower bound of RNA tpm", default=1, type=int)
    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    ## process the file if the input files exist
    if (args.i!=None) & (args.u!=None):
        ## Load the content of the table
        sys.stderr.write("[status]\tstatsmodel version:\t" + str(statsmodels.__version__) + "\n")
        sys.stderr.write("[status]\tReading the ribosome counts: " + str(args.i) + "\n")
        sys.stderr.write("[status]\tReading the un-mappable regions: " + str(args.u) + "\n")
        df_fn = args.i
        unmap_fn = args.u
        tpm_lb = args.l
        ## start model fitting
        sys.stderr.write("[execute]\tStart the modelling of TE" + "\n")
        mod = glmTE(df_fn, unmap_fn, tpm_lb)
        sys.stderr.write("[execute]\tCalculate lengths of chromosomes and filter the df" + "\n")
        mod.loadDat()
        mod.filterDf()
        mod.varScaling()
        sys.stderr.write("[execute]\tFitting the GLM" + "\n")
        #mod.checkrpy()
        #mod.rpyFit()
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
