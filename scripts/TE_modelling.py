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

class GLM_TE:
    ''' class to perform glm for TE estimation
    '''
    def __init__(self, fn):
        self.fn = fn
        self.df = pd.read_table(self.fn,  header=0)

    def nb_glm(self):
        ## convert TPM to log(TPM)
        self.df['TPM'] = np.log(self.df['TPM'])

        ## construct feature array and outcome array
        formula = 'ribosome_count ~ gene_name + codon + TPM + pair_prob' 
        # y, X = dmatrices('ribosome_count ~ gene_name + codon + TPM + pair_prob', data=self.df, return_type='dataframe')

        ## model fitting
        mod = smf.glm(formula, self.df, family=sm.families.NegativeBinomial())
        #for sovler in ["newton", "nm", "bfgs", "lbfgs", "cg", "ncg"]:
        sovler = "bfgs"
        print("[status]\tsolver: " + sovler, flush=True)
        res = mod.fit(method=sovler)
        print ( res.summary() )


## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input the data-frame for modelling")
    parser.add_argument("-o", help="output result")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.o!=None):

        ## Load the content of the table
        print("[status]\tstatsmodel version:\t"+ str(statsmodels.__version__), flush=True)
        print("[status]\tReading the file: " + str(args.i), flush=True)
        df_fn = args.i
        print("[execute]\tStart the modelling of TE", flush=True)
        model_hdl = GLM_TE(df_fn)
        print("[execute]\tFitting the GLM", flush=True)
        model_hdl.nb_glm()


    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
