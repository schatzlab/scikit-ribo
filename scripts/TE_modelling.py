#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for normalizing read counts
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

from __future__ import print_function
import sys
import argparse
import numpy as np
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
        y, X = dmatrices('ribosome_count ~ gene_name + codon + TPM + pair_prob', data=self.df, return_type='dataframe')
        print(y)
        print(X)

        ## model fitting
        #mod = smf.glm(y, X, family=sm.families.NegativeBinomial())
        #mod.fit(skip_hessian=True)
        #print ( mod.summary() )


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
        df_fn = args.i
        model_hdl = GLM_TE(df_fn)
        model_hdl.nb_glm()

        # cnt_array=pd.io.parsers.read_table(args.i,header=None,names=('count', 'gene', 'codon', 'rna_reads', 'rpf_reads') )
        # cnt_array = np.genfromtxt(args.i, delimiter='\t', dtype=None, names=('count', 'gene', 'codon', 'rna_reads', 'rpf_reads','len','tpm','pair_prob'))

        ## Calculate the normalization factor for each sample
        #nb_glm(cnt_array)

    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
