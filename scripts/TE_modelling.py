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
#from memory_profiler import profile


class glmTE(object):
    ''' class to perform glm for TE estimation
    '''

    def __init__(self, fn, unMappableFn):
        self.fn = fn
        self.unMappableFn = unMappableFn

    def getLen(self):
        self.df = pd.read_table(self.fn,  header=0)
        tmp = self.df[["gene_name", "codon_idx"]]
        self.geneLen = tmp.groupby('gene_name').max().reset_index()
        self.geneLen.columns = ['gene_name', "gene_length"]
        self.geneLen["gene_length"] = self.geneLen["gene_length"] - 8

    def filterDf(self):
        ## get a summary of the starting df
        sys.stderr.write("[status]\tStarting number of observations:\t" + str(self.df.shape[0]) + "\n")
        sys.stderr.write("[status]\tNumber of variables:\t" + str(self.df.shape[1]-1) + "\n")
        ## rolling mean of 2' structure pairing probabilities
        tmp = self.df[["gene_name", "codon_idx", "pair_prob"]]
        rollingAvg = tmp.groupby("gene_name").rolling(window=11, center=True).mean().reset_index(0, drop=True)
        rollingAvg.columns = ["gene_name", "codon_idx", "avg_prob"]
        self.df = pd.merge(self.df, rollingAvg, how='left', on = ["gene_name", 'codon_idx'])
        self.df = self.df[["chrom", "asite_start", "asite_end", "gene_name", "codon_idx", "gene_strand", "codon", "TPM", "avg_prob", "ribosome_count"]]
        ## filter out unmappable regions
        self.unMappable = pbt.BedTool(self.unMappableFn)
        colNames = list(self.df.columns.values)
        bt = pbt.BedTool.from_dataframe(self.df)
        btMapped = bt.intersect(self.unMappable, v=True)
        self.df = btMapped.to_dataframe(names=colNames)
        ## filter 5/3 utr codons
        self.df = pd.merge(self.df, self.geneLen, how = "left", on = ["gene_name"])
        self.df = self.df[(self.df.codon_idx <= self.df.gene_length) & (self.df.codon_idx >=1)]
        ## remove genes that have too few codons remained (<10)
        tmp = self.df.groupby('gene_name').size().reset_index(name="remained")
        tmp["keep"] = np.where(tmp.remained>10, True, False)
        tmp = tmp[(tmp.keep == True)][['gene_name']]
        self.df = pd.merge(self.df, tmp, how = "inner", on = ["gene_name"])
        ## codons to exclude , stop: TAG, TAA, TGA
        self.df = self.df[((self.df.codon != 'TAG') & (self.df.codon != 'TAA') & (self.df.codon != 'TGA'))]
        ## define TPM lower bound, filter df, require TPM > 10, ribosome count > 100
        tpmLB = 10
        self.df = self.df[(self.df.TPM >= tpmLB)]
        numGenes = np.unique(self.df.gene_name).shape[0]
        ## group by gene names
        # grouped = self.df.groupby("gene_name")
        # ribosome_count_sum = np.array(grouped.aggregate(np.sum).TPM)
        ## calculate log(TPM)
        logTPM = np.log(self.df['TPM'])
        self.df['logTPM'] = logTPM
        tmp = self.df[["gene_name", "logTPM"]]
        tmp = tmp.drop_duplicates()
        tmp["logTPM_scaled"] = (tmp.logTPM - np.mean(tmp.logTPM)) / np.std(tmp.logTPM)
        tmp = tmp[["gene_name", "logTPM_scaled"]]
        self.df = pd.merge(self.df, tmp, how='left', on=["gene_name"])
        ## variable scaling
        self.df["avg_prob_scaled"] = (self.df.avg_prob - np.mean(self.df.avg_prob)) / np.std(self.df.avg_prob)
        self.df = self.df[["chrom", "asite_start", "asite_end", "gene_name", "codon_idx", "gene_strand", "codon", "logTPM_scaled", "avg_prob_scaled", "ribosome_count"]]
        ## remove missing values
        self.df = self.df.dropna()
        self.df.to_csv(path_or_buf='filtered.txt', sep='\t', header=True, index=False)
        ## print status
        sys.stderr.write("[status]\tTPM lower bound:\t" + str(tpmLB) + "\n")
        sys.stderr.write("[status]\tNumber of observations after filtering:\t" + str(self.df.shape[0]) + "\n")
        sys.stderr.write("[status]\tNumber of genes after filtering:\t" + str(numGenes) + "\n")

    def nbGlm(self):
        ## define model formula
        self.df = self.df[['ribosome_count', 'gene_name', 'codon', 'avg_prob_scaled', 'logTPM_scaled']]
        formula = 'ribosome_count ~ C(gene_name) + C(codon) + avg_prob_scaled'
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
        sys.stderr.write("[status]\tstatsmodel version:\t" + str(statsmodels.__version__) + "\n")
        sys.stderr.write("[status]\tReading the ribosome counts: " + str(args.i) + "\n")
        sys.stderr.write("[status]\tReading the un-mappable regions: " + str(args.u) + "\n")
        df_fn = args.i
        unmap_fn = args.u
        ## start model fitting
        sys.stderr.write("[execute]\tStart the modelling of TE" + "\n")
        mod = glmTE(df_fn, unmap_fn)
        sys.stderr.write("[execute]\tCalculate lengths of chromosomes and filter the df" + "\n")
        mod.getLen()
        mod.filterDf()
        sys.stderr.write("[execute]\tFitting the GLM" + "\n")
        mod.nbGlm()
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
