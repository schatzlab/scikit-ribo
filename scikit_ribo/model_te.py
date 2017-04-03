#!/usr/bin/env python

# ----------------------------------------
# scikit-ribo
# ----------------------------------------
# a module for fitting the GLM of TE
# ----------------------------------------
# author: Han Fang
# contact: hanfang.cshl@gmail.com
# website: hanfang.github.io
# date: 3/31/2016
# ----------------------------------------

from __future__ import print_function
import sys, os
import argparse
import scipy
import pandas as pd
import numpy as np
import pybedtools as pbt
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glmnet_py.glmnet as glmnet
import glmnet_py.dataprocess as dataprocess
import glmnet_py.glmnetCoef as glmnetCoef
import glmnet_py.cvglmnet as cvglmnet
import glmnet_py.cvglmnetCoef as cvglmnetCoef
import glmnet_py.cvglmnetPlot as cvglmnetPlot
import glmnet_py.cvglmnetPredict as cvglmnetPredict


class ModelTE(object):
    """ perform glm for TE estimation
    """
    def __init__(self, df=None, unMappableFn=None, output=None, tpmLB=1):
        """
        :param fn:
        :param unMappableFn:
        :param tpmLB:
        """
        self.df = df
        self.unMappableFn = unMappableFn
        self.output = output
        self.tpmLB = tpmLB
        self.unMappable = None

    def loadDat(self):
        """
        loading data: self.df, self.unMappable
        :return: None
        """
        # read the df and un-mappable regions
        # self.df = pd.read_table(self.fn,  header=0)
        if self.unMappableFn:
            self.unMappable = pbt.BedTool(self.unMappableFn)
        # compute the codon length
        codons = self.df[["gene", "codon_idx"]]
        codonLen = codons.groupby('gene').max().reset_index()
        codonLen.columns = ['gene', "numCodons"]
        codonLen["numCodons"] -= 8  # removed the appended 3 utr codons
        self.df = pd.merge(self.df, codonLen, how="left", on=["gene"])
        self.df.rename(columns={'pair_prob': 'downstreamSL'}, inplace=True)
        self.df = self.df.dropna()
        self.df = self.df[(self.df["ribosome_count"] > 0)]

    def filterDf(self):
        """
        filter dataframe, self.df
        :return: None
        """
        # get a summary of the starting df
        sys.stderr.write("[status]\tStarting number of observations:\t" + str(self.df.shape[0]) + "\n")
        sys.stderr.write("[status]\tTPM lower bound:\t" + str(self.tpmLB) + "\n")
        # 1. exclude 5/3 utr codons, stop codons: TAG, TAA, TGA
        self.df = self.df[(self.df["codon_idx"] <= self.df["numCodons"]) & (self.df["codon_idx"] > 0)]
        self.df = self.df[((self.df["codon"] != 'TAG') & (self.df["codon"] != 'TAA') & (self.df["codon"] != 'TGA'))]
        # 2. exclude genes w/ low Ribo TPM
        subCols = ["gene", "codon_idx", "ribosome_count"]
        genes = self.df[subCols].groupby("gene").sum().ribosome_count.reset_index(name="riboCnt")
        codonLen = self.df[["gene", "numCodons"]].drop_duplicates(subset=["gene"])
        genes = pd.merge(genes, codonLen, on="gene")
        genes["rpkm"] = genes["riboCnt"] * 1000000000 / (np.sum(genes["riboCnt"]) * genes["numCodons"] * 3)
        genes["tpm"] = genes["rpkm"] * 1000000 / np.sum(genes["rpkm"])
        genesToKeep = genes[(genes["tpm"] >= self.tpmLB)][["gene"]]
        self.df = pd.merge(self.df, genesToKeep, on="gene")
        # 3. exclude genes w/ low RNA TPM, ribosome count > 100
        self.df = self.df[(self.df["TPM"] >= self.tpmLB)]
        # 4. filter out unmappable regions
        if self.unMappable:
            colNames = list(self.df.columns.values)
            bt = pbt.BedTool.from_dataframe(self.df).intersect(self.unMappable, v=True)
            self.df = bt.to_dataframe(names=colNames)
        # 5. remove genes that have too few codons remained (<10)
        geneFilter = self.df.groupby('gene').size().reset_index(name="remained")
        geneFilter["keep"] = np.where(geneFilter.remained > 10, True, False)
        geneFilter = geneFilter[(geneFilter["keep"] == True)][['gene']]
        self.df = pd.merge(self.df, geneFilter, how="inner", on=["gene"])
        # summary statistics of the data-frame
        numGenesRemained = np.unique(self.df.gene).shape[0]
        sys.stderr.write("[status]\tNumber of observations after filtering:\t" + str(self.df.shape[0]) + "\n")
        sys.stderr.write("[status]\tNumber of genes after filtering:\t" + str(numGenesRemained) + "\n")

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
        # remove unnecessary cols and NAs
        self.df = self.df.drop(["TPM", "numCodons"], axis=1).dropna()

    def sparseDfToCsc(self, df):
        """
        convert a pandas sparse dataframe to a scipy sparse array
        :param df: pandas dataframe
        :return: sparse array
        """
        columns = df.columns
        dat, rows = map(list, zip(
            *[(df[col].sp_values - df[col].fill_value, df[col].sp_index.to_int_index().indices) for col in columns]))
        cols = [np.ones_like(a) * i for (i, a) in enumerate(dat)]
        datF, rowsF, colsF = np.concatenate(dat), np.concatenate(rows), np.concatenate(cols)
        arr = sparse.coo_matrix((datF, (rowsF, colsF)), df.shape, dtype=np.float64)
        return arr.tocsc()

    def glmnetArr(self):
        """

        :return:
        """
        # get col names
        numGenesToInfer, numCodonsToInfer = len(self.df.gene.unique()), len(self.df.codon.unique())
        cateVars = pd.get_dummies(self.df[["gene", "codon"]], sparse=True)
        cateVarsNames = list(cateVars.columns.values)
        varsNamesLst = cateVarsNames + ["secondary_structure"]
        # prepare input arrays
        cateVars = dataprocess().sparseDf(cateVars)
        downstreamSL = np.array(self.df["downstreamSL"].as_matrix().reshape(int(len(self.df["downstreamSL"])), 1),
                                dtype=np.float64)
        sparseX = sparse.hstack([cateVars, downstreamSL], format="csc", dtype=np.float64)
        sparseY = np.array(self.df["ribosome_count"], dtype=np.float64)
        offsetsArr = np.array(self.df["logTPM_scaled"], dtype=np.float64).reshape(int(len(self.df["logTPM_scaled"])), 1)
        return sparseX, sparseY, offsetsArr, numCodonsToInfer, numGenesToInfer, varsNamesLst

    def glmnetFit(self, X, y, offsets, numCodons, numGenes, varsNames, lambda_min):
        """

        :param X:
        :param y:
        :param offsets:
        :param numCodons:
        :param numGenes:
        :param varsNames:
        :param lambda_min:
        :return:
        """
        # fit the model
        if not lambda_min:
            fit = cvglmnet(x=X.copy(), y=y.copy(), family='poisson',
                           offset=offsets, alpha=0, parallel=True, lambda_min=np.array([0]))
            coefs = cvglmnetCoef(fit, s=fit['lambda_min'])  # lambda_min lambda_1se
        else:
            fit = glmnet(x=X.copy(), y=y.copy(), family='poisson',
                         offset=offsets, alpha=0, lambda_min=np.array([0]))
            coefs = glmnetCoef(fit, s=scipy.float64([lambda_min]))
        # parse and scale coefficients
        intercept = coefs[0][0]
        geneBetas = pd.DataFrame([[varsNames[i-1].split("_")[1], coefs[i][0]]
                                  for i in range(1, numGenes+1)], columns=["gene", "beta"])
        geneBetas["log2_TE"] = (geneBetas["beta"] - np.median(geneBetas["beta"])) / np.log(2)
        geneBetas.drop(["beta"], inplace=True, axis=1)
        codonBetas = pd.DataFrame([[varsNames[i-1].split("_")[1], coefs[i][0]]
                                   for i in range(numGenes+1, numGenes + numCodons + 1)], columns=["codon", "beta"])
        codonBetas["log_codon_dwell_time"] = (codonBetas["beta"] - np.median(codonBetas["beta"]))
        codonBetas["codon_dwell_time"] = np.exp(codonBetas["log_codon_dwell_time"])
        codonBetas.drop(["beta", "log_codon_dwell_time"], inplace=True, axis=1)
        downstreamSLBeta = coefs[numGenes + numCodons + 1][0]
        #  export to local
        geneBetas.to_csv(path_or_buf=self.output + '/genesTE.csv', sep='\t',
                         header=True, index=False, float_format='%.4f')
        codonBetas.to_csv(path_or_buf=self.output + '/codons.csv', sep='\t',
                          header=True, index=False, float_format='%.4f')
        # print results
        if lambda_min:
            sys.stderr.write("[results]\tpre-defined lambda: " + str(lambda_min) + "\n")
        else:
            sys.stderr.write("[results]\tlambda that gives minimum mean cv error: " + str(fit['lambda_min']) + "\n")
            sys.stderr.write("[results]\tlambda 1 se away: " + str(fit['lambda_1se']) + "\n")
        sys.stderr.write("[results]\tintercept: " + str(intercept) + "\n")
        sys.stderr.write("[results]\tbetas for 2' structure windows: " + str(downstreamSLBeta) + "\n")
        # plot
        if not lambda_min:
            plt.figure()
            cvglmnetPlot(fit)
            plt.gcf()
            plt.savefig(self.output + "/" + "lambda_cv.pdf")
            plt.clf()

    def nbGlm(self):
        """

        :return:
        """
        import statsmodels
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
        sys.stderr.write("[status]\tstatsmodels version:\t" + str(statsmodels.__version__) + "\n")
        # define model formula
        self.df = self.df[['ribosome_count', 'gene', 'codon', 'avgProb_scaled', 'logTPM_scaled']]
        formula = 'ribosome_count ~ C(gene) + C(codon) + avgProb_scaled'
        sys.stderr.write("[status]\tFormula: " + str(formula) + "\n")
        # define model fitting options
        solver = "IRLS"  # "lbfgs"
        tolerance = 1e-4
        numIter = 100
        sys.stderr.write("[status]\tSolver: " + solver + "\n")
        sys.stderr.write("[status]\tConvergence tolerance: " + str(tolerance) + "\n")
        sys.stderr.write("[status]\tMaxiter: " + str(numIter) + "\n")
        # model fitting NegativeBinomial GLM
        sys.stderr.write("[status]\tModel: smf.glm(formula, self.df, family=sm.families.NegativeBinomial(), "
                         "offset=self.df['logTPM_scaled']" + "\n")
        model = smf.glm(formula, self.df, family=sm.families.NegativeBinomial(), offset=self.df['logTPM_scaled'])
        res = model.fit(method=solver, tol=tolerance, maxiter=numIter)
        # print model output
        print(res.summary(), flush=True)


# main
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input data-frame for modelling")
    parser.add_argument("-u", help="un-mappable regions", default=None)
    parser.add_argument("-l", help="lower bound of RNA tpm", default=1, type=int)
    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # process the file if the input files exist
    if args.i is not None:
        # Load the content of the table
        sys.stderr.write("[status]\tReading the ribosome counts: " + str(args.i) + "\n")
        if args.u:
            sys.stderr.write("[status]\tReading the un-mappable regions: " + str(args.u) + "\n")
        else:
            sys.stderr.write("[status]\tThe list of un-mappable regions was not provided""\n")
        df_fn = args.i
        unmap_fn = args.u
        tpm_lb = args.l
        lambda_min = None
        # start model fitting
        sys.stderr.write("[execute]\tStart the modelling of TE" + "\n")
        mod = ModelTE(df_fn, unmap_fn, tpm_lb)
        sys.stderr.write("[execute]\tLoading data" + "\n")
        mod.loadDat()
        sys.stderr.write("[execute]\tFiltering the df" + "\n")
        mod.filterDf()
        sys.stderr.write("[execute]\tScaling the variables" + "\n")
        mod.varScaling()
        sys.stderr.write("[execute]\tFitting the GLM" + "\n")
        X, y, offsets, numCodons, numGenes, varsNames = mod.glmnetArr()
        mod.glmnetFit(X, y, offsets, numCodons, numGenes, varsNames, lambda_min)
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
