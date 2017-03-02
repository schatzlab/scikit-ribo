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
import sys
import argparse
import scipy
import pandas as pd
import numpy as np
from patsy import dmatrices
import pybedtools as pbt
from statsmodels.base.model import GenericLikelihoodModel
from scipy import sparse
import glmnet_python.glmnet as glmnet
import glmnet_python.dataprocess as dataprocess
import glmnet_python.glmnetCoef as glmnetCoef
import glmnet_python.cvglmnet as cvglmnet
import glmnet_python.cvglmnetCoef as cvglmnetCoef
import glmnet_python.cvglmnetPlot as cvglmnetPlot
import glmnet_python.cvglmnetPredict as cvglmnetPredict
#import statsmodels
#import statsmodels.api as sm
#import statsmodels.formula.api as smf


class modelTE(object):
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
        self.unMappable = None
        self.df = None

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
        codonLen["numCodons"] -= 8 # removed the appended 3 utr files
        self.df = pd.merge(self.df, codonLen, how="left", on=["gene"])
        # rolling mean of 2' structure pairing probabilities
        paring = self.df[["gene", "codon_idx", "pair_prob"]]
        rollingAvg = paring.groupby("gene").rolling(window=11, center=True).mean().reset_index(0, drop=True)
        rollingAvg.columns = ["gene", "codon_idx", "avg_prob"]
        # merge 2' prob
        self.df = pd.merge(self.df, rollingAvg, how='left', on=["gene", 'codon_idx']).drop(["pair_prob"], axis=1)
        self.df = self.df.dropna()

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
        gene["tpm"] = gene["rpkm"] * 1000000 / np.sum(gene["rpkm"])
        gene = gene[(gene["tpm"] >= self.tpmLB)][["gene"]]
        self.df = pd.merge(self.df, gene, on="gene")
        # 3. exclude genes w/ low RNA TPM, ribosome count > 100
        self.df = self.df[(self.df["TPM"] >= self.tpmLB)]
        # 4. filter out unmappable regions
        colNames = list(self.df.columns.values)
        bt = pbt.BedTool.from_dataframe(self.df).intersect(self.unMappable, v=True)
        self.df = bt.to_dataframe(names=colNames)
        # 5. remove genes that have too few codons remained (<10)
        geneFilter = self.df.groupby('gene').size().reset_index(name="remained")
        geneFilter["keep"] = np.where(geneFilter.remained > 10, True, False)
        geneFilter = geneFilter[(geneFilter["keep"] == True)][['gene']]
        self.df = pd.merge(self.df, geneFilter, how="inner", on=["gene"])
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

    def nbGlm(self):
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

    def sparseDfToCsc(self, df):
        columns = df.columns
        dat, rows = map(list, zip(
            *[(df[col].sp_values - df[col].fill_value, df[col].sp_index.to_int_index().indices) for col in columns]))
        cols = [np.ones_like(a) * i for (i, a) in enumerate(dat)]
        datF, rowsF, colsF = np.concatenate(dat), np.concatenate(rows), np.concatenate(cols)
        arr = sparse.coo_matrix((datF, (rowsF, colsF)), df.shape, dtype=np.float64)
        return arr.tocsc()

    def glmnetFit(self):
        # get col names
        self.df = pd.read_table("filtered.txt",  header=0)
        numGenes, numCodons = len(self.df.gene.unique()), len(self.df.codon.unique())
        cateVars = pd.get_dummies(self.df[["gene", "codon"]], sparse=True)
        cateVarsNames = list(cateVars.columns.values)
        varsNames = cateVarsNames + ["secondary_structure"]
        # prepare input arrays
        avgProb = self.df["avgProb_scaled"]
        avgProbArr = np.array(avgProb.as_matrix().reshape(len(avgProb), 1), dtype=np.float64)
        cateVars = dataprocess().sparseDf(cateVars)
        X = sparse.hstack([cateVars, avgProbArr], format="csc", dtype=np.float64)
        y = np.array(self.df["ribosome_count"], dtype=np.float64)
        offsets = np.array(self.df["logTPM_scaled"], dtype=np.float64).reshape(len(self.df["logTPM_scaled"]), 1)
        # fit the model
        cvfit = cvglmnet(x=X.copy(), y=y.copy(), family='poisson', offset=offsets, alpha=0, parallel=True, lambda_min=np.array([0]))
        cvcoefs = cvglmnetCoef(cvfit, s=cvfit['lambda_min'])
        # parse and scale coefficients
        intercept = cvcoefs[0][0]
        geneBetas = pd.DataFrame([[varsNames[i-1].split("_")[1], cvcoefs[i][0]] for i in range(1, numGenes)], columns=["gene", "beta"])
        geneBetas["log2_TE"] = (geneBetas["beta"] - np.median(geneBetas["beta"])) / np.log(2)
        geneBetas.drop(["beta"], inplace=True, axis=1)
        codonBetas = pd.DataFrame([[varsNames[i-1].split("_")[1], cvcoefs[i][0]] for i in range(numGenes, numGenes + numCodons)], columns=["codon", "beta"])
        codonBetas["log_codon_elongation_rate"] = (codonBetas["beta"] - np.median(codonBetas["beta"]))
        codonBetas["codon_elongation_rate"] = np.exp(codonBetas["log_codon_elongation_rate"])
        codonBetas.drop(["beta", "log_codon_elongation_rate"], inplace=True, axis=1)
        pairProbBeta = cvcoefs[-1][0]
        # export to local
        geneBetas.to_csv(path_or_buf='genesTE.csv', sep='\t', header=True, index=False, float_format='%.4f')
        codonBetas.to_csv(path_or_buf='codons.csv', sep='\t', header=True, index=False, float_format='%.4f')
        # parse prefix
        base, dir = os.path.basename(self.fn), os.path.dirname(self.fn)
        prefix = dir + "/" + os.path.splitext(base)[0]
        # plot
        plt.figure()
        cvglmnetPlot(cvfit)
        plt.gcf()
        plt.savefig(prefix + ".lambda_cv.pdf")
        plt.clf()
        # print results
        sys.stderr.write("[results]\tlambda that gives minimum mean cv error: " + str(cvfit['lambda_min']) + "\n")
        sys.stderr.write("[results]\tintercept: " + str(intercept) + "\n")
        sys.stderr.write("[results]\tbeta for secondary structure pairing probability: " + str(pairProbBeta) + "\n")
        # lambdas = np.concatenate((np.arange(1600, 100, -400), np.arange(100, 1, -1), np.arange(1, 0, -0.01), np.array([0])))
        # fit = glmnet(x=X.copy(), y=y.copy(), family='poisson', offset=offsets, alpha=0, lambda_min=np.array([0]))# lambdau=lambdas)
        # coefs = glmnetCoef(fit, s=scipy.float64([0]))

# main
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input data-frame for modelling")
    parser.add_argument("-u", help="un-mappable regions")
    parser.add_argument("-l", help="lower bound of RNA tpm", default=1, type=int)
    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # process the file if the input files exist
    if (args.i is not None) & (args.u is not None):
        # Load the content of the table
        sys.stderr.write("[status]\tReading the ribosome counts: " + str(args.i) + "\n")
        sys.stderr.write("[status]\tReading the un-mappable regions: " + str(args.u) + "\n")
        df_fn = args.i
        unmap_fn = args.u
        tpm_lb = args.l
        # start model fitting
        sys.stderr.write("[execute]\tStart the modelling of TE" + "\n")
        mod = modelTE(df_fn, unmap_fn, tpm_lb)
        sys.stderr.write("[execute]\tLoading data" + "\n")
        #mod.loadDat()
        sys.stderr.write("[execute]\tFiltering the df" + "\n")
        #mod.filterDf()
        sys.stderr.write("[execute]\tScaling the variables" + "\n")
        #mod.varScaling()
        sys.stderr.write("[execute]\tFitting the GLM" + "\n")
        mod.glmnetFit()
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
