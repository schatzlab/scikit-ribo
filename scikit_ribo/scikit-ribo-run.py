#!/usr/bin/env python

# ----------------------------------------
# scikit-ribo
# ----------------------------------------
# main wrapper
# ----------------------------------------
# author: Han Fang
# contact: hanfang.cshl@gmail.com
# website: hanfang.github.io
# date: 3/28/2017
# ----------------------------------------

from __future__ import print_function, division
import os
import sys
import argparse
import pybedtools as pbt
import pysam
import pandas as pd
import numpy as np
import csv
import errno
from datetime import datetime
import scikit_ribo
from bam_process import BamProcess
from asite_predict import VisualizeAsite
from asite_predict import PredictAsite
from model_te import ModelTE
from colorama import init
from termcolor import cprint
from pyfiglet import figlet_format
init(strip=not sys.stdout.isatty())


def warn(*args, **kwargs):
    """
    skip warnings
    :param args:
    :param kwargs:
    :return:
    """
    pass


def log_status(bam, directory, prefix, out, cv):
    """
    Load data, log status and start
    :param bam: str, input bam file
    :param directory: str, folder to pre-built files
    :param prefix: str, prefix for pre-built files
    :param out: str, output folder, usually sample id
    :param cv: bool, cross-validation flag
    :return: str * 4
    """
    print("[status]\tStarted scikit-ribo at " + str(datetime.now()), file=sys.stderr)
    dir_pre = directory + "/" + prefix
    start, cds = dir_pre + ".start.bed", dir_pre + ".cds.bed"
    posIdx, nt = dir_pre + ".pos_ranges.txt", dir_pre + '.nt_table.txt'
    lambda_min = 0.13 if not cv else None
    # create folder
    cmd = 'mkdir -p ' + out
    os.system(cmd)
    print("[status]\tProcessing the input bam file: " + bam, file=sys.stderr)
    print("[setting]\tOutput path: " + out, file=sys.stderr)
    print("[setting]\tOutput path: " + out, flush=True)
    print("[status]\tReading the start codon BED file: " + start, file=sys.stderr)
    print("[status]\tReading the open reading frame codon BED file: " + cds, file=sys.stderr)
    print("[status]\tReading the position-phase file: " + posIdx, file=sys.stderr)
    print("[status]\tReading the nt table file: " + nt, file=sys.stderr)
    if lambda_min:
        print("[setting]\tpre-defined lambda: " + str(lambda_min), file=sys.stderr)
    else:
        print("[setting]\tCross-validation on glmnet enabled", file=sys.stderr)
    sys.stderr.flush()
    return start, cds, posIdx, nt, lambda_min


def aln_module(bam, mapQual, minReadLen, maxReadLen, rele, out, start, cds, posIdx, nt):
    """
    Module of alignment process
    :param bam: str, bam
    :param mapQual: int, 20
    :param minReadLen: int, 15
    :param maxReadLen: int, 35
    :param rele: bool, False
    :param out: str, output directory
    :param start:
    :param cds:
    :param posIdx:
    :param nt:
    :return: pandas DataFrame * 2
    """
    print("[setting]\tKeep reads with length range: [" + str(minReadLen) + "," + str(maxReadLen) + "]", file=sys.stderr)
    print("[setting]\tMinimum mapq allowed: " + str(mapq), file=sys.stderr)
    aln = BamProcess(bam, mapQual, minReadLen, maxReadLen, rele, out, start, cds, posIdx, nt)
    print("[execute]\tFiltering out overlapping regions", file=sys.stderr)
    aln.filterRegion()
    print("[execute]\tImporting the position ranges", file=sys.stderr)
    aln.posIndex()
    print("[execute]\tFiltering out un-reliable alignment from bam files", file=sys.stderr)
    aln.filterBam()
    print("[execute]\tCreate training dataframe at " + str(datetime.now()), file=sys.stderr)
    trainingData = aln.makeTrainingData()
    print("[execute]\tCreate CDS dataframe at " + str(datetime.now()), file=sys.stderr)
    cdsData = aln.makeCdsData()
    print("[status]\tBam processing finished at " + str(datetime.now()), file=sys.stderr)
    sys.stderr.flush()
    return trainingData, cdsData


def asite_module(trainingData, cdsData, rele, prefix, out, directory):
    """
    Module of predicting a-site
    :param trainingData: pandas DataFrame
    :param cdsData: pandas DataFrame
    :param rele: bool, False
    :param prefix: str
    :param out: str, output folder
    :param directory: str, input folder
    :return: dataFrame
    """
    # silent
    import warnings
    warnings.warn = warn
    # plot
    print("[execute]\tPlotting the a-site location distribution", file=sys.stderr)
    trainingDist = VisualizeAsite(trainingData, rele, out)
    trainingDist.plot()
    # predict a-site
    print("[execute]\tstart the process of a-site prediction", file=sys.stderr)
    classifier = "rf"
    model = PredictAsite(trainingData, cdsData, classifier, rele, prefix, out, directory)
    print("[execute]\tperform model training and cross validation on the training data", file=sys.stderr)
    model.rfFit()
    print("[execute]\tplotting the bar plot of the feature importance", file=sys.stderr)
    model.rfImportance()
    print("[execute]\tplot roc curve based on cross validation", file=sys.stderr)
    model.rocCurve()
    print("[execute]\tpredicting the a-site from the cds regions", file=sys.stderr)
    model.rfPredict()
    print("[execute]\tlocalize the a-site codon and create coverage df", file=sys.stderr)
    dataFrame = model.recoverAsite()
    sys.stderr.flush()
    return dataFrame


def te_module(dataFrame, unmap, out, lambda_min):
    """
    Module of inferring translation efficiency
    :param dataFrame: dataFrame
    :param unmap: str, un-mappable region
    :param out: str
    :return: None
    """
    # silent
    import warnings
    warnings.warn = warn
    if unmap:
        print("[status]\tReading the un-mappable regions: " + str(args.u), file=sys.stderr)
    else:
        print("[status]\tThe list of un-mappable regions was not provided", file=sys.stderr)
    tpmLB = 1
    # start model fitting
    print("[execute]\tStart the modelling of translation efficiency (TE)", file=sys.stderr)
    mod = ModelTE(dataFrame, unmap, out, tpmLB)
    print("[execute]\tLoading data", file=sys.stderr)
    mod.loadDat()
    print("[execute]\tFiltering the df", file=sys.stderr)
    mod.filterDf()
    print("[execute]\tScaling the variables", file=sys.stderr)
    mod.varScaling()
    print("[execute]\tFitting the GLM", file=sys.stderr)
    X, y, offsets, numCodons, numGenes, varsNames = mod.glmnetArr()
    mod.glmnetFit(X, y, offsets, numCodons, numGenes, varsNames, lambda_min)
    sys.stderr.flush()


def scikit_ribo(bam, directory, prefix, mapQual, minReadLen, maxReadLen, rele, out, unmap, cv):
    """
    Wrapper master function
    :param bam: str, bam file location
    :param directory: str, input directory
    :param prefix: str, prefix
    :param mapQual: int, 20
    :param minReadLen: int, 15
    :param maxReadLen: int, 35
    :param rele: bool, False
    :param out: str, output folder
    :param unmap: str, un-mappable file
    :return: None
    """
    # log status
    start, cds, posIdx, nt, lambda_min = log_status(bam, directory, prefix, out, cv)
    # load data
    trainingData, cdsData = aln_module(bam, mapQual, minReadLen, maxReadLen, rele, out, start, cds, posIdx, nt)
    # predict a-site
    dataFrame = asite_module(trainingData, cdsData, rele, prefix, out, directory)
    # model te
    te_module(dataFrame, unmap, out, lambda_min)
    # Finish
    print("[status]\tScikit-ribo finished at " + str(datetime.now()), file=sys.stderr)
    sys.stderr.flush()

# ----------------------------------------
# parse input arguments
# ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input bam file")
    parser.add_argument("-f", help="path to the Folder of BED/index files generated by the pre-processing module")
    parser.add_argument("-p", help="Prefix for BED/index files")
    parser.add_argument("-q", help="minimum mapQ allowed, Default: 20", default=20, type=int)
    parser.add_argument("-s", help="Shortest read length allowed, Default: 10", default=10, type=int)
    parser.add_argument("-l", help="Longest read length allowed, Default: 35", default=35, type=int)
    parser.add_argument("-c", help="enable cross validation for glmnet", action='store_true')
    parser.add_argument("-r", help="setting this flag will enable the RelE mode", action='store_true')
    parser.add_argument("-u", help="Un-mappable regions", default=None)
    parser.add_argument("-o", help="Output path, recommend using the sample id")
    # check if there is any argument
    if len(sys.argv[1:]) == 0:
        cprint(figlet_format('scikit-ribo', font='ogre'), 'green', attrs=['bold'])
        parser.print_help()
        parser.exit()
    if len(sys.argv) <= 1:
        cprint(figlet_format('scikit-ribo', font='ogre'), 'green', attrs=['bold'])
        parser.print_usage()
        sys.exit(1)
    else:
        cprint(figlet_format('scikit-ribo', font='ogre'), 'green', attrs=['bold'])
        args = parser.parse_args()
    # process the file if the input files exist
    if args.i != None and args.f != None and args.p != None and args.o != None:
        # specify inputs
        input_bam = args.i
        folder = args.f
        pre = args.p
        mapq = args.q
        minRL, maxRL = args.s, args.l
        cv = args.c
        RelE = args.r
        output = args.o
        unmap_fn = args.u
        scikit_ribo(input_bam, folder, pre, mapq, minRL, maxRL, RelE, output, unmap_fn, cv)
    else:
        print("[error]\tmissing argument", file=sys.stderr)
        parser.print_usage()
