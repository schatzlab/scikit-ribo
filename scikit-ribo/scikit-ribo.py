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
from bam_process import BamProcess
from asite_predict import VisualizeAsite
from asite_predict import PredictAsite
from model_te import ModelTE
from colorama import init
init(strip=not sys.stdout.isatty())
from termcolor import cprint
from pyfiglet import figlet_format


def log_status(inBam, folder, pre, output):
    """
    load data, log status and start
    :param inBam:
    :param output:
    :param pre:
    :return:
    """
    dir_pre = folder + "/" + pre
    start, cds = dir_pre + ".start.bed", dir_pre + ".cds.bed"
    posIdx, nt = dir_pre + ".pos_ranges.txt", dir_pre + '.nt_table.txt'
    # create folder
    cmd = 'mkdir -p ' + output
    os.system(cmd)
    print("[status]\tProcessing the input bam file: " + inBam, file=sys.stderr)
    print("[status]\tOutput path: " + output, file=sys.stderr)
    print("[status]\tReading the start codon BED file: " + start, file=sys.stderr)
    print("[status]\tReading the open reading frame codon BED file: " + cds, file=sys.stderr)
    print("[status]\tReading the position-phase file: " + posIdx, file=sys.stderr)
    print("[status]\tReading the nt table file: " + nt, file=sys.stderr)
    sys.stderr.flush()
    return start, cds, posIdx, nt

def aln_module(inBam, mapq, minRL, maxRL, RelE, output, start, cds, posIdx, nt):
    """

    :param inBam:
    :param mapq:
    :param output:
    :param minRL:
    :param maxRL:
    :param start:
    :param cds:
    :param posIdx:
    :param RelE:
    :param nt:
    :return:
    """
    print("[execute]\tKeep reads with length range: [" + str(minRL) + "," + str(maxRL) + "]", file=sys.stderr)
    print("[execute]\tMinimum mapq allowed: " + str(mapq), file=sys.stderr)
    aln = BamProcess(inBam, mapq, minRL, maxRL, RelE, output, start, cds, posIdx, nt)
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

def asite_module(trainingData, cdsData, RelE, pre, output):
    """

    :param trainingData:
    :param cdsData:
    :param RelE:
    :param pre:
    :param output:
    :return:
    """
    # silent
    def warn(*args, **kwargs):
        pass
    import warnings
    warnings.warn = warn
    # plot
    print("[execute]\tPlotting the a-site location distribution", file=sys.stderr)
    trainingDist = VisualizeAsite(trainingData, RelE, output)
    trainingDist.plot()
    # predict a-site
    print("[execute]\tstart the process of a-site prediction", file=sys.stderr)
    classifier = "rf"
    model = PredictAsite(trainingData, cdsData, classifier, RelE, pre, output)
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

def te_module(dataFrame, unmap_fn, output):
    """

    :param dataFrame:
    :param unmap_fn:
    :param output:
    :return:
    """
    # silent
    def warn(*args, **kwargs):
        pass
    import warnings
    warnings.warn = warn
    if unmap_fn:
        print("[status]\tReading the un-mappable regions: " + str(args.u), file=sys.stderr)
    else:
        print("[status]\tThe list of un-mappable regions was not provided", file=sys.stderr)
    lambda_min = None
    tpmLB = 1
    # start model fitting
    print("[execute]\tStart the modelling of translation efficiency (TE)", file=sys.stderr)
    mod = ModelTE(dataFrame, unmap_fn, output, tpmLB)
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

def scikit_ribo(inBam, folder, pre, mapq, minRL, maxRL, RelE, output, unmap_fn):
    """

    :param inBam:
    :param pre:
    :param mapq:
    :param minRL:
    :param maxRL:
    :param RelE:
    :param output:
    :param unmap_fn:
    :return:
    """
    print("[status]\tStarted scikit-ribo at " + str(datetime.now()), file=sys.stderr)
    sys.stderr.flush()
    # log status
    start, cds, posIdx, nt = log_status(inBam, folder, pre, output)
    # load data
    trainingData, cdsData = aln_module(inBam, mapq, minRL, maxRL, RelE, output, start, cds, posIdx, nt)
    # predict a-site
    dataFrame = asite_module(trainingData, cdsData, RelE, pre, output)
    # model te
    te_module(dataFrame, unmap_fn, output)
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
    if (args.i != None and args.f != None and args.p != None and args.o != None):
        # specify inputs
        inBam = args.i
        folder = args.f
        pre = args.p
        mapq = args.q
        minRL, maxRL = args.s, args.l
        RelE = args.r
        output = args.o
        unmap_fn = args.u
        scikit_ribo(inBam, folder, pre, mapq, minRL, maxRL, RelE, output, unmap_fn)
    else:
        print("[error]\tmissing argument", file=sys.stderr)
        parser.print_usage()
