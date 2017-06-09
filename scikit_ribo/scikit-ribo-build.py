#!/usr/bin/env python

# ----------------------------------------
# scikit-ribo
# ----------------------------------------
# pre-processing module
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
from gtf_preprocess import GtfPreProcess
from process_rnafold import ProcessRnafold
from merge_df import MergeDF


def log_status(gtf_fn, ref_fn, prefix, rnafold_fn, tpm_fn, out_dir):
    """
    Logging the status
    :param gtf_fn: str, gtf file
    :param ref_fn: str, ref fasta file
    :param prefix: str, prefix
    :param rnafold_fn: str, rnafold file path
    :param tpm_fn: str, tpm file path
    :param out_dir: str, output directory
    :return: None
    """
    # create output folder
    cmd = 'mkdir -p ' + out_dir
    os.system(cmd)
    print("[status]\tStarted the pre-processing module", file=sys.stderr)
    print("[status]\tImport the gtf file: " + gtf_fn, file=sys.stderr)
    print("[status]\tImport the ref genome fasta file: " + ref_fn, file=sys.stderr)
    print("[status]\tImport RNAfold file: " + rnafold_fn, file=sys.stderr)
    print("[status]\tImport TPM file of RNAseq sample: " + tpm_fn, file=sys.stderr)
    print("[setting]\tPrefix to use: " + prefix, file=sys.stderr)
    print("[setting]\tOutput path: " + out_dir, file=sys.stderr)
    sys.stderr.flush()


def module_gtf(gtf_fn, ref_fn, prefix, out_dir):
    """
    Module for processing gtf and ref
    :param gtf_fn: str
    :param ref_fn: str
    :param prefix: str
    :param out_dir: str
    :return: None
    """
    worker = GtfPreProcess(gtf_fn, ref_fn, prefix, out_dir)
    print("[execute]\tLoading the the gtf file in to sql db", file=sys.stderr)
    worker.convertGtf()
    print("[execute]\tCalculating the length of each chromosome", file=sys.stderr)
    worker.getChrLen()
    print("[execute]\tExtracting the start codons' positions from the gtf db", file=sys.stderr)
    worker.getStartCodon()
    print("[execute]\tExtracting the sequences for each gene", file=sys.stderr)
    worker.getSeq()
    print("[execute]\tBuilding the index for each position at the codon level", file=sys.stderr)
    worker.getCodons()
    worker.getNts()
    print("[execute]\tCreating the codon table for the coding region", file=sys.stderr)
    worker.createCodonTable()
    print("[status]\tGtf processing module finished", file=sys.stderr)
    sys.stderr.flush()


def module_merge(prefix, tpm_fn, rnafold_fn, out_dir):
    """
    merge data
    :param prefix: prefix for the files
    :param tpm_fn: tmp file name
    :param rnafold_fn: rnafold file path
    :param out_dir: output directory
    :return:
    """
    bed = out_dir + "/" + prefix + ".codons.bed"
    # execute
    dat = MergeDF(bed, rnafold_fn, tpm_fn, out_dir)
    print("[execute]\tTransforming the dataframe of RNA secondary structure pairing probabilities", file=sys.stderr)
    dat.transformPairProb()
    print("[execute]\tLoading tpm", file=sys.stderr)
    dat.loadTpm()
    print("[execute]\tMerging all the df together", file=sys.stderr)
    dat.mergeDf()
    sys.stderr.flush()


def scikit_ribo_build(gtf_fn, ref_fn, prefix, rnafold_fn, tpm_fn, out):
    """

    :param gtf_fn:
    :param ref_fn:
    :param prefix:
    :param rnafold_fn:
    :param tpm_fn:
    :param out:
    :return: None
    """
    log_status(gtf_fn, ref_fn, prefix, rnafold_fn, tpm_fn, out)
    module_gtf(gtf_fn, ref_fn, prefix, out)
    module_merge(prefix, tpm_fn, rnafold_fn, out)
    print("[status]\tPre-processing module finished", file=sys.stderr)
    sys.stderr.flush()


# ----------------------------------------
# parse input arguments
# ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="Gtf file, required")
    parser.add_argument("-f", help="Fasta file, required")
    parser.add_argument("-p", help="Prefix to use, required")
    parser.add_argument("-r", help="Path to the Rnafold file, required")
    parser.add_argument("-t", help="TPM of RNAseq sample, required")
    parser.add_argument("-o", help="Output path of the built indexes, required")
    # check if there is any argument
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # process the file if the input files exist
    if args.g != None and args.f != None and args.p != None and args.r != None and args.t != None and args.o != None:
        sys.stderr.write("[status]\tReading the input file: " + args.g + "\n")
        gtf = args.g
        fasta = args.f
        pre = args.p
        rnafold = args.r
        tpm = args.t
        output = args.o
        scikit_ribo_build(gtf, fasta, pre, rnafold, tpm, output)
    else:
        print("[error]\tmissing argument", file=sys.stderr)
        parser.print_usage()
