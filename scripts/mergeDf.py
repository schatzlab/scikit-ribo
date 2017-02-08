#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for preprocessing gtf/bed files
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 10/28/2016
## ----------------------------------------

from __future__ import print_function, division
import os
import sys
import csv
import argparse
import pybedtools as pbt
import pandas as pd
import numpy as np
from pybedtools.featurefuncs import gff2bed
from itertools import groupby


class mergeDf(object):
    ''' class to sort and get start codon from a gtf file
    '''

    def __init__(self, fn=None, pairprobFn=None, tpmFn=None, output=None):
        self.fn = fn
        self.pairprobFn = pairprobFn
        self.tpmFn = tpmFn
        self.output = output

    def transformPairProb(self):
        ## read the pairing prob arrays then convert it to a df
        pairProb = []
        with open(self.pairprobFn, 'r') as fl:
            for line in fl:
                row = line.split("\t")
                codonIdx = -8
                geneName = row[0]
                probs = row[1].split(" ")
                for prob in probs:
                    pairProb.append([geneName, codonIdx, float(prob)])
                    codonIdx += 1
        self.pairProb = pd.DataFrame(pairProb, columns=["gene", "codon_idx", "pair_prob"])

    def loadTpm(self):
        self.tpm = pd.read_table(self.tpmFn, header=0)
        tpmColNames = set(list(self.tpm.columns.values))
        if 'TPM' in tpmColNames:
            tool = 'Salmon'
            self.tpm = self.tpm[["Name", "TPM"]]
            self.tpm.columns = ["gene", "TPM"]
        elif 'tpm' in tpmColNames:
            tool = 'Kallisto'
            self.tpm = self.tpm[["target_id", "tpm"]]
            self.tpm.columns = ["gene", "TPM"]
        else:
            exit("Check file format, only support Salmon or Kallisto")
        sys.stderr.write("[status]\tTPM input: " + str(tool) + "\n")

    def mergeDf(self):
        # import codon df
        codons = pd.read_table(self.fn, header=0)
        codons['codon_idx'].astype(int)
        ## import the salmon df, rna secondary structure, and merge with cds df
        codons = pd.merge(codons, self.tpm, how="left", on="gene")
        codons = pd.merge(codons, self.pairProb, how="left", on=["gene", "codon_idx"]).fillna('NA')
        codons = codons[["chrom", "start", "end", "gene", "codon_idx", "gene_strand", "codon", "TPM", "pair_prob"]]
        # parse prefix, export file
        base = os.path.basename(self.fn)
        prefix = self.output + "/" + os.path.splitext(base)[0]
        codons.to_csv(path_or_buf=prefix + '.bed', sep='\t', header=True, index=False)


## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="dataframe file of codon table, required")
    parser.add_argument("-s", help="arrays of RNA secondary structure pairing probabilities, required")
    parser.add_argument("-t", help="pre-computed tpm salmon data-frame from RNA-seq data, required")
    parser.add_argument("-o", help="output path, required")
    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    ## process the file if the input files exist
    if (args.f != None) & (args.s != None) & (args.t != None) & (args.o != None):
        sys.stderr.write("[status]\tReading the input file: " + args.f + "\n")
        fn = args.f
        pairprob = args.s
        tpm = args.t
        output = args.o
        # create output folder
        cmd = 'mkdir -p ' + output
        os.system(cmd)
        ## execute
        sys.stderr.write("[execute]\tStarting the pre-processing module" + "\n")
        dat = mergeDf(fn, pairprob, tpm, output)
        sys.stderr.write("[execute]\tTransforming the dataframe of RNA 2' structure pairing probabilities" + "\n")
        dat.transformPairProb()
        sys.stderr.write("[execute]\tLoading tpm" + "\n")
        dat.loadTpm()
        sys.stderr.write("[execute]\tMerging all the df together" + "\n")
        dat.mergeDf()
        ## finish
        sys.stderr.write("[status]\tData merging module finished" + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
