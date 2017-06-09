#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## rna-fold output processing
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 12/28/2016
## ----------------------------------------

from __future__ import division, with_statement
import os
import argparse
import sys
import functools
import csv
import multiprocessing
from collections import defaultdict
import pandas as pd

class ProcessRnafold(object):
    ''' process rnafold outputs
    '''
    def __init__(self, fa, input, prefix, output):
        self.fa = fa
        self.input = input
        self.prefix = prefix
        self.output = output
        self.lenDf = None
        self.header = None
        self.probDic = multiprocessing.Manager().dict()

    def loadFa(self):
        with open(self.fa, 'r') as f:
            firstLine = f.readline()
            if firstLine[0] != ">":
                exit("check format")
            self.header = "split" if "|" in firstLine else "gene"
        # parse the file
        lst = []
        with open(self.fa, 'r') as f:
            dat = f.read().split("\n")
            n = len(dat)
            i = 0
            if self.header == "gene":
                while i < n-1:
                    line = dat[i]
                    if line and line[0] == ">":
                        geneName = line[1:]
                        tmp = ""
                    else:
                        tmp += line.rstrip()
                        if i < n-2 and dat[i+1][0] == ">" or i == n-2:
                            lst.append([geneName, len(tmp)])
                    i += 1
            elif self.header == "split":
                while i < n-1:
                    line = dat[i]
                    if line and line[0] == ">":
                        contig, pos = line.split("|")[3], line.split("|")[4].split("-")[0]
                        tmp = ""
                    else:
                        tmp += line.rstrip()
                        if i < n-2 and dat[i+1][0] == ">" or i == n-2:
                            lst.append([contig, pos, len(tmp)])
                    i += 1
        # convert to a dataframe
        if self.header == "split":
            self.lenDf = pd.DataFrame(lst, columns=["contig", "pos", "length"])
        elif self.header == "gene":
            self.lenDf = pd.DataFrame(lst, columns=["geneName", "length"])
            self.geneSet = set(self.lenDf['geneName'].tolist())

    def loadDpps(self, fn):
        filePath = self.input + "/" + fn
        # parse file names
        if self.header == "split":
            contig = fn.split("|")[3]
            pos = fn.split("|")[4].split("-")[0]
            geneLength = self.lenDf[((self.lenDf['contig'] == contig) & (self.lenDf['pos'] == pos))]['length'].values[0]
        elif self.header == "gene":
            geneName = fn.split("_")[0]
            # find intersection of chr length df and 2' structure probabilities
            if geneName not in self.geneSet:
                return
            geneLength = self.lenDf[(self.lenDf['geneName'] == geneName)]['length'].values[0]
        # parse file
        lbox, ubox = [], []
        with open(filePath, 'r') as hdl:
            for line in hdl:
                tmp = line.strip().split(" ")
                if len(tmp) == 4 and tmp[3] == 'ubox':
                    pos, target, prob = int(tmp[0]), int(tmp[1]), float(tmp[2]) ** 2
                    ubox.append([pos, target, prob])
                    ubox.append([target, pos, prob])
                if len(tmp) == 4 and tmp[3] == 'lbox':
                    pos, target, prob = int(tmp[0]), int(tmp[1]), float(tmp[2])
                    lbox.append([pos, target, prob])
        # lbox
        lbox_probs = pd.DataFrame(lbox, columns=["pos", "target", "probability"])
        lbox_grouped = lbox_probs.groupby('pos').max().reset_index()
        lbox_probs = lbox_grouped[["pos", "probability"]].fillna(value=0)
        # ubox
        ubox_probs = pd.DataFrame(ubox, columns=["pos", "target", "probability"])
        ubox_grouped = ubox_probs.groupby('pos').max().reset_index()
        ubox_probs = ubox_grouped[["pos", "probability"]].fillna(value=0)
        # adding zeros to missing values
        lst = [[i] for i in range(geneLength)]
        fullDf = pd.DataFrame(lst, columns=["pos"])
        lbox_df = pd.merge(fullDf, lbox_probs, how="left").fillna(0)
        lbox_df = lbox_df.sort_values("pos")
        # ubox
        ubox_df = pd.merge(fullDf, ubox_probs, how="left").fillna(0)
        ubox_df = ubox_df.sort_values("pos")
        # save to dic
        gene = None
        if self.header == "gene":
            gene = geneName
        elif self.header == "split":
            gene = contig + "|" + pos
        # save a list to dic
        self.probDic[gene] = [lbox_df["probability"].tolist(), ubox_df["probability"].tolist()]

    def loadAll(self):
        # save all file names to a list
        fileNames = []
        for file in os.listdir(self.input):
            if file.endswith("dp.ps"):
                fileNames.append(file)
        # pool and run in parallel
        cpus = 10 if multiprocessing.cpu_count() > 10 else max(1, multiprocessing.cpu_count() - 1)
        pool = multiprocessing.Pool(cpus, maxtasksperchild=2)
        pool.map(self.loadDpps, fileNames)
        pool.close()
        pool.join()
        sys.stderr.write("[status]\tFinished loading rnafold results for all" + "\n")

    def mergeAll(self):
        # write lbox
        lboxFile = open(self.output + '/' + self.prefix + '.rnafold_lbox.txt', 'w')
        for k, v in self.probDic.items():
            lbox = v[0]
            lboxFile.write(k + "\t" + " ".join(str(i) for i in lbox) + "\n")
        lboxFile.close()
        # write ubox
        uboxFile = open(self.output + '/' + self.prefix + '.rnafold_ubox.txt', 'w')
        for k, v in self.probDic.items():
            ubox = v[1]
            uboxFile.write(k + "\t" + " ".join(str(i) for i in ubox) + "\n")
        uboxFile.close()

# the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="transcript fasta file, required")
    parser.add_argument("-i", help="input folder of rnafold files, required")
    parser.add_argument("-o", help="output folder, required")

    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.f!=None):
        sys.stderr.write("[status]\tReading the input fasta file: " + args.f + "\n")
        sys.stderr.write("[status]\tReading the rnafold files from " + args.i + "\n")
        fasta = args.f
        input = args.i
        output = args.o
        #
        rna = ProcessRnafold(fasta, input, output)
        sys.stderr.write("[status]\tParsing fasta file" + "\n")
        rna.loadFa()
        sys.stderr.write("[status]\tParsing the pairing probability file" + "\n")
        rna.loadAll()
        sys.stderr.write("[status]\tMerging the pairing probabilities into one file" + "\n")
        rna.mergeAll()
        sys.stderr.write("[status]\tMerging finished" + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()