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

class rnafold(object):
    ''' process rnafold outputs
    '''
    def __init__(self, fa, folder, output):
        self.fa = fa
        self.folder = folder
        self.output = output
        self.lengDf = None
        self.header = None
        self.probDic = multiprocessing.Manager().dict()

    def loadFa(self):
        with open(self.fa, 'r') as f:
            firstLine = f.readline()
            if firstLine[0] != ">":
                exit("check format")
            self.header = "split" if "|" in firstLine else "gene"
        ## parse the file
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
        # self.lenDf.to_csv(path_or_buf= 'debug.txt', sep='\t', index=False)

    def loadDpps(self, fn):
        filePath = self.folder + "/" + fn
        # parse file names
        if self.header == "split":
            contig = fn.split("|")[3]
            pos = fn.split("|")[4].split("-")[0]
            geneLength = self.lenDf[((self.lenDf['contig'] == contig) & (self.lenDf['pos'] == pos))]['length'].values[0]
        elif self.header == "gene":
            geneName = fn.split("_")[0]
            geneLength = self.lenDf[(self.lenDf['geneName'] == geneName)]['length'].values[0]
        # parse file
        dp = []
        with open(filePath, 'r') as hdl:
            for line in hdl:
                tmp = line.strip().split(" ")
                if len(tmp) == 4 and tmp[3] == 'ubox':
                    pos, target, prob = int(tmp[0]), int(tmp[1]), float(tmp[2])
                    dp.append([pos, target, prob])
        probs = pd.DataFrame(dp, columns=["pos", "target", "probability"])
        grouped = probs.groupby('pos').max().reset_index()
        probs = grouped[["pos", "probability"]].fillna(value=0)
        # adding zeros to missing values
        lst = [[i] for i in range(geneLength)]
        fullDf = pd.DataFrame(lst, columns=["pos"])
        df = pd.merge(fullDf, probs, how = "left").fillna(0)
        # df.to_csv(path_or_buf= self.output + geneName + '.pair_prob.txt', sep='\t', header=False, index=False)
        # save to dic
        if self.header == "gene":
            gene = geneName
        elif self.header == "split":
            gene = contig + "|" + pos
        self.probDic[gene] = df["probability"].tolist()

    def loadAll(self):
        fileNames = []
        for file in os.listdir(self.folder):
            if file.endswith("dp.ps"):
                fileNames.append(file)
        pool = multiprocessing.Pool(16)
        pool.map(self.loadDpps, fileNames)
        print("[status]\tFinished loading rnafold results for all.", flush=True)

    def mergeAll(self):
        csvFile = open(self.output + '.txt', 'w')
        for k, v in self.probDic.items():
            csvFile.write(k + "\t" + " ".join(str(i) for i in v) + "\n")
        csvFile.close()


## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", help="transcript fasta file, required")
    parser.add_argument("-f", help="folder of rnafold outputs, required")
    parser.add_argument("-o", help="output file prefix, required")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.f!=None):
        print ("[status]\tReading the input fasta file: " + args.r, flush=True)
        print ("[status]\tReading the rnafold files from " + args.f, flush=True)
        fasta = args.r
        folder = args.f
        output = args.o
        #
        rna = rnafold(fasta, folder, output)
        print ("[status]\tParsing fasta file", flush=True)
        rna.loadFa()
        print ("[status]\tParsing the pairing probability file", flush=True)
        rna.loadAll()
        print ("[status]\tMerging the pairing probabilities into one file", flush=True)
        rna.mergeAll()
        print ("[status]\tFinished.", flush=True)
    else:
        print ("[error]\tmissing argument", flush=True)
        parser.print_usage()