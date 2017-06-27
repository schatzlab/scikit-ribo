#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## Run RNAfold in parallel
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
from itertools import groupby
from process_rnafold import ProcessRnafold


class CallRnafold(object):
    ''' Run RNAfold in parallel
    '''
    def __init__(self, fastaFn, rnafold, prefix, processes, output_dir):
        self.fastaFn = fastaFn
        self.rnafold = rnafold
        self.prefix = prefix
        self.processes = processes
        self.output_dir = output_dir
        self.fastaDic = defaultdict(str)

    def fastaIter(self):
        fh = open(self.fastaFn)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            geneName = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            self.fastaDic[geneName] = seq

    def splitFa(self):
        cmd = "mkdir -p fastaFiles; mkdir -p tmpFiles"
        os.system(cmd)
        for geneName in self.fastaDic.keys():
            fa = open("fastaFiles/" + geneName + ".fasta", 'w')
            fa.write(">" + geneName + "\n" + self.fastaDic[geneName] + "\n")
            fa.close()
    
    def runRnafold(self, geneName):
        cmd = self.rnafold + ' -i ' + 'fastaFiles/' + geneName + '.fasta' + ' -p --outfile tmpFiles/' + geneName
        os.system(cmd)
        cmd = "mv " + geneName + "*ps tmpFiles/"
        os.system(cmd)

    def executeInParallel(self):
        geneNames = self.fastaDic.keys()
        cpus = self.processes if self.processes > multiprocessing.cpu_count() else max(1, multiprocessing.cpu_count()-1)
        pool = multiprocessing.Pool(cpus, maxtasksperchild=2)
        pool.map(self.runRnafold, geneNames)
        pool.close()
        pool.join()
        print("[status]\tFinished calling rnafold for all.", file=sys.stderr)

    def parseData(self):
        """
        parse RNAfold outputs and merge into one file
        :return: None
        """
        cmd = "mkdir -p " + self.output_dir
        os.system(cmd)
        rna = ProcessRnafold(self.fastaFn, "tmpFiles", self.prefix, self.output_dir)
        print("[status]\tParsing fasta file", file=sys.stderr)
        rna.loadFa()
        print("[status]\tParsing the pairing probability file", file=sys.stderr)
        rna.loadAll()
        print("[status]\tMerging the pairing probabilities into one file", file=sys.stderr)
        rna.mergeAll()
        print("[status]\tRNAfold processing module finished", file=sys.stderr)
        sys.stderr.flush()

    def rmTmpFiles(self):
        cmd = 'rm -rf fastaFiles/; rm -rf tmpFiles/'
        os.system(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="transcript Fasta file, required")
    parser.add_argument("-r", help="path to Rnafold binary, required")
    parser.add_argument("-p", help="Prefix for output file name, required")
    parser.add_argument("-n", help="Number of parallel processes, Default: 4", default=4, type=int)
    parser.add_argument("-o", help="Output folder, required")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.f!=None) and (args.r!=None) and (args.p!=None) and (args.o!=None):
        sys.stderr.write("[status]\tReading the input fasta file: " + args.f + "\n")
        sys.stderr.write("[status]\tPath to rnafold binary:" + args.r + "\n")
        sys.stderr.write("[status]\tPrefix for output file name:" + args.p + "\n")
        sys.stderr.write("[status]\tNumber of parallel processes:" + str(args.n) + "\n")
        sys.stderr.write("[status]\tOutput folder name:" + args.o + "\n")
        # parse arguments
        fasta = args.f
        rnafold = args.r
        prefix = args.p
        processes = args.n
        output = args.o
        # run
        mod = CallRnafold(fasta, rnafold, prefix, processes, output)
        sys.stderr.write("[status]\tIterate through the fasta file" + "\n")
        mod.fastaIter()
        sys.stderr.write("[status]\tSplitting the fasta file" + "\n")
        mod.splitFa()
        sys.stderr.write("[status]\tExecute RNAfold on all sub fasta files" + "\n")
        mod.executeInParallel()
        sys.stderr.write("[status]\tCombine data into one file" + "\n")
        mod.parseData()
        sys.stderr.write("[status]\tRemove temp files" + "\n")
        mod.rmTmpFiles()
        sys.stderr.write("[status]\tRNAfold processing module finished" + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()