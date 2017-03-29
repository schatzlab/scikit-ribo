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


class CallRnafold(object):
    ''' Run RNAfold in parallel
    '''
    def __init__(self, fastaFn, rnafold):
        self.fastaFn = fastaFn
        self.rnafold = rnafold
        self.fastaDic = defaultdict(str)

    def fastaIter(self):
        fh = open(self.fastaFn)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            geneName = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            self.fastaDic[geneName] = seq

    def splitFa(self):
        cmd = "mkdir -p fastaFiles"
        os.system(cmd)
        for geneName in self.fastaDic.keys():
            fa = open("fastaFiles/" + geneName + ".fasta", 'w')
            fa.write(">" + geneName + "\n" + self.fastaDic[geneName] + "\n")
            fa.close()
    
    def callRnafold(self, geneName):
        cmd = self.rnafold + ' -i ' + 'fastaFiles/' + geneName + '.fasta' + ' -p --outfile ' + geneName
        os.system(cmd)

    def runAll(self):
        geneNames = self.fastaDic.keys()
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(self.callRnafold, geneNames)
        print("[status]\tFinished calling rnafold for all.", flush=True)
    
    def rmTmpFa(self):
        cmd = 'rm -rf fastaFiles/'
        os.system(cmd)

if __name__ == '__main__':
    fasta = sys.argv[1]
    rnafold = sys.argv[2]
    worker = CallRnafold(fasta, rnafold)
    worker.fastaIter()
    worker.splitFa()
    worker.runAll()
    worker.rmTmpFa()
