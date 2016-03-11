#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for preprocessing bam files
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

from __future__ import print_function, division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import argparse

# import pybedtools as pbt
# import pysam
# import pandas as pd
# import numpy as np
# from collections import defaultdict
# from re import match
# from gtf_preprocess import *
# import multiprocessing

class VisualizeAsite:
    ''' define a class for plotting a-site location distribution
    '''
    def __init__(self, asite_fn):
        self.asite_fn = asite_fn
        self.asite_df = pd.read_table(self.asite_fn,  header=0)
        print(self.asite_df)
    def plot(self):
        sns.set(font_scale=2)
        g0 = sns.FacetGrid(self.asite_df, row="offset", col="read_length", margin_titles=True , col_order= list(range(25,36)), row_order= list(range(3))) # col_order= [26,28,32]
        bins = np.linspace(12, 21, 10)
        g0.map(plt.hist, "a_site", color="steelblue", bins=bins, lw=0,normed=True)
        g0.set(xticks=[12,15,18,21])
        plt.gcf()
        plt.savefig( asite_fn +'.asite.pdf')

## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input a-site file")
    #parser.add_argument("-b", help="a BED file with start codons")
    #parser.add_argument("-q", help="minimum mapq allowed", default=10, type=int)
    #parser.add_argument("-o", help="output bam file")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i != None):
        print("[status]\tprocessing the input file: " + args.i)
        asite_fn = args.i

        print("[execute]\tplotting the a-site location distribution from " + str(asite_fn))
        asite_loc = VisualizeAsite(asite_fn)
        asite_loc.plot()

    else:
        print("[error]\tmissing argument")
        parser.print_usage()

    '''
    def get_chr(self):
        hdl = list(dict())
        for read in self.in_bam_fn.fetch():
            hdl.append({'chr': self.in_bam_fn.getrname(read.reference_id) ,
                        'start': read.reference_start,
                        'end': read.reference_end,
                        'len': read.query_length,
                        'cigar': read.cigar,
                        'seq': read.query_sequence})

        hdl_df = pd.DataFrame(hdl, columns = ['chr', 'start', 'end', 'len','cigar','seq'])

        print(hdl_df)
    '''