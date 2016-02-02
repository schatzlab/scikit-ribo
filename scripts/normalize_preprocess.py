#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for normalizing read counts
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

import sys
import argparse
import numpy as np
from scipy.stats.mstats import gmean
np.set_printoptions(precision=3)

class NormGeo:
    ''' define a class for calculating the normalization factor
    '''
    def __init__(self, input):
        self.in_array = input
        
    def normalize(self):
        input_array = self.in_array
        
        if np.any(input_array < 0, axis=None):
            sys.stderr.write('Error: negative value(s) detected.\n')
            sys.exit(1)

        ## add one to the connt and calculate geometric mean
        plus1_array = input_array.copy() + 1
        geo_means = gmean(plus1_array, axis=1 )
        
        ## initiate a array for normailization factors
        norm_factor = np.zeros(input_array.shape[1])

        for i in range(input_array.shape[1]):
            idx = input_array[:, i] > 0
            norm_factor[i] = np.median(input_array[idx, i] / geo_means[idx])
        return norm_factor

## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input the summarize count table")
    parser.add_argument("-o", help="output file")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.o!=None):
        print ("[status]\tProcessing the input file: " + args.i)

        ## Load in the file and extract the header
        with open(args.i, 'r') as fl:
            fl_header = np.array(fl.readline().strip().split('\t'))
        
        ## Load the content of the table
        cnt_array = np.loadtxt(args.i, dtype=int, delimiter='\t', skiprows=1, usecols=range(1, fl_header.size))

        ## Summary statistics of the table
        print ('[summary]\tTotal number of samples: %i' % cnt_array[0, :].size)
        print ('[summary]\tTotal number of genes: %i' % cnt_array[:, 0].size)

        ## Calculate the normalization factor for each sample
        norm_obj = NormGeo(cnt_array)
        norm_factor_list = norm_obj.normalize()

        print ("[summary]\t" + "\t".join(fl_header[1:]))
        print ("[summary]\tNormalization factors:", norm_factor_list)

    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
