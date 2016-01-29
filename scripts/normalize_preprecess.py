#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for normalizing
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

import sys
import argparse
import numpy as np

def lib_size(countNdarray):
    """ Calculate normalization factor    
    @args countNdarray: read count data
    @type countNdarray: numpy array
    """

    if np.any(countNdarray < 0, axis=None):
        sys.stderr.write('Error: read count smaller than one is detected, please check input file.\n')
        sys.exit()

    countNdarrayTmp = countNdarray.copy()
    countNdarrayTmp = countNdarrayTmp + 1
    geoMeans = np.exp(np.mean(np.log(countNdarrayTmp), axis=1))

    librarySizes = np.zeros(countNdarray.shape[1])

    for i in range(countNdarray.shape[1]):
        idx = countNdarray[:, i] > 0
        librarySizes[i] = np.median(countNdarray[idx, i] / geoMeans[idx])

    return librarySizes

## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input count file")
    parser.add_argument("-o", help="output file")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.o!=None):  # & (args.sort!=None) & (args.start!=None):
        print ("[status]\tprocessing the input file: " + args.i)

        with open(args.i, 'r') as FileIn:
            header = np.array(FileIn.readline().strip().split('\t'))

        count = np.loadtxt(args.i, dtype=int, delimiter='\t', skiprows=1, usecols=range(1, header.size))

        print ('[summary]\tTotal number of records: %i' % count[:, 0].size)
        libSizes = lib_size(count)
        np.set_printoptions(precision=3)
        print (header[1:])
        print ("[summary]\tNormalization factors:", libSizes)

    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
