#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for preprocessing gtf files
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

import os
import argparse
import sys
import pandas as pd
import functools


def __mergeCount__(in_list,out_fn):
    output_fn = out_fn 
    fileList = in_list # sorted(iglob('*.skip'))
    fileDict = {}

    for fl in fileList:
        fileDict[fl] = pd.io.parsers.read_table(fl, sep='\t')

    first_col_id = str((list(fileDict.values())[0].columns.values[0]))
    func_merge = functools.partial(pd.merge, on = first_col_id, how='inner')
    merge_df = functools.reduce(func_merge, fileDict.values())
    merge_df.to_csv(output_fn, sep='\t',index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input files", type=argparse.FileType('r'), nargs='+')
    parser.add_argument("-o", help="output file")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.o!=None):
        print ("[status]\tprocessing the input file: " + str(args.i) )
        __mergeCount__(args.i,args.o)

    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
