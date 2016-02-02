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

class MergeCount:
    ''' define a class for merging read count tables
    '''
    def __init__(self, in_list, out_fn):
        self.in_list = in_list
        self.out_fn = out_fn
        
    def merge(self):
        fileDict = {}
        for fl in self.in_list:
            fileDict[fl] = pd.io.parsers.read_table(fl, sep='\t')

        first_col_id = str((list(fileDict.values())[0].columns.values[0]))
        func_merge = functools.partial(pd.merge, on = first_col_id, how='inner')
        merge_df = functools.reduce(func_merge, fileDict.values())
        merge_df.to_csv(self.out_fn, sep='\t',index=False)

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
        print ("[status]\tProcessing the input file..."  )
        in_list = args.i
        out_fn = args.o
        table_obj = MergeCount(in_list, out_fn)
        table_obj.merge()
        print ("[status]\tFinshed."  )

    else:
        print ("[error]\tmissing argument")
        parser.print_usage()