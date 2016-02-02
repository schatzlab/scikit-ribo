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
import sys
import cmd
import argparse
import pybedtools as pbt
from pybedtools.featurefuncs import gff2bed

class Gtf2Bed:
    ''' class to sort and get start codon from a gtf file '''
    def __init__(self, in_gtf): 
        self.in_gtf = in_gtf
        self.bedtool = None
        self.start = None
        self.out_bed = None

    def convert(self):
        ''' create a bedtool object '''
        self.bedtool = pbt.BedTool(self.in_gtf)
        
        ''' extract start codon entries '''
        self.start = self.bedtool.filter(lambda x: x[2] == "start_codon")
        
        ''' sort by coordinates '''
        self.start_sort = self.start.sort()

        ''' save the sorted start codon records to a bed file'''
        # self.fn = self.in_gtf.strip( '.gtf' )
        # self.out_bed = self.start_sort.each(gff2bed).saveas(self.fn + '.sort.start.bed') # , trackline='track name=test')

## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input gtf file")
    
    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    
    ## process the file if the input files exist
    if (args.i!=None):  # & (args.sort!=None) & (args.start!=None):
        print ("[status]\tprocessing the input file: " + args.i)
        input_fl = args.i
        bt_obj = Gtf2Bed(input_fl)
        bt_obj.convert()
        
    else:
        print ("[error]\tmissing argument")
        parser.print_usage() 