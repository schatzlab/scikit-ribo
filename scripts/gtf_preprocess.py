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
        self.bedtool = pbt.BedTool(self.in_gtf)
        self.fn = os.path.splitext(self.in_gtf)[0]

    def convert(self):
        ## sort by coordinates
        self.bedtool_sort = self.bedtool.sort()

        ## extract start codon entries, save the records to a bed file
        self.start = self.bedtool_sort.filter(lambda x: x[2] == "start_codon")
        self.out_bed = self.start.each(gff2bed).saveas(self.fn + '.sort.start.bed')

        ## extract CDS entries, save the records to a bed file
        self.start = self.bedtool_sort.filter(lambda x: x[2] == "CDS")
        self.out_bed = self.start.each(gff2bed).saveas(self.fn + '.sort.CDS.bed')

        #self.in_gtf.strip( '.gtf' )
        #self.out_bed = self.start.each(gff2bed).saveas(self.fn + '.sort.start.bed')

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