#!/usr/bin/env python

## ----------------------------------------
## Photon 
## ----------------------------------------
## a module for preprocessing bam files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 7/24/2015
## ----------------------------------------

## usage: python bam_preprocess.py bam_file

import os
import sys
import cmd
import argparse

## define global variables
global samtools
samtools="../tools/samtools-1.2/samtools"

## define a class for extracting alignment based on MAPQ (>10)
class extract_uniqe_alignment:
    
   def __init__(self, in_bam, out_bam):
      self.in_bam  = in_bam
      self.out_bam = out_bam

   def samtools_filter_mapq(self):
      cmd = "{0} view -bq 10 {1} > {2}".format( samtools, self.in_bam, self.out_bam )
      print "[execute]\tcommand: " + cmd      
      os.system( cmd )

## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-o", help="output bam file")
    
    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    
    ## process the file if the input files exist
    if (args.i!=None) & (args.o!=None):
        print "[status]\tprocessing the input bam file: " + args.i
        unique_mappers = extract_uniqe_alignment( args.i, args.o)
        unique_mappers.samtools_filter_mapq()
    else:
        print "[error]\tmissing argument"
        parser.print_usage() 