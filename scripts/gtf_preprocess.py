#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for preprocessing gtf files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 7/24/2015
## ----------------------------------------

import os
import sys
import cmd
import argparse

## sort a gtf file
class sort_gtf:
   def __init__(self, in_gtf, out_gtf):
      self.in_gtf  = in_gtf
      self.out_gtf = out_gtf

   def sort_gtf_by_coordinate(self):
      cmd = "sort -k1,1V -k 4,4g -k 5,5g {0} > {1}".format( self.in_gtf, self.out_gtf )
      print "[execute]\t" + cmd      
      os.system( cmd )
## extract the gtf records for start codon location
class extract_startcodon_location:
   def __init__(self, in_gtf, out_gtf):
      self.in_gtf  = in_gtf
      self.out_gtf = out_gtf
    
   def grep_startcodon(self):
      cmd = "grep 'start_codon' {0} > {1}".format( self.in_gtf, self.out_gtf )
      print "[execute]\t" + cmd      
      os.system( cmd )

## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input gtf file")
    parser.add_argument("-sort", help="output sorted gtf file")
    parser.add_argument("-start", help="output start codon location in gtf format")
    
    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    
    ## process the file if the input files exist
    if (args.i!=None) & (args.sort!=None) & (args.start!=None):
        print "[status]\tprocessing the input file: " + args.i
        
        sorted_gtf = sort_gtf( args.i, args.sort )
        sorted_gtf.sort_gtf_by_coordinate()
        startcodon = extract_startcodon_location( args.sort, args.start )
        startcodon.grep_startcodon()
        
    else:
        print "[error]\tmissing argument"
        parser.print_usage() 
