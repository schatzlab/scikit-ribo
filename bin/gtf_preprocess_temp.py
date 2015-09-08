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
# import cmd
import argparse

import gzip
import pandas as pd
import re
from collections import defaultdict

## parse and sort a gtf file
class parse_gtf:
    def __init__(self, infile, outfile):
        self.infile  = infile
        self.outfile = outfile
    
    
    ## define gtf fields and regular expression patterns
    gtf_fileds   = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    re_semicolon = re.compile(r'\s*;\s*')
    re_comma     = re.compile(r'\s*,\s*')
    re_keyvalue  = re.compile(r'(\s+|\s*=\s*)')
    
    
    def pd_dataframe(filename):
        """Open a GTF file and return a pandas.DataFrame.
        """
        # Each column is a list stored as a value in this dict.
        result = defaultdict(list)

        for i, line in enumerate(lines(filename)):
            for key in line.keys():
                # This key has not been seen yet, so set it to None for all
                # previous lines.
                if key not in result:
                    result[key] = [None] * i

            # Ensure this row has some value for each column.
            for key in result.keys():
                result[key].append(line.get(key, None))

        return pd.DataFrame(result)
    
    def readline_gtf(filename):
        """ Open a GTF file (gzipped or not) and generate a dict for each line.
        """
        if filename.endswith('.gz'):
            open_fn = gzip.open
        else:
            open_fn = open

        with open_fn(filename) as lines:
            for line in lines:
                if line.startswith('#'):
                    continue
                else:
                    yield parse(line)


    def open_df(filename):
        """Open an optionally gzipped GTF file and return a pandas.DataFrame.
        """
        
class sort_gtf:
   def __init__(self, in_gtf, out_gtf):
      self.in_gtf  = in_gtf
      self.out_gtf = out_gtf

   def sort_gtf_by_coordinate(self):
      cmd = "sort -k1,1V -k 4,4g -k 5,5g {0} > {1}".format( self.in_gtf, self.out_gtf )
      print ("[execute]\t" + cmd)
      os.system( cmd )
## extract the gtf records for start codon location
class extract_startcodon_location:
   def __init__(self, in_gtf, out_gtf):
      self.in_gtf  = in_gtf
      self.out_gtf = out_gtf
    
   def grep_startcodon(self):
      cmd = "grep 'start_codon' {0} > {1}".format( self.in_gtf, self.out_gtf )
      print ("[execute]\t" + cmd)
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
        print ("[status]\tprocessing the input file: " + args.i)
        
        sorted_gtf = sort_gtf( args.i, args.sort )
        sorted_gtf.sort_gtf_by_coordinate()
        startcodon = extract_startcodon_location( args.sort, args.start )
        startcodon.grep_startcodon()
        
    else:
        print ("[error]\tmissing argument")
        parser.print_usage() 
