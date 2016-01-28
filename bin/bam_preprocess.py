#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for preprocessing bam files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 7/24/2015
## ----------------------------------------

import os
import sys
import argparse
import pybedtools as pbt
import pysam
import pandas as pd
from collections import defaultdict

## ----------------------------------------
## class
## ----------------------------------------
class Get_Uniq_Aln:
    ''' define a class for extracting alignment based on MAPQ (>10)
    '''
    def __init__(self, in_bam_fn, user_mapq, out_bam_fn):
        self.in_bam_fn = in_bam_fn
        self.out_bam_fn = out_bam_fn
        self.user_mapq = int(user_mapq)

    def filter_by_mapq(self):
        self.pysam_hdl = pysam.AlignmentFile(self.in_bam_fn, "rb")
        self.pysam_ftd = pysam.AlignmentFile(self.out_bam_fn, "wb", template=self.pysam_hdl)

        if self.user_mapq > 10:
            mapq = self.user_mapq
        else:
            mapq = 10

        for read in self.pysam_hdl.fetch():
            if read.mapping_quality > mapq:
                self.pysam_ftd.write(read)
        
        return (self.pysam_ftd)
        print ("[execute]\t" + "filter by mapq of " + str(mapq) )
        pysam_ftd.close()
        pysam_hdl.close()
        
class Get_Aln_Info:
    ''' prepare a summary pandas table for alignment
    '''
    def __init__(self, pysam_hdl, bed): # , aln_info):
        self.pysam_hdl = pysam_hdl
        self.bed = bed
        #self.aln_info = aln_info

    def make_bedtool(self):
         self.bt_hdl = pbt.BedTool(self.pysam_hdl)
         print(self.pysam_hdl)
         print(self.bt_hdl)
         
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

        print(len(hdl_df.cigar))
        return(hdl_df)
 
## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help = "input bam file")
    parser.add_argument("-b", help = "a BED file with start codons")
    parser.add_argument("-q", help = "minimum mapq allowed")
    parser.add_argument("-o", help = "output bam file")

    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None and 
        args.b!=None and 
        args.q!=None and 
        args.o!=None) :
        
        print ("[status]\tprocessing the input bam file: " + args.i)
        fetch_aln = Get_Uniq_Aln(args.i, args.q, args.o)
        fetch_aln.filter_by_mapq()

        gen_aln = Get_Aln_Info(fetch_aln.pysam_hdl, args.b)
        gen_aln.make_bedtool()
        
    else:
        print ("[error]\tmissing argument")
        parser.print_usage() 