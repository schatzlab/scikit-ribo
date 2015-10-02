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
import pysam
import pandas as pd
from collections import defaultdict

class get_uniq_aln:
    ''' define a class for extracting alignment based on MAPQ (>10)
    '''
    def __init__(self, in_bam, user_mapq, out_bam):
        self.in_bam    = in_bam
        self.out_bam   = out_bam
        self.user_mapq = int(user_mapq)

    def filter_by_mapq(self):
        self.samfile = pysam.AlignmentFile(self.in_bam, "rb")
        self.uniq_mappers = pysam.AlignmentFile(self.out_bam, "wb", template=self.samfile)

        if self.user_mapq > 10:
            mapq = self.user_mapq
        else:
            mapq = 10

        for read in self.samfile.fetch():
            if read.mapping_quality > mapq:
                self.uniq_mappers.write(read)
        
        return (self.samfile)
        print ("[execute]\t" + "filter by mapq of " + str(mapq) )
        unique_mappers.close()
        samfile.close()
        
class get_aln_info:
    ''' prepare a summary pandas table for alignment
    '''
    def __init__(self, in_bam, bed): # , aln_info):
        self.in_bam  = in_bam
        self.bed     = bed
        #self.aln_info = aln_info

    def get_chr(self):
        hdl = list(dict())
        for read in self.in_bam.fetch():    
            hdl.append({'chr': self.in_bam.getrname(read.reference_id) , 'start': read.reference_start, 'end': read.reference_end, 'len':read.query_length, 'seq':read.query_sequence})

        hdl_df = pd.DataFrame(hdl,columns=['chr', 'start', 'end', 'len','seq'])
        print (hdl_df)
        return (hdl_df)
        # bed = pysam.asBed(self.bed)
        #i=0
        #while i < fasta.nreferences:
        #    print(self.in_bam.getrname(i))
        #    i+=1   

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-b", help="a BED file with start codons")
    parser.add_argument("-q", help="minimum mapq allowed")
    parser.add_argument("-o", help="output bam file")

    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.b!=None) & (args.q!=None) & (args.o!=None) :
        print ("[status]\tprocessing the input bam file: " + args.i)
        unique_mappers = get_uniq_aln( args.i, args.q, args.o)
        unique_mappers.filter_by_mapq()

        get_aln = get_aln_info(unique_mappers.samfile,args.b)
        get_aln.get_chr()
    else:
        print ("[error]\tmissing argument")
        parser.print_usage() 