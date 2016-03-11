#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for preprocessing bam files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

from __future__ import print_function, division
import os
import sys
import argparse
import pybedtools as pbt
import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
from re import match
from gtf_preprocess import *
# import multiprocessing

class FilterAln:
    ''' define a class for extracting alignment based on MAPQ (>10) and read length meets [25,35]
    '''

    def __init__(self, in_bam_fn, user_mapq, out_bam_fn):
        self.in_bam_fn = in_bam_fn
        self.user_mapq = user_mapq
        self.out_bam_fn = out_bam_fn
        self.read_list = list()

    def filter(self):
        pysam_hdl = pysam.AlignmentFile(self.in_bam_fn, "rb")
        pysam_ftd = pysam.AlignmentFile(self.out_bam_fn, "wb", template = pysam_hdl)

        for read in pysam_hdl.fetch():
            if read.mapping_quality > self.user_mapq and \
            read.query_length >= 25 and read.query_length <= 35 and \
            match(r'\d\dM$', read.cigarstring) :
                pysam_ftd.write(read)
                self.read_list.append ([ read.query_name, read.query_length, read.query_sequence[0:2], read.query_sequence[-2:][::-1] ])
        return(self.read_list)
        pysam_ftd.close()
        pysam_hdl.close()

class MakeTrainingSet():
    ''' prepare a pandas df of training data from the alignment
    '''

    def __init__(self, object, bed, bam):
        self.bed = bed
        self.bam = bam
        self.bt_hdl = pbt.BedTool(self.bam)
        self.read_list = object.read_list

    def overlap_start(self):
        ## create pandas dataframes
        read_df = pd.DataFrame(self.read_list, columns=['name', 'read_length', 'start_seq', 'end_seq'] )
        bam_start = self.bt_hdl.intersect(self.bed, bed=True, wa=True, wb=True, sorted=True, split=True)
        bam_start_df = bam_start.to_dataframe( names = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                                                       'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                                                       'sc_chrom', 'sc_start', 'sc_end', 'sc_name', 'sc_score', 'sc_strand']                                    )

        ## retrieve the read length and seq information from the bam file
        bam_start_df_read = pd.merge(bam_start_df, read_df, on='name')
        bam_start_df_read['a_site'] = np.where(bam_start_df_read['sc_strand']=='+',
                                             bam_start_df_read['sc_start'] - bam_start_df_read['start'] + 3 ,
                                             bam_start_df_read['end'] - bam_start_df_read['sc_end'] + 3 )
        bam_start_df_read['offset'] = bam_start_df_read['a_site'] % 3

        ## filter a read by whether it has a-site that satisfies [12,18]
        bam_start_df_read_filter = bam_start_df_read[ ( (bam_start_df_read['a_site'] >= 12) &
                                                    (bam_start_df_read['a_site'] <= 18) ) ]

        ## slice the dataframe to the variables needed for training data
        bam_start_df_read_filter_out = bam_start_df_read_filter[["read_length", "offset", "a_site", "start_seq", "end_seq" ]]
        bam_start_df_read_filter_out.to_csv(path_or_buf=self.bam+'.asite.txt', sep='\t', header=True, index=False)

## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-b", help="a BED file with start codons")
    parser.add_argument("-q", help="minimum mapq allowed", default=20, type=int)
    parser.add_argument("-o", help="output bam file")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i != None and
                args.b != None and
                args.o != None):

        print("[status]\tprocessing the input bam file: " + args.i)
        in_bam_fn = args.i
        user_mapq = args.q
        out_bam_fn = args.o
        bed_fl = args.b

        print("[execute]\tfilter by read length between [25,35] and mapq of " + str(user_mapq))
        fetch_aln = FilterAln(in_bam_fn, user_mapq, out_bam_fn)
        fetch_aln.filter()

        print("[execute]\tconstruct pandas dataframe of traning data")
        training_set = MakeTrainingSet(fetch_aln, bed_fl, out_bam_fn)
        training_set.overlap_start()
        # gen_aln.get_chr()

    else:
        print("[error]\tmissing argument")
        parser.print_usage()

    '''
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

        print(hdl_df)
    '''