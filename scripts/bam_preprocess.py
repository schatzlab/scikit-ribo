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
from gtf_preprocess import *


class FilterAln:
    ''' define a class for extracting alignment based on MAPQ (>10) and read length meets [25,35]
    '''

    def __init__(self, in_bam_fn, user_mapq, out_bam_fn):
        self.in_bam_fn = in_bam_fn
        self.user_mapq = user_mapq
        self.out_bam_fn = out_bam_fn
        self.read_list = list()

    ## TODO: add a function to see whether there are too many soft-clipped alignment
    def filter(self):
        pysam_hdl = pysam.AlignmentFile(self.in_bam_fn, "rb")
        pysam_ftd = pysam.AlignmentFile(self.out_bam_fn, "wb", template=pysam_hdl)
        
        ## read a bam file and extract info
        for read in pysam_hdl.fetch():
            cigar_to_exclude = ('I','D','S','H')
            if read.mapping_quality > self.user_mapq and \
                            read.query_length >= 25 and read.query_length <= 35 and \
                            not any(c in read.cigarstring for c in cigar_to_exclude):
                pysam_ftd.write(read)
                self.read_list.append([read.query_name, read.query_length, read.query_sequence[0:2],
                                       read.query_sequence[-2:][::-1] ]) ## remove read
        return (self.read_list)
        pysam_ftd.close()
        pysam_hdl.close()


class ConvertBam():
    ''' prepare a pandas df of training data from the alignment
    '''

    def __init__(self, object, sc_bed, cds_bed, bam):
        self.sc_bed = pbt.BedTool(sc_bed)
        self.cds_bed = pbt.BedTool(cds_bed)
        self.bam = bam
        self.bt_hdl = pbt.BedTool(self.bam).bam_to_bed(stream=True, bed12=True).sort()
        self.read_list = object.read_list
        self.read_df = pd.DataFrame(self.read_list, columns=['name', 'read_length', 'start_seq', 'end_seq']) ## remove read


    def training_data(self):
        ## create pandas dataframes
        bam_start = self.bt_hdl.intersect(self.sc_bed, wa=True, wb=True, sorted=True)
        bam_start_df = bam_start.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart',
                                                     'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                                                     'sc_chrom', 'sc_start', 'sc_end', 'gene_name', 'sc_score',
                                                     'gene_strand'])

        ## retrieve the read length and seq information from the bam file
        bam_start_df_read = pd.merge(bam_start_df, self.read_df, on='name')
        bam_start_df_read['asite'] = np.where(bam_start_df_read['gene_strand'] == '+',
                                               bam_start_df_read['sc_start'] - bam_start_df_read['start'] + 3,
                                               bam_start_df_read['end'] - bam_start_df_read['sc_end'] + 3 )
        bam_start_df_read['offset'] = bam_start_df_read['asite'] % 3

        ## filter a read by whether it has a-site that satisfies [12,18]
        bam_start_df_read_filter = bam_start_df_read[((bam_start_df_read['asite'] >= 12) &
                                                      (bam_start_df_read['asite'] <= 18))]

        ## slice the dataframe to the variables needed for training data
        bam_start_df_read_filter_out = bam_start_df_read_filter[
            ["asite", "read_length", "offset", "start_seq", "end_seq", "gene_strand"]]
        bam_start_df_read_filter_out.to_csv(path_or_buf=self.bam + '.asite.txt', sep='\t', header=True, index=False)

    def testing_data(self):
        ## create pandas dataframes from bedtools intersect, /*important*/ to use split=True here
        bam_cds = self.bt_hdl.intersect(self.cds_bed, wa=True, wb=True, split=True) #sorted=True,
        bam_cds_df = bam_cds.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart',
                                                 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                                                 'gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score',
                                                 'gene_strand', 'gene_thickStart', 'gene_thickEnd', 'gene_itemRgb',
                                                 'gene_blockCount', 'gene_blockSizes', 'gene_blockStarts'],
                                          dtype={'blockSizes':'object','blockStarts':'object', 
                                                 'gene_blockSizes':'object','gene_blockStarts':'object'})

        ## retrieve the read length and seq information from the bam file
        bam_cds_df_read = pd.merge(bam_cds_df, self.read_df, on='name')
        bam_cds_df_read['distance'] = np.where(bam_cds_df_read['gene_strand'] == '+',
                                               bam_cds_df_read['gene_start'] - bam_cds_df_read['start'] ,
                                               bam_cds_df_read['end'] - bam_cds_df_read['gene_end'] )
        bam_cds_df_read['offset'] = bam_cds_df_read['distance'] % 3

        ## slice the dataframe to the variables needed for training data
        
        bam_cds_df_read_out = bam_cds_df_read[[ "read_length", "offset", "start_seq", "end_seq", "gene_chrom", 
                                                "gene_start", "gene_end", "gene_name", "gene_strand", "chrom", 
                                                "start", "end", "name", "score", "strand"]] ## remove read
        bam_cds_df_read_out.to_csv(path_or_buf=self.bam + '.cds.txt', sep='\t', header=True, index=False)
        

## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-b", help="a BED file with start codons")
    parser.add_argument("-g", help="a BED file with CDS (gene) entries")
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
                args.b != None and
                args.o != None):

        print("[status]\tprocessing the input bam file: " + args.i, flush=True)
        in_bam_fn = args.i
        user_mapq = args.q
        out_bam_fn = args.o
        start_bed_fl = args.b
        cds_bed_fl = args.g

        print("[execute]\tfilter by read length between [25,35] and mapq of " + str(user_mapq), flush=True)
        fetch_aln = FilterAln(in_bam_fn, user_mapq, out_bam_fn)
        fetch_aln.filter()

        print("[execute]\tstart converting the filtered bam file", flush=True)
        converted_bam = ConvertBam(fetch_aln, start_bed_fl, cds_bed_fl, out_bam_fn)
        print("[execute]\tconstruct pandas dataframe of training data", flush=True)
        converted_bam.training_data()
        print("[execute]\tconstruct pandas dataframe of testing data", flush=True)
        converted_bam.testing_data()

    else:
        print("[error]\tmissing argument", flush=True)
        parser.print_usage()

