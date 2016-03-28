#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for preprocessing gtf/bed files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

from __future__ import print_function, division
import os
import sys
import csv
import argparse
import pybedtools as pbt
from pybedtools.featurefuncs import gff2bed
import numpy as np
import gffutils


class Gtf2Bed:
    ''' class to sort and get start codon from a gtf file
    '''
    def __init__(self, gtf, fasta):
        self.gtf = gtf
        self.fasta = fasta
        self.start = None
        self.bedtool = pbt.BedTool(self.gtf)
        self.prefix = os.path.splitext(self.gtf)[0]
        self.gene_bed12s = list()

    def convert_gtf(self):
        ## create a db from the gtf file
        self.db = gffutils.create_db(self.gtf, ":memory:", force=True)

        ## retrieve a list of gene names with type "CDS" from db
        gene_names = list()
        for gene in self.db.features_of_type("CDS", order_by=("seqid","start") ):
            gene_name = gene['gene_id'][0]
            if gene_name not in gene_names:
                gene_names.append(gene_name)

        ## convert a gtf/gff3 file to bed12 and save to a nested list
        for gene_name in gene_names:
            gene_bed12 = self.db.bed12(gene_name, name_field='gene_id')
            self.gene_bed12s.append(gene_bed12.split("\t"))

        ## extract CDS entries, write to a bed12 file
        with open( self.prefix + '.sort.CDS.bed', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(self.gene_bed12s)

    def get_startcodon(self):
        ## return a set of feature types
        feature_types = set(list(self.db.featuretypes()))

        ## sort by coordinates, extract start codon entries, save the records to a bed file
        if "start_codon" in feature_types:
            self.bedtool_sort = self.bedtool.sort()
            self.start = self.bedtool_sort.filter(lambda x: x[2] == "start_codon")
            self.start_bed = self.start.each(gff2bed).saveas(self.prefix + '.sort.start.bed')
        else:
            gene_bed6s = list()
            for bed12 in self.gene_bed12s:
                bed6 = bed12[0:6]
                bed6[2] = int(bed6[1]) + 3
                gene_bed6s.append(bed6)
            with open( self.prefix + '.sort.start.bed', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerows(gene_bed6s)

    def get_fasta(self):
        ## extract cds sequence from the reference genome
        CDS_bedtool = pbt.BedTool(self.prefix + '.sort.CDS.bed')
        CDS_bedtool.sequence(fi=self.fasta, fo= self.prefix + '.fasta', s=True, name=True, tab=True, split=True)

        CDS_bedtool_df = CDS_bedtool.to_dataframe(names=['gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score',
                                                         'gene_strand', 'gene_thickStart', 'gene_thickEnd', 'gene_reserved',
                                                         'gene_blockCount','gene_blockSizes', 'gene_blockStarts'])

        ##  construct the df from the bed12 file
        cds_codonidxs = list()
        for gene_name in CDS_bedtool_df.gene_name:
            gene_start = CDS_bedtool_df.loc[ CDS_bedtool_df["gene_name"]==gene_name ]["gene_start"].values[0]
            gene_end = CDS_bedtool_df.loc[ CDS_bedtool_df["gene_name"]==gene_name ]["gene_end"].values[0]

            chrom = CDS_bedtool_df.loc[ CDS_bedtool_df["gene_name"]==gene_name ]["gene_chrom"].values[0]
            gene_strand = CDS_bedtool_df.loc[ CDS_bedtool_df["gene_name"]==gene_name ]["gene_strand"].values[0]

            exon_starts = CDS_bedtool_df.loc[ CDS_bedtool_df["gene_name"]==gene_name ]["gene_blockStarts"].values[0].split(",")
            exon_starts = list( map(int, exon_starts))
            exon_sizes = CDS_bedtool_df.loc[ CDS_bedtool_df["gene_name"]==gene_name ]["gene_blockSizes"].values[0].split(",")
            exon_sizes = list(map(int, exon_sizes))

            cds_intervals = list()
            exon_idx = 0
            for exon_start in exon_starts:
                cds_intervals.extend( list(range( gene_start + exon_start, gene_start + exon_start + exon_sizes[exon_idx] )) )
                exon_idx += 1

            if gene_strand == "+":
                cds_idx = 0
                for start in cds_intervals:
                    if cds_idx % 3 == 0:
                        codon_idx = cds_idx // 3
                        cds_codonidxs.append( [chrom, start, start+3, gene_name, codon_idx, gene_strand])
                    cds_idx += 1

            elif gene_strand == "-":
                cds_idx = len(cds_intervals)
                for start in cds_intervals:
                    if cds_idx % 3 == 0:
                        codon_idx = cds_idx // 3 - 1
                        cds_codonidxs.append( [chrom, start, start+3, gene_name, codon_idx, gene_strand])
                    cds_idx -= 1

        ## export the df to a bed3 file plus two fields
        with open( self.prefix + '.cds_codonidxs.bed', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(cds_codonidxs)


## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="input gtf file")
    parser.add_argument("-r", help="input fasta file")

    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    
    ## process the file if the input files exist
    if (args.g!=None) & (args.r!=None):  # & (args.sort!=None) & (args.start!=None):
        print ("[status]\tprocessing the input file: " + args.g, flush=True)
        input_gtf = args.g
        input_ref = args.r
        gtf_hdl = Gtf2Bed(input_gtf, input_ref)
        gtf_hdl.convert_gtf()
        gtf_hdl.get_startcodon()
        gtf_hdl.get_fasta()

        print ("[status]\tFinished.", flush=True)

    else:
        print ("[error]\tmissing argument", flush=True)
        parser.print_usage() 