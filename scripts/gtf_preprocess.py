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

from __future__ import print_function
import os
import sys
import csv
import argparse
import pybedtools as pbt
from pybedtools.featurefuncs import gff2bed
import numpy as np
import gffutils


class Gtf2Bed:
    ''' class to sort and get start codon from a gtf file '''
    def __init__(self, gtf, fasta):
        self.gtf = gtf
        self.fasta = fasta
        self.bedtool = None
        self.start = None
        self.out_bed = None
        self.bedtool = pbt.BedTool(self.gtf)
        self.fn = os.path.splitext(self.gtf)[0]

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
        self.gene_bed12s = list()
        for gene_name in gene_names:
            gene_bed12 = self.db.bed12(gene_name, name_field='gene_id')
            self.gene_bed12s.append(gene_bed12.split("\t"))

        ## extract CDS entries, write to a bed12 file
        with open( self.fn + '.sort.CDS.bed', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(self.gene_bed12s)

    def get_startcodon(self):
        ## return a set of feature types
        feature_types = set(list(self.db.featuretypes()))

        ## sort by coordinates, extract start codon entries, save the records to a bed file
        if "start_codon" in feature_types:
            self.bedtool_sort = self.bedtool.sort()
            self.start = self.bedtool_sort.filter(lambda x: x[2] == "start_codon")
            self.start_bed = self.start.each(gff2bed).saveas(self.fn + '.sort.start.bed')
        else:
            gene_bed6s = list()
            for bed12 in self.gene_bed12s:
                bed6 = bed12[0:6]
                bed6[2] = int(bed6[1]) + 3
                gene_bed6s.append(bed6)
            with open( self.fn + '.sort.start.bed', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerows(gene_bed6s)

    def get_fasta(self):
        ## extract cds sequence from the reference genome
        cds_fa = self.bedtool_sort.sequence(fi=self.fasta, s=True, name=True, tab=True, split=True)
        # print(open(cds_fa.seqfn).read())

        ## calculate the ribosome position, +/- strand is different
        self.bedtool_sort_df = self.bedtool_sort.to_dataframe(names=['gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score',
                                                 'gene_strand'])

        ## calculate the transcript length in codon space
        self.bedtool_sort_df["length_in_codons"] = (self.bedtool_sort_df["gene_start"] - self.bedtool_sort_df["gene_end"]) / 3
        self.bedtool_sort_df["length_in_codons"] = self.bedtool_sort_df["length_in_codons"].abs().astype(int)

        ## create an empty df for the entire set of transcripts
        gene_size_df =  self.bedtool_sort_df[["gene_name","length_in_codons"]].drop_duplicates()
        gene_list = list()
        for i in gene_size_df.gene_name:
            gene_length=gene_size_df.loc[ gene_size_df["gene_name"]==i ]["length_in_codons"].values[0]
            for j in range(0, gene_length):
                               gene_list.append([i,j])


        '''

        ## create an empty df for the entire set of transcripts
        gene_size_df =  self.testing_df_out[["gene_name","length_in_codons"]].drop_duplicates()
        gene_list = list()
        for i in gene_size_df.gene_name:
            gene_length=gene_size_df.loc[ gene_size_df["gene_name"]==i ]["length_in_codons"].values[0]
            for j in range(0, gene_length):
                               gene_list.append([i,j])
        '''

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
        # gtf_hdl.get_fasta()

        print ("[status]\tFinished.", flush=True)

    else:
        print ("[error]\tmissing argument", flush=True)
        parser.print_usage() 