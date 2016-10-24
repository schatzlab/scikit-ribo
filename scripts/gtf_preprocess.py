#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for preprocessing gtf/bed files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 10/28/2016
## ----------------------------------------

from __future__ import print_function, division
import os
import sys
import csv
import argparse
import pybedtools as pbt
import pandas as pd
import numpy as np
import gffutils
from pybedtools.featurefuncs import gff2bed
from itertools import groupby


class gtf2Bed:
    ''' class to sort and get start codon from a gtf file
    '''
    def __init__(self, gtf, fasta, pairprob, tpm):
        self.gtf = gtf
        self.fasta = fasta
        self.pairprob = pairprob
        self.tpm = tpm
        self.bedtool = pbt.BedTool(self.gtf)
        self.base=os.path.basename(self.gtf)
        self.prefix = os.path.splitext(self.base)[0]
        self.gene_bed12s = []

    def convertGtf(self):
        ## create a db from the gtf file
        self.db = gffutils.create_db(self.gtf, ":memory:", force=True)

        ## retrieve a list of gene names with type "CDS" from db
        gene_names = []
        for gene in self.db.features_of_type("CDS", order_by=("seqid","start") ):
            gene_names.append(gene['gene_id'][0])
        gene_names = list(set(gene_names))

        ## convert a gtf/gff3 file to bed12 and save to a nested list
        for gene_name in gene_names:
            gene_bed12 = self.db.bed12(gene_name, name_field='gene_id')
            self.gene_bed12s.append(gene_bed12.split("\t"))

        ## extract CDS entries, write to a bed12 file
        with open( self.prefix + '.sort.CDS.bed', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(self.gene_bed12s)

    def getStartCodon(self):
        ## return a set of feature types
        feature_types = set(list(self.db.featuretypes()))

        ## sort by coordinates, extract start codon entries, save the records to a bed file
        if "start_codon" in feature_types:
            self.bedtool_sort = self.bedtool.sort() # sort the entries
            self.start_codon = self.bedtool_sort.filter(lambda x: x[2] == "start_codon")
            self.start_codon_bed = self.start_codon.each(gff2bed).saveas(self.prefix + '.sort.start.bed')
        else:
            gene_bed6s = []
            for bed12 in self.gene_bed12s:
                bed6 = bed12[0:6]
                bed6[2] = int(bed6[1]) + 3
                gene_bed6s.append(bed6)
            with open( self.prefix + '.sort.start.bed', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerows(gene_bed6s)

    def transformPairProb(self):
        ## read the pairing prob arrays then convert it to a df
        pairprob = []
        with open(self.pairprob, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                codon_idx = 0
                gene_name = row[0]
                for prob in row[1].split(" "):
                    pairprob.append([gene_name, codon_idx, float(prob)])
                    codon_idx += 1
        self.pairprob = pd.DataFrame(pairprob, columns= ["gene_name","codon_idx","pair_prob"])

    def fastaIter(self, fasta_fn):
        fh = open(fasta_fn)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        fasta_dict = {}
        for header in faiter:
            gene_name = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            seq_list = [seq[i:i+3] for i in range(0, len(seq), 3)]
            fasta_dict[gene_name] = seq_list
        return fasta_dict

    def getCodons(self):
        ## extract cds sequence from the ref genome, iterate the transcript sequence and yield codons
        cds_bt = pbt.BedTool(self.prefix + '.sort.CDS.bed')
        cds_bt.sequence(fi=self.fasta, fo= self.prefix + '.fasta', s=True, name=True,  split=True)
        self.cds_bt_df = cds_bt.to_dataframe(names=['chrom', 'start', 'end', 'gene_name', 'score', 'strand',
                                                    'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes',
                                                    'blockStarts'])
        self.codons_dict = self.fastaIter(self.prefix + '.fasta') # parse a fasta file to a list of codons

    def createCodonTable(self):
        ## construct the gene level df from the bed12 file
        codons = []
        pos_ranges_writer = open(self.prefix + "_pos_ranges.txt", "w")
        pos_ranges_writer.write("#gene_name\tstrand\tchrom\tpos_ranges\n")
        for gene_name in self.cds_bt_df.gene_name:
            ## get the info from df
            curr = self.cds_bt_df.loc[self.cds_bt_df["gene_name"] == gene_name]
            gene_start = curr["start"].values[0]
            chrom = curr["chrom"].values[0]
            gene_strand = curr["strand"].values[0]
            exon_starts = list(map(int, curr["blockStarts"].values[0].split(",")))
            exon_sizes = list(map(int, curr["blockSizes"].values[0].split(",")))

            ## extract the cds nt index, the order of exons, and save to a nested list
            pos_ranges = [(gene_start + exon_starts[i], gene_start + exon_starts[i] + exon_sizes[i])
                          for i in range(len(exon_starts))]

            ## add 10 codons upstream of a gene 
            upstream_counter = 10
            while upstream_counter > 0:
                if gene_strand == "+":
                    pos_ranges_upstream = [(pos_ranges[0][0]-30, pos_ranges[0][0])] + pos_ranges
                    for start in range(pos_ranges[0][0]-30, pos_ranges[0][0], 3):
                        codons.append([chrom, start, start+3, gene_name, -upstream_counter, gene_strand, "NA", 0])
                        upstream_counter -= 1
                else:
                    pos_ranges_upstream = [(pos_ranges[0][1]+30, pos_ranges[0][1])] + pos_ranges[::-1]
                    for end in range(pos_ranges[0][1]+30, pos_ranges[0][1], -3):
                        codons.append([chrom, end-3, end, gene_name, -upstream_counter, gene_strand, "NA", 0])
                        upstream_counter -= 1

            ## create position ranges file for look ups
            pos_ranges_writer.write(gene_name + "\t" + gene_strand + "\t" + chrom + "\t" +
                                    str(pos_ranges_upstream).strip('[]') + "\n")

            ## coding region of a gene
            if gene_strand == "+":
                starts = [start for start_range in pos_ranges for start in range(start_range[0], start_range[1])]
                for nt_idx in range(0, len(starts), 3):
                    codon_idx = int(nt_idx/3)
                    codon = self.codons_dict[gene_name][codon_idx]
                    pos = starts[nt_idx]
                    codons.append([chrom, pos, pos+3, gene_name, codon_idx, gene_strand, codon, 0])
            else:
                ends = [end for start_range in pos_ranges for end in range(start_range[0], start_range[1])]
                for nt_idx in range(0, len(ends), 3):
                    codon_idx = int(nt_idx/3)
                    codon = self.codons_dict[gene_name][codon_idx]
                    pos = ends[::-1][nt_idx] # make sure to reverse
                    codons.append([chrom, pos-3, pos, gene_name, codon_idx, gene_strand, codon, 0])

        ## convert nested list to df
        self.codons_df = pd.DataFrame(codons, columns=["chrom", "asite_start", "asite_end", "gene_name", "codon_idx",
                                                       "gene_strand", "codon", "reading_frame"])

    def mergeDf(self):
        ## import the salmon df, rna secondary structure, and merge with cds df
        tpm_df = pd.read_table(self.tpm,  header=0)
        codons_tpm = pd.merge(self.codons_df, tpm_df, how = "inner", left_on = ["gene_name"], right_on="Name")
        codons_tpm_pairprob = pd.merge(codons_tpm, self.pairprob, how = "inner", on = ["gene_name", "codon_idx"])
        codons_tpm_pairprob_out = codons_tpm_pairprob[["chrom", "asite_start", "asite_end", "gene_name", "codon_idx",
                                                       "gene_strand", "codon", "TPM", "pair_prob"]]
        codons_tpm_pairprob_out.to_csv(path_or_buf=self.prefix + '.cds_codons.bed', sep='\t', header=True, index=False)

## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="gtf file, required")
    parser.add_argument("-r", help="fasta file, required")
    parser.add_argument("-s", help="arrays of RNA secondary structure pairing probabilities, required")
    parser.add_argument("-t", help="pre-computed tpm salmon data-frame from RNA-seq data, required")

    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    
    ## process the file if the input files exist
    if (args.g!=None) & (args.r!=None) & (args.s!=None) & (args.t!=None):
        print ("[status]\tReading the input file: " + args.g, flush=True)
        input_gtf = args.g
        input_ref = args.r
        input_pairprob = args.s
        input_tpm = args.t
        ## execute
        print("[execute]\tStarting the pre-processing module", flush=True)
        gtf_hdl = gtf2Bed(input_gtf, input_ref, input_pairprob, input_tpm)
        print("[execute]\tConverting the the gtf file in to sql db", flush=True)
        gtf_hdl.convertGtf()
        print("[execute]\tExtracting the start codon from the gtf db", flush=True)
        gtf_hdl.getStartCodon()
        print("[execute]\tTransforming and preparing the df for RNA secondary structure pairing probabilities data", flush=True)
        gtf_hdl.transformPairProb()
        print("[execute]\tBuilding the index for each position at the codon level", flush=True)
        gtf_hdl.getCodons()
        print("[execute]\tCreating the codon table for the coding region", flush=True)
        gtf_hdl.createCodonTable()
        print("[execute]\tMerging all the df together", flush=True)
        gtf_hdl.mergeDf()
        ## finish
        print ("[status]\tPre-processing module finished.", flush=True)
    else:
        print ("[error]\tmissing argument", flush=True)
        parser.print_usage() 
