#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for processing bam files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 10/28/2016
## ----------------------------------------

from __future__ import print_function, division
import os
import sys
import argparse
import pybedtools as pbt
import pysam
import pandas as pd
import numpy as np
import csv
import errno
from datetime import datetime


class BamProcess(object):
    '''extracting alignment, prepare training/testing data
    '''
    def __init__(self, bam=None, mapq=20, minRL=10, maxRL=35, RelE=None, output=None, startCodons=None, cds=None,
                 posRanges=None, nts=None):
        self.bam = bam
        self.mapq = mapq
        self.output = output
        self.minRL = minRL
        self.maxRL = maxRL
        self.startCodons = startCodons
        self.cds = cds
        self.RelE = RelE
        self.posRanges = posRanges
        self.ntsFn = nts

    def filterRegion(self):
        # create bedtool, filter and sort
        self.startCodons = pbt.BedTool(self.startCodons).sort()
        self.cds = pbt.BedTool(self.cds).bed6().sort()
        # find overlapping regions
        distinctStartCodons = self.startCodons.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctCDSs = self.cds.merge(c=    "1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctStartCodons = distinctStartCodons.intersect(distinctCDSs, wa=True, sorted=True)
        # filter start codon
        startCodonHash = set([(i[0], i[1], i[2]) for i in distinctStartCodons])
        self.startCodons = self.startCodons.filter(lambda x: (x[0], x[1], x[2]) in startCodonHash)

    # TODO: add a function to see whether there are too many soft-clipped alignment
    def filterBam(self):
        # create a template/header from the input bam file
        inBam = pysam.AlignmentFile(self.bam, "rb")
        self.prefix = self.output + "/riboseq"
        outBam = pysam.AlignmentFile(self.prefix + ".bam", "wb", template=inBam)
        # read a bam file and extract info
        cigar_to_exclude = set([1,2,3,4,5]) #set(['I','D','S','H'])
        for read in inBam.fetch():
            cigars = set([c[0] for c in read.cigartuples])
            if read.mapping_quality > self.mapq and self.minRL <= read.query_length <= self.maxRL and \
            not cigars.intersection(cigar_to_exclude):
            #len(cigars.intersection(cigar_to_exclude)) == 0: # read.reference_id != 'chrM': # and 'N' not in edges:
                read.mapping_quality = read.query_length # change to bamtobed
                outBam.write(read)
        inBam.close()
        outBam.close()
        sys.stderr.write("[status]\tFinished filtering the bam file" + "\n")
        # save the bedtool to local
        self.bedtool = pbt.BedTool(self.prefix + ".bam")
        self.bedtool = self.bedtool.bam_to_bed(bed12=False)
        self.bedtool.saveas(self.prefix + '.bed')

    def posIndex(self):
        # create a dict to store the position read-frame and index info
        self.posOffsets, self.negOffsets = [], []
        with open(self.posRanges, 'r') as f:
            next(f) # skip header
            for line in f:
                gene, chr, strand, ranges = line.rstrip("\n").split("\t")
                boxes = [(int(i[0]), int(i[1]), int(i[2])) for i in [j.split(",") for j in ranges.split("|") ]]
                if strand == "+":
                    self.posOffsets.extend([[chr, pos, (abs(pos-(box[0]-15)) + box[2])%3] for box in boxes for pos in range(box[0]-15, box[1]+12)])
                else:
                    boxes = boxes[::-1] # flip the order
                    self.negOffsets.extend([[chr, pos, (abs(pos-(box[1]+15)) + box[2])%3] for box in boxes for pos in range(box[1]+15, box[0]-12, -1)])
        # convert to dataframe
        self.posOffsets = pd.DataFrame(self.posOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(subset=["chrom", "pos"])
        self.negOffsets = pd.DataFrame(self.negOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(subset=["chrom", "pos"])

    def sortBam(self):
        self.bedtool = pbt.BedTool(self.prefix + ".bam")
        self.bedtool = self.bedtool.bam_to_bed(bed12=False)
        self.bedtool = self.bedtool.sort()

    def makeTrainingData(self):
        # intersect with start codons
        self.bedtool = pbt.BedTool(self.prefix + '.bed')
        training = self.bedtool.intersect(self.startCodons, wa=True, wb=True, sorted=True, s=True)
        time = str(datetime.now())
        sys.stderr.write("[status]\tFinished intersecting the bedtool with start codons: " + time + "\n")
        # convert bedtool to df
        colNames = ['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'sc_chrom', 'sc_start', 'sc_end', 'gene',
                    'sc_score', 'gene_strand']
        training = training.to_dataframe(names=colNames)
        # a-site
        if not self.RelE:
            training['asite'] = np.where(training['gene_strand'] == '+', training['sc_start'] - training['start'] + 3,
                                         training['end'] - training['sc_end'] + 3 )
        else:
            training['asite'] = np.where(training['gene_strand'] == '+', training['end'] - training['sc_start'] - 3,
                                         training['sc_end'] - training['start'] - 3 )
        # phasing 5'
        tmpA = pd.merge(training, self.posOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
        tmpB = pd.merge(training, self.negOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
        training = pd.concat([tmpA, tmpB])
        training.rename(columns={'offset':'5_offset'}, inplace=True)
        # phasing 3'
        tmpA = pd.merge(training, self.posOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
        tmpB = pd.merge(training, self.negOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
        training = pd.concat([tmpA, tmpB])
        training.rename(columns={'offset':'3_offset'}, inplace=True)
        # filter a read by whether it has a-site that satisfies [9,18] or [1,8]
        if not self.RelE:
            training = training[((training['asite'] >= 9) & (training['asite'] <= 18))]
            training = training[(training['asite'] >= training['read_length'] / 2 - 1)]
        else:
            training = training[((training['asite'] >= 1) & (training['asite'] <= 5))]
        # get nts
        self.nts = pd.read_table(self.ntsFn, header=0)
        training['pos_-1'] = np.where(training['gene_strand'] == '+', training['start']-1,  training['end'])
        training['pos_0'] = np.where(training['gene_strand'] == '+', training['start'], training['end']-1)
        training['pos_n-1'] = np.where(training['gene_strand'] == '+', training['end']-1, training['start'])
        training['pos_n'] = np.where(training['gene_strand'] == '+', training['end'], training['start']-1)
        # merge training with nts
        training = pd.merge(training, self.nts, left_on=["chrom", "pos_-1"], right_on=["chrom", "pos"])
        training.drop(["pos_-1", "pos"], axis=1)
        training.rename(columns={'nt': 'nt_-1'}, inplace=True)
        training = pd.merge(training, self.nts, left_on=["chrom", "pos_0"], right_on=["chrom", "pos"])
        training.drop(["pos_0", "pos"], axis=1)
        training.rename(columns={'nt': 'nt_0'}, inplace=True)
        training = pd.merge(training, self.nts, left_on=["chrom", "pos_n-1"], right_on=["chrom", "pos"])
        training.drop(["pos_n-1", "pos"], axis=1)
        training.rename(columns={'nt': 'nt_n-1'}, inplace=True)
        training = pd.merge(training, self.nts, left_on=["chrom", "pos_n"], right_on=["chrom", "pos"])
        training.drop(["pos_n", "pos"], axis=1)
        training.rename(columns={'nt': 'nt_n'}, inplace=True)
        # slice the dataframe to the variables needed for training data, removed "start_nt", "end_nt"
        training = training[["asite", "read_length", "5_offset", "3_offset", "gene_strand", "nt_-1", "nt_0", "nt_n-1", "nt_n"]]
        return training

    def makeCdsData(self):
        # create pandas df from bedtools intersect
        self.bedtool = pbt.BedTool(self.prefix + '.bed')
        cds = self.bedtool.intersect(self.cds, wa=True, wb=True, sorted=True, s=True)
        time = str(datetime.now())
        sys.stderr.write("[status]\tFinished intersecting the bedtool within cds: " + time + "\n")
        colNames = ['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'gene_chrom', 'gene_start', 'gene_end',
                    'gene', 'gene_score', 'gene_strand']
        cds = cds.to_dataframe(names=colNames)
        # phasing 5'
        tmpA = pd.merge(cds, self.posOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        tmpB = pd.merge(cds, self.negOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        cds = pd.concat([tmpA, tmpB])
        cds.rename(columns={'offset':'5_offset'}, inplace=True)
        # phasing 3'
        tmpA = pd.merge(cds, self.posOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        tmpB = pd.merge(cds, self.negOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        cds = pd.concat([tmpA, tmpB])
        cds.rename(columns={'offset':'3_offset'}, inplace=True)
        # get nts
        cds['pos_-1'] = np.where(cds['gene_strand'] == '+', cds['start']-1,  cds['end'])
        cds['pos_0'] = np.where(cds['gene_strand'] == '+', cds['start'], cds['end']-1)
        cds['pos_n-1'] = np.where(cds['gene_strand'] == '+', cds['end']-1, cds['start'])
        cds['pos_n'] = np.where(cds['gene_strand'] == '+', cds['end'], cds['start']-1)
        # merge cds with nt
        cds = pd.merge(cds, self.nts, left_on=["chrom", "pos_-1"], right_on=["chrom", "pos"])
        cds.drop(["pos_-1", "pos"], axis=1, inplace=True)
        cds.rename(columns={'nt': 'nt_-1'}, inplace=True)
        cds = pd.merge(cds, self.nts, left_on=["chrom", "pos_0"], right_on=["chrom", "pos"])
        cds.drop(["pos_0", "pos"], axis=1, inplace=True)
        cds.rename(columns={'nt': 'nt_0'}, inplace=True)
        cds = pd.merge(cds, self.nts, left_on=["chrom", "pos_n-1"], right_on=["chrom", "pos"])
        cds.drop(["pos_n-1", "pos"], axis=1, inplace=True)
        cds.rename(columns={'nt': 'nt_n-1'}, inplace=True)
        cds = pd.merge(cds, self.nts, left_on=["chrom", "pos_n"], right_on=["chrom", "pos"])
        cds.drop(["pos_n", "pos"], axis=1, inplace=True)
        cds.rename(columns={'nt': 'nt_n'}, inplace=True)
        # slice the dataframe to the variables needed for training data
        cds = cds[["read_length", "5_offset", "3_offset", "gene_strand", "chrom", "start", "end", "nt_-1", "nt_0", "nt_n-1", "nt_n", 'strand']]
        return cds


# ----------------------------------------
#  main
# ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-p", help="prefix for BED/index files")
    parser.add_argument("-q", help="minimum mapq allowed, Default: 20", default=20, type=int)
    parser.add_argument("-l", help="shortest read length allowed, Default: 10", default=10, type=int)
    parser.add_argument("-u", help="longest read length allowed, Default: 35", default=35, type=int)
    parser.add_argument("-e", help="whether the sample involved RelE, Default: F", default='F', type=str)
    parser.add_argument("-o", help="output path")
    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # process the file if the input files exist
    if (args.i != None and args.p != None and args.o != None):
        # specify inputs
        inBam = args.i
        pre = args.o + "/" + args.p
        start, cds, posIdx, nt = pre + ".start.bed", pre + ".cds.bed", pre + ".pos_ranges.txt", pre + '.nt_table.txt'
        mapq = args.q
        minRL, maxRL = args.l, args.u
        RelE = False if args.e == 'F' else True
        output = args.o
        time = str(datetime.now())
        sys.stderr.write("[status]\tStart the module at " + time + "\n")
        sys.stderr.write("[status]\tProcessing the input bam file: " + inBam + "\n")
        sys.stderr.write("[status]\tOutput path: " + output + "\n")
        sys.stderr.write("[status]\tReading the start codon BED file: " + start + "\n")
        sys.stderr.write("[status]\tReading the open reading frame codon BED file: " + cds + "\n")
        sys.stderr.write("[status]\tReading the position-phase file: " + posIdx + "\n")
        # filter alignment
        sys.stderr.write("[execute]\tKeep reads with length [" + str(minRL) + ","+ str(maxRL) + "] and mapq >= " +
                         str(mapq) + "\n")
        aln = BamProcess(inBam, mapq, output, minRL, maxRL, start, cds, posIdx, RelE, nt)
        sys.stderr.write("[execute]\tFilter out overlapping regions" + "\n")
        aln.filterRegion()
        sys.stderr.write("[execute]\tImport the position ranges" + "\n")
        aln.posIndex()
        sys.stderr.write("[execute]\tFilter out un-reliable alignment from bam files" + "\n")
        aln.filterBam()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tCreate training dataframe at " + time + "\n")
        aln.makeTrainingData()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tCreate CDS dataframe at " + time + "\n")
        aln.makeCdsData()
        time = str(datetime.now())
        sys.stderr.write("[status]\tBam processing finished at " + time + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()

