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
from gtf_preprocess import gtf2Bed


class processAln(object):
    '''extracting alignment, prepare training/testing data
    '''
    def __init__(self, inBam, mapq, outBam, minRL, maxRL, startCodons, orf, posRanges, RelE):
        self.inBam = inBam
        self.mapq = mapq
        self.outBam = outBam
        self.minRL = minRL
        self.maxRL = maxRL
        self.startCodons = startCodons
        self.orf = orf
        self.RelE = RelE
        self.posRanges = posRanges

    def filterRegion(self):
        # create bedtool, filter and sort
        self.startCodons = pbt.BedTool(self.startCodons).filter(lambda x: x.chrom != 'chrM')
        self.startCodons = self.startCodons.sort()
        self.orf = pbt.BedTool(self.orf).filter(lambda x: x.chrom != 'chrM')
        self.orf = self.orf.sort()
        # find overlapping regions
        distinctStartCodons = self.startCodons.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctOrfs = self.orf.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctStartCodons = distinctStartCodons.intersect(distinctOrfs, wa=True, sorted=True)
        # filter start codon
        startCodonHash = set([(i[0], i[1], i[2]) for i in distinctStartCodons])
        self.startCodons = self.startCodons.filter(lambda x: (x[0], x[1], x[2]) in startCodonHash)

    ## TODO: add a function to see whether there are too many soft-clipped alignment
    def filterBam(self):
        # create a template/header from the input bam file
        inBamHdl = pysam.AlignmentFile(self.inBam, "rb")
        outBamHdl = pysam.AlignmentFile(self.outBam, "wb", template=inBamHdl)
        ## read a bam file and extract info
        cigar_to_exclude = set([1,2,3,4,5]) #set(['I','D','S','H'])
        for read in inBamHdl.fetch():
            cigars = set([c[0] for c in read.cigartuples])
            if read.mapping_quality > self.mapq and \
            self.minRL <= read.query_length <= self.maxRL and \
            len(cigars.intersection(cigar_to_exclude)) == 0 and \
            read.reference_id != 'chrM': # and 'N' not in edges:
                read.mapping_quality = read.query_length
                outBamHdl.write(read)
        inBamHdl.close()
        outBamHdl.close()
        sys.stderr.write("[status]\tFinished filtering the bam file")
        # save the bedtool to local
        self.bedtool = pbt.BedTool(self.outBam)
        self.bedtool = self.bedtool.bam_to_bed(bed12=True)
        self.bedtool.saveas(self.outBam + '.bed')

    def posIndex(self):
        ## create a dict to store the position read-frame and index info
        self.posOffsets = []
        self.negOffsets = []
        with open(self.posRanges, 'r') as f:
            next(f)
            for line in f:
                gene, chr, strand, ranges = line.rstrip("\n").split("\t")
                boxes = [(int(i[0]), int(i[1]), int(i[2])) for i in [j.split(",") for j in ranges.split("|") ]]
                if strand == "+":
                    self.posOffsets.extend([[chr, pos, (abs(pos-(box[0]-15)) + box[2])%3] for box in boxes for pos in range(box[0]-15, box[1]+12)])
                else:
                    boxes = boxes[::-1] # flip the order
                    self.negOffsets.extend([[chr, pos, (abs(pos-(box[1]+15)) + box[2])%3] for box in boxes for pos in range(box[1]+15, box[0]-12, -1)])
        ## convert to dataframe
        self.posOffsets = pd.DataFrame(self.posOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(subset=["chrom", "pos"])
        self.negOffsets = pd.DataFrame(self.negOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(subset=["chrom", "pos"])

    def sortBam(self):
        self.bedtool = pbt.BedTool(self.outBam)
        self.bedtool = self.bedtool.bam_to_bed(bed12=True)
        self.bedtool = self.bedtool.sort()

    def makeTrainingData(self):
        # intersect with start codons
        self.bedtool = pbt.BedTool(self.outBam + '.bed')
        trainingData = self.bedtool.intersect(self.startCodons, wa=True, wb=True, sorted=True)
        time = str(datetime.now())
        sys.stderr.write("[status]\tFinished intersecting the bedtool with start codons: " + time)
        # convert bedtool to df
        colNames = ['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
                    'blockCount', 'blockSizes', 'blockStarts', 'sc_chrom', 'sc_start', 'sc_end', 'gene_name',
                    'sc_score', 'gene_strand']
        trainingDf = trainingData.to_dataframe(names=colNames)
        ### a-site
        if not self.RelE:
            trainingDf['asite'] = np.where(trainingDf['gene_strand'] == '+',
                                           trainingDf['sc_start'] - trainingDf['start'] + 3,
                                           trainingDf['end'] - trainingDf['sc_end'] + 3 )
        else:
            trainingDf['asite'] = np.where(trainingDf['gene_strand'] == '+',
                                           trainingDf['end'] - trainingDf['sc_start'] - 3,
                                           trainingDf['sc_end'] - trainingDf['start'] - 3 )
        ## phasing 5'
        trainingA = pd.merge(trainingDf, self.posOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        trainingB = pd.merge(trainingDf, self.negOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        trainingDf = pd.concat([trainingA, trainingB])
        trainingDf.rename(columns={'offset':'five_offset'}, inplace=True)
        ## phasing 3'
        trainingA = pd.merge(trainingDf, self.posOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        trainingB = pd.merge(trainingDf, self.negOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        trainingDf = pd.concat([trainingA, trainingB])
        trainingDf.rename(columns={'offset':'three_offset'}, inplace=True)
        ## filter a read by whether it has a-site that satisfies [9,18] or [1,8]
        if not self.RelE:
            trainingDf = trainingDf[((trainingDf['asite'] >= 9) & (trainingDf['asite'] <= 18))]
            trainingDf = trainingDf[(trainingDf['asite'] >= trainingDf['read_length'] / 2 - 1)]
        else:
            trainingDf = trainingDf[((trainingDf['asite'] >= 1) & (trainingDf['asite'] <= 8))]
        ## slice the dataframe to the variables needed for training data, removed "start_nt", "end_nt"
        trainingDf = trainingDf[["asite", "read_length", "five_offset", "three_offset", "gene_strand"]]
        trainingDf.to_csv(path_or_buf=self.outBam + '.training.txt', sep='\t', header=True, index=False)

    def makeTestingData(self):
        ## create pandas df from bedtools intersect
        self.bedtool = pbt.BedTool(self.outBam + '.bed')
        testingData = self.bedtool.intersect(self.orf, wa=True, wb=True, sorted=True)
        time = str(datetime.now())
        sys.stderr.write("[status]\tFinished intersecting the bedtool with ORFs: " + time)
        colNames = ['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
                    'blockCount', 'blockSizes', 'blockStarts', 'gene_chrom', 'gene_start', 'gene_end', 'gene_name',
                    'gene_score', 'gene_strand', 'gene_thickStart', 'gene_thickEnd', 'gene_itemRgb', 'gene_blockCount',
                    'gene_blockSizes', 'gene_blockStarts']
        testingDf = testingData.to_dataframe(names=colNames,
                                             dtype={'blockSizes':'object','blockStarts':'object',
                                                    'gene_blockSizes':'object','gene_blockStarts':'object'})
        ## phasing 5'
        testingA = pd.merge(testingDf, self.posOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        testingB = pd.merge(testingDf, self.negOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        testingDf = pd.concat([testingA, testingB])
        testingDf.rename(columns={'offset':'five_offset'}, inplace=True)
        ## phasing 3'
        testingA = pd.merge(testingDf, self.posOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        testingB = pd.merge(testingDf, self.negOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        testingDf = pd.concat([testingA, testingB])
        testingDf.rename(columns={'offset':'three_offset'}, inplace=True)
        ## slice the dataframe to the variables needed for training data
        testingDf = testingDf[["read_length", "five_offset", "three_offset", "gene_strand", "chrom", "start", "end"]]
        testingDf.to_csv(path_or_buf=self.outBam + '.testing.txt', sep='\t', header=True, index=False)

## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-p", help="prefix for BED/index files")
    parser.add_argument("-q", help="minimum mapq allowed, Default: 20", default=20, type=int)
    parser.add_argument("-l", help="shortest read length allowed, Default: 20", default=20, type=int)
    parser.add_argument("-u", help="longest read length allowed, Default: 35", default=35, type=int)
    parser.add_argument("-e", help="whether the sample involved RelE, Default: F", default='F', type=str)
    parser.add_argument("-o", help="output filtered bam file")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i != None and args.p != None and args.o != None):
        ## specify inputs
        inBam = args.i
        mapq = args.q
        outBam = args.o
        RelE = False if args.e == 'F' else True
        minRL = args.l
        maxRL = args.u
        startCodons = args.p + ".sort.start.bed"
        orf = args.p + ".sort.CDS.bed"
        posRanges = args.p + ".pos_ranges.txt"
        time = str(datetime.now())
        sys.stderr.write("[status]\tStart the module at " + time)
        sys.stderr.write("[status]\tProcessing the input bam file: " + inBam)
        sys.stderr.write("[status]\tOutput bam file name: " + outBam)
        sys.stderr.write("[status]\tReading the start codon BED file: " + startCodons)
        sys.stderr.write("[status]\tReading the open reading frame codon BED file: " + orf)
        sys.stderr.write("[status]\tReading the position-phase file: " + posRanges)
        ## filter alignment
        sys.stderr.write("[execute]\tKeep reads with length [" + str(minRL) + ","+ str(maxRL) + "] and mapq >= " + str(mapq))
        aln = processAln(inBam, mapq, outBam, minRL, maxRL, startCodons, orf, posRanges, RelE)
        sys.stderr.write("[execute]\tFilter out overlapping regions")
        aln.filterRegion()
        sys.stderr.write("[execute]\tImport the position ranges")
        aln.posIndex()
        sys.stderr.write("[execute]\tFilter out un-reliable alignment from bam files")
        aln.filterBam()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tCreate dataframe - training data at " + time)
        aln.makeTrainingData()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tCreate dataframe - testing data at " + time)
        aln.makeTestingData()
        time = str(datetime.now())
        sys.stderr.write("[status]\tBam pre-processing finished at " + time)

    else:
        sys.stderr.write("[error]\tmissing argument")
        parser.print_usage()

