#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for preprocessing bam files
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
import multiprocessing
import errno
import pyximport; pyximport.install()
import skrCython as skrc
from timeit import default_timer as timer
from collections import defaultdict


class processAln(object):
    '''extracting alignment based on MAPQ and read length, prepare training/testing data
    '''

    def __init__(self, inBam, mapq, outBam, minRL, maxRL, startCodons, orf, posRanges):
        self.inBam = inBam
        self.mapq = mapq
        self.outBam = outBam
        self.minRL = minRL
        self.maxRL = maxRL
        self.reads = []
        self.startCodons = pbt.BedTool(startCodons).filter(lambda x: x.chrom != 'chrM')
        self.startCodons = self.startCodons.sort()
        self.orf = pbt.BedTool(orf).filter(lambda x: x.chrom != 'chrM')
        self.orf = self.orf.sort()
        self.posRanges = posRanges
        self.posDic = defaultdict(list)
        ## temp
        #with open(self.outBam + '.reads.txt', 'r') as f:
        #    for line in f:
        #        self.reads.append(line.rstrip("\n").split("\t"))

    def filterRegion(self):
        # find overlapping regions
        distinctStartCodons = self.startCodons.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctOrfs = self.orf.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctStartCodons = distinctStartCodons.intersect(distinctOrfs, wa=True, sorted=True)
        # filter start codon
        startCodonHash = set([(i[0], i[1], i[2]) for i in distinctStartCodons])
        self.startCodons = self.startCodons.filter(lambda x: (x[0], x[1], x[2]) in startCodonHash)
        # filter orf
        # distinctOrfs = self.orf.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
        #orfHash = set([(i[0], i[1], i[2]) for i in distinctOrfs])
        #self.orf = self.orf.filter(lambda x: (x[0], x[1], x[2]) in orfHash)

    ## TODO: add a function to see whether there are too many soft-clipped alignment
    def filterBam(self):
        # create a template/header from the input bam file
        inBamHdl = pysam.AlignmentFile(self.inBam, "rb")
        outBamHdl = pysam.AlignmentFile(self.outBam, "wb", template=inBamHdl)
        ## read a bam file and extract info
        for read in inBamHdl.fetch():
            cigar_to_exclude = set([1,2,3,4,5]) #set(['I','D','S','H'])
            cigars = [c[0] for c in read.cigartuples]
            # start_nt, end_nt = read.query_sequence[:2], read.query_sequence[-2:][::-1]
            # edges = set(list(start_nt + end_nt))
            # filter the bam file
            if read.mapping_quality > self.mapq and \
            self.minRL <= read.query_length <= self.maxRL and \
            not any(c in cigar_to_exclude for c in cigars): # and 'N' not in edges:
                outBamHdl.write(read)
                self.reads.append([read.query_name, read.query_length])
                # self.reads.append([read.query_name, read.query_length, start_nt, end_nt])
        inBamHdl.close()
        outBamHdl.close()
        # save the bedtool to local
        self.bedtool = pbt.BedTool(self.outBam)
        self.bedtool = self.bedtool.bam_to_bed(bed12=True).filter(lambda x: x.chrom != 'chrM')
        self.bedtool.saveas(self.outBam + '.bed')
        # export reads to local
        with open(self.outBam + '.reads.txt', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(self.reads)

    def posIndex(self):
        ## create a dict to store the position read-frame and index info
        with open(self.posRanges, 'r') as f:
            next(f)
            for line in f:
                gene, ranges = line.rstrip("\n").split("\t")
                self.posDic[gene] = [(int(i[0]), int(i[1]), int(i[2])) for i in [j.split(",") for j in ranges.split("|") ]]

    def sortBam(self):
        self.bedtool = pbt.BedTool(self.outBam)
        self.bedtool = self.bedtool.bam_to_bed(bed12=True)
        self.bedtool = self.bedtool.sort()
        #tmpfile = self.bedtool._tmp()
        #os.system('sort {0} -k1,1 -nk4,4 > {1}'.format(self.bedtool.fn, tmpfile))
        #self.bedtool = pbt.BedTool(tmpfile)

    def trainingData(self):
        # intersect with start codons
        self.bedtool = pbt.BedTool(self.outBam + '.bed')
        trainingData = self.bedtool.intersect(self.startCodons, wa=True, wb=True, sorted=True)

        trainingDf = trainingData.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart',
                                                      'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                                                      'sc_chrom', 'sc_start', 'sc_end', 'gene_name', 'sc_score',
                                                      'gene_strand'],
                                               dtype={'start':'int','end':'int', 'score':'int','sc_start':'int',
                                                      'sc_end':'int','sc_score':'int'} )
        self.readsDf = pd.DataFrame(self.reads, columns=['name', 'read_length'])
        #self.readsDf = pd.DataFrame(self.reads, columns=['name', 'read_length', 'start_nt', 'end_nt'])

        ## retrieve the read length and seq information from the bam file
        trainingDf = pd.merge(trainingDf, self.readsDf, on='name')
        ### a-site
        trainingDf['asite'] = np.where(trainingDf['gene_strand'] == '+',
                                       trainingDf['sc_start'] - trainingDf['start'] + 3,
                                       trainingDf['end'] - trainingDf['sc_end'] + 3 )
        ## phasing
        trainingDf['five_offset'] = trainingDf.apply(lambda x: self.getDistance(x['gene_name'], x['start'], x['end'],
                                                               x['gene_strand'], self.posDic, 5), axis=1)

        trainingDf['three_offset'] = trainingDf.apply(lambda x: self.getDistance(x['gene_name'], x['start'], x['end'],
                                                                x['gene_strand'], self.posDic, 3), axis=1)
        ## filter a read by whether it has a-site that satisfies [12,18]
        trainingDf = trainingDf[((trainingDf['asite'] >= 12) & (trainingDf['asite'] <= 18))]
        ## slice the dataframe to the variables needed for training data
        trainingDataOut = trainingDf[["asite", "read_length", "five_offset", "three_offset", "gene_strand"]]
        # trainingDataOut = trainingDf[["asite", "read_length", "five_offset", "three_offset", "start_nt", "end_nt", "gene_strand"]]
        trainingDataOut.to_csv(path_or_buf=self.outBam + '.asite.txt', sep='\t', header=True, index=False)

    def testingData(self):
        ## create pandas df from bedtools intersect
        self.bedtool = pbt.BedTool(self.outBam + '.bed')
        testingData = self.bedtool.intersect(self.orf, wa=True, wb=True, sorted=True)
        testingDf = testingData.to_dataframe(names=['chrom', 'start', 'end', 'name', 'mapq', 'strand', 'thickStart',
                                                    'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                                                    'gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score',
                                                    'gene_strand', 'gene_thickStart', 'gene_thickEnd', 'gene_itemRgb',
                                                    'gene_blockCount', 'gene_blockSizes', 'gene_blockStarts'],
                                             dtype={'start':'int','end':'int', 'mapq':'int', 'gene_start':'int',
                                                    'gene_end':'int', 'gene_score':'int',
                                                    'blockSizes':'object','blockStarts':'object',
                                                    'gene_blockSizes':'object','gene_blockStarts':'object'})
        ## compute offset
        #start = timer()
        testingDf['five_offset'] = testingDf.apply(lambda x: self.getDistance(x['gene_name'], x['start'], x['end'],
                                                             x['gene_strand'], self.posDic, 5), axis=1)
        testingDf['three_offset'] = testingDf.apply(lambda x: self.getDistance(x['gene_name'], x['start'], x['end'],
                                                              x['gene_strand'], self.posDic, 3), axis=1)
        #end = timer()
        #print(end-start, flush=True)
        testingDf = testingDf[(testingDf['five_offset'] != -1) & (testingDf['three_offset'] != -1)]
        testingDf = pd.merge(testingDf, self.readsDf, on='name', how='inner')
        ## slice the dataframe to the variables needed for training data
        testingDataOut = testingDf[["read_length", "five_offset", "three_offset", "gene_strand", "chrom", "start", "end"]]
        # testingDataOut = testingDf[["read_length", "five_offset", "three_offset", "start_nt", "end_nt", "gene_strand"]
        #                            "gene_name", "chrom", "start", "end", "name", "strand"]]
        testingDataOut.to_csv(path_or_buf=self.outBam + '.cds.txt', sep='\t', header=True, index=False)

    def getDistance(self, gene_name, start, end, gene_strand, dic, mode):
        pos_ranges = dic[gene_name]
        if gene_strand == "-": pos_ranges = pos_ranges[::-1]
        num_regions = len(pos_ranges)
        if mode == 5:
            pos = start + 15 if gene_strand == "+" else end - 15
            return self.offsetHelper(num_regions, pos_ranges, gene_strand, pos)
        elif mode == 3:
            pos = end - 12 if gene_strand == "+" else start + 12
            return self.offsetHelper(num_regions, pos_ranges, gene_strand, pos)
        else:
            exit(errno.EINVAL)

    def offsetHelper(self, numRegions, posRanges, geneStrand, pos):
            i = 0
            while i < numRegions:
                left, right, phase = posRanges[i]
                if left <= pos <= right:
                    up = left if geneStrand == "+" else right
                    offset = ( abs(pos - up) + phase) % 3
                    return offset
                else:
                    i += 1
            return -1


## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-p", help="prefix for BED/index files")
    parser.add_argument("-q", help="minimum mapq allowed, Default: 20", default=20, type=int)
    parser.add_argument("-l", help="shortest read length allowed, Default: 25", default=20, type=int)
    parser.add_argument("-u", help="longest read length allowed, Default: 35", default=35, type=int)
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
        in_bam_fn = args.i
        mapq = args.q
        out_bam_fn = args.o
        min_read_length = args.l
        max_read_length = args.u
        start_codons = args.p + ".sort.start.bed"
        orf = args.p + ".sort.CDS.bed"
        pos_ranges = args.p + ".pos_ranges.txt"
        print("[status]\tprocessing the input bam file: " + in_bam_fn, flush=True)
        print("[status]\toutput bam file name: " + out_bam_fn, flush=True)
        print("[status]\treading the start codon BED file: " + start_codons, flush=True)
        print("[status]\treading the open reading frame codon BED file: " + orf, flush=True)
        print("[status]\treading the position-phase file: " + pos_ranges, flush=True)

        ## filter alignment
        print("[execute]\tkeep only reads with length between [" + str(min_read_length) + "," + str(max_read_length) +
              "] and mapq of " + str(mapq), flush=True)
        aln = processAln(in_bam_fn, mapq, out_bam_fn, min_read_length, max_read_length, start_codons, orf, pos_ranges)
        print("[execute]\tfilter out overlapping regions", flush=True)
        aln.filterRegion()
        print("[execute]\timport the position ranges", flush=True)
        aln.posIndex()
        print("[execute]\tfilter out un-reliable alignment from bam files", flush=True)
        aln.filterBam()
        print("[execute]\tcreate dataframe - training data", flush=True)
        aln.trainingData()
        print("[execute]\tcreate dataframe - testing data", flush=True)
        aln.testingData()
        print("[status]\tBam pre-processing finished", flush=True)

    else:
        print("[error]\tmissing argument", flush=True)
        parser.print_usage()

