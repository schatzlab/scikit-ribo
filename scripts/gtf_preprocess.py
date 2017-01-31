#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for processing gtf/bed files
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


class gtf2Bed(object):
    ''' class to sort and get start codon from a gtf file
    '''
    def __init__(self, gtf=None, fasta=None, output=None):
        self.gtf = gtf
        self.fasta = fasta
        self.output = output
        self.geneBed12s = []
        self.geneNames = []
        self.bedtool = None
        
    def convertGtf(self):
        # parse prefix
        self.base = os.path.basename(self.gtf)
        self.prefix = self.output + "/" + os.path.splitext(self.base)[0]
        ## remove
        gtfDf = pd.read_table(self.gtf, comment="#", header=None)
        gtfDf = gtfDf[(gtfDf.iloc[:, 2] != 'gene') & (gtfDf.iloc[:, 2] != 'transcript')]
        gtfDedup = self.prefix + '.dedup.gtf'
        gtfDf.to_csv(path_or_buf=gtfDedup, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
        ## create bedtool from gtf
        self.bedtool = pbt.BedTool(gtfDedup)
        ## create a db from the gtf file
        self.db = gffutils.create_db(gtfDedup, ":memory:", force=True)
        ## retrieve a list of gene names with type "CDS" from db
        for gene in self.db.features_of_type("CDS", order_by=("seqid","start")):
            self.geneNames.append(gene['gene_id'][0])
        self.geneNames = list(set(self.geneNames))
        ## convert a gtf/gff3 file to bed12 and save to a nested list
        for geneName in self.geneNames:
            geneBed12 = self.db.bed12(geneName, name_field='gene_id')
            row = geneBed12.split("\t")
            self.geneBed12s.append([row[0], int(row[1]), int(row[2])] + row[3:])
        ## sort the file
        self.geneBed12s = sorted(self.geneBed12s, key = lambda x: (x[0], x[1], x[2]))
        ## extract CDS entries, write to a bed12 file
        with open(self.prefix + '.sort.CDS.bed', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(self.geneBed12s)

    def getChrLen(self):
        fh = open(self.fasta)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        self.chrDic = {}
        for header in faiter:
            chr = header.__next__()[1:].strip().split(" ")[0]
            seq = "".join(s.strip() for s in faiter.__next__())
            self.chrDic[chr] = len(seq)

    def getStartCodon(self):
        ## return a set of feature types
        featureTypes = set(list(self.db.featuretypes()))
        ## sort by coordinates, extract start codon entries, save the records to a bed file
        if "start_codon" in featureTypes:
            self.bedtool = self.bedtool.sort() # sort the entries
            self.startCodon = self.bedtool.filter(lambda x: x[2] == "start_codon")
            self.startCodonBed = self.startCodon.each(gff2bed).saveas(self.prefix + '.sort.start.bed')
        else:
            geneBed6s = []
            for bed12 in self.geneBed12s:
                bed6 = bed12[0:6]
                bed6[2] = int(bed6[1]) + 3
                geneBed6s.append(bed6)
            with open(self.prefix + '.sort.start.bed', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerows(geneBed6s)

    def fastaIter(self, fastaFn, mode):
        fh = open(fastaFn)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        codonDic, fastaDic, ntDic = {}, {}, {}
        # iterate
        for header in faiter:
            id = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
            nts = list(seq)
            codonDic[id] = codons
            fastaDic[id] = seq
            ntDic[id] = nts
        # return wrt to mode
        if mode == "codon":
            return codonDic
        elif mode == "seq":
            return fastaDic
        elif mode == "nt":
            return ntDic

    def fiveUtrBed(self, feature):
        fiveUtr = feature[:6]
        fiveUtr[1], fiveUtr[2] = int(fiveUtr[1]), int(fiveUtr[2])
        if feature.strand == "+":
            fiveUtr[1], fiveUtr[2] = max(0, fiveUtr[1]-24), fiveUtr[1]
        elif feature.strand == "-":
            fiveUtr[1], fiveUtr[2] = fiveUtr[2], min(fiveUtr[2]+24, self.chrDic[feature.chrom])
        return pbt.create_interval_from_list(fiveUtr)

    def threeUtrBed(self, feature):
        threeUtr = feature[:6]
        threeUtr[1], threeUtr[2] = int(threeUtr[1]), int(threeUtr[2])
        if feature.strand == "+":
            threeUtr[1], threeUtr[2] = threeUtr[2], min(threeUtr[2]+24, self.chrDic[feature.chrom])
        elif feature.strand == "-":
            threeUtr[1], threeUtr[2] = max(0, threeUtr[1]-24), threeUtr[1]
        return pbt.create_interval_from_list(threeUtr)

    def getSeq(self):
        cdsBt = pbt.BedTool(self.prefix + '.sort.CDS.bed')
        cdsBt.sequence(fi=self.fasta, fo= self.prefix + '.cds.fasta', s=True, name=True,  split=True)
        self.fastaDic = self.fastaIter(self.prefix + '.cds.fasta', "seq")
        # utr
        fiveUtrBt = cdsBt.each(self.fiveUtrBed)
        threeUtrBt = cdsBt.each(self.threeUtrBed)
        fiveUtrBt.sequence(fi=self.fasta, fo=self.prefix + '.5utr.fasta', s=True, name=True,  split=True)
        threeUtrBt.sequence(fi=self.fasta, fo=self.prefix + '.3utr.fasta', s=True, name=True,  split=True)
        self.fiveUtrDic = self.fastaIter(self.prefix + '.5utr.fasta', "seq")
        self.threeUtrDic = self.fastaIter(self.prefix + '.3utr.fasta', "seq")
        ## write expanded cds fasta to local
        expandedFasta = open(self.prefix + '.expandCDS.fasta', 'w')
        for geneName in self.geneNames:
            expandedFasta.write(">" + geneName + "\n" +
                                self.fiveUtrDic[geneName] + self.fastaDic[geneName] + self.threeUtrDic[geneName] + "\n")
        expandedFasta.close()

    def getCodons(self):
        ## extract cds sequence from the ref genome, iterate the transcript sequence and yield codons
        cdsBt = pbt.BedTool(self.prefix + '.sort.CDS.bed')
        ## create df
        colNames = ['chrom', 'start', 'end', 'gene_name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved',
                    'blockCount', 'blockSizes', 'blockStarts']
        self.cdsDf = cdsBt.to_dataframe(names=colNames)
        self.codonsDic = self.fastaIter(self.prefix + '.expandCDS.fasta', "codon") # parse a fasta file to a list of codons

    def createCodonTable(self):
        ## construct the gene level df from the bed12 file
        codons = []
        posRangesWriter = open(self.prefix + ".pos_ranges.txt", "w")
        posRangesWriter.write("#gene_name\tchr\tstrand\tpos_ranges\n")
        ## iterate over each gene
        for geneName in self.cdsDf.gene_name:
            ## get the info from df
            row = self.cdsDf.loc[self.cdsDf["gene_name"] == geneName]
            geneStart = row["start"].values[0]
            chrom = row["chrom"].values[0]
            geneStrand = row["strand"].values[0]
            exonStarts = list(map(int, row["blockStarts"].values[0].split(",")))
            exonSizes = list(map(int, row["blockSizes"].values[0].split(",")))
            ## extract the cds nt index, the order of exons, and save to a nested list
            posRanges = []
            if geneStrand == "+":
                phase, distance = 0, 0
                for i in range(len(exonStarts)):
                    posRanges.append((geneStart + exonStarts[i], geneStart + exonStarts[i] + exonSizes[i], phase))
                    distance += exonSizes[i]
                    phase = distance % 3
            elif geneStrand == "-":
                distance = 0
                for i in range(len(exonStarts)):
                    distance += exonSizes[i]
                    phase = (3 - (distance % 3)) % 3
                    posRanges.append((geneStart + exonStarts[i], geneStart + exonStarts[i] + exonSizes[i], phase))
            # add 8 codons upstream and downstream
            posRanges = [(posRanges[0][0]-24, posRanges[0][0], 0)] + posRanges + [(posRanges[-1][1], posRanges[-1][1]+24, 0)]
            ## create position ranges file for look ups
            posRangesStr = "|".join([ str(i[0]) + "," + str(i[1]) + "," + str(i[2]) for i in posRanges])
            posRangesWriter.write(geneName + "\t" + chrom + "\t" + geneStrand + "\t" + posRangesStr + "\n")
            ## upstream + coding region + downstream of a gene
            if geneStrand == "+":
                starts = [start for posRange in posRanges for start in range(posRange[0], posRange[1])]
                for ntIdx in range(0, len(starts), 3):
                    expandIdx = int(ntIdx/3)
                    codonIdx = expandIdx-8
                    codon = self.codonsDic[geneName][expandIdx]
                    pos = starts[ntIdx]
                    codons.append([chrom, pos, pos+3, geneName, codonIdx, geneStrand, codon])
            else:
                ends = [end+1 for posRange in posRanges for end in range(posRange[0], posRange[1])][::-1] # reverse
                for ntIdx in range(0, len(ends), 3):
                    expandIdx = int(ntIdx/3)
                    codon = self.codonsDic[geneName][expandIdx]
                    codonIdx = expandIdx-8
                    pos = ends[ntIdx]
                    codons.append([chrom, pos-3, pos, geneName, codonIdx, geneStrand, codon])
        ## convert nested list to df
        self.codonsDf = pd.DataFrame(codons, columns=["chrom", "asite_start", "asite_end", "gene_name", "codon_idx",
                                                      "gene_strand", "codon"])
        self.codonsDf.to_csv(path_or_buf=self.prefix + '.codonsDf.txt', sep='\t', header=True, index=False)


## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="gtf file, required")
    parser.add_argument("-r", help="fasta file, required")
    parser.add_argument("-o", help="output path, required")
    ## check if there is any argument
    if len(sys.argv) <= 1: 
        parser.print_usage() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    ## process the file if the input files exist
    if (args.g!=None) & (args.r!=None) & (args.o!=None):
        print ("[status]\tReading the input file: " + args.g, flush=True)
        input_gtf = args.g
        input_ref = args.r
        output = args.o
        # create output folder
        cmd = 'mkdir -p ' + output
        os.system(cmd)
        ## execute
        print("[execute]\tStarting the pre-processing module", flush=True)
        gtf_hdl = gtf2Bed(input_gtf, input_ref, output)
        print("[execute]\tLoading the the gtf file in to sql db", flush=True)
        gtf_hdl.convertGtf()
        print("[execute]\tCalculating the length of each chromosome", flush=True)
        gtf_hdl.getChrLen()
        print("[execute]\tExtracting the start codons' positions from the gtf db", flush=True)
        gtf_hdl.getStartCodon()
        print("[execute]\tExtracting the sequences for each gene", flush=True)
        gtf_hdl.getSeq()
        print("[execute]\tBuilding the index for each position at the codon level", flush=True)
        gtf_hdl.getCodons()
        print("[execute]\tCreating the codon table for the coding region", flush=True)
        gtf_hdl.createCodonTable()
        print ("[status]\tGtf pre-processing module finished.", flush=True)
    else:
        print ("[error]\tmissing argument", flush=True)
        parser.print_usage() 
