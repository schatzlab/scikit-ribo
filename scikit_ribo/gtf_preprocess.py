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


class GtfPreProcess(object):
    ''' class to sort and get start codon from a gtf file
    '''
    def __init__(self, gtf=None, fasta=None, prefix=None, output=None):
        self.gtf = gtf
        self.fasta = fasta
        self.prefix = prefix
        self.output = output
        self.geneBed12s = []
        self.geneNames = []
        self.bedtool = None
        
    def convertGtf(self):
        # remove gene/transcript entries, remove chrM
        gtfDf = pd.read_table(self.gtf, comment="#", header=None)
        gtfDf = gtfDf[(gtfDf.iloc[:, 2] != 'gene') & (gtfDf.iloc[:, 2] != 'transcript')]
        gtfDf = gtfDf[(gtfDf.iloc[:, 0] != 'chrM') & (gtfDf.iloc[:, 0] != 'chrMT') & (gtfDf.iloc[:, 0] != 'MT')]
        gtfDedup = self.output + "/" + self.prefix + '.dedup.gtf'
        gtfDf.to_csv(path_or_buf=gtfDedup, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
        # create bedtool from gtf and a db from the gtf file
        self.bedtool = pbt.BedTool(gtfDedup)
        self.db = gffutils.create_db(gtfDedup, ":memory:", force=True, verbose=False)
        # retrieve a list of gene names with type "CDS" from db
        for gene in self.db.features_of_type("CDS", order_by=("seqid","start")):
            self.geneNames.append(gene['gene_id'][0])
        self.geneNames = list(set(self.geneNames))
        # convert a gtf/gff3 file to bed12 and save to a nested list
        for geneName in self.geneNames:
            geneBed12 = self.db.bed12(geneName, name_field='gene_id')
            row = geneBed12.split("\t")
            self.geneBed12s.append([row[0], int(row[1]), int(row[2])] + row[3:])
        # sort the file
        self.geneBed12s = sorted(self.geneBed12s, key = lambda x: (x[0], x[1], x[2]))
        # extract CDS entries, write to a bed12 file
        with open(self.output + "/" + self.prefix + '.cds.bed', 'w', newline='') as csvfile:
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
        # return a set of feature types
        featureTypes = set(list(self.db.featuretypes()))
        # sort by coordinates, extract start codon entries, save the records to a bed file
        if "start_codon" in featureTypes:
            self.bedtool = self.bedtool.sort() # sort the entries
            self.startCodon = self.bedtool.filter(lambda x: x[2] == "start_codon")
            self.startCodonBed = self.startCodon.each(gff2bed).saveas(self.output + "/" + self.prefix + '.start.bed')
        else:
            geneBed6s = []
            for bed12 in self.geneBed12s:
                bed6 = bed12[0:6]
                bed6[2] = int(bed6[1]) + 3
                geneBed6s.append(bed6)
            with open(self.output + "/" + self.prefix + '.start.bed', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerows(geneBed6s)

    def fastaIter(self, fastaFn, mode):
        fh = open(fastaFn)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        codonDic, fastaDic, ntDic = {}, {}, {}
        # iterate
        for header in faiter:
            key = header.__next__()[1:].strip().split(" ")[0]
            seq = "".join(s.strip() for s in faiter.__next__())
            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
            nts = list(seq)
            codonDic[key] = codons
            fastaDic[key] = seq
            ntDic[key] = nts
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
        cdsBt = pbt.BedTool(self.output + "/" + self.prefix + '.cds.bed')
        cdsBt.sequence(fi=self.fasta, fo= self.output + "/" + self.prefix + '.cds.fasta', s=True, name=True,  split=True)
        self.fastaDic = self.fastaIter(self.output + "/" + self.prefix + '.cds.fasta', "seq")
        # utr
        fiveUtrBt = cdsBt.each(self.fiveUtrBed)
        threeUtrBt = cdsBt.each(self.threeUtrBed)
        fiveUtrBt.sequence(fi=self.fasta, fo=self.output + "/" + self.prefix + '.5utr.fasta', s=True, name=True,  split=True)
        threeUtrBt.sequence(fi=self.fasta, fo=self.output + "/" + self.prefix + '.3utr.fasta', s=True, name=True,  split=True)
        self.fiveUtrDic = self.fastaIter(self.output + "/" + self.prefix + '.5utr.fasta', "seq")
        self.threeUtrDic = self.fastaIter(self.output + "/" + self.prefix + '.3utr.fasta', "seq")
        # write expanded cds fasta to local
        expandedFasta = open(self.output + "/" + self.prefix + '.expandCDS.fasta', 'w')
        for geneName in self.geneNames:
            expandedFasta.write(">" + geneName + "\n" +
                                self.fiveUtrDic[geneName] + self.fastaDic[geneName] + self.threeUtrDic[geneName] + "\n")
        expandedFasta.close()

    def getCodons(self):
        # extract cds sequence from the ref genome, iterate the transcript sequence and yield codons
        cdsBt = pbt.BedTool(self.output + "/" + self.prefix + '.cds.bed')
        # create df
        colNames = ['chrom', 'start', 'end', 'gene', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved',
                    'blockCount', 'blockSizes', 'blockStarts']
        self.cdsDf = cdsBt.to_dataframe(names=colNames)
        self.codonsDic = self.fastaIter(self.output + "/" + self.prefix + '.expandCDS.fasta', "codon") # parse a fasta file to [codons]

    def createCodonTable(self):
        # construct the gene level df from the bed12 file
        codons = []
        posRangesWriter = open(self.output + "/" + self.prefix + ".pos_ranges.txt", "w")
        posRangesWriter.write("#gene\tchr\tstrand\tpos_ranges\n")
        # iterate over each gene
        for geneName in self.cdsDf.gene:
            # get the info from df
            row = self.cdsDf.loc[self.cdsDf["gene"] == geneName]
            geneStart = row["start"].values[0]
            chrom = row["chrom"].values[0]
            geneStrand = row["strand"].values[0]
            exonStarts = list(map(int, row["blockStarts"].values[0].split(",")))
            exonSizes = list(map(int, row["blockSizes"].values[0].split(",")))
            # extract the cds nt index, the order of exons, and save to a nested list
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
            # create position ranges file for look ups
            posRangesStr = "|".join([ str(i[0]) + "," + str(i[1]) + "," + str(i[2]) for i in posRanges])
            posRangesWriter.write(geneName + "\t" + chrom + "\t" + geneStrand + "\t" + posRangesStr + "\n")
            # upstream + coding region + downstream of a gene
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
        # convert nested list to df
        self.codonsDf = pd.DataFrame(codons, columns=["chrom", "start", "end", "gene", "codon_idx", "gene_strand", "codon"])
        self.codonsDf.to_csv(path_or_buf=self.output + "/" + self.prefix + '.codons.bed',
                             sep='\t', header=True, index=False)

    def getNts(self):
        # parse nts
        self.ntDic = self.fastaIter(self.fasta, "nt") # parse a fasta file to [nts]
        lst = []
        for chr in sorted(self.ntDic.keys()):
            nts = self.ntDic[chr]
            for i in range(len(nts)):
                nt = nts[i]
                lst.append([chr, i, i+1, nt])
        df = pd.DataFrame(lst, columns=["chrom", "start", "end", "nt"])
        bt = pbt.BedTool.from_dataframe(df).sort()
        # gene bedtool
        geneBt = pbt.BedTool(self.output + "/" + self.prefix + '.cds.bed')
        colNames = ['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
                    'blockCount', 'blockSizes', 'blockStarts']
        geneDf = geneBt.to_dataframe(names=colNames)
        cdsDf = geneDf[['chrom', 'start', 'end']]
        geneDf["start"] = cdsDf["start"] - 30
        geneDf["end"] = cdsDf["end"] + 30
        geneBt = pbt.BedTool.from_dataframe(geneDf)
        bt = bt.intersect(geneBt, wa=True, u=True, sorted=True)
        df = bt.to_dataframe(names=["chrom", "start", "end", "nt"])
        df = df[["chrom", "start", "nt"]].drop_duplicates(subset=["chrom", "start"])
        df.rename(columns={'start': 'pos'}, inplace=True)
        # export
        df.to_csv(path_or_buf=self.output + "/" + self.prefix + '.nt_table.txt', sep='\t', header=True, index=False)


# the main process
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
        sys.stderr.write("[status]\tReading the input file: " + args.g + "\n")
        gtf = args.g
        ref = args.r
        output = args.o
        # create output folder
        cmd = 'mkdir -p ' + output
        os.system(cmd)
        ## execute
        sys.stderr.write("[execute]\tStarting the pre-processing module" + "\n")
        worker = GtfPreProcess(gtf, ref, output)
        sys.stderr.write("[execute]\tLoading the the gtf file in to sql db" + "\n")
        worker.convertGtf()
        sys.stderr.write("[execute]\tCalculating the length of each chromosome" + "\n")
        worker.getChrLen()
        sys.stderr.write("[execute]\tExtracting the start codons' positions from the gtf db" + "\n")
        worker.getStartCodon()
        sys.stderr.write("[execute]\tExtracting the sequences for each gene" + "\n")
        worker.getSeq()
        sys.stderr.write("[execute]\tBuilding the index for each position at the codon level" + "\n")
        worker.getCodons()
        worker.getNts()
        sys.stderr.write("[execute]\tCreating the codon table for the coding region" + "\n")
        worker.createCodonTable()
        sys.stderr.write("[status]\tGtf processing module finished" + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage() 
