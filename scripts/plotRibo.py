#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for visualization riboseq data
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 10/28/2016
## ----------------------------------------

from __future__ import print_function, division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import pybedtools as pbt
import argparse
import multiprocessing


class figures(object):
    ''' define a class for plotting figures
    '''
    def __init__(self, fn=None, output=None):
        self.fn = fn
        self.output = output

    def loadDat(self):
        self.riboDf = pd.read_table(self.fn,  header=0)

    def plotCoverageOnGene(self, geneName):
        # check if the gene exist
        if geneName not in set(self.riboDf["gene"]):
            sys.stderr.write("[error]\tthe database does not have gene: " + str(geneName) + "\n")
            return
        # retrieve the data for this gene, construct numpy array
        geneDf = self.riboDf[self.riboDf["gene"] == geneName]
        riboCnt = np.array(geneDf["ribosome_count"])
        pairProb = np.array(geneDf["pair_prob"])
        # reverse the array if the strand is negative
        if geneDf["gene_strand"].values[0] == "-":
            riboCnt, pairProb = riboCnt[::-1], pairProb[::-1]
        # rolling mean of 2' structure pairing probabilities
        pairing = geneDf[["gene", "codon_idx", "pair_prob"]]
        rollingAvg = pairing.groupby("gene").rolling(window=11, center=True).mean().reset_index(0, drop=True)
        rollingAvg.columns = ["gene", "codon_idx", "avg_prob"]
        # merge 2' prob
        geneDf = pd.merge(geneDf, rollingAvg, how='left', on=["gene", 'codon_idx']).fillna(value=0)
        avgProb = geneDf["avg_prob"]
        # plot the ribosome count along a transcript
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        f.suptitle(geneName)
        ax1.plot(riboCnt, sns.xkcd_rgb["denim blue"], lw = 2)
        # indicate start/stop codons
        startCodonPos, stopCodonPos = 8, riboCnt.size-8
        ax1.axvline(startCodonPos, color="#999999",dashes=[3,2],zorder=-1)
        ax1.axvline(stopCodonPos ,color="#999999",dashes=[3,2],zorder=-1)
        ax1.set_ylabel("Ribosome counts")
        ax2 = plt.subplot(2, 1, 2)
        ax2.plot(pairProb, sns.xkcd_rgb["medium green"], label="One Codon")
        ax2.plot(avgProb, sns.xkcd_rgb["amber"], label="10-codon average")
        ax2.set_ylim([0,1])
        ax2.set_ylabel("2' Pairing probability")
        ax2.set_xlabel("Position in gene (5' -> 3')")
        ax2.axvline(startCodonPos, color="#999999",dashes=[3,2],zorder=-1)
        ax2.axvline(stopCodonPos,color="#999999",dashes=[3,2],zorder=-1)
        ax2.legend()
        plt.gcf()
        plt.savefig(self.output + "/plots/" + str(geneName) + ".riboCount.pdf")
        plt.clf()
        plt.cla()
        plt.close(f)

    def plotAllGenes(self):
        # loop over all genes and plot
        geneNames = set(self.riboDf["gene"])
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(self.plotCoverageOnGene, geneNames)
        sys.stderr.write("[status]\tFinished plotting all genes" + "\n")


## ----------------------------------------
## main
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input data frame, required")
    parser.add_argument("-g", help="gene of interest, or type [all] - automatically plot all genes, required")
    parser.add_argument("-o", help="output path, required")
    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    ## process the file if the input files exist
    if (args.i != None and args.g != None and args.o != None):
        sys.stderr.write("[status]\tprocessing the input file: " + args.i + "\n")
        df_fn = args.i
        gene = args.g
        output = args.o
        # create folders
        cmd = "mkdir -p " + output + "/plots"
        os.system(cmd)
        # all genes or one gene
        fig = figures(df_fn, output)
        fig.loadDat()
        if gene != "all":
            sys.stderr.write("[execute]\tplotting ribosome coverage for gene: " + str(gene) + "\n")
            fig.plotCoverageOnGene(gene)
        else:
            sys.stderr.write("[execute]\tplotting ribosome coverage for each gene" + "\n")
            fig.plotAllGenes()
        ## end
        sys.stderr.write("[status]\tPlotting module finished" + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
