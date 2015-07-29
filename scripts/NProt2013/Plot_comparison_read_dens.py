"""
Supplementary Note 11: Plot comparison of read density 

Author: Annemarie Becker

inputFileP: 
comparison of RPM-normalized read density files for protein coding regions on plus strand from samples 1 and 2 (Supplementary Note 9)
    col0: RPM-normalized read density of sample 1
    col1: RPM-normalized read density of sample 2

inputFileM: 
comparison of RPM-normalized read density files for protein coding regions on minus strand from samples 1 and 2 (Supplementary Note 9)
    col0: RPM-normalized read density of sample 1
    col1: RPM-normalized read density of sample 2

outputFile: 
graph

"""


import matplotlib.pyplot as plt
import numpy as np


def graph(inputFileP, inputFileM):

    x_series = []
    y_series = []

    inFileP = open(inputFileP, 'r')
    line = inFileP.readline()
    while line != '':
        fields = line.split()
        x_series.append(float(fields[0]))
        y_series.append(float(fields[1]))
        line = inFileP.readline()
        
    inFileM = open(inputFileM, 'r')
    line = inFileM.readline()
    while line != '':
        fields = line.split()
        x_series.append(float(fields[0]))
        y_series.append(float(fields[1]))
        line = inFileM.readline()

    fig = plt.figure(figsize=(8, 8))
    plt.plot(x_series, y_series, ',')   #use 'o' for larger, '.' for medium, ',' for smaller data points
    
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.xlabel("")
    plt.ylabel("")
    plt.title("")

    plt.savefig("")     #add folderpath and name; possible formats: png, pdf, ps, eps, svg

    


if __name__ == '__main__':
    inputFileP = ''
    inputFileM = ''

    graph(inputFileP, inputFileM)
