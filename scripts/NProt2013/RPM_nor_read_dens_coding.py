"""
Supplementary Note 8: RPM-normalized read densities in protein coding regions

Author: Annemarie Becker

inputFileP:
complete RPM-normalized read density file for plus strand (Supplementary Note 7)
    col0: position along genome
    col1: RPM-normalized read density at that position

inputFileM:
complete RPM-normalized read density file for minus strand (Supplementary Note 7)
    col0: position along genome
    col1: RPM-normalized read density at that position

inputListP:
E. coli MC4100 gene list of the plus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene

inputListM
E. coli MC4100 gene list of the minus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene

outputFileP: 
RPM-normalized read densities in protein coding regions on plus strand
    col0: position along genome
    col1: RPM-normalized read density at that position

outputFileM: 
RPM-normalized read densities in protein coding regions on minus strand
    col0: position along genome
    col1: RPM-normalized read density at that position

"""


def RPMinGenes(inputFileP, inputFileM, inputListP, inputListM, outputFileP, outputFileM):
    
### PLUS STRAND ###

  # Upload plus strand data as a dictionary

    DictP = {}

    inFile = open(inputFileP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        DictP[col0] = col1
        line = inFile.readline()

  # Upload plus strand gene list as a dictionary
    
    geneDictP = {}         #create dictionary with col0=gene name; col1=read number
    
    inFile = open(inputListP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = str(fields[0])           #plus strand gene list: gene name
        col1 = int(fields[1])           #plus strand gene list: start
        col2 = int(fields[2])           #plus strand gene list: stop
    
  # Select only positions in genes

        for Z in range(col1, col2 + 1):
            geneDictP[Z] = DictP[Z]
        line = inFile.readline()

    ListP = geneDictP.items()
    ListP.sort()

    outFile = open(outputFileP, 'w')
    for J in ListP:
        outFile.write(str(J[0]) + '\t' + str(J[1]) + '\n')


    
### MINUS STRAND ###

  # Upload minus strand data as a dictionary

    DictM = {}

    inFile = open(inputFileM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        DictM[col0] = col1
        line = inFile.readline()

  # Upload minus strand gene list as a dictionary
    
    geneDictM = {}         #create dictionary with col0=gene name; col1=read number
    
    inFile = open(inputListM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = str(fields[0])           #minus strand gene list: gene name
        col1 = int(fields[1])           #minus strand gene list: start
        col2 = int(fields[2])           #minus strand gene list: stop

  # Select positions in genes
    
        for Z in range(col2, col1 + 1):
            geneDictM[Z] = DictM[Z]
        line = inFile.readline()

    ListM = geneDictM.items()
    ListM.sort()

    outFile = open(outputFileM, 'w')
    for J in ListM:
        outFile.write(str(J[0]) + '\t' + str(J[1]) + '\n')


            
if __name__=='__main__':
    inputFileP = ''
    inputFileM = ''
    inputListP = ''
    inputListM = ''
    outputFileP = ''
    outputFileM = ''

    RPMinGenes(inputFileP, inputFileM, inputListP, inputListM, outputFileP, outputFileM)
    
