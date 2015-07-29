"""
Supplementary Note 3: Total reads

Author: Annemarie Becker

inputFileP:
read density file for plus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

inputFileM:
read density file for minus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

outputFile: 
total read number as float

"""


def countReads(inputFileP,inputFileM, outputFile):

    inFileP = open(inputFileP, 'r')
    inFileM = open(inputFileM, 'r')
    outFile = open(outputFile, 'w')

    line = inFileP.readline()
    i = 0
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        i = i+col1		#count reads on plus strand
        line = inFileP.readline()
    totReadsP = i

    line = inFileM.readline()
    j = 0
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        j = j+col1		#count reads on minus strand
        line = inFileM.readline()
    totReadsM = j

    totalReads = i + j
    outFile.write(str(totalReads))
        
        
if __name__ == '__main__':
    inputFileP = ''
    inputFileM = ''
    outputFile = ''

    countReads(inputFileP,inputFileM, outputFile)

