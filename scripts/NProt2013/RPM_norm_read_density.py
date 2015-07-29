"""
Supplementary Note 6: RPM-normalized read densities

Author: Annemarie Becker

inputFileP:
read density file for plus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

inputFileM:
read density file for minus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

inputNumber:
total read number as float (Supplementary Note 3)

outputFileP: 
RPM-normalized read density file for plus strand
    col0: position along genome
    col1: RPM-normalized read density at that position

outputFileM: 
RPM-normalized read density file for minus strand
    col0: position along genome
    col1: RPM-normalized read density at that position

"""


def norm(inputFileP, inputFileM, inputNumber, outputFileP, outputFileM):

### PLUS STRAND ###

    inFile = open(inputFileP, 'r')
    inNumber = open(inputNumber, 'r') 
    outFile = open(outputFileP, 'w')

    line = inFile.readline()
    number = inNumber.readline()
    totalReads = int(float(number))

    i = 0
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        RPM = col1 / totalReads * 1000000
        outFile.write(str(col0) + '\t' + str(RPM) + '\n')
        line = inFile.readline()


### MINUS STRAND ###

    inFile = open(inputFileM, 'r')
    inNumber = open(inputNumber, 'r') 
    outFile = open(outputFileM, 'w')

    line = inFile.readline()
    number = inNumber.readline()
    totalReads = int(float(number))

    i = 0
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        RPM = col1 / totalReads * 1000000
        outFile.write(str(col0) + '\t' + str(RPM) + '\n')
        line = inFile.readline()


            
if __name__=='__main__':
    inputFileP = ''
    inputFileM = ''
    inputNumber = ''
    outputFileP = ''
    outputFileM = ''

    norm(inputFileP, inputFileM, inputNumber, outputFileP, outputFileM)
