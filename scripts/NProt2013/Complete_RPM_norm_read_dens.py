"""
Supplementary Note 7: Complete RPM-normalized read densities

Author: Annemarie Becker

inputFileP:
RPM-normalized read density file for plus strand (Supplementary Note 6)
    col0: position along genome
    col1: RPM-normalized read density at that position

inputFileM:
RPM-normalized read density file for minus strand (Supplementary Note 6)
    col0: position along genome
    col1: RPM-normalized read density at that position

outputFileP: 
complete RPM-normalized read density file for all positions along the genome on plus strand
    col0: position along genome
    col1: RPM-normalized read density at that position

outputFileM: 
complete RPM-normalized read density file for all positions along the genome on minus strand
    col0: position along genome
    col1: RPM-normalized read density at that position

"""


def RPMcomplete(inputFileP, inputFileM, outputFileP, outputFileM):
    
### PLUS STRAND###

    DictP = {}

    inFile = open(inputFileP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        DictP[col0] = col1
        line = inFile.readline()

    outFile = open(outputFileP, 'w')
    
    for elem in range(1,4578160):          #genome length, change for different genome
        if elem not in DictP:
            DictP[elem] = 0.0
        
    ListP = DictP.items()
    ListP.sort()
            
    for J in ListP:
        outFile.write(str(J[0]) + '\t' + str(J[1]) + '\n')
       
 
### MINUS STRAND###

    DictM = {}

    inFile = open(inputFileM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        DictM[col0] = col1
        line = inFile.readline()

    outFile = open(outputFileM, 'w')
    
    for elem in range(1,4578160):          #genome length, change for different genome
        if elem not in DictM:
            DictM[elem] = 0.0
        
    ListM = DictM.items()
    ListM.sort()
            
    for J in ListM:
        outFile.write(str(J[0]) + '\t' + str(J[1]) + '\n')


            
if __name__=='__main__':
    inputFileP = ''
    inputFileM = ''
    outputFileP = ''
    outputFileM = ''

    RPMcomplete(inputFileP, inputFileM, outputFileP, outputFileM)
    
