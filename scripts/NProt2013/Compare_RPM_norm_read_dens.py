"""
Supplementary Note 9: Compare RPM-normalized read densities 

Author: Annemarie Becker

inputFile1: 
RPM-normalized read densities along the whole genome or in protein coding regions on plus or minus strand from sample 1 (Supplementary Note 7 or 8)
    col0: position along genome
    col1: RPM-normalized read density at that position

inputFile2: 
RPM-normalized read densities along the whole genome or in protein coding regions on plus or minus strand from sample 2 (Supplementary Note 7 or 8)
    col0: position along genome
    col1: RPM-normalized read density at that position

outputFile: 
comparison of RPM-normalized read density files for protein coding regions on plus or minus strand from samples 1 and 2
    col0: RPM-normalized read density of sample 1
    col1: RPM-normalized read density of sample 2
    
"""


def matchRPM(inputFile1, inputFile2, outputFile):

  # Upload list of sample 1
    
    list1 = []
    
    inFile1 = open(inputFile1, 'r')
    line = inFile1.readline()
    while line != '':
        fields = line.split()
        list1.append(fields)
        line = inFile1.readline()
        
    
  # Upload list of sample 2
    
    list2 = []
    
    inFile2 = open(inputFile2, 'r')
    line = inFile2.readline()
    while line != '':
        fields = line.split()
        list2.append(fields)
        line = inFile2.readline()
        

  # Compile both lists
        
    listA = zip(list1, list2)

  # Output files

    outFile = open(outputFile, 'w')

    for Z in listA:
        position = int(Z[0][0])
        read1 = float(Z[0][1])
        read2 = float(Z[1][1])
        outFile.write(str(read1) + '\t' + str(read2) + '\n')

           
        
if __name__ == '__main__':
    inputFile1 = ''
    inputFile2 = ''
    outputFile = ''

    matchRPM(inputFile1, inputFile2, outputFile)           
