"""
Supplementary Note 14: Enrichment efficiency 

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
ratio of RPM-normalized read densities for protein coding regions on plus or minus strand from samples 1 and 2
    col0: position along genome
    col1: ratio of RPM-normalized read densities from samples 1 and 2
       
"""

def ratio(inputFile1, inputFile2, outputFile):

  # Upload input1

    Dict1 = {}
    inFile = open(inputFile1, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        Dict1[col0] = col1
        line = inFile.readline()
            
  # Upload input2
    
    Dict2 = {}
    inFile = open(inputFile2, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        Dict2[col0] = col1
        line = inFile.readline()


  # calculate window +-20

    outFile = open(outputFile, 'w')

    sum1 = 0
    sum2 = 0
    start = 1
    end = 4578159                  #change this value dependent on the genome
    start_sum = start + 20                  #start_sum = 21                        
    end_sum = end - 20                      #end_sum = 4578139
    

    for J in range(start, start_sum + 1):       #lines 1-21
        for X in range(start, J+20+1):          #sum 1-(J+20)
            sum1 += float(Dict1[X])
            sum2 += float(Dict2[X])
        if sum2 != 0:
            ratio = sum1 / sum2    
        else:
            ratio = 0.0
        outFile.write(str(J) + '\t' + str(ratio) + '\n')
        sum1 = 0
        sum2 = 0

    for K in range(start_sum + 1, end_sum + 1):         #lines 22-4578139
        for Y in range(K-20, K+20+1):                   #sum (K-20)-(K+20)
            sum1 += float(Dict1[Y])
            sum2 += float(Dict2[Y])
        if sum2 != 0:
            ratio = sum1 / sum2    
        else:
            ratio = 0.0
        outFile.write(str(K) + '\t' + str(ratio) + '\n')
        sum1 = 0
        sum2 = 0

    for L in range(end_sum + 1, end + 1):               #lines 4578140-4578159
        for Z in range(L-20, end + 1):
            sum1 += float(Dict1[Z])
            sum2 += float(Dict2[Z])
        if sum2 != 0:
            ratio = sum1 / sum2    
        else:
            ratio = 0.0
        outFile.write(str(L) + '\t' + str(ratio) + '\n')
        sum1 = 0
        sum2 = 0

        
        
if __name__ == '__main__':
    inputFile1 = ''
    inputFile2 = ''
    outputFile = ''

    ratio(inputFile1, inputFile2, outputFile)
