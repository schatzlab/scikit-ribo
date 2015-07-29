"""
Supplementary Note 5: Compare gene expression levels

Authors: Eugene Oh, Annemarie Becker

inputFileP1:
read densities per gene on plus strand for sample 1 (Supplementary Note 4)
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: RPKM-normalized sum of read densities

inputFileM1:
read densities per gene on minus strand for sample 1 (Supplementary Note 4)
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: RPKM-normalized sum of read densities

inputFileP2:
read densities per gene on plus strand for sample 2 (Supplementary Note 4)
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: RPKM-normalized sum of read densities

inputFileM2:
read densities per gene on minus strand for sample 2 (Supplementary Note 4)
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: RPKM-normalized sum of read densities

inputNumber1:
total read number as float (Supplementary Note 3)

inputNumber2:
total read number as float (Supplementary Note 3)

outputFile: 
RPKM-normalized comparison of gene expression levels from samples 1 and 2
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: RPKM-normalized gene expression levels for sample 1
    col4: RPKM-normalized gene expression levels for sample 2
    
"""


def expression(inputFileP1, inputFileM1, inputFileP2, inputFileM2, inputNumber1, inputNumber2, outputFile):

  # Upload input files for sample 1
    
    list1 = []
    
    inFile = open(inputFileP1, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        list1.append(fields)
        line = inFile.readline()
        
    inFile = open(inputFileM1, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        list1.append(fields)	#add minus strand input file to the same list
        line = inFile.readline()
    
  # Upload input files for sample 2
    
    list2 = []
    
    inFile = open(inputFileP2, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        list2.append(fields)
        line = inFile.readline()
        
    inFile = open(inputFileM2, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        list2.append(fields)
        line = inFile.readline()

  # Compile both lists
        
    listA = zip(list1, list2)
        
  # Normalize read densities per gene by RPKM metric using a comparison cutoff value
  
    listB = []
  
    import math
    
    inNumber1 = open(inputNumber1, 'r')
    number1 = inNumber1.readline()
    totalReads1 = int(float(number1))

    inNumber2 = open(inputNumber2, 'r')
    number2 = inNumber2.readline()
    totalReads2 = int(float(number2))

    for Z in listA:			#Z = line
        val1 = float(Z[0][3])		#0th dimension (list1), 3rd column
        val2 = float(Z[1][3])		#1th dimension (list2), 3rd column
        gene_sum = val1 + val2		#sum of reads per gene from both experiments
        if gene_sum > 100:		#take only genes with at least 100 reads in both experiments --> this value has to be determined by variability analysis
            col0 = str(Z[0][0])		#gene name (of list1)
            col1 = int(Z[0][1])		#start position (of list1)
            col2 = int(Z[0][2])		#stop position (of list1)
            length = abs(col1 - col2) + 1   #gene length
            col3 = val1 / totalReads1 * 1000000 / length * 1000
            col4 = val2 / totalReads2 * 1000000 / length * 1000
            fields = (col0, col1, col2, col3, col4)
            listB.append(fields)
            
  # Output file
    
    outFile = open(outputFile, 'w')
    for J in listB:
        for K in range(5):
            outFile.write(str(J[K]))
            if K < 4:
                outFile.write('\t')
        outFile.write('\n')
        


if __name__ == '__main__':
    inputFileP1 = ''
    inputFileM1 = ''
    inputFileP2 = ''
    inputFileM2 = ''
    inputNumber1 = ''
    inputNumber2 = ''
    outputFile = ''
    
    expression(inputFileP1, inputFileM1, inputFileP2, inputFileM2, inputNumber1, inputNumber2, outputFile)
            
