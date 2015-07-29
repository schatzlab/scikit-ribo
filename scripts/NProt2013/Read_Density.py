"""
Supplementary Note 4: Read density per gene

Authors: Eugene Oh

inputFileP:
read density file for plus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

inputFileM:
read density file for minus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

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
read densities per gene on plus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: sum of read densities

outputFileM: 
read densities per gene on minus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: sum of read densities

"""


def expression(inputFileP, inputFileM, inputListP, inputListM, outputFileP, outputFileM):

### PLUS STRAND ###

  # Upload read density file from plus strand as a dictionary
    
    inFileDictP = {}         #create dictionary that looks like input file
    
    inFile = open(inputFileP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        inFileDictP[col0] = col1
        line = inFile.readline()
        
  # Upload plus strand gene list as a dictionary and list
    
    geneDictP = {}         #create dictionary with col0=gene name; col1=read number
    geneListP = []         #create list that looks like input gene list
    
    inFile = open(inputListP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        geneListP.append(fields)        #add an item to the end of the list
        col0 = str(fields[0])           #plus strand gene list: gene name
        col1 = int(fields[1])           #plus strand gene list: start
        col2 = int(fields[2])           #plus strand gene list: stop

  # Sum up and write read densities per protein coding region in dictionary

        for Z in range(col1, col2 + 1):
            if Z in inFileDictP and col0 in geneDictP:
                geneDictP[col0] += inFileDictP[Z]
            elif Z in inFileDictP:
                geneDictP[col0] = inFileDictP[Z]
        line = inFile.readline()
    
  # Assign gene expression levels to all genes
            
    tupledlistP = geneDictP.items()        
    for J in geneListP:
        match = 0
        for K in tupledlistP:
            if J[0] == K[0]:        
                match = 1
                J.append(K[1])
        if match == 0:		#list genes that don't have any reads
            J.append(0)
            
  # Output file for plus strand
    
    outFile = open(outputFileP, 'w')
    for J in geneListP:
        outFile.write(str(J[0]) + '\t' + str(J[1]) + '\t' + str(J[2]) + '\t' + str(J[3]) + '\n')

        
### MINUS STRAND ###     
    
  # Upload read density file from minus strand as a dictionary
    
    inFileDictM = {}
    
    inFile = open(inputFileM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        inFileDictM[col0] = col1
        line = inFile.readline()
        
  # Upload minus strand gene list as a dictionary and list
    
    geneDictM = {}	#create dictionary with col0=gene name; col1=read number
    geneListM = []	#create list that looks like input gene list
    
    inFile = open(inputListM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        geneListM.append(fields)	#add an item to the end of the list
        col0 = str(fields[0])           #minus strand gene list: gene name
        col1 = int(fields[1])           #minus strand gene list: start
        col2 = int(fields[2])           #minus strand gene list: stop

  # Sum up and write read densities per protein coding region in dictionary

        for Z in range(col2, col1 + 1):
            if Z in inFileDictM and col0 in geneDictM:
                geneDictM[col0] += inFileDictM[Z]
            elif Z in inFileDictM:
                geneDictM[col0] = inFileDictM[Z]
        line = inFile.readline()
    
  # Assign gene expression levels to all genes
    
    tupledlistM = geneDictM.items()
    for J in geneListM:
        match = 0
        for K in tupledlistM:
            if J[0] == K[0]:
                match = 1
                J.append(K[1])
        if match == 0:
            J.append(0)
            
  # Output file for minus strand
    
    outFile = open(outputFileM, 'w')
    for J in geneListM:
        outFile.write(str(J[0]) + '\t' + str(J[1]) + '\t' + str(J[2]) + '\t' + str(J[3]) + '\n')

    
if __name__ == '__main__':
    inputFileP = ''
    inputFileM = ''
    inputListP = ''
    inputListM = ''
    outputFileP = ''
    outputFileM = ''

    expression(inputFileP, inputFileM, inputListP, inputListM, outputFileP, outputFileM)

    
