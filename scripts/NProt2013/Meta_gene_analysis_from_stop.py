"""
Supplementary Note 13: Meta-gene analysis from stop codon

Authors: Eugene Oh, Annemarie Becker

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

outputFile: 
average read density along an average transcript
    col0: position along the average transcript, counted from the stop codon
    col1: mean read density at that position

"""


def metagene(inputFileP, inputFileM, inputListP, inputListM, outputFile):
    
    mylist = range(-1500, 51)
    rangeDict = dict([i, 0] for i in mylist)
    countDict = dict([i, 0] for i in mylist)

    
  ### Plus Strand ###

  # Upload plus strand data in dictionary
    
    DictP = {}
    inFile = open(inputFileP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        DictP[col0] = col1
        line = inFile.readline()

  # Upload plus strand gene list
        
    inFile = open(inputListP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = str(fields[0])       #gene name
        col1 = int(fields[1])       #start of gene
        col2 = int(fields[2])       #stop of gene
        length = abs(col1 - col2) + 1

  # Select genes
        
        if length > 400:                    #use genes longer than 400 nt; this value can be changed
            normVal = 0                     
            for Z in range(col1, col2 + 1):
                if Z in DictP:
                    normVal += DictP[Z]     #determine expression level
            if normVal > 50:                #continue if sum of reads is greater than 50; this can be changed if desired
                meanNorm = normVal / length    #calculate read density (average reads per base)
                for K in range(-1500, 51):
                    elem = col2 + K
                    if elem in DictP and abs(K) <= length:
                        reads = DictP[elem] / meanNorm      #divide reads at one position by average number of reads per base for this gene                                   
                        rangeDict[K] += reads               #sum up reads at that position from all genes
                for K in countDict:                     #how often was this position counted in the calculation
                    if abs(K) <= length:
                        countDict[K] += 1
        elif length <= 400:
            pass
        line = inFile.readline()


  ### Minus Strand ###

  ## Upload minus strand data in dictionary

    DictM = {}
    inFile = open(inputFileM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = int(fields[0])
        col1 = float(fields[1])
        DictM[col0] = col1
        line = inFile.readline()

  # Upload minus strand gene list

    inFile = open(inputListM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = str(fields[0])
        col1 = int(fields[1])
        col2 = int(fields[2])
        length = abs(col1 - col2) + 1

  # Select genes
        
        if length > 400:
            normVal = 0
            for Z in range(col2, col1 + 1):
                if Z in DictM:
                    normVal += DictM[Z]
            if normVal > 50:
                meanNorm = normVal / length
                for K in range(-1500, 51):
                    elem = col2 - K
                    if elem in DictM and abs(K) <= length:
                        reads = DictM[elem] / meanNorm
                        rangeDict[K] += reads
                for K in countDict:
                    if abs(K) <= length:
                        countDict[K] += 1
        elif length <= 400:
            pass
        line = inFile.readline()
        
### Output data ###

    tupledlist1 = rangeDict.items()
    tupledlist1.sort()
    tupledlist2 = countDict.items()
    tupledlist2.sort()
    
    fullDict = {}
    zippedlist = zip(tupledlist1, tupledlist2)
    for elem in zippedlist:
        col0 = elem[0][0]       #list0 col0 = position (K)
        col1 = elem[0][1]       #list0 col1 = norm read number 
        col2 = elem[1][1]       #list1 col1 = how often was position counted
        fullDict[col0] = col1 / col2        #normalization2

        
  # Finish output

    tupledlist = fullDict.items()
    tupledlist.sort()
    
    outFile = open(outputFile, 'w')
    
    for J in tupledlist:
        for K in range(2):
            outFile.write(str(J[K]))
            if K < 1:
                outFile.write('\t')
        outFile.write('\n')




if __name__ == '__main__':
    inputFileP = ''
    inputFileM = ''
    inputListP = ''
    inputListM = ''
    outputFile = ''

    metagene(inputFileP, inputFileM, inputListP, inputListM, outputFile)

