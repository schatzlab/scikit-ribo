"""
Supplementary Note 10: Pearson correlation coefficient 

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
Pearson correlation coefficient as float, 2-tailed p-value as float

"""


from scipy.stats.stats import pearsonr


def Pearson(inputFileP, inputFileM, outputFile):

    List1 = []
    List2 = []

    inFile = open(inputFileP, 'r')

    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = float(fields[0])
        col1 = float(fields[1])
        List1.append(col0)
        List2.append(col1)
        line = inFile.readline()

    inFile = open(inputFileM, 'r')
    
    line = inFile.readline()
    while line != '':
        fields = line.split()
        col0 = float(fields[0])
        col1 = float(fields[1])
        List1.append(col0)
        List2.append(col1)
        line = inFile.readline()

    correlation = pearsonr(List1, List2)   

    outFile = open(outputFile, 'w')
    outFile.write(str(correlation))

    

        
if __name__ == '__main__':
    inputFileP = ''
    inputFileM = ''
    outputFile = ''

    Pearson(inputFileP, inputFileM, outputFile)
