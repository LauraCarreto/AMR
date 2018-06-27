#!/usr/bin/env python
'''
This script parses data from Kmer Resistance and derives a matrix to compare to phenotypic data.

Requires:
1 - Kmer Resistance output files.
2 - Rules file, setting a correspondence between gene name and AMR.
    Rules file format: column 1 = Gene name; column 2 = list of antibiotics with elements separated by a comma.
3 - Phenotypic matix file, with the reference AMR.
    Phenotypic matrix file format: Line 1 = headers; other lines = samples; column 1 header = samples;
    other columns' headers = antibiotic names, corresponding to names in rules (check spelling!).

'''
import sys
import os
import re
import argparse
import glob
import fileinput
import numpy as np
import matplotlib.pyplot as plt
import itertools

__version__ = '0.1'
__date__ = 'April2017'
__author__ = 'laura'

OUTDIRNAME = "KmerMatrixOut"

# -------------------------------------------------------------------------------------------------

def parse_args():
    """
    Parse arguments

    Parameters
    ----------
    no inputs

    Returns
    -------
    oArgs: obj
        arguments object
    """

    sDescription = 'version %s, date %s, author %s' %(__version__, __date__, __author__)

    parser = argparse.ArgumentParser(description=sDescription)

    parser.add_argument("--input",
                        "-i",
                        type=str,
                        metavar="STRING",
                        dest="input",
                        required=True,
                        help="The folder where the kmer input files are found.")

    parser.add_argument("--ref",
                        "-r",
                        type=str,
                        metavar="STRING",
                        dest="ref",
                        required=True,
                        help="The folder where the phenotypic matrix file and the rules file are found.")

    parser.add_argument("--output",
                        "-o",
                        type=str,
                        metavar="STRING",
                        dest="output",
                        required=True,
                        help="The folder where the output matrix file will be stored.")


    oArgs = parser.parse_args()
    return oArgs

# -------------------------------------------------------------------------------------------------


def main():
    '''
    Main funtion

    Parameters
    ----------
    no inputs

    Returns
    -------
    0 on success else 1
    '''

    oArgs = parse_args()

    bDirExisted = False
    outputFolder = os.path.join(oArgs.output, OUTDIRNAME)
    try:
        os.mkdir(outputFolder)
    except OSError:
        bDirExisted = True

    aStrainFileNames = glob.glob(oArgs.input + "/*_ResOutput.out")
    print "Number of samples:", len(aStrainFileNames)

    rulesFile = oArgs.ref + "/rules_for_salmonella.txt"
    #check if rules file exists and exit script if it doesn't
    try:
        with open(rulesFile) as file:
            print 'AMR rules file', rulesFile
            pass
    except IOError:
        print 'Unable to open rules file; check if rules_for_salmonella.txt exists in input directory.'
        return 0

    phenoMatrixFile = oArgs.ref + "/phenotypic_matrix.txt"
    #phenoMatrixFile = oArgs.ref + "/kmer_res_matrix_ali.txt" #to test against other kmer matrix output

    #check if phenotypic matrix file exists and exit script if it doesn't
    try:
        with open(phenoMatrixFile) as file:
            print 'Phenotypic matrix file', phenoMatrixFile
            pass
    except IOError:
        print 'Unable to open phenotypic matrix file; check if my_phenotypic_matrix.txt exists in input directory.'
        return 0

    amList = ['Amikacin',
            'Amoxicillin',
            'Ampicillin',
            'Apramycin',
            'Cefotaxime',
            'Ceftazidime',
            'Chloramphenicol',
            'Ciprofloxacin',
            'Furazolidone',
            'Gentamicin',
            'Nalidixic-acid',
            'Neomycin',
            'Streptomycin',
            'Sulphonamide-compounds',
            'Sulphamethoxazole-trimethoprim',
            'Tetracycline']


    get_genes_in_sample(aStrainFileNames)

    get_rules(rulesFile)

    make_matrix(amList, dSample2Genes, dGene2AMR, outputFolder)

    compare_matrices(phenoMatrixFile, dResult, outputFolder)

    confusion_matrix(dPhenoData, dResult, outputFolder)

    return 0

#--------------------------------------------------------------------------------
def get_genes_in_sample(aStrainFileNames):

    #open each file and make a dictionary of sample:genes
    global dSample2Genes
    dSample2Genes ={}
    for filename in aStrainFileNames:
        #get name of sample on separator '_'
        sample_id =[x.strip() for x in os.path.basename(filename).split("_")][0]
        #print sample_id

        genes=[]
        #open file, look for line that starts with sample id and get genes into dictionary
        lines = iter(fileinput.input([filename]))
        for line in lines:
            if line.startswith(('#', 'gi')):
                next
            else:
                gene = line.split ("\t")[0]
                clean_gene = re.sub(r"\s+", "", gene, flags=re.UNICODE) #strip spaces from gene name
                #clean_gene = re.sub('[()]', '', clean_gene1) #strip parentheses from gene name
                genes.append(clean_gene)

        dSample2Genes[sample_id] = genes

    #for key in  dSample2Genes.keys():
        #if key == 'L01667-06':
            #print dSample2Genes[key]

    return

#---------------------------------------------------------------------------------
def get_rules(rulesFile):

    #assert file is present and contains genes, if not exit with message to user

    lines = iter(fileinput.input([rulesFile]))

    global dGene2AMR
    dGene2AMR = {}
    #mdrList = []
    for line in lines:
        details = line.split("\t")
        gene = details[0] #original gene name
        #gene = gene.translate(None, "()'-") #to remove characters "()'-" fom string, if required
        amr = line.split("\t")[1]
        amr = amr.rstrip()
        dGene2AMR[gene] = amr
        #print gene, amr

    #print dGene2AMR['tet(A)']
    return
#-----------------------------------------------------------------------------------
def make_matrix(amList, dSample2Genes, dGene2AMR, outputFolder):

    noRules_list = [] #for genes that are not in the rules
    generalNames = ['Qnr', 'qnr', 'cml'] #list of general names in rules

    global dResult
    dResult = {} #dictionary linking sample to AMR

    for sample in sorted(dSample2Genes.keys()):#for each sample:genes
        resistanceList = [] #create list of AMR for sample
        for gene in dSample2Genes[sample]: #for each gene in sample
            if gene in sorted(dGene2AMR.keys()):#if gene in rules
                #print gene
                res = dGene2AMR[gene] #get resistance related to gene
                #print gene, res
            else:
                #if gene ends in numbers, remove the numbers and compare the remaining string to list generalNames
                if re.search(r'\d+$', gene) is not None: #if end of string has numbers
                    split_gene = re.split('\d+$', gene)[0] #remove numbers

                    if split_gene in sorted(dGene2AMR.keys()): #and check again if in dGene2AMR rules
                        res = dGene2AMR[split_gene]

                    elif split_gene.startswith(tuple(generalNames)): #if split_gene starts with string in list generalNames
                        res = dGene2AMR[split_gene[:3]] #get first 3 characters and check again if in dGene2AMR rules

                    else:
                        if gene in noRules_list:
                            next
                        else:
                            noRules_list.append(gene)

            resFromGene = res.split(', ') #split string into list

            #append resFromGene to dictionary of sample: list of AMR
            for each in resFromGene:
                if each in resistanceList:
                    next
                else:
                    resistanceList.append(each)

        #print sample, resistanceList

        #for each antimicrobial in amList, check if matches to element in resistanceList: if yes, append 1 to dResult[sample], else, append 0
        dResult[sample]= [] # initialise list for sample in dictionary dSample2Resistance
        for am in amList: #for each coumpound in amList
            if am in resistanceList:
                result = '1'
            else:
                result = '0'

            dResult[sample].append(result)

    #print dResult

    print 'List of genes not in rules:', noRules_list

    file_matrix(dResult, amList, outputFolder)

    plot_matrix(dResult, amList, outputFolder, 'Greys',
                title = 'Kmer_matrix',
                xlabel = 'Antimicrobial resistance (in black)')


    return
#----------------------------------------------------------------------------------
def file_matrix(dResult, amList, outputFolder):

    matrixFile = os.path.join(outputFolder, 'Kmer_matrix.txt') # output file
    print 'Resistance matrix file', matrixFile

    #open output file
    openFile = open(matrixFile, 'w')
    #write header line
    amList.insert(0, 'Sample')
    print >> openFile, '\t'.join(amList)
    #append a line for each sample with sample name and binary code for presence or absence of resistance
    for sample in sorted(dResult.keys()):
        print >> openFile, sample + '\t'+'\t'.join(dResult[sample])

    openFile.close()

    return
#----------------------------------------------------------------------------------

def plot_matrix(dResult, amList, outputFolder, cmap, title, xlabel):

    #create a 2D array for plotting - need to convert '1' and '0' strings to int using map(int(myList))
    array = np.array([map(int, dResult[key]) for key in sorted(dResult.keys())])
    #print array

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #different color range if cmap ='Greys' or cmap = 'bwr'
    if cmap == 'Greys':
        ax.matshow(array, cmap = cmap, aspect = 'auto')
    elif cmap == 'bwr': #set limits for color range in difference matrix plot
        ax.matshow(array, cmap = cmap, vmin = -1, vmax = 1, aspect = 'auto')

    plt.xlabel(xlabel, fontsize = 8)
    plt.ylabel('Sample', fontsize = 8)

    #Use MultipleLocator ticker class to place ticks on each column
    from matplotlib.ticker import MultipleLocator
    majorLocator = MultipleLocator(1)
    ax.xaxis.set_major_locator(majorLocator)
    #format ticks and adjust space
    xtickNames = plt.setp(ax, xticklabels= amList)
    plt.setp(xtickNames, rotation = 30, ha='left', fontsize = 6)
    #general layout
    plt.suptitle(title)
    plt.tight_layout()
    #plt.show()

    #print to outputFolder
    from matplotlib.backends.backend_pdf import PdfPages
    pd = PdfPages('%s/%s.pdf' % (outputFolder, title))
    pd.savefig()
    pd.close()

    return

#---------------------------------------------------------------------------------
def compare_matrices(phenoMatrixFile, dResult, outputFolder):

    #NOTE: both files need to have the same order of antimicrobials

    #open file with phenotypic matrix and make a dictionary similar to that of the kmer results
    global dPhenoData
    dPhenoData = {} #dictionary linking sample to AMR in phenotypic data
    phenoAMList = []# list of antimicrobials in phenoMatrixFile
    dDifferences = {} #dictionary linking sample to differences in AMR between results matrix and phenotypic data matrix

    lines = iter(fileinput.input([phenoMatrixFile]))
    #get antimicrobial list
    for line in lines:
        if line.startswith('Sample'):
            line = line.rstrip('\n')
            ams = line.split("\t")
            for am in ams [0:]:
                #print am
                phenoAMList.append(am)

        #get sample names and make a dictionary with sample: AMR result
        for sample in sorted(dResult.keys()):
        #print sample
            #populate dPhenoData
            if line.startswith(sample):
                #print sample
                dPhenoData[sample] = []
                dDifferences[sample] = []

                line = line.rstrip('\n')
                details = line.split("\t")
                for detail in details [1:]:
                    dPhenoData[sample].append(detail)

    for sample in sorted(dResult.keys()):
        if sample in sorted(dPhenoData.keys()):
            #convert list of strings to list of integers
            result = map(int, dResult[sample])
            phenoData = map(int, dPhenoData[sample])
            #compare with resistance in dPhenoData[sample]
            difference = list(np.array(result)- np.array(phenoData))
            diffString = []
            for d in difference:
                diffString.append(str(d))
            #print sample, diffString
            dDifferences[sample] = diffString

    #print phenoAMList
    #print dPhenoData
    #print dDifferences
#title, xlabel
    plot_matrix(dPhenoData, phenoAMList, outputFolder, 'Greys',
                title = 'Phenotypic_matrix',
                xlabel = 'Antimicrobial resistance (in black)')

    plot_matrix(dDifferences, phenoAMList, outputFolder, 'bwr',
                title = 'Difference_matrix',
                xlabel = 'Differences in antimicrobial resistance\nRed = only predicted by software\nBlue = only phenotypic resistance\nWhite = no difference')

    return
#----------------------------------------------------------------------------------
def confusion_matrix(dActual, dPredicted, outputFolder):

    #NOTE: both dActual and dPredicted need to have the same order of antimicrobials

    #initialise variables for confusion matrix
    TP = 0
    TN = 0
    FP = 0
    FN = 0

    #loop through dictionaries with actual (phenotypic) and predicted (software) data
    for sample in sorted(dActual.keys()): #phenotypic data
        if sample in sorted(dPredicted.keys()): #predicted by software (e.g., kmer resistance)
            #compare data lists for each sample
            actualData = dActual[sample]
            predictedData = dPredicted[sample]
            for i in range (0, len(actualData)):
                #compare lists actualData and predictedData
                if actualData[i]=='1' and predictedData[i]=='1': #True positives (pheno + kmer = 1)
                    TP += 1
                if actualData[i]=='0' and predictedData[i]=='0': #True negatives (pheno + kmer = 0)
                    TN += 1
                if actualData[i]=='0' and predictedData[i]=='1': #False positives (pheno = 0; kmer = 1)
                    FP += 1
                if actualData[i]=='1' and predictedData[i]=='0': #False negatives (pheno = 1; kmer = 0)
                    FN += 1

    print 'TP, TN, FP, FN'
    print TP, TN, FP, FN

    #create confusion matrix
    mConfusion = [[TN, FP], [FN,TP]]
    #print mConfusion

    #print confusion matrix to pdf
    plt.figure()
    plot_confusion_matrix(outputFolder, mConfusion,
                            classes = ['Negative', 'Positive'],
                            title = 'Confusion_matrix',
                            cmap = plt.cm.Blues)

    #calculations from confusion matrix
    #initialise dictionary to store confusion matrix values and derived calculations; calc number:[symbol, value, description]
    dConfusion = {'TP':['TP', TP, 'True positive'],
                    'TN':['TN', TN, 'True negative'],
                    'FP':['FP', FP, 'False positive'],
                    'FN':['FN', FN, 'False negative']}

    #convert values to float before calculations
    fTP = float(TP)
    fTN = float(TN)
    fFP = float(FP)
    fFN = float(FN)

    #Specificity (SPC) or true negative rate
    SPC = round(fTN/(fTN + fFP), 3)
    dConfusion['calc1']=['SPC', SPC, 'Specificity (true negative rate)']

    #Sensitivity or true positive rate (TPR)
    TPR = round(fTP/(fTP + fFN), 3)
    dConfusion['calc2']=['TPR', TPR, 'Sensitivity (true positive rate)']

    #Precision or positive predictive value (PPV)
    PPV = round(fTP/(fTP + fFP), 3)
    dConfusion['calc3']=['PPV', PPV, 'Positive predictive value (precision)']

    #Negative predictive value (NPV)
    NPV = round(fTN/(fTN + fFN), 3)
    dConfusion['calc4']=['NPV', NPV,'Negative predictive value']

    #Accuracy (ACC)
    ACC = round((fTP + fTN)/(fTP + fFP + fTN + fFN), 3)
    dConfusion['calc5']=['ACC', ACC, 'Accuracy']

    #False positive rate (FPR)
    FPR = 1-SPC
    dConfusion['calc6']=['FPR', FPR, 'False positive rate']

    #False negative rate (FNR)
    FNR = 1-TPR
    dConfusion['calc7']=['FNR', FNR, 'False negative rate']

    #False discovery rate (FDR)
    FDR = 1-PPV
    dConfusion['calc8']=['FDR', FDR, 'False discovery rate']

    #Matthews correlation coeficient (MCC) - returns value between -1 and +1:
    #+1 perfect prediction; 0 no better than random prediction;-1 disagreement between prediction and observation
    import math
    MCC =round( (fTP*fTN - fFP*fFN)/math.sqrt((fTP+fFP)*(fTP+fFN)*(fTN+fFP)*(fTN+fFN)), 3)
    dConfusion['calc9']=['MCC', MCC,'Matthews Correlation Coefficient - values between -1 (worse) and +1 (best)']

    for key in sorted(dConfusion.keys()):
        if key not in ['TP', 'TN', 'FP', 'FN']:
            print dConfusion[key][1], dConfusion[key][2]

    #print calculations from confusion matrix to file
    confusionMatrix_to_file(dConfusion, outputFolder, fileName= 'Confusion_matrix_calculations.txt')

    return

#----------------------------------------------------------------------------------
def plot_confusion_matrix(outputFolder, cm, classes, title, cmap = plt.cm.Blues):

    plt.imshow(cm, interpolation = 'nearest', cmap = cmap)
    plt.title(title)
    #plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes) #, rotation = 45)
    plt.yticks(tick_marks, classes)

    for i, j in itertools.product(range(0,2), range(0,2)):#for i in range (0,2) for j in range (0,2)
        plt.text(j,i,cm[i][j],
        horizontalalignment = 'center',
        color = 'black',
        size = 16)

    plt.ylabel('Actual (phenotypic)')
    plt.xlabel('Predicted (software)')

    #print to outputFolder
    from matplotlib.backends.backend_pdf import PdfPages
    pd = PdfPages('%s/%s.pdf' % (outputFolder, title))
    pd.savefig()
    pd.close()

    return

#----------------------------------------------------------------------------------
def confusionMatrix_to_file(dictionary, outputFolder, fileName):

    outFile = os.path.join(outputFolder, fileName) # output file
    print 'Confusion Matrix file', outFile

    #open output file
    openFile = open(outFile, 'w')
    #write header line
    headers = ['Symbol', 'Value', 'Description']
    print >> openFile, '\t'.join(headers)
    #append a line for each symbol-value-description from confusion matrix
    for key in sorted(dictionary.keys()):
        print >> openFile, '\t'.join(str(x) for x in dictionary[key])

    openFile.close()

    return
#----------------------------------------------------------------------------------

#command:
#python AMR_parseKmer.py -i /home/laura/Desktop/PHE/AMR_data/Kmer_resistance_output -o /home/laura/Desktop/PHE/AMR -r /home/laura/Desktop/PHE/AMR_data
#------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
