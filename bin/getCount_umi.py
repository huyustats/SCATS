#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
import math, sys, os, re, pysam, time
import my_functions as my

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)


###############################################################################
###  ARGUMENT SETTINGS
###############################################################################

# checking whether argument is valid or not
validArgList = ["-bam", "-ref", "-out", "-gpinfo", "-umitag"]
addAbsPath = [1,1,3,1,0]
warnMessage = "-bam, -ref, -out, -gpinfo, -umitag"        
inputFile = my.parse_argument(validArgList, addAbsPath, warnMessage)
bamFile = inputFile[0]
refGeneFile = inputFile[1]
outFile = inputFile[2]
gpinfoFile = inputFile[3]
umitag = inputFile[4]


# load gene information
geneStructureInformation = auto_dict()
geneLineCount = auto_dict()

with open(refGeneFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        tmpinf = line.split("\t")
        gene = tmpinf[0]
        
        if not bool(geneStructureInformation[gene]):
            geneLineCount[gene] = 0
            geneStructureInformation[gene][geneLineCount[gene]] = line
        else:
            geneLineCount[gene] += 1
            geneStructureInformation[gene][geneLineCount[gene]] = line

# load group information

groupInformation = auto_dict()
geneLineCount1 = auto_dict()
with open(gpinfoFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        tmpinf = line.split("\t")
        tmpinf[5] = tmpinf[5].strip(",")
        gene = tmpinf[0]
        
        groupInformation[gene][tmpinf[1]][tmpinf[3]] = tmpinf[5]




#####################################
## Using pysam to read in bam file !!
#####################################
bamFilePysam = pysam.Samfile(bamFile,"rb")


## RESULTS FILE
OUT = open(outFile, 'w')


###########################################################################################################################
###  START TO ANALYZE DATA FOR EACH GENE ###
##########################################################################################################################

geneCount = 0

startTime = time.time()

umiSet = auto_dict()

#OUT.write("GeneName\tIsoformName\tNumberOfReads\tRelativeAbundance\n") ## Header of Results

for gene in geneStructureInformation:

    countResults = auto_dict()

    geneCount += 1
    tmpTime = (time.time() - startTime)/60.0
    
    
    sameReadCount = auto_dict()
    readStart = auto_dict()
    readEnd = auto_dict()
    readCigar = auto_dict()

    numofExons = geneLineCount[gene]
    tmpgeneinf = geneStructureInformation[gene][0].split("\t")
    geneChr = tmpgeneinf[1]
    geneStart = int(tmpgeneinf[3])
    geneEnd = int(tmpgeneinf[4])
    readCount = 0
    if bamFilePysam.get_tid(geneChr) == -1:
        continue

    ## load all reads information which were mapped to the specific gene within this loop using pysam
    for read in bamFilePysam.fetch(geneChr, geneStart, geneEnd):
        line = str(read)
        tmpinf = line.split("\t")
        tmpReadName = tmpinf[0]
        try:
            tmpUMI = read.get_tag(umitag)
        except:
            continue
        
    
        tmpReadChr = geneChr
        tmpReadStart = int(tmpinf[3]) + 1
        tmpReadCigar = ""
        
        ## Adjust to different Pysam Version!! ##

        if ")]" in tmpinf[5]: ## vector format
            
            tmpinf[5] = tmpinf[5].rstrip(")]")
            tmpinf[5] = tmpinf[5].lstrip("[(")
            tmpinfcigar = tmpinf[5].split("), (")
            for cc in tmpinfcigar:
                ttcc = cc.split(", ")
                if ttcc[0] == "3":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "N"
                if ttcc[0] == "2":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "D"
                if ttcc[0] == "1":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "I"
                if ttcc[0] == "0":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "M"
                if not (ttcc[0] == "3" or ttcc[0] == "2" or ttcc[0] == "1" or ttcc[0] == "0"):
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "X"
        else:      ## 100M10N100M format
            tmpReadCigar = tmpinf[5]
        #print(tmpReadCigar)                            
        if not bool(sameReadCount[tmpReadName]):
            sameReadCount[tmpReadName] = 1
            umiSet[tmpReadName] = tmpUMI
        else:
            sameReadCount[tmpReadName] += 1
                                        
        readStart[tmpReadName][sameReadCount[tmpReadName]] = tmpReadStart
        readCigar[tmpReadName][sameReadCount[tmpReadName]] = tmpReadCigar


    ## load structure information of the specific gene within this loop                    
                        
    tmpgeneinf[5] = tmpgeneinf[5].rstrip(",")
    isoformNames = tmpgeneinf[5].split(",")
    exonStarts = [None] * numofExons
    exonEnds = [None] * numofExons
    exonIndicators = auto_dict()
    
    for i in range(1,numofExons+1):
        tmpinf = geneStructureInformation[gene][i].split("\t")
        exonStarts[i-1] = int(tmpinf[3])+1
        exonEnds[i-1] = int(tmpinf[4])
        tmpinf[5] = tmpinf[5].rstrip(",")
        tmpExonIndicators = tmpinf[5].split(",")

        for j in range(len(tmpExonIndicators)):
            exonIndicators[isoformNames[j]][i-1] = int(tmpExonIndicators[j])

    lociIndicators = auto_dict()
    for i in range(len(isoformNames)):
        for j in range(len(exonStarts)):
            if exonIndicators[isoformNames[i]][j] == 1:
                for k in range(exonStarts[j], exonEnds[j]+1):
                    lociIndicators[isoformNames[i]][k] = 1

    #########################################################################################################################################
    ## START TO ANALYZE EACH READ 
    ##################################################################################################################################################

    qualifiedRead = auto_dict()
    readSet = []
    fragmentStart = auto_dict()
    fragmentEnd = auto_dict()
    CompatibleMatrix = auto_dict()
    tmpCompatibleMatrix = auto_dict()
    
    for readName in sameReadCount:

        # load CIGAR information
        cigarNumberRead1 = auto_dict()
        cigarNumberRead2 = auto_dict()
        cigarMatchRead1 = auto_dict()
        cigarMatchRead2 = auto_dict()
        cigarInfCountRead1 = 0
        cigarInfCountRead2 = 0
        cigarInfCountRead1tmp = 0
        cigarInfCountRead2tmp = 0
        
        tmp1 = re.split("([A-Z])",readCigar[readName][1])
        for i in range(len(tmp1)-1):
            if tmp1[i].isalpha():
                cigarMatchRead1[cigarInfCountRead1] = tmp1[i]
                cigarInfCountRead1 += 1
            else:
                cigarNumberRead1[cigarInfCountRead1] = int(tmp1[i])
                cigarInfCountRead1tmp += 1
                
        if sameReadCount[readName] == 2:
            tmp2 = re.split("([A-Z])",readCigar[readName][2])
            for i in range(len(tmp2)-1):
                if tmp2[i].isalpha():
                    cigarMatchRead2[cigarInfCountRead2] = tmp2[i]
                    cigarInfCountRead2 += 1
                else:
                    cigarNumberRead2[cigarInfCountRead2] = int(tmp2[i])
                    cigarInfCountRead2tmp += 1
                    
        # calculate read end positions
        readEnd[readName][1] = readStart[readName][1]
        for i in range(cigarInfCountRead1):
            readEnd[readName][1] += cigarNumberRead1[i]
            
        if sameReadCount[readName] == 2:
            readEnd[readName][2] = readStart[readName][2]
            for i in range(cigarInfCountRead2):
                readEnd[readName][2] += cigarNumberRead2[i]

        # calculate fragment START and END positions
        if sameReadCount[readName] == 2:
            fragmentStart[readName] = readStart[readName][2] if readStart[readName][1] >= readStart[readName][2] else readStart[readName][1]
            fragmentEnd[readName] = readEnd[readName][1] if readEnd[readName][1] >= readEnd[readName][2] else readEnd[readName][2]

        if sameReadCount[readName] == 1:
            fragmentStart[readName] = readStart[readName][1]
            fragmentEnd[readName] = readEnd[readName][1]
            
        ##################################################################################################################################    
        ## Obtain compatible matrix of isoforms with respect to reads
        #################################################################################################################################
        if (readStart[readName][1] >= geneStart and readStart[readName][1] <= geneEnd):
        #if (readStart[readName][1] >= geneStart and readStart[readName][1] <= geneEnd) or (readStart[readName][2] >= geneStart and readStart[readName][2] <= geneEnd and sameReadCount[readName]==2) :
            if cigarInfCountRead1 == cigarInfCountRead1tmp and cigarInfCountRead2 == cigarInfCountRead2tmp:
                base1 = readStart[readName][1] - 1
                exonIndicatorRead1 = [0] * numofExons
                if sameReadCount[readName] == 2:
                    base2 = readStart[readName][2] - 1
                    exonIndicatorRead2 = [0] * numofExons
                compatibleVector = [1] * len(isoformNames)

                ##############################################################################################################################################
                ### SET TUP COMPATIBLE INDICATOR VECTOR ###############
                ###############################################################################################################################################
                ## READ 1 ##
                # find exons where read 1 mapped to
                for i in range(cigarInfCountRead1):
                    
                    if cigarMatchRead1[i] == "M" or cigarMatchRead1[i] == "I": ## matched CIGAR

                        for j in range(1,cigarNumberRead1[i]+1):
                            tmpbase = base1 + j
                            for k in range(len(exonStarts)):
                                if exonIndicatorRead1[k] == 1: continue
                                if tmpbase >= exonStarts[k] and tmpbase <= exonEnds[k]: exonIndicatorRead1[k] = 1 ## confirm that the read covers this exon
        
                        base1 += cigarNumberRead1[i] # jump to next match information

                    if cigarMatchRead1[i] == "N": ## skipping area
                        base1 += cigarNumberRead1[i] # jump to next match information directly

                # set up indicator vector
                tmpcount1 = 0
                tmpcount11 = 0 ## these two variable are used to rule out skipping exons
                for i in range(len(exonIndicatorRead1)):
                    if exonIndicatorRead1[i] == 1: tmpcount1 += 1
                for i in range(len(exonIndicatorRead1)):

                    if exonIndicatorRead1[i] == 1:
                        tmpcount11 += 1
                        for j in range(len(isoformNames)):
                            if exonIndicators[isoformNames[j]][i] == 0: compatibleVector[j] = 0 ## rule out isoform j if reads covers skipping area of isoform j

                    if exonIndicatorRead1[i] == 0: #aim to rule out isforms which includes exons which skipped by read
                        if tmpcount1 > 1 and tmpcount11 >= 1 and tmpcount11 < tmpcount1: ## confirm the exon i is skipped by read!!
                            for j in range(len(isoformNames)):
                                if exonIndicators[isoformNames[j]][i] == 1: compatibleVector[j] = 0

                    
                ## READ 2 ## SAME AS READ 1
                tmpcount2 = 0
                if sameReadCount[readName] == 2: ## ONLY WHEN THE READ IS PAIRED-END READ!!!
                    # find exons where read 2 mapped to
                    for i in range(cigarInfCountRead2):
                        
                        if cigarMatchRead2[i] == "M" or cigarMatchRead2[i] == "I": ## matched CIGAR
                            
                            for j in range(1,cigarNumberRead2[i]+1):
                                tmpbase = base2 + j
                                for k in range(len(exonStarts)):
                                    if exonIndicatorRead2[k] == 1: continue
                                    if tmpbase >= exonStarts[k] and tmpbase <= exonEnds[k]: exonIndicatorRead2[k] = 1 ## confirm that the read covers this exon
                                    
                            base2 += cigarNumberRead2[i] # jump to next match information
                                    
                        if cigarMatchRead2[i] == "N": ## skipping area
                            base2 += cigarNumberRead2[i] # jump to next match information directly
                                        
                    # set up indicator vector
                    tmpcount2 = 0
                    tmpcount22 = 0 ## these two variable are used to rule out skipping exons
                    for i in range(len(exonIndicatorRead2)):
                        if exonIndicatorRead2[i] == 1: tmpcount2 += 1
                    for i in range(len(exonIndicatorRead2)):
                            
                        if exonIndicatorRead2[i] == 1:
                            tmpcount22 += 1
                            for j in range(len(isoformNames)):
                                if exonIndicators[isoformNames[j]][i] == 0: compatibleVector[j] = 0 ## rule out isoform j if reads covers skipping area of isoform j
                                                    
                        if exonIndicatorRead2[i] == 0: #aim to rule out isforms which includes exons which skipped by read
                            if tmpcount2 > 1 and tmpcount22 >= 1 and tmpcount22 < tmpcount2: ## confirm the exon i is skipped by read!!
                                for j in range(len(isoformNames)):
                                    if exonIndicators[isoformNames[j]][i] == 1: compatibleVector[j] = 0
                                                                
                ##################################################################################################################################################
                ## fill in compatible matrix ##
                if tmpcount1 > 0 or (tmpcount2 > 0 and sameReadCount[readName] == 2):
                    #umibarcode = readName.split("_")
                    #umibarcode = umibarcode[len(umibarcode)-1]
                    umibarcode = umiSet[readName]
                    readSet.append(umibarcode)
                    qualifiedRead[readName] = 1
                    readCount += 1
                    for i in range(len(isoformNames)):
                        CompatibleMatrix[readName][isoformNames[i]] = compatibleVector[i]
                        tmpCompatibleMatrix[readName][isoformNames[i]] = compatibleVector[i]



    ### COMPATIBLE MATRIX OBTAINED !!!
    ###############################################################################################################

    #readCount = len(set(readSet))
    if readCount == 0: continue
    print(gene+"\t"+str(readCount)+" reads detected...")
    #print(umibarcode)
    for weight in groupInformation[gene]:
        countResults[weight]["+"] = []
        countResults[weight]["-"] = []
        
        isosetplus = groupInformation[gene][weight]["+"].split(",")
        isosetminus = groupInformation[gene][weight]["-"].split(",")
        
        for readName in qualifiedRead:
            umibarcode = readName.split("_")                                                                                                                                                
            umibarcode = umibarcode[len(umibarcode)-1]  
            
            if qualifiedRead[readName] == 0: continue
            #print(umibarcode)
            sumindexplus = 0
            for index in isosetplus:
                if CompatibleMatrix[readName][isoformNames[int(index)]] == 1: sumindexplus += 1
            sumindexminus = 0
            for index in isosetminus:
                if CompatibleMatrix[readName][isoformNames[int(index)]] == 1: sumindexminus += 1
            if sumindexplus == 0:
                countResults[weight]["+"].append(umibarcode)
            if sumindexminus == 0:
                countResults[weight]["-"].append(umibarcode)
        count_plus = len(set(countResults[weight]["+"]))
        count_minus = len(set(countResults[weight]["-"]))
        OUT.write(gene+"\t"+str(readCount)+"\t"+weight+"\t"+"+"+"\t"+str(count_plus)+"\n")
        OUT.write(gene+"\t"+str(readCount)+"\t"+weight+"\t"+"-"+"\t"+str(count_minus)+"\n")

OUT.close()
            
                            


                    
