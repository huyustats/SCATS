#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import math, sys, os, re, time

# check whether meta file is qualified
def check_meta(metaFile, umiRun, onebam):
    cdtset = []
    with open (metaFile, "r") as FP:
        for line in FP:
            line = line.strip("\n")
            tmpinf = line.split("\t")
            cellbc = tmpinf[0]
            condition = tmpinf[1]
            bamfile = tmpinf[2]
            
            cdtset.append(condition)
            exists = os.path.isfile(bamfile)
            if not exists:
                print(bamfile+" does not exist!")
                sys.exit()
            exists = os.path.isfile(bamfile+".bai")
            if not exists:
                print(bamfile+".bai does not exist! Please index BAM file.")
                sys.exit()

            if umiRun == "yes":
                umitag = tmpinf[3]
                if not umitag:
                    print("Please specify UMI tag name for each cell at 4th column of meta file!")
                    sys.exit()
                if onebam == "yes":
                    celltag = tmpinf[4]
                    if not celltag:
                        print("Please specify cell tag name for each cell at 5th column of meta file!")
                        sys.exit()
            if umiRun == "no":
                if onebam == "yes":
                    celltag = tmpinf[3]
                    if not celltag:
                        print("Please specify cell tag name for each cell at 4th column of meta file!")
                        sys.exit()

    # check number of conditions
    #cdtset = len(set(cdtset))
    #if cdtset < 2:
        #print("Please specify 2 conditions at 2nd column of meta file!")
        #sys.exit()
    #if cdtset > 2:
        #print("Please specify only 2 conditions at 2nd column of meta file!")
        #sys.exit()
    return

# check count file
def check_count_file(metaFile, tmpDir):
    check_meta(metaFile, "no", "no")
    cdtset = []
    with open (metaFile, "r") as FP:
        for line in FP:
            line = line.strip("\n")
            tmpinf = line.split("\t")
            countFile = tmpDir + "/count_" + tmpinf[0] + ".out"
            check = os.path.exists(countFile)
            if not check:
                print(countFile+" does not exist! Please run SCATS.py -task count to obtain read count files.")
                sys.exit()
            condition = tmpinf[1]
            cdtset.append(condition)
    
    cdtset = list(set(cdtset))
    cdtset.sort()
    return cdtset


## write count sh file to tmp directory
def write_count_sh(fileAbsPath, umiRun, onebam, metaFile, tmpDir, refgeneFile, gpinfoFile):
    check_meta(metaFile, umiRun, onebam)
    if umiRun == "yes" and onebam == "yes":
        with open (metaFile, "r") as FP:
            for line in FP:
                line = line.strip("\n")
                tmpinf = line.split("\t")
                cellbc = tmpinf[0]
                bamfile = tmpinf[2]
                umitag = tmpinf[3]
                celltag = tmpinf[4]
                outFile = tmpDir + "/count_" + cellbc + ".sh"
                OUT = open(outFile, "w")
                outwrite = "python " + fileAbsPath + "/bin/getCount_umi_cellid.py -bam " + bamfile + " -ref " + refgeneFile + " -gpinfo " + gpinfoFile + " -out " + tmpDir + "/count_" + cellbc + ".out"
                outwrite += " -cellid " + cellbc + " -celltag " + celltag + " -umitag " + umitag + "\n"
                OUT.write(outwrite)
                OUT.close()
    if umiRun == "yes" and onebam == "no":
        with open (metaFile, "r") as FP:
            for line in FP:
                line = line.strip("\n")
                tmpinf = line.split("\t")
                cellbc = tmpinf[0]
                bamfile = tmpinf[2]
                umitag = tmpinf[3]
                #celltag = tmpinf[4]
                outFile = tmpDir + "/count_" + cellbc + ".sh"
                OUT = open(outFile, "w")
                outwrite = "python " + fileAbsPath + "/bin/getCount_umi.py -bam " + bamfile + " -ref " + refgeneFile + " -gpinfo " + gpinfoFile + " -\
out " + tmpDir + "/count_" + cellbc + ".out"
                outwrite += " -umitag " + umitag + "\n"
                OUT.write(outwrite)
                OUT.close()
    if umiRun == "no" and onebam == "yes":
        with open (metaFile, "r") as FP:
            for line in FP:
                line = line.strip("\n")
                tmpinf = line.split("\t")
                cellbc = tmpinf[0]
                bamfile = tmpinf[2]
                #umitag = tmpinf[3]
                celltag = tmpinf[3]
                outFile = tmpDir + "/count_" + cellbc + ".sh"
                OUT = open(outFile, "w")
                outwrite = "python " + fileAbsPath + "/bin/getCount_cellid.py -bam " + bamfile + " -ref " + refgeneFile + " -gpinfo " + gpinfoFile + " -\
out " + tmpDir + "/count_" + cellbc + ".out"
                outwrite += " -cellid " + cellbc + " -celltag " + celltag  + "\n"
                OUT.write(outwrite)
                OUT.close()
    if umiRun == "no" and onebam == "no":
        with open (metaFile, "r") as FP:
            for line in FP:
                line = line.strip("\n")
                tmpinf = line.split("\t")
                cellbc = tmpinf[0]
                bamfile = tmpinf[2]
                outFile = tmpDir + "/count_" + cellbc + ".sh"
                OUT = open(outFile, "w")
                outwrite = "python " + fileAbsPath + "/bin/getCount.py -bam " + bamfile + " -ref " + refgeneFile + " -gpinfo " + gpinfoFile + " -\
out " + tmpDir + "/count_" + cellbc + ".out"
                OUT.write(outwrite)
                OUT.close()

    return
