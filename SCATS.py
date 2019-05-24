#!/usr/bin/python

from bin import my_functions as my
from bin import scats_functions as sc
import os, sys
fileAbsPath = os.path.abspath(os.path.dirname(__file__))
crtAbsPath = os.getcwd()

task = ""
taskList = ["refgene", "group", "count", "gene", "das", "sum"]
for i in range(1,len(sys.argv)):
    if sys.argv[i] == "-task" and len(sys.argv)!=i+1:
        task = sys.argv[i+1]
if (task not in taskList):
    print("\nPlease specify task (SCATS.py -task <task>):\n")
    print("\trefgene:   preprocess reference file\n")
    print("\tgroup:   group alternative splicing exon\n")
    print("\tcount:   count informative reads from indexed BAM file\n")
    #print("\tabkt:   calculate technical parameters (alpha beta kappa tau)")
    print("\tgene:    estimate mean gene expression for each single cell condition\n")
    print("\tdas:    detect differential alternative splicing (DAS) for each exon group between conditions\n")
    print("\tsum:    summarize DAS test results\n")

if task == "refgene":
    validArgList = ["-task", "-ref", "-out"]
    addAbsPath = [0, 1, 3]
    message = "SCATS.py -task refgene -ref <reference_file> -out <output_file>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    refFile = inputs[1]
    outFile = inputs[2]
    myCommand = "perl " + fileAbsPath + "/bin/PreProcess.pl -r " + refFile + " -o " + outFile
    os.system(myCommand)

if task == "group":
    validArgList = ["-task", "-refgene", "-out"]
    addAbsPath = [0, 1, 3]
    message = "SCATS.py -task group -refgene <refgene_file> -out <output_file>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    refFile = inputs[1]
    outFile = inputs[2]
    myCommand = "perl " + fileAbsPath + "/bin/getgroupinfo.pl " + refFile + " > " + outFile
    os.system(myCommand)

if task == "count":
    umiRun = ""
    onebam = ""
    for i in range(1,len(sys.argv)):
        if sys.argv[i] == "-umi" and len(sys.argv)!=i+1:
            umiRun = sys.argv[i+1]
        if sys.argv[i] == "-onebam" and len(sys.argv)!=i+1:
            onebam = sys.argv[i+1]
    if (umiRun not in ["yes", "no"]) or (onebam not in ["yes", "no"]):
        print("\nPlease specify umi and onebam option (SCATS.py -task count -umi <yes/no> -onebam <yes/no>):\n")
        print("\tumi:   collect UMI count or not (if yes umitag is required to be specified)\n")
        print("\tonebam:   whether all aligned reads are merged in one BAM file (if yes celltag and cellbc are required to be specified)\n")
        sys.exit()

    validArgList = ["-task", "-umi", "-onebam", "-meta", "-refgene", "-gpinfo"]
    addAbsPath = [0, 0, 0, 1, 1, 1]
    message = "SCATS.py -task count -umi yes -onebam -yes -meta <meta_file> -refgene <refgene_file> -gpinfo <gpinfo_file>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    metaFile = inputs[3]
    tmpDir = crtAbsPath + "/tmp"
    my.mk_dir(tmpDir)
    tmpDir = tmpDir + "/count_script"
    my.mk_dir(tmpDir)
    refgeneFile = inputs[4]
    gpinfoFile = inputs[5]

    # generate sh files for read counting process
    sc.write_count_sh(fileAbsPath, umiRun, onebam, metaFile, tmpDir, refgeneFile, gpinfoFile)
    print("\nPlease run all scripts (count_\*.sh files) under directory: " + tmpDir + "\n")

if task == "gene":
    
    validArgList = ["-task", "-ncore", "-meta"]
    addAbsPath = [0, 0, 1]
    message = "SCATS.py -task gene -ncore <# cores> -meta <metaFile>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    ncore = inputs[1]
    metaFile = inputs[2]
    tmpDir = crtAbsPath + "/tmp/gene_script"
    my.mk_dir(tmpDir)
    tmpDir = crtAbsPath + "/tmp/count_script"
    cdtList = sc.check_count_file(metaFile, tmpDir)

    outFile = crtAbsPath + "/tmp/celltypes"
    OUT = open(outFile, "w") # create celltype file
    for i in range(0, len(cdtList)):
        OUT.write(cdtList[i]+"\n")
    OUT.close()
    
    # estimate alpha
    my.mk_dir(crtAbsPath+"/tmp/abkt")
    myCommand = "perl " + fileAbsPath + "/bin/getalpha.pl " + metaFile + " " + tmpDir + " " + crtAbsPath+"/tmp/abkt/abkt_umi"
    os.system(myCommand)

    # estimate gene expression
    tmpDir = crtAbsPath + "/tmp/gene_script"
    my.mk_dir(tmpDir+"/data")
    
    outFile = crtAbsPath + "/tmp/comparegroup"
    OUT = open(outFile, "w") # create compare group file    
    for i in range(0,1):
        for j in range(i+1, len(cdtList)):
            OUT.write(cdtList[i]+"\t"+cdtList[j]+"\n")
            myCommand = "perl " + fileAbsPath + "/bin/getgenelevelcount.pl " + cdtList[i] + " " + cdtList[j]
            myCommand += " " + crtAbsPath + "/tmp/abkt/abkt_umi " + metaFile + " " + crtAbsPath + "/tmp/count_script " + crtAbsPath + "/tmp/gene_script/data"
            os.system(myCommand)
            myCommand = "perl " + fileAbsPath + "/bin/gettascdata.pl " + cdtList[i] + " " + cdtList[j]
            myCommand += " " + crtAbsPath + "/tmp/gene_script/data " + crtAbsPath + "/tmp/gene_script/data";
            os.system(myCommand)
            # generate sh files for gene expression estimation
            tmpDir = crtAbsPath + "/tmp/gene_script/data"
            mywrite = "mpirun -n " + ncore + " --bind-to none python " +  fileAbsPath + "/bin/model_selection_das_umi.py -y " + tmpDir + "/tascdata_" + cdtList[i] + "_" + cdtList[j]
            mywrite += " -k " + tmpDir + "/abktfile_" + cdtList[i] + "_" + cdtList[j] + " -x " + tmpDir + "/condition_" + cdtList[i] + "_" + cdtList[j]
            mywrite += " -t 4 -o " + tmpDir + "/outgene_" + cdtList[i] + "_" + cdtList[j] + "\n"
            myoutsh = crtAbsPath + "/tmp/gene_script" + "/gene_" + cdtList[i] + "_" + cdtList[j] + ".sh"
            os.system("echo \"" + mywrite +"\" > " + myoutsh)

    OUT.close()
    tmpDir = crtAbsPath + "/tmp/gene_script"
    print("\nPlease run all scripts (gene_\*.sh files) under directory: " + tmpDir + "\n")

if task == "das":
    ############### read count need to be filtered #### check getexonlevelcount.pl file to specify #######################
    validArgList = ["-task", "-ncore", "-meta", "-gpinfo"] 
    addAbsPath = [0, 0, 1, 1]
    message = "SCATS.py -task das -ncore <# cores> -meta <metaFile> -gpinfo <gpinfo_file>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    ncore = inputs[1]
    metaFile = inputs[2]
    gpinfoFile = inputs[3]
    tmpDir = crtAbsPath + "/tmp/das_script"
    my.mk_dir(tmpDir)
    my.mk_dir(tmpDir+"/data")
    cdtList = sc.check_count_file(metaFile, crtAbsPath + "/tmp/count_script")
    
    #collect gene expression and bursting rate
    myCommand = "perl " + fileAbsPath + "/bin/getgeneleveltheta_umi.pl " + crtAbsPath + "/tmp"
    os.system(myCommand)
    #collect informative read counts
    for i in range(0, len(cdtList)-1):
        for j in range(i+1, len(cdtList)):
            myCommand = "perl " + fileAbsPath + "/bin/getexonlevelcount_umi.pl " + cdtList[i] + " " + cdtList[j]
            myCommand += " " + crtAbsPath+"/tmp " + metaFile + " " + gpinfoFile
            os.system(myCommand)
            # generate sh files
            tmpDir = crtAbsPath + "/tmp/das_script/data"
            mywrite = "mpirun -n " + ncore + " --bind-to none python " +  fileAbsPath + "/bin/model_selection_das_umi.py -y " + tmpDir + "/countdata_" + cdtList[i] + "_" + cdtList[j]
            mywrite += " -k " + tmpDir + "/abktfile_" + cdtList[i] + "_" + cdtList[j] + " -x " + tmpDir + "/condition_" + cdtList[i] + "_" + cdtList[j]
            mywrite += " -t 6 -o " + tmpDir + "/out_" + cdtList[i] + "_" + cdtList[j] + "\n"
            myoutsh = crtAbsPath + "/tmp/das_script" + "/das_" + cdtList[i] + "_" + cdtList[j] + ".sh"
            os.system("echo \"" + mywrite +"\" > " + myoutsh)

    tmpDir = crtAbsPath + "/tmp/das_script"
    print("\nPlease run all scripts (das_\*.sh files) under directory: " + tmpDir + "\n")


if task == "sum":
    validArgList = ["-task", "-gpinfo"]
    addAbsPath = [0, 1]
    message = "SCATS.py -task das -gpinfo <gpinfo_file>"
    inputs = my.parse_argument(validArgList, addAbsPath, message)
    gpinfoFile = inputs[1]
    tmpDir = crtAbsPath + "/summary"
    my.mk_dir(tmpDir)
    outFile = tmpDir + "/DAS_results"
    tmpDir = crtAbsPath + "/tmp"
    compareFile = crtAbsPath + "/tmp/comparegroup"
    my.check_file(compareFile,"Please run SCATS.py -task gene.")
    
    myCommand = "perl " + fileAbsPath + "/bin/summarizedas.pl " + tmpDir + " " + gpinfoFile + " " + outFile
    os.system(myCommand)
