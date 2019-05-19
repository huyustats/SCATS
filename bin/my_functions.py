#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import math, sys, os, re, time

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)

# make a directory
def mk_dir(path):
    check = os.path.isdir(path)
    if not check:
        os.system("mkdir " + path)
    return

# parse arguments
def parse_argument(validArgList, addAbsPath, warnMessage):    
    for argIndex in range(1,len(sys.argv)):
        if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
            print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
            sys.exit()

    # assign arguments to a list
    outList = []
    for i in range(0, len(validArgList)):
        for argIndex in range(1,len(sys.argv)):
            if sys.argv[argIndex] == validArgList[i]:
                argIndex += 1
                if "~" in sys.argv[argIndex]:
                    sys.argv[argIndex] = os.path.expanduser(sys.argv[argIndex])
                fileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
                fileTmp = sys.argv[argIndex].split("/")
                if addAbsPath[i] == 1: # target file
                    fileTmp = fileAbsPath + "/" + fileTmp[len(fileTmp)-1]
                    check = os.path.exists(fileTmp)
                    if not check:
                        print(fileTmp+" does not exist!")
                        sys.exit()
                if addAbsPath[i] == 3: # create target file
                    fileTmp = fileAbsPath + "/" + fileTmp[len(fileTmp)-1]
                if addAbsPath[i] == 0: # value
                    fileTmp = fileTmp[len(fileTmp)-1]
                if addAbsPath[i] == 2: # target directory
                    fileTmp = os.path.abspath(sys.argv[argIndex])
                    check = os.path.isdir(fileTmp)
                    if not check:
                        print(fileTmp+" does not exist!")
                        sys.exit()
                
                outList.append(fileTmp)
    
    if len(outList) != len(validArgList):
        print(warnMessage)
        sys.exit()
    return outList


