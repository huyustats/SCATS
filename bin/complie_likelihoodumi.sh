#!/bin/bash
rm ./likelihoodumi.c
rm ./likelihoodumi.html
rm ./likelihoodumi.so
cython -a ./likelihoodumi.pyx
gcc -shared -pthread -fPIC `python-config --cflags` -o likelihoodumi.so likelihoodumi.c
