## Prerequisites:

#### Make sure you have the following libraries and packages installed on your system.
```
    cmake >= 2.6.4
    gcc >= 4.4
    Python 2.7
    python packages:
    	pysam
	numpy
	scipy
	cython
    OpenMPI
    SAMTOOLS
    Perl
```
## Installation

#### Download SCATS from github.
```
git clone https://github.com/huyustats/SCATS.git
```

#### Complie C functions using `Cython` and `gcc`
```
cd SCATS/bin/
bash complie_likelihoodumi.sh
```

#### Check programs and python packages installed or not
```
python check_software.py
```


