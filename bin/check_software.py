import my_functions as my

list = ['pysam', "numpy", "scipy", "cython"]
#for x in list:
    #my.check_module(x)

try:
    import pysam
    print("Module 'pysam' is installed.")
except:
    print("Module 'pysam' is NOT installed!")

try:
    import numpy
    print("Module 'numpy' is installed.")
except:
    print("Module 'numpy' is NOT installed!")

try:
    import scipy
    print("Module 'scipy' is installed.")
except:
    print("Module 'scipy' is NOT installed!")

try:
    import cython
    print("Module 'cython' is installed.")
except:
    print("Module 'cython' is NOT installed!")



list = ["mpirun", "samtools"]
for x in list:
    my.check_program(x)
