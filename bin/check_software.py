import my_functions as my

list = ['pysam', "numpy", "scipy", "cython"]
for x in list:
    my.check_module(x)


list = ["mpirun", "samtools"]
for x in list:
    my.check_program(x)
