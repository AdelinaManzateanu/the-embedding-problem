# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements the arguments of Section 3.4.1 and Section 5.3 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step3xJustBdp.sage CMfieldnumber
# e.g. sage ./Find_solutions_Step3xJustBdp.sage 34 runs Step 3 for CMfield 34
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 3: Find the bound for primes p up to which we may find solutions [x, dOVERn, a, cOVERn, b, gamma, n].
# This step works for fields in Section 6.3.
# The bound for primes p is as described in Section 5.3.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# A file "CMfield"+CMfieldnumber+"Step3xJustBdp.txt" is created which contains the bound for the prime p for each solution [x, dOVERn, a, cOVERn, b, gamma, n] coming from Step 2 and also the maximum of these bounds.
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv

######################################################
CMfieldnumber = sage_eval(sys.argv[1])
filelocation = "./"
T = time()

load('CM_Fields.sage')

    
print ("CM field number ",CMfieldnumber,": ","A = ",A,", B = ",B,", C = ",C)
#print version()
print ("----------------------------------------------------------------------")
sys.stdout.flush()


print ("reading file with solutions from Step 2...")
sys.stdout.flush()
filename = "CMfield"+CMfieldnumber.str()+"Step2.csv"
reader=csv.DictReader(open(filelocation+filename))


#write the results in a file
filename2 = "CMfield"+CMfieldnumber.str()+"Step3JustBdp.txt"
file2 = open(filelocation+filename2, 'w')


@parallel
def find_bound_for_p(CMfieldnumber): 
    #tm = time()
    all_bounds = []
    for row in reader:
        x = int(row["x"])
        dOVERn = int(row["dOVERn"])
        a = int(row["a"])
        cOVERn = int(row["cOVERn"])
        b = int(row["b"])
        gamma = int(row["gamma"])
        n = int(row["n"])
        bdp = 4*a*gamma*x^2
        print ("The bound for ", [x, dOVERn, a, cOVERn, b, gamma, n], "is ", bdp)
        file2.write(str("The bound for "))
        file2.write(str([x, dOVERn, a, cOVERn, b, gamma, n]))
        file2.write(str(" is "))
        file2.write(str(bdp))
        file2.write(str(".\n"))
        all_bounds.append(bdp)
        sys.stdout.flush()
    print ("----------------------------------------------------------------------")
    print (all_bounds)
    all_bounds = sorted(all_bounds)
    L = len(all_bounds)
    print ("There are ", L, "solutions.")
    file2.write(str("There are "))
    file2.write(str(L))
    file2.write(str(" solutions.\n"))
    print ("The largest bound for p is", all_bounds[L-1],".")
    file2.write(str("The largest bound for p is "))
    file2.write(str(all_bounds[L-1]))
    file2.write(str("."))
    sys.stdout.flush()
    return all_bounds


bound_for_p = find_bound_for_p(CMfieldnumber)

file2.close()


print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")