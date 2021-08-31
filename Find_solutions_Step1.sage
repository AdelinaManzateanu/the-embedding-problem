# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This corresponds to Algorithm 2 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step1.sage CMfieldnumber 
# e.g. sage ./Find_solutions_Step1.sage 1 runs Step 1 for CM field 1.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 1: Find [x, dOVERn, a, cOVERn, b, gamma, n] for each field using bounds in Section 5.1.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# One file "CMfield"+CMfieldnumber+"Step1.csv" contains the list of [x, dOVERn, a, cOVERn, b, gamma, n] obtained
# One file "CMfield"+CMfieldnumber.str()+"Step1_Primes.txt" contains the list of primes dividing any n occuring in the previous file. This list includes 1 if there are solutions with n=1.
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


filename = "CMfield"+CMfieldnumber.str()+"Step1.csv"
csvfile = open(filelocation+filename, 'w')
csvwriter = csv.writer(csvfile)
csvwriter.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n'))

filename2 = "CMfield"+CMfieldnumber.str()+"Step1_Primes.txt"
file2 = open(filelocation+filename2, 'w')

@parallel
def find_all_integer_potential_solutions(CMfieldnumber):
    primes = []
    OUTPUT = []
    upx = floor(A^2-2*B)
    for x in [1..upx]:
        dOVERn = A-x
        upa = floor((A^2-2*B-x^2)/2)
        for a in [1..upa]:
            cOVERn = dOVERn*x - a - B
            b=cOVERn*x + dOVERn*a + C
            gamma=cOVERn*a+dOVERn*b
            if gamma > b^2/a: #includes gamma >0 and n = a*gamma - b^2 >0
                n = a*gamma - b^2
                if n >1:
                    F = list(factor(n))
                    for f in F:
                        if ZZ(f[0]) not in primes:
                            primes += [ZZ(f[0])]
                elif n == 1:
                    if 1 not in primes:
                        primes += [1]
                print ('[x, dOVERn, a, cOVERn, b, gamma, n] = ', [x, dOVERn, a, cOVERn, b, gamma, n])
                csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n])
                csvfile.flush()
                OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n]])
                sys.stdout.flush()
    print ("----------------------------------------------------------------------")
    print (len(OUTPUT)," solutions [x, dOVERn, a, cOVERn, b, gamma, n] found")
    print (sorted(primes)," are the primes dividing n. This list includes 1 if there are solutions with n=1.")
    file2.write(str(sorted(primes)))
    file2.write(" are the primes dividing n")
    sys.stdout.flush()
    return OUTPUT

potsolutions = find_all_integer_potential_solutions(CMfieldnumber)

csvfile.close()
file2.close()

print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")
