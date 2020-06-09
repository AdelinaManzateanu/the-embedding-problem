# run as sage ./Find_solutions_Step1.sage CMfieldnumber 
# e.g. sage ./Find_solutions_Step1.sage 1 runs Step 1 for CM field 1.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 1: Find [x, dOVERn, a, cOVERn, b, gamma, n] for each field using bounds.
#-------------------------------------------------------------------------------------------------------------------------------------------
from sage.all_cmdline import *  
from time import time
from random import randint
import csv

#---------------------------------------
CMfieldnumber = sage_eval(sys.argv[1])
#p = sage_eval(sys.argv[2])
filelocation = "./"
T = time()

if CMfieldnumber == 1:
    A = 13; B = 50; C = 49
#2a
elif CMfieldnumber == 2:
    A = 6; B = 9; C = 1
#2b
elif CMfieldnumber == 22:
    A = 9; B = 6; C = 1
#3a
elif CMfieldnumber == 3:
    A = 5; B = 6; C = 1 
#3b
elif CMfieldnumber == 33:
    A = 6; B = 5; C = 1 
elif CMfieldnumber == 4: 
    A = 7; B = 14; C = 7
#5a
elif CMfieldnumber == 5:
    A = 42; B = 441; C = 847
#5b
elif CMfieldnumber == 55:
    A = 63; B = 126; C = 63
#6a
elif CMfieldnumber == 6:
    A = 29; B = 180; C = 64
#6b
elif CMfieldnumber == 66:
    A = 30; B = 257; C = 484
#6c
elif CMfieldnumber == 666:
    A = 37; B = 356; C = 1024
elif CMfieldnumber == 7:
    A = 21; B = 116; C = 64
#8a
elif CMfieldnumber == 8:
    A = 42; B = 441; C = 784
#8b
elif CMfieldnumber == 88:
    A = 45; B = 612; C = 2304
#X9
elif CMfieldnumber == 9:
    A = 18; B = 56; C = 8 
#CM plane quartic curve X1 in https://arxiv.org/pdf/1701.06489.pdf (p 21)
elif CMfieldnumber == 10:
    A = 63; B = 686; C = 343 
else:
    print ("unknown CM field number")


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
    #tm = time()
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
                    print ('[x, dOVERn, a, cOVERn, b, gamma, n] = ', [x, dOVERn, a, cOVERn, b, gamma, n])
                    csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n])
                    csvfile.flush()
                    OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n]])
                    sys.stdout.flush()
    print ("----------------------------------------------------------------------")
    print (len(OUTPUT)," solutions [x, dOVERn, a, cOVERn, b, gamma, n] found")
    print (sorted(primes)," are the primes dividing n")
    file2.write(str(sorted(primes)))
    file2.write(" are the primes dividing n")
    sys.stdout.flush()
    return OUTPUT

potsolutions = find_all_integer_potential_solutions(CMfieldnumber)

csvfile.close()
file2 .close()

print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")
