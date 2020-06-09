# run as sage ./Find_solutions_Step3.sage CMfieldnumber 
# e.g. sage ./Find_solutions_Step3.sage 1 runs Step 2 Part 2 for CM field 1.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 2 (Part 2): Find out which primes divide n for [x, dOVERn, a, cOVERn, b, gamma, n] that lift.
#-------------------------------------------------------------------------------------------------------------------------------------------
from sage.all_cmdline import * 
from time import time
from random import randint
import csv

#------------------------------------------
CMfieldnumber = sage_eval(sys.argv[1])
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


print ("reading file with solutions from Step 2 Part 1...")
sys.stdout.flush()
filename = "CMfield"+CMfieldnumber.str()+"Step2.csv"
reader=csv.DictReader(open(filelocation+filename))

filename2 = "CMfield"+CMfieldnumber.str()+"Step2_Primes.txt"
file2 = open(filelocation+filename2, 'w')

@parallel
def find_primes_dividing_n(CMfieldnumber): 
    #tm = time()
    primes = []
    for row in reader:
        n = int(row["n"])
        if n >1:
            F = list(factor(n))
            for f in F:
                if ZZ(f[0]) not in primes:
                    primes += [ZZ(f[0])]
            sys.stdout.flush()
    print ("----------------------------------------------------------------------")
    print (sorted(primes)," are the primes dividing n")
    file2.write(str(sorted(primes)))
    file2.write(" are the primes dividing n")
    sys.stdout.flush()
    return primes

potsolutions = find_primes_dividing_n(CMfieldnumber)

file2 .close()

print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")