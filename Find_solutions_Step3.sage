# run as sage ./Find_solutions_Step3.sage CMfieldnumber p
# e.g. sage ./Find_solutions_Step3.sage 1 runs Step 3 for CMfield 1 and prime p = 7
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 3: Separate solutions [x, dOVERn, a, cOVERn, b, gamma, n] from Step 2 by prime.
# This step does nothing when performed on CM field 4, 5, 55, 9 or 10. 
# This step affects only the CM fields with extra automorphisms which lead to the condition that p must divide n for the solution to be valid.
#-------------------------------------------------------------------------------------------------------------------------------------------
from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv

#---------------------------------------
CMfieldnumber = sage_eval(sys.argv[1])
p = sage_eval(sys.argv[2])
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
#5c
elif CMfieldnumber == 555:
    A = 210; B = 2541; C = 4375
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


print ("reading file with solutions from Step 2...")
sys.stdout.flush()
filename = "CMfield"+CMfieldnumber.str()+"Step2.csv"
reader=csv.DictReader(open(filelocation+filename))


#once the thing below is done write all potsolutions in csv files
filename2 = "CMfield"+CMfieldnumber.str()+"Step3"+"p"+p.str()+".csv"
csvfile = open(filelocation+filename2, 'w')
csvwriter = csv.writer(csvfile)
csvwriter.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n'))


@parallel
def separate_potential_solutions(CMfieldnumber,p): 
    #tm = time()
    OUTPUT = []
    for row in reader:
        solutions = []
        x = int(row["x"])
        dOVERn = int(row["dOVERn"])
        a = int(row["a"])
        cOVERn = int(row["cOVERn"])
        b = int(row["b"])
        gamma = int(row["gamma"])
        n = int(row["n"])
        if n > 1 and CMfieldnumber in [1,2,22,3,33,6,66,666,7,8,88]:
            if n/p in ZZ: #Need p | n due to argument at end of Section: Extra automorphisms
                print ('[x, dOVERn, a, cOVERn, b, gamma, n] = ', [x, dOVERn, a, cOVERn, b, gamma, n])
                csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n])
                csvfile.flush()
                OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n]])
                sys.stdout.flush()
        else: #if CMfieldnumber in [4,5,55,9,10] keep all solutions for each prime
            print ('[x, dOVERn, a, cOVERn, b, gamma, n] = ', [x, dOVERn, a, cOVERn, b, gamma, n])
            csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n])
            csvfile.flush()
            OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n]])
            sys.stdout.flush()
    print ("----------------------------------------------------------------------")
    print (len(OUTPUT)," solutions [x, dOVERn, a, cOVERn, b, gamma, n] found for the prime ", p)
    sys.stdout.flush()
    return OUTPUT


potsolutions = separate_potential_solutions(CMfieldnumber,p)

csvfile.close()

print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")
