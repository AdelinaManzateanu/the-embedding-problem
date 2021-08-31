# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This corresponds to Section 5.7.3 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step2.sage CMfieldnumber 
# e.g. sage ./Find_solutions_Step2.sage 1 runs Step 2 for CM field 1.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 2 
#1) Which [x, dOVERn, a, cOVERn, b, gamma, n] lift?
#2) Find out which primes divide n for [x, dOVERn, a, cOVERn, b, gamma, n] that lift.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# One file "CMfield"+CMfieldnumber+"Step2.csv" contains the list of [x, dOVERn, a, cOVERn, b, gamma, n] obtained
# One file "CMfield"+CMfieldnumber.str()+"Step2_Primes.txt" contains the list of primes dividing any n occuring in the previous file. This list includes 1 if there are solutions with n=1.
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


print ("reading file with solutions from Step 1...")
sys.stdout.flush()
filename = "CMfield"+CMfieldnumber.str()+"Step1.csv"
reader=csv.DictReader(open(filelocation+filename))

filename = "CMfield"+CMfieldnumber.str()+"Step2.csv"
csvfile = open(filelocation+filename, 'w')
csvwriter = csv.writer(csvfile)
csvwriter.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n'))

filename2 = "CMfield"+CMfieldnumber.str()+"Step2_Primes.txt"
file2 = open(filelocation+filename2, 'w')

@parallel
def find_primes_dividing_n_and_lift(CMfieldnumber): 
    #tm = time()
    primes = []
    OUTPUT = []
    for row in reader:
        x = int(row["x"])
        dOVERn = int(row["dOVERn"])
        a = int(row["a"])
        cOVERn = int(row["cOVERn"])
        b = int(row["b"])
        gamma = int(row["gamma"])
        n = int(row["n"])
        U = matrix(QQ,3,3,[x,a,b,1,0,cOVERn,0,1,dOVERn])
        I = matrix(QQ,3,3,[1,0,0,0,1,0,0,0,1])
        if CMfieldnumber == 1:
            N = U^2
        elif CMfieldnumber == 2:
            N = U^2
        elif CMfieldnumber == 22:
            N = U^2
        elif CMfieldnumber == 3:
            N = U^2
        elif CMfieldnumber == 33:
            N = U^2
        elif CMfieldnumber == 4:
            N = U^2
        elif CMfieldnumber == 5:
            N = 1/55*(U^2 + 45*U + 11*I) 
        elif CMfieldnumber == 55:
            N = 1/15*(U^2 + 6) 
        elif CMfieldnumber == 555:
            N = 1/4983*(U^2 + 2185*U + 3466*I) 
        elif CMfieldnumber == 5555:
            N = U^2
        elif CMfieldnumber == 6:
            N =  1/44*(U^2 + 7*U + 36*I) 
        elif CMfieldnumber == 66:
            N = 1/8*(U^2 - 5*U - 4*I) 
        elif CMfieldnumber == 666:
            N = 1/4*(U^2 - U) 
        elif CMfieldnumber == 6666:
            N = 1/4*(U^2 - U - 2)
        elif CMfieldnumber == 66666:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber == 7:
            N = 1/4*(U^2 - U) 
        elif CMfieldnumber == 77:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber == 8:
            N = 1/56*(U^2 + 7*U) 
        elif CMfieldnumber == 88:
            N = 1/12*(U^2 + 3*U)
        elif CMfieldnumber == 888:
            N = 1/124*(U^2 + 35*U + 100*I)
        elif CMfieldnumber == 8888:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber == 88888:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber == 9:
            N = 1/20*(U^2 + 16*I)
        elif CMfieldnumber == 99:
            N = 1/4*U^2
        elif CMfieldnumber == 999:
            N = U^2
        elif CMfieldnumber == 10:
            N = 1/245*(U^2 + 196*I)
        elif CMfieldnumber == 1010:
            N = 1/49*U^2
        elif CMfieldnumber == 1011:
            N = U^2
        elif CMfieldnumber ==  221:
            N = 1/49*U^2
        elif CMfieldnumber ==  2211:
            N = U^2
        elif CMfieldnumber ==  23:
            N = 1/63945*(U^2 + 2240*U + 53361)
        elif CMfieldnumber ==  233:
            N = 1/225*(U^2 + 28*U + 25*I)
        elif CMfieldnumber ==  2333:
            N = 1/3*(U^2 + 2*U)
        elif CMfieldnumber ==  25:
            N = 1/8379*(U^2 + 581*U + 7056*I)
        elif CMfieldnumber ==  255:
            N = 1/27*(U^2 - 5*U + 4)
        elif CMfieldnumber ==  2555:
            N = 1/297*(U^2 + 182*U + 225*I)
        elif CMfieldnumber ==  25555:
            N = 1/1341*(U^2 + 449*U + 315*I)
        elif CMfieldnumber ==  255555:
            N = 1/3*(U^2 + U)
        elif CMfieldnumber ==  26:
            N = 1/26901*(U^2 + 2177*U + 21609*I)
        elif CMfieldnumber ==  266:
            N = 1/207*(U^2 + 89*U + 18*I)
        elif CMfieldnumber ==  2666:
            N = 1/3*(U^2 + 2*U)
        elif CMfieldnumber ==  26666:
            N = 1/3*(U^2 + U)
        elif CMfieldnumber ==  27:
            N = 1/539*(U^2 + 14*U + 441*I)
        elif CMfieldnumber ==  277:
            N = 1/49*U^2
        elif CMfieldnumber ==  2777:
            N = U^2
        elif CMfieldnumber ==  28:
            N = 1/441*(U^2 + 35*U)
        elif CMfieldnumber ==  288:
            N = 1/3*(U^2 + 2*I)
        elif CMfieldnumber ==  2888:
            N = 1/3*(U^2 + U)
        elif CMfieldnumber ==  28888:
            N = 1/3*(U^2 + U)
        elif CMfieldnumber ==  29:
            N = 1/320*(U^2 + 256*I)
        elif CMfieldnumber ==  299:
            N = 1/20*(U^2 + 16*I)
        elif CMfieldnumber ==  2999:
            N = 1/4*U^2
        elif CMfieldnumber ==  29999:
            N = U^2
        elif CMfieldnumber ==  210:
            N = 1/64*U^2
        elif CMfieldnumber ==  2102:
            N = 1/4*U^2
        elif CMfieldnumber ==  2103:
            N = 1/4*U^2
        elif CMfieldnumber ==  2104:
            N = 1/52*(U^2 + 16*I)
        elif CMfieldnumber ==  2105:
            N = 1/116*(U^2 + 20*U + 96*I)
        elif CMfieldnumber ==  2106:
            N = 1/52*(U^2 + 40*I)
        elif CMfieldnumber ==  2107:
            N = 1/4*U^2
        elif CMfieldnumber ==  2108:
            N = 1/64*(U^2 - 4*U + 4*I)
        elif CMfieldnumber ==  2109:
            N = 1/364*(U^2 + 112*U + 168*I)
        elif CMfieldnumber ==  21010:
            N = 1/64*U^2
        elif CMfieldnumber ==  21011:
            N = U^2
        elif CMfieldnumber ==  211:
            N = 1/256*(U^2 - 8*U)
        elif CMfieldnumber ==  2112:
            N = 1/16*(U^2 - 2*U)
        elif CMfieldnumber ==  2113:
            N = 1/4*(U^2 - U - 2*I)
        elif CMfieldnumber ==  2114:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber ==  212:
            N = 1/121*U^2
        elif CMfieldnumber ==  2122:
            N = U^2
        elif CMfieldnumber ==  213:
            N = 1/5324*(U^2 + 77*U + 4356)
        elif CMfieldnumber ==  2133:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber ==  214:
            N = 1/1331*(U^2 - 99*U - 242)
        elif CMfieldnumber ==  2144:
            N = U^2
        elif CMfieldnumber ==  215:
            N = 1/7*(U^2 + 5*U + 4)
        elif CMfieldnumber ==  2155:
            N = U^2
        elif CMfieldnumber ==  216:
            N = 1/361*U^2
        elif CMfieldnumber ==  2166:
            N = U^2
        elif CMfieldnumber ==  217:
            N = 1/200716*(U^2 + 3173*U + 174724*I)
        elif CMfieldnumber ==  2177:
            N = 1/416*(U^2 + 177*U + 260*I)
        elif CMfieldnumber ==  21777:
            N = 1/192*(U^2 + 129*U + 132*I)
        elif CMfieldnumber ==  217777:
            N = 1/1024*(U^2 - 1003*U - 960*I)
        elif CMfieldnumber ==  2177777:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber ==  218:
            N = 1/532836*(U^2 + 21185*U + 12996*I)
        elif CMfieldnumber ==  2188:
            N = 1/6*(U^2 + 3*U + 2*I)
        elif CMfieldnumber ==  21888:
            N = 1/6*(U^2 + U)
        elif CMfieldnumber ==  218888:
            N = 1/4*(U^2 + U)
        elif CMfieldnumber ==  2188888:
            N = 1/4*(U^2 + U - 2*I)
        elif CMfieldnumber ==  219:
            N = 1/81356*(U^2 + 301*U + 66564*I)
        elif CMfieldnumber ==  2199:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber ==  220:
            N = 1/202005*(U^2 + 2345*U + 40401*I)
        elif CMfieldnumber ==  2200:
            N = 1/3*(U^2 + 2*U)
        elif CMfieldnumber ==  2201:
            N = 1/5*(U^2 + U)
        elif CMfieldnumber ==  31:
            N = 1/4*(U^2 - U - 2*I)
        elif CMfieldnumber ==  311:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber ==  32:
            N = U^2
        elif CMfieldnumber ==  34:
            N = 1/3*(U^2 + U)
        elif CMfieldnumber ==  35:
            N = U^2
        elif CMfieldnumber ==  36:
            N = 1/3*(U^2 + 2*I)
        elif CMfieldnumber ==  37:
            N = U^2
        elif CMfieldnumber ==  38:
            N = 1/8*(U^2 - 3*U - 4*I)
        elif CMfieldnumber ==  388:
            N = 1/4*(U^2 + U*I)
        elif CMfieldnumber ==  39:
            N = U^2
        elif CMfieldnumber ==  40:
            N = U^2
        elif CMfieldnumber ==  41:
            N = U^2
        elif CMfieldnumber ==  42:
            N = 1/4*(U^2 - U - 2*I)
        elif CMfieldnumber ==  422:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber ==  43:
            N = U^2
        elif CMfieldnumber ==  45:
            N = U^2
        elif CMfieldnumber ==  46:
            N = U^2
        elif CMfieldnumber ==  47:
            N = 1/2*(U^2 + U)
        elif CMfieldnumber ==  48:
            N = U^2
        elif CMfieldnumber ==  488:
            N = U^2
        elif CMfieldnumber ==  49:
            N = 1/2*(U^2 + I)
        elif CMfieldnumber ==  50:
            N = U^2
        elif CMfieldnumber ==  51:
            N = U^2
        elif CMfieldnumber ==  52:
            N = U^2
        elif CMfieldnumber ==  53:
            N = U^2
        elif CMfieldnumber ==  54:
            N = U^2
        elif CMfieldnumber ==  56:
            N = U^2
        elif CMfieldnumber ==  57:
            N = 1/2*(U^2 + U)
        if (N[0,0] in ZZ) and (n*N[1,0] in ZZ) and (n*N[2,0] in ZZ) and (N[0,1] in ZZ) and (N[0,2] in ZZ) and (n*N[1,1] in ZZ)  and (n*N[1,2] in ZZ) and (n*N[2,1] in ZZ) and (n*N[2,2] in ZZ):
            if n > 1:
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
    return primes

potsolutions = find_primes_dividing_n_and_lift(CMfieldnumber)

csvfile.close()
file2.close()

print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")

#-----------------------------------------------------------------------------------------------------------------------
#Explanation: How to find N?
#Magma
#A:=42; B:=441; C:=847; // modify this accordingly
#R<x>:=PolynomialRing(Rationals());
#f:=x^3-A*x^2+B*x-C; 
#K<U>:=NumberField(f);
#OK := MaximalOrder(K); O:= EquationOrder(K);
#//Index(OK,O);
#//T:=TransformationMatrix(OK,O); T;
#Basis(OK, NumberField(OK)); // this gives N
