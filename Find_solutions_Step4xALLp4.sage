# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements Algorithm 3 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step4xALLp4.sage CMfieldnumber minbdp maxbdp 
# e.g. sage ./Find_solutions_Step4xALLp4.sage 34 1 1000 runs Step 4 for CMfield 34 for primes 1 < p < 1000
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 4: For each prime p between min_bdp and  max_bdp, read all solutions [x, dOVERn, a, cOVERn, b, gamma, n] from the corresponding file "CMfield"+CMfieldnumber+"Step3"+"p"+p.csv" created in Step 3. 
# For each [x, dOVERn, a, cOVERn, b, gamma, n] from Step 3, find solutions [x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3] that satisfy the following:
# 1) the 3 equations in Proposition 3.5 and the implied bounds in Lemma 5.5.
# 2) Nx1, Nx2, Nx3, N(x1+x2), N(x1+x3), N(x2+x3), N(x1+x2+x3) are achievable norms of elements with trace 0 in the maximal order (i.e. satisfy certain congruence conditions mod p as in Section 5.7.4)
# 3) |Tr(xi*xj)| \leq floor(2*sqrt(Nxi*Nxj)) for all i,j=1,2,3 (see Lemma Lemma 5.3)
# 4) Dxixj = 4*Nxi*Nxj - Trxixj^2 are integers divisible by p for all i,j=1,2,3 (see Lemma 5.14 in Section 5.6)
# 5) Dx1x2x3 = 4*Nx1*Nx2*Nx3 - Nx1*Trx2x3^2 - Nx2*Trx1x2^2 - Nx3*Trx1x3^2 - Trx1x2*Trx1x3*Trx2x3 is an integer divisible by p and is a square (see Lemma 5.14 in Section 5.6)
# 6) if Nxi = 0 then Trxixj = Trxjxi = 0 for all i,j=1,2,3.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Brief explanation for 2) above:
# Know: x1,x2,x3 in max order O in quaternion algebra QA with Tr(x1)=Tr(x2)=Tr(x3)=0, N(x1)=Nx1, N(x2)=Nx2, N(x3)=Nx3, Tr(x1*x2)=Trx1x2, Tr(x1*x3)=Trx1x3, Tr(x2*x3)=Trx2x3.
# Know: Tr(x1*x2)= N(x1)+N(x2)-N(x1+x2) and Tr(x1*x2) + Tr(x1*x3) + Tr(x2*x3) = N(x1) + N(x2) + N(x3) - N(x1+x2+x3)
# Know: Tr(x1*x2*x3) = Tr(x2*x3*x1) = Tr(x3*x1*x2) = - Tr(x1*x3*x2) = - Tr(x3*x2*x1) = - Tr(x2*x1*x3).
# Check: N(x1)+N(x2)-Tr(x1*x2), N(x1)+N(x3)-Tr(x1*x3), N(x2)+N(x3)-Tr(x2*x3), N(x1) + N(x2) + N(x3) - Tr(x1*x2) - Tr(x1*x3) - Tr(x2*x3) in GoodNorms.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# For each prime p between min_bdp and  max_bdp, one file "CMfield"+CMfieldnumber+"Step4"+"p"+p.csv" is created which contains contains solutions [x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3].
#-------------------------------------------------------------------------------------------------------------------------------------------
# Remark: max_bdp can be chosen to be the actual maximum bound for primes found in Step 3.
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import *   # import sage library
from time import time
from random import randint

import os.path
from os import path

import csv

################################################################################
T = time()

filelocation = "./"
CMfieldnumber = sage_eval(sys.argv[1])
min_bdp = sage_eval(sys.argv[2])
max_bdp = sage_eval(sys.argv[3])


load('CM_Fields.sage')


print ("CMfieldnumber = ",CMfieldnumber,", max_bdp = ",max_bdp)
print (sys.version)
print ("----------------------------------------------------------------------")
sys.stdout.flush()


load('GoodNorms.sage')


@parallel
def find_potential_solutions_new_and_elim_bad2(reader, p):
    OUTPUT = []
    R = Integers(p)
    GoodNorms = good_norms(p)
    for row in reader: #x, dOVERn, a, cOVERn, b, gamma, n
        x = QQ(row["x"])
        dOVERn = QQ(row["dOVERn"])
        a = QQ(row["a"])
        cOVERn = QQ(row["cOVERn"])
        b = QQ(row["b"])
        gamma = QQ(row["gamma"])
        n = QQ(row["n"])
        gOVERn = QQ(gamma/n)
        aOVERn = QQ(a/n)
        bOVERn = QQ(b/n)
        k1 = QQ(aOVERn*dOVERn - aOVERn*x - bOVERn)
        k2 = QQ(aOVERn*cOVERn + bOVERn*x)
        k3 = QQ(k1*cOVERn - k2*dOVERn)
        k4 = QQ(k2*cOVERn - bOVERn*b)
        k5 = QQ(aOVERn*a - k1*dOVERn - k2)
        k6 = QQ(aOVERn*b + k1*cOVERn)
        l1 = a*cOVERn + b*x
        l2 = a*dOVERn - a*x - b
        l3 = dOVERn*x - x^2 - a + cOVERn
        for Nx1 in [0..x-1]: # if Nx1 = x, then the bounds imply Nx2 = 0 and Nx3 = 0. Then the second equation would give a = 0 which contradicts a > 0.
            if Nx1 == 0: # see Explanation for case Nx1 = 0 below
                Trx1x2 = 0
                Trx1x3 = 0
                Nx2 = -dOVERn*x^2 + x^3 - a*dOVERn + 2*a*x - cOVERn*x + b
                Nx3 = cOVERn*dOVERn*x^2 - cOVERn*x^3 + a*cOVERn*dOVERn - a*cOVERn*x + cOVERn^2*x + b*x^2 + a*b - b*cOVERn
                Trx2x3 = dOVERn^2*x^2 - dOVERn*x^3 + a*dOVERn^2 - a*dOVERn*x + cOVERn*dOVERn*x - a*x^2 - a^2 - b*dOVERn - b*x
                if Nx2 >=0 and Nx2 in ZZ and Nx3>=0 and Nx3 in ZZ and R(Nx3) in GoodNorms and R(Nx2) in GoodNorms:
                    bd23 = floor(2*sqrt(Nx2*Nx3))
                    Dx2x3 = 4*Nx2*Nx3 - Trx2x3^2 
                    Dx1x2 = 0
                    Dx1x3 = 0
                    Dx1x2x3 = 0
                    if abs(Trx2x3) <= bd23 and Dx2x3/p in ZZ and R(Nx2+Nx3-Trx2x3) in GoodNorms and Dx2x3 != 0:
                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3], [Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nx1, Nx2, Nx2, Trx1x2, Trx1x3, Trx2x3], [Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]])
                        OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
                        sys.stdout.flush()
            else:  # Nx1 != 0
                if R(Nx1) in GoodNorms:
                    upNx2 = a*(x-Nx1)
                    l12 = l1*(x - Nx1) + a*b
                    l22 = l2*(x - Nx1) - a^2
                    l32 = l3*(x - Nx1) + a*dOVERn - b - a*x
                    for Nx2 in [0..upNx2]:
                        if Nx2 == 0: # see Explanation for case Nx1 != 0 and Nx2 = 0:
                            Trx1x2 = 0
                            Trx2x3 = 0
                            if k1 != 0 and (x-Nx1)*n/a in ZZ and (x-Nx1)*n/a == a/k1:
                                Nx3 = a/k1 
                                if R(Nx3) in GoodNorms:
                                    Trx1x3 = - b - (a*aOVERn - k1*dOVERn - k2)*Nx3
                                    bd13 = floor(2*sqrt(Nx1*Nx3))
                                    Dx1x2 = 0
                                    Dx2x3 = 0
                                    Dx1x3 = 4*Nx1*Nx3 - Trx1x3^2 
                                    Dx1x2x3 = 0
                                    if abs(Trx1x3) <= bd13 and Dx1x3/p in ZZ and Dx1x3 != 0 and R(Nx1+Nx3-Trx1x3) in GoodNorms:
                                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3], [Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nx1, Nx2, Nx2, Trx1x2, Trx1x3, Trx2x3], [Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]])
                                        OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
                                        sys.stdout.flush()
                        else: # Nx1 !=0 and Nx2 !=0
                            if R(Nx2) in GoodNorms:
                                l123 = l12 - cOVERn*Nx2
                                l223 = l22 - dOVERn*Nx2
                                l323 = l32 + Nx2
                                bd12 = floor(2*sqrt(Nx1*Nx2))
                                for Trx1x2 in [-bd12..bd12]: # see Explanation for Nx1 !=0 and Nx2 !=0
                                    Dx1x2 = 4*Nx1*Nx2 - Trx1x2^2
                                    if Dx1x2/p in ZZ:
                                        Nx3 = l123 + b*Trx1x2 # Nx3 = ((a*cOVERn + b*x)*x + a*b) - (a*cOVERn + b*x)*Nx1 - cOVERn*Nx2 + b*Trx1x2
                                        if Nx3 in ZZ and R(Nx3) in GoodNorms and Nx3 >=0 and R(Nx1+Nx2-Trx1x2) in GoodNorms:   
                                            upNx3 = min(floor(2/a*(n*(x-Nx1)+(b^2-n)*Nx2/a)),gamma*(x-Nx1))
                                            if Nx3 <= upNx3:
                                                Trx1x3 = l323 + (dOVERn - x)*Trx1x2 #Trx1x3 = (dOVERn*x - x^2 - 2*a + cOVERn)*x + a*dOVERn - b) - (dOVERn*x - x^2 - a + cOVERn)*Nx1 + Nx2 + (dOVERn - x)*Trx1x2
                                                if (Nx3 == 0 and Trx1x3 == 0) or (Nx3 != 0):
                                                    Dx1x3 = 4*Nx1*Nx3 - Trx1x3^2
                                                    bd13 = floor(2*sqrt(Nx1*Nx3))
                                                    upNx2 = floor(2/gamma*(n*(x-Nx1)+(b^2-n)*Nx3/gamma))
                                                    if Trx1x3 in ZZ and abs(Trx1x3) <= bd13 and Dx1x3/p in ZZ and R(Nx1+Nx3-Trx1x3) in GoodNorms and Nx2 <= upNx2:
                                                        Trx2x3 = l223 - a*Trx1x2 #Trx2x3 = ((a*dOVERn - a*x - b)*x - a^2) - (a*dOVERn - a*x - b)*Nx1 - dOVERn*Nx2 - a*Trx1x2
                                                        if ((Nx3 == 0 and Trx2x3 == 0) or (Nx3 != 0)) and R(Nx1+Nx2+Nx3-Trx1x2-Trx2x3-Trx1x3) in GoodNorms:
                                                            Dx2x3 = 4*Nx2*Nx3 - Trx2x3^2
                                                            bd23 = floor(2*sqrt(Nx2*Nx3))
                                                            if Trx2x3 in ZZ and abs(Trx2x3) <= bd23 and Dx2x3/p in ZZ and R(Nx2+Nx3-Trx2x3) in GoodNorms:
                                                                Dx1x2x3 = 4*Nx1*Nx2*Nx3 - Nx1*Trx2x3^2 - Nx2*Trx1x3^2 - Nx3*Trx1x2^2 - Trx1x2*Trx1x3*Trx2x3 # = (Trx1x2x3)^2
                                                                if Dx1x2x3.is_square() and Dx1x2x3/p in ZZ:
                                                                    if  (Dx1x2 == 0 and Dx1x3 !=0 and Dx2x3 !=0 ) or (Dx1x3 == 0 and Dx1x2 !=0 and Dx2x3 !=0 ) or (Dx2x3 == 0 and Dx1x2 !=0 and Dx1x3 !=0 ) or (Dx1x2 == 0 and Dx1x3 ==0 and Dx2x3 !=0 ) or (Dx1x2 == 0 and Dx1x3 !=0 and Dx2x3 ==0 ) or (Dx1x2 != 0 and Dx1x3 ==0 and Dx2x3 ==0 ) or (Dx1x2 != 0 and Dx1x3 !=0 and Dx2x3 !=0 ):
                                                                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3], [Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nx1, Nx2, Nx2, Trx1x2, Trx1x3, Trx2x3], [Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]])
                                                                        OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
                                                                    sys.stdout.flush()
    if len(OUTPUT) > 0:
        print ("For p = ", p, " there are ", len(OUTPUT)," solutions.")
        print ("----------------------------------------------------------------------")
        filename2 = "CMfield"+CMfieldnumber.str()+"Step4"+"p"+p.str()+".csv"
        csvfile2 = open(filelocation+filename2, 'w')
        csvwriter2 = csv.writer(csvfile2)
        csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'Nx1', 'Nx2', 'Nx3', 'Trx1x2', 'Trx1x3', 'Trx2x3', 'Dx1x2', 'Dx1x3', 'Dx2x3', 'Dx1x2x3'))
        for S in OUTPUT:
            csvwriter2.writerow([S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], S[9], S[10], S[11], S[12], S[13], S[14], S[15], S[16]])
            csvfile2.flush()
            sys.stdout.flush()
        csvfile2.close()
        sys.stdout.flush()
    sys.stdout.flush()
    return OUTPUT


for p in list(primes(min_bdp, max_bdp+1)):
    print ("p = ", p)
    if path.exists("CMfield"+CMfieldnumber.str()+"Step3"+"p"+p.str()+".csv") == True:
        filename3 = "CMfield"+CMfieldnumber.str()+"Step3"+"p"+p.str()+".csv"
        reader3 = csv.DictReader(open(filelocation+filename3, 'r'))
        potsolutions = find_potential_solutions_new_and_elim_bad2(reader3,p)
        sys.stdout.flush()


print ("All done in ", time()- T)