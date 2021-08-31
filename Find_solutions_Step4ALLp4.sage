# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements Algorithm 3 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step4ALLp4.sage CMfieldnumber min_bdp max_bdp
# e.g. sage ./Find_solutions_Step4ALLp4.sage 1 2 11 runs Step 4 for CMfield 1 for primes 1 < p < 12
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 4: For each prime p between min_bdp and  max_bdp, read all solutions [x, dOVERn, a, cOVERn, b, gamma, n] from the corresponding file "CMfield"+CMfieldnumber+"Step3"+"p"+p.csv" created in Step 3. 
# For each [x, dOVERn, a, cOVERn, b, gamma, n] from Step 3, find solutions [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3] that satisfy the following:
# 1) the 3 equations in Proposition 3.8 and the implied bounds in Lemma 5.6.
# 2) Nd1, Nd2, Nd3, N(d1+d2), N(d1+d3), N(d2+d3), N(d1+d2+d3) are achievable norms of elements with trace 0 in the maximal order (i.e. satisfy certain congruence conditions mod p as in Section 5.7.4)
# 3) |Tr(di*dj)| \leq floor(2*sqrt(Ndi*Ndj)) for all i,j=1,2,3 (see Lemma Lemma 5.3)
# 4) Ddidj = 4*Ndi*Ndj - Trdidj^2 are integers divisible by p for all i,j=1,2,3 (see Lemma 5.14 in Section 5.6)
# 5) Dd1d2d3 = 4*Nd1*Nd2*Nd3 - Nd1*Trd2d3^2 - Nd2*Trd1d2^2 - Nd3*Trd1d3^2 - Trd1d2*Trd1d3*Trd2d3 is an integer divisible by p and is a square (see Lemma 5.14 in Section 5.6)
# 6) if Ndi = 0 then Trdidj = Trdjdi = 0 for all i,j=1,2,3.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Brief explanation for 2) above:
# Know: d1,d2,d3 in max order O in quaternion algebra QA with Tr(d1)=Tr(d2)=Tr(d3)=0, N(d1)=Nd1, N(d2)=Nd2, N(d3)=Nd3, Tr(d1*d2)=Trd1d2, Tr(d1*d3)=Trd1d3, Tr(d2*d3)=Trd2d3.
# Know: Tr(d1*d2)= N(d1)+N(d2)-N(d1+d2) and Tr(d1*d2) + Tr(d1*d3) + Tr(d2*d3) = N(d1) + N(d2) + N(d3) - N(d1+d2+d3)
# Know: Tr(d1*d2*d3) = Tr(d2*d3*d1) = Tr(d3*d1*d2) = - Tr(d1*d3*d2) = - Tr(d3*d2*d1) = - Tr(d2*d1*d3).
# Check: N(d1)+N(d2)-Tr(d1*d2), N(d1)+N(d3)-Tr(d1*d3), N(d2)+N(d3)-Tr(d2*d3), N(d1) + N(d2) + N(d3) - Tr(d1*d2) - Tr(d1*d3) - Tr(d2*d3) in GoodNorms.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# For each prime p between min_bdp and  max_bdp, one file "CMfield"+CMfieldnumber+"Step4"+"p"+p.csv" is created which contains contains solutions [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3].
#-------------------------------------------------------------------------------------------------------------------------------------------
# Remark: max_bdp can be chosen to be the actual maximum bound for primes found in Step 3.
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import *   # import sage library
from time import time
from random import randint

import os.path
from os import path

import csv

######################################################
T = time()

filelocation = "./"
CMfieldnumber = sage_eval(sys.argv[1])
min_bdp = sage_eval(sys.argv[2])
max_bdp = sage_eval(sys.argv[3])


load('CM_Fields.sage')
    
        
print ("CMfieldnumber = ",CMfieldnumber,", D = ",D,", min_bdp = ",min_bdp,", max_bdp = ",max_bdp)
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
        l = dOVERn*x - x^2 + cOVERn
        l1 = a*cOVERn + b*x
        l2 = a*dOVERn - a*x - b
        l3 = l - a
        for Nd1 in [0..D]: # 
            if Nd1 == 0: # see Explanation for case Nd1 = 0 below
                Trd1d2 = 0
                Trd1d3 = 0
                Nd2 = -(l-a)*D #Nd2 = -(dOVERn*x - x^2 + cOVERn - a)*D
                Nd3 = (cOVERn*l + b*x)*D #Nd3 = (cOVERn*(dOVERn*x - x^2 + cOVERn) + b*x)*D
                Trd2d3 = (dOVERn*l - a*x - b)*D #Trd2d3 = (dOVERn*(dOVERn*x - x^2 + cOVERn) - a*x - b)*D
                if Nd2 >=0 and Nd2 in ZZ and Nd3>=0 and Nd3 in ZZ and R(Nd3) in GoodNorms and R(Nd2) in GoodNorms:
                    bd23 = floor(2*sqrt(Nd2*Nd3))
                    Dd2d3 = 4*Nd2*Nd3 - Trd2d3^2 
                    Dd1d2 = 0
                    Dd1d3 = 0
                    Dd1d2d3 = 0
                    if abs(Trd2d3) <= bd23 and Dd2d3/p in ZZ and R(Nd2+Nd3-Trd2d3) in GoodNorms:
                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                        OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                        sys.stdout.flush()
            else:  # Nd1 != 0
                if R(Nd1) in GoodNorms:
                    upNd2 = a*(D - Nd1)
                    l12 = l1*(D - Nd1)
                    l22 = l2*(D - Nd1)
                    l32 = l3*(D - Nd1) 
                    for Nd2 in [0..upNd2]:
                        if Nd2 == 0: # see Explanation for case Nd1 != 0 and Nd2 = 0:
                            Trd1d2 = 0
                            Trd2d3 = 0
                            if k1 != 0 and Nd1 == D and D*n/p in ZZ: #D*n/p in ZZ is due to Prop 4.3
                                Nd3 = 0
                                Trd1d3 = 0
                                Dd1d2 = 0
                                Dd2d3 = 0
                                Dd1d3 = 0
                                Dd1d2d3 = 0
                                print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                                sys.stdout.flush()
                            elif k1 == 0:
                                Nd3 = (D-Nd1)*n/a
                                if R(Nd3) in GoodNorms:
                                    Trd1d3 = - k5*Nd3
                                    bd13 = floor(2*sqrt(Nd1*Nd3))
                                    Dd1d2 = 0
                                    Dd2d3 = 0
                                    Dd1d3 = 4*Nd1*Nd3 - Trd1d3^2 
                                    Dd1d2d3 = 0
                                    if abs(Trd1d3) <= bd13 and Dd1d3/p in ZZ and R(Nd1+Nd3-Trd1d3) in GoodNorms:
                                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                        OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                                        sys.stdout.flush()
                        else: # Nd1 !=0 and Nd2 !=0
                            if R(Nd2) in GoodNorms:
                                l123 = l12 - cOVERn*Nd2
                                l223 = l22 - dOVERn*Nd2
                                l323 = l32 + Nd2
                                bd12 = floor(2*sqrt(Nd1*Nd2))
                                lbd12 = max(-bd12, floor(-l123/b))
                                for Trd1d2 in [lbd12..bd12]: # see Explanation for Nd1 !=0 and Nd2 !=0
                                    Dd1d2 = 4*Nd1*Nd2 - Trd1d2^2
                                    if Dd1d2/p in ZZ: 
                                        Nd3 = l123 + b*Trd1d2 # Nd3 = (a*cOVERn + b*x)*(D-Nd1) + Trd1d2*b - cOVERn*Nd2
                                        if Nd3 in ZZ and R(Nd3) in GoodNorms and Nd3 >=0 and R(Nd1+Nd2-Trd1d2) in GoodNorms:   
                                            upNd3 = min(floor(2/a*(n*(D-Nd1)+(b^2-n)*Nd2/a)),gamma*(D-Nd1))
                                            if Nd3 <= upNd3:
                                                Trd1d3 = l323 + (dOVERn - x)*Trd1d2 #Trd1d3 = (dOVERn*x - x^2 - a + cOVERn)*(D - Nd1) + Nd2 + (dOVERn - x)*Trd1d2
                                                if (Nd3 == 0 and Trd1d3 == 0) or (Nd3 != 0):
                                                    Dd1d3 = 4*Nd1*Nd3 - Trd1d3^2
                                                    bd13 = floor(2*sqrt(Nd1*Nd3))
                                                    upNd2 = floor(2/gamma*(n*(x-Nd1)+(b^2-n)*Nd3/gamma))
                                                    if Trd1d3 in ZZ and abs(Trd1d3) <= bd13 and Dd1d3/p in ZZ and R(Nd1+Nd3-Trd1d3) in GoodNorms and Nd2 <= upNd2:
                                                        Trd2d3 = l223 - a*Trd1d2 # Trd2d3 = (a*dOVERn - a*x - b)*(D-Nd1) - dOVERn*Nd2 - Trd1d2*a
                                                        if (Nd3 == 0 and Trd2d3 == 0) or (Nd3 != 0):
                                                            Dd2d3 = 4*Nd2*Nd3 - Trd2d3^2
                                                            bd23 = floor(2*sqrt(Nd2*Nd3))        
                                                            if Trd2d3 in ZZ and abs(Trd2d3) <= bd23 and Dd2d3/p in ZZ and R(Nd2+Nd3-Trd2d3) in GoodNorms and R(Nd1+Nd2+Nd3-Trd1d2-Trd2d3-Trd1d3) in GoodNorms:
                                                                Dd1d2d3 = 4*Nd1*Nd2*Nd3 - Nd1*Trd2d3^2 - Nd2*Trd1d3^2 - Nd3*Trd1d2^2 - Trd1d2*Trd1d3*Trd2d3 # = (Trd1d2d3)^2
                                                                if Dd1d2d3.is_square() and Dd1d2d3/p in ZZ:
                                                                    if (Dd1d2 != 0) or (Dd1d3 != 0) or (Dd1d2d3 != 0) or (Dd1d2 == 0 and Dd1d3 == 0 and Dd2d3 == 0 and Dd1d2d3 == 0):
                                                                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                                                        OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                                                                        sys.stdout.flush()
    if len(OUTPUT) > 0:
        print ("For p = ", p, " there are ", len(OUTPUT)," solutions.")
        print ("----------------------------------------------------------------------")
        filename2 = "CMfield"+CMfieldnumber.str()+"Step4"+"p"+p.str()+".csv"
        csvfile2 = open(filelocation+filename2, 'w')
        csvwriter2 = csv.writer(csvfile2)
        csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'Nd1', 'Nd2', 'Nd3', 'Trd1d2', 'Trd1d3', 'Trd2d3', 'Dd1d2', 'Dd1d3', 'Dd2d3', 'Dd1d2d3'))
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


#------------------------------------------------------------------------------------------------------------------
# Explanation for case Nd1 = 0:
# x, dOVERn,a,cOVERn,b,gamma,n,D = var('x, dOVERn,a,cOVERn,b,gamma,n,D')
# aOVERn = a/n
# bOVERn = b/n
# gammaOVERn = gamma/n
# gamma = a*cOVERn + b*dOVERn
# k1 = aOVERn*dOVERn - aOVERn*x - bOVERn
# k2 = aOVERn*cOVERn + bOVERn*x
# A = matrix(SR,Integer(3),Integer(3),[-gammaOVERn,-aOVERn,-bOVERn, k1*cOVERn-k2*dOVERn,k1,-k2, k2*cOVERn-b*bOVERn, -a*aOVERn+k1*dOVERn+k2, -a*bOVERn-k1*cOVERn])
# B = matrix(SR,3,1,[-D,0,0])
# A.determinant().expand().factor()
# = (a*cOVERn*dOVERn*x + b*dOVERn^2*x - a*cOVERn*x^2 - b*dOVERn*x^2 + a*cOVERn^2 + b*cOVERn*dOVERn - dOVERn*gamma*x + gamma*x^2 - b^2 + a*gamma - cOVERn*gamma)*(a^2*cOVERn + a*b*dOVERn - b^2)/n^3
# (a*cOVERn*dOVERn*x + b*dOVERn^2*x - a*cOVERn*x^2 - b*dOVERn*x^2 + a*cOVERn^2 + b*cOVERn*dOVERn - dOVERn*gamma*x + gamma*x^2 - b^2 + a*gamma - cOVERn*gamma).expand()
# = a^2*cOVERn + a*b*dOVERn - b^2
# Then det(A) = (a^2*cOVERn + a*b*dOVERn - b^2)^2/n^3.
# Use b^2= a*gamma - n = a^2*cOVERn + a*b*dOVERn - n. Thus, n = a^2*cOVERn + a*b*dOVERn - b^2.
# Then det(A)=1/n. So A is always invertible.
# S=(A.inverse()*B)
# Nd2=S[0]
# Nd3=S[1]
# Trd2d3=S[2] 
# S[0]
# = (-D*(((cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*n/gamma + (dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*(dOVERn*(a*dOVERn/n - a*x/n - b/n) - (cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*a/gamma - a^2/n + a*cOVERn/n + b*x/n)*n/(((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n)*gamma))*(b/gamma - a*((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*b/gamma - a*cOVERn/n - b*x/n)/(((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n)*gamma))/(cOVERn*(a*dOVERn/n - a*x/n - b/n) + (dOVERn*(a*dOVERn/n - a*x/n - b/n) - (cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*a/gamma - a^2/n + a*cOVERn/n + b*x/n)*((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*b/gamma - a*cOVERn/n - b*x/n)/((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n) + (cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*b/gamma + a*b/n) + (dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a*n/(((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n)*gamma^2) - n/gamma))
# (copy what S[0] displays).expand().factor()
# = -(dOVERn*x - x^2 - a + cOVERn)*D*n/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Nd2 = -(dOVERn*x - x^2 + cOVERn - a)*D
# Similarly, 
# S[1] gives (cOVERn*dOVERn*x - cOVERn*x^2 + cOVERn^2 + b*x)*D*n/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Nd3 = (cOVERn*dOVERn*x - cOVERn*x^2 + cOVERn^2 + b*x)*D
# Thus, Nd3 = (cOVERn*(dOVERn*x - x^2 + cOVERn) + b*x)*D = cOVERn*(dOVERn*x - x^2 + cOVERn)*D + b*x*D = cOVERn*(a*D-Nd2) + b*x*D = -cOVERn*Nd2 + (cOVERn*a + b*x)*D
# and
# S[2] gives (dOVERn^2*x - dOVERn*x^2 + cOVERn*dOVERn - a*x - b)*D*n/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Trd2d3 = (dOVERn^2*x - dOVERn*x^2 + cOVERn*dOVERn - a*x - b)*D
# Thus, Trd2d3 = (dOVERn*(dOVERn*x - x^2 + cOVERn) - a*x - b)*D = dOVERn*(a*D-Nd2) - (a*x + b)*D = -dOVERn*Nd2 + (dOVERn*a - a*x - b)*D

#------------------------------------------------------------------------------------------------------------------
# Explanation for case Nd1 != 0 and Nd2 = 0:
# Then Trd1d2 = Trd2d3 = 0
# The 3 equations become: D = Nd1 + (a/n)*Nd3, 0 = k1*Nd3, 0 = - (a^2/n - k1*d/n - k2)*Nd3 - Trd1d3 = -k5*Nd3 - Trd1d3.
# If k1 = 0 then Nd3 = (D-Nd1)*n/a and Trd1d3 = - k5*Nd3.
# Otherwise, if k1 !=0 (which is actually always the case), then Nd3 = 0, D = Nd1, Trd1d3 = 0.

#------------------------------------------------------------------------------------------------------------------
# Explanation for Nd1 !=0 and Nd2 !=0
# x, dOVERn,a,cOVERn,b,gamma,n,D = var('x, dOVERn,a,cOVERn,b,gamma,n,D')
# aOVERn = a/n
# bOVERn = b/n
# gammaOVERn = gamma/n
# gamma = a*cOVERn + b*dOVERn
# k1 = aOVERn*dOVERn - aOVERn*x - bOVERn
# k2 = aOVERn*cOVERn + bOVERn*x
# A = matrix(SR,3,3,[aOVERn,bOVERn,0, k1,-k2,0, -a*aOVERn+k1*dOVERn+k2, -a*bOVERn-k1*cOVERn,-1])
# Nd1,Nd2,Trd1d2 = var('Nd1,Nd2,Trd1d2')
# B = matrix(SR,3,1,[D-Nd1-gammaOVERn*Nd2, -(k1*cOVERn-k2*dOVERn)*Nd2 + Trd1d2, -(k2*cOVERn-b*bOVERn)*Nd2])
# A.determinant().expand().factor()
# = (a^2*cOVERn + a*b*dOVERn - b^2)/n^2
# Use b^2= a*gamma - n = a^2*cOVERn + a*b*dOVERn - n. Thus, n = a^2*cOVERn + a*b*dOVERn - b^2.
# Then det(A)= 1/n. So A is always invertible.
# S = (A.inverse()*B)
# Nd3 = S[0]
# Trd2d3 = S[1]
# Trd1d3 = S[2] = k4*Nd2 - k5*Nd3 - k6*Trd2d3 - b

# S[0]
# = ((D - Nd1 - Nd2*gamma/n)*(n/a - b*(a*dOVERn/n - a*x/n - b/n)*n/(a^2*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))) + ((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*Nd2 + Trd1d2)*b/(a*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n)))
# (copy what S[0] displays).expand().factor()
# = -(Nd2*a^2*cOVERn^2 + Nd2*a*b*cOVERn*dOVERn - Nd2*b^2*cOVERn - D*a*cOVERn*n + Nd1*a*cOVERn*n - D*b*n*x + Nd1*b*n*x - Trd1d2*b*n)/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Nd3 = (D*a*cOVERn*n + D*b*n*x + Trd1d2*b*n - (a*cOVERn*n + b*n*x)*Nd1 - (a^2*cOVERn^2 + a*b*cOVERn*dOVERn - b^2*cOVERn)*Nd2)/n
# = ((a*cOVERn + b*x)*D*n + Trd1d2*b*n - (a*cOVERn + b*x)*n*Nd1 - (a^2*cOVERn + a*b*dOVERn - b^2)*cOVERn*Nd2)/n
# = ((a*cOVERn + b*x)*D*n + Trd1d2*b*n - (a*cOVERn + b*x)*n*Nd1 - n*cOVERn*Nd2)/n
# = (a*cOVERn + b*x)*D + Trd1d2*b - (a*cOVERn + b*x)*Nd1 - cOVERn*Nd2
# Thus, Nd3 = (a*cOVERn + b*x)*(D - Nd1) - cOVERn*Nd2 + Trd1d2*b
# Put l1 = a*cOVERn + b*x
# Thus, Nd3 = l1*(D-Nd1) - cOVERn*Nd2 + Trd1d2*b
# Put l12 = l1*(D-Nd1)
# Thus, Nd3 = l12 - cOVERn*Nd2 + Trd1d2*b
# Put l123 = l12 - cOVERn*Nd2
# Thus, Nd3 = l123 + b*Trd1d2

# S[1]
# = ((D - Nd1 - Nd2*gamma/n)*(a*dOVERn/n - a*x/n - b/n)*n/(a*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n)) - ((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*Nd2 + Trd1d2)/(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))
# (copy what S[1] displays).expand().factor()
# = -(Nd2*a^2*cOVERn*dOVERn + Nd2*a*b*dOVERn^2 - Nd2*b^2*dOVERn - D*a*dOVERn*n + Nd1*a*dOVERn*n + D*a*n*x - Nd1*a*n*x + Trd1d2*a*n + D*b*n - Nd1*b*n)/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Trd2d3 = -(Nd2*a^2*cOVERn*dOVERn + Nd2*a*b*dOVERn^2 - Nd2*b^2*dOVERn - D*a*dOVERn*n + Nd1*a*dOVERn*n + D*a*n*x - Nd1*a*n*x + Trd1d2*a*n + D*b*n - Nd1*b*n)/n
# = ((a*dOVERn - a*x - b)*D*n - (a*dOVERn - a*x - b)*n*Nd1 - (a^2*cOVERn + a*b*dOVERn - b^2)*dOVERn*Nd2 - Trd1d2*a*n)/n
# = ((a*dOVERn - a*x - b)*D*n - (a*dOVERn - a*x - b)*n*Nd1 - n*dOVERn*Nd2 - Trd1d2*a*n)/n
# = (a*dOVERn - a*x - b)*D - (a*dOVERn - a*x - b)*Nd1 - dOVERn*Nd2 - Trd1d2*a
# Thus, Trd2d3 = (a*dOVERn - a*x - b)*(D - Nd1) - dOVERn*Nd2 - Trd1d2*a
# Put l2 = a*dOVERn - a*x - b
# Thus, Trd2d3 = l2*(D-Nd1) - dOVERn*Nd2 - Trd1d2*a
# Put l22 = l2*(D - Nd1)
# Thus, Trd2d3 = l22 - dOVERn*Nd2 - a*Trd1d2 
# Put l223 = l22 - dOVERn*Nd2
# Thus, Trd2d3 = l223 - a*Trd1d2 

# S[2]
# = ((cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*Nd2 + (D - Nd1 - Nd2*gamma/n)*((dOVERn*(a*dOVERn/n - a*x/n - b/n) - a^2/n + a*cOVERn/n + b*x/n)*n/a - (cOVERn*(a*dOVERn/n - a*x/n - b/n) + (dOVERn*(a*dOVERn/n - a*x/n - b/n) - a^2/n + a*cOVERn/n + b*x/n)*b/a + a*b/n)*(a*dOVERn/n - a*x/n - b/n)*n/(a*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))) + ((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*Nd2 + Trd1d2)*(cOVERn*(a*dOVERn/n - a*x/n - b/n) + (dOVERn*(a*dOVERn/n - a*x/n - b/n) - a^2/n + a*cOVERn/n + b*x/n)*b/a + a*b/n)/(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))
# (copy what S[2] displays).expand().factor()
# = (Nd2*a^2*cOVERn + Nd2*a*b*dOVERn + D*dOVERn*n*x - Nd1*dOVERn*n*x - D*n*x^2 + Nd1*n*x^2 - Nd2*b^2 - D*a*n + Nd1*a*n + D*cOVERn*n - Nd1*cOVERn*n + Trd1d2*dOVERn*n - Trd1d2*n*x)/n
# Thus, Trd1d3 = (Nd2*a^2*cOVERn + Nd2*a*b*dOVERn + D*dOVERn*n*x - Nd1*dOVERn*n*x - D*n*x^2 + Nd1*n*x^2 - Nd2*b^2 - D*a*n + Nd1*a*n + D*cOVERn*n - Nd1*cOVERn*n + Trd1d2*dOVERn*n - Trd1d2*n*x)/n
# = ((dOVERn*x - x^2 - a + cOVERn)*D*n - (dOVERn*x - x^2 - a + cOVERn)*n*Nd1 + (a^2*cOVERn + a*b*dOVERn - b^2)*Nd2 + (dOVERn*n - n*x)*Trd1d2)/n
# = ((dOVERn*x - x^2 - a + cOVERn)*D*n - (dOVERn*x - x^2 - a + cOVERn)*n*Nd1 + n*Nd2 + (dOVERn - x)*n*Trd1d2)/n
# Thus, Trd1d3 = (dOVERn*x - x^2 - a + cOVERn)*(D - Nd1) + Nd2 + (dOVERn - x)*Trd1d2
# Put l3 = dOVERn*x - x^2 - a + cOVERn = l - a
# Thus, Trd1d3 = l3*(D - Nd1) + Nd2 + (dOVERn - x)*Trd1d2
# Put l32 = l3*(D - Nd1) 
# Thus, Trd1d3 = l32 + Nd2 + (dOVERn - x)*Trd1d2
# Put l323 = l32 + Nd2
# Thus, Trd1d3 = l323 + (dOVERn - x)*Trd1d2
#------------------------------------------------------------------------------------------------------------------
