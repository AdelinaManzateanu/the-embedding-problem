# run as sage ./From_Lambda3_to_Lambda2.sage p max_order CMfieldnumber D
# e.g. sage ./From_Lambda3_to_Lambda2.sage 3 1 0 0 0 0 1 0 0 0 1/2 0 1/2 1/2 0 1/2 0 2 1 runs p=3, O=1 with basis [1, i, i/2 + k/2, 1/2 + j/2], for CM field 2, D = 1.
#-------------------------------------------------------------------------------------------------------------------------------------------
# 1) Read solutions Lambda1, Lambda3 in the form [x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3].
# 2) Find expression of Lambda2 in terms of those elements and Lambda1 using Lambda2^2 = - Lambda1.
# 3) Compute Lambda2
# 4) Express solution in the form [x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3]
# 5) Check it lifts to the maximal order and the 3 equations are satisfied.
# 6) Write solutions [x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3] in a file.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Explanation on how to find expresion of II = Lambda3 (Magma):
# A:=13; B:=50; C:=49; //Modify accordingly
# R<x>:=PolynomialRing(Rationals());
# f:=x^6+A*x^4+B*x^2+C;
# K<V>:=NumberField(f); //V is Lambda2
# OK := MaximalOrder(K); O:= EquationOrder(K);

# L<II>:=NumberField(x^2+1);  //So II^2 = has -D in diagonal and here D = 1. Modify accordingly.
# OL := MaximalOrder(L); OO:= EquationOrder(L);
# // What is II=Lambda3 in terms of Lambda2?
# Roots(PolynomialRing(K)!DefiningPolynomial(L));
# Thus can choose II = 1/7*(V^5 + 6*V^3 + V) or  II = 1/7*(-V^5 - 6*V^3 - V)

# 2) Find expression of Lambda2=V in terms of those elements and Lambda1=U.
# I = matrix(QQ,3,3, [1,0,0, 0,1,0, 0,0,1])
# II = matrix(QA,3,3, [QA(elem with norm D and trace 0),0,0, 0,QA(elem with norm D and trace 0),0, 0,0,QA(elem with norm D and trace 0)])
# CMfieldnumber = 1, D = 1
# II = (V+6*V^3+V^5)/7 = V*(I+6*V^2+V^4)/7 = V*(I-6*U+U^2)/7 => V:= 7*II*(I-6*U+U^2)^(-1)
# II = (-V^5 - 6*V^3 - V)/7 = V*(-V^4 - 6*V^2 - I)/7 = V*(-U^2 + 6*U - I)/7 => V:= 7*II*(-U^2 + 6*U - I)^(-1)
# CMfieldnumber = 2, D = 1
# II = 3*V+V^3 = V*(3*I+V^2) = V*(3*I - U) => V:= II*(3*I - U)^(-1);
# II = - 3*V-V^3 = V*(-3*I-V^2) = V*(-3*I + U) => V:= II*(-3*I + U)^(-1);
# CMfieldnumber = 22, D = 1
# II = 9*V + 26*V^3 + 3*V^5 = V*(9*I + 26*V^2 + 3*V^4) = V*(9*I - 26*U + 3*U^2) => V := II*(9*I - 26*U + 3*U^2)^(-1)
# CMfieldnumber = 3, D = 1
# II = 3*V + 4*V^3 + V^5 = V*(3*I + 4*V^2 + V^4) = V*(3*I - 4*U + U^2) => V := II*(3*I - 4*U + U^2)^(-1);
# CMfieldnumber = 33, D = 1
# II = 5*V + 11*V^3 + 2*V^5 = V*(5*I + 11*V^2 + 2*V^4) = V*(5*I - 11*U + 2*U^2) => V := II*(5*I - 11*U + 2*U^2)^(-1);
# CMfieldnumber = 4, D = 7
# II = 7*V + 6*V^3 + V^5 = V*(7*I + 6*V^2 + V^4) = V*(7*I - 6*U + U^2) => V := II*(7*I - 6*U + U^2)^(-1);
# CMfieldnumber = 5, D = 7
# II = (21*V + V^3)/11 = V*(21*I + V^2)/11 = V*(21*I - U)/11 => V := 11*II*(21*I - U)^(-1);
# CMfieldnumber = 55, D = 7
# II = (63*V + 62*V^3 + V^5)/3 = V*(63*I + 62*V^2 + V^4)/3 = V*(63*I - 62*U + U^2)/3 => V := 3*II*(63*I - 62*U + U^2)^(-1);
# CMfieldnumber = 6, D = 1
# II = (300*V + 37*V^3 + V^5)/176 = V*(300*I + 37*V^2 + V^4)/176 = V*(300*I - 37*U + U^2)/176 => V := 176*II*(300*I - 37*U + U^2)^(-1);
# CMfieldnumber = 66, D = 1
# II = (70*V + 19*V^3 + V^5)/44 = V*(70*I + 19*V^2 + V^4)/44 = V*(70*I - 19*VU + U^2)/44 => V := 44*II*(70*I - 19*U + U^2)^(-1);
# CMfieldnumber = 666, D = 1
# II = (788*V + 153*V^3 + 5*V^5)/64 = V*(788*I + 153*V^2 + 5*V^4)/64 = V*(788*I - 153*U + 5*U^2)/64 => V := 64*II*(788*I - 153*U + 5*U^2)^(-1);
# CMfieldnumber = 7, D = 1
# II = (28*V + 13*V^3 + V^5)/16 = V*(28*I + 13*V^2 + V^4)/16 = V*(28*I - 13*U + U^2)/16 => V := 16*II*(28*I - 13*U + U^2)^(-1);
# CMfieldnumber = 8, D = 1
# II = (21*V + V^3)/28 = V*(21*I + V^2)/28 = V*(21*I - U)/28 => V := 28*II*(21*I - U)^(-1);
# CMfieldnumber = 88, D = 1
# II = (180*V+29*V^3 +V^5)/96 = V*(180*I + 29*V^2 + V^4)/96 = V*(180*I - 29*U + U^2)/96 => V:= 96*II*(180*I - 29*U + U^2)^(-1);
# CMfieldnumber = 9, D = 2
# II = (76*V + 20*V^3 + V^5)/20 = V*(76*I + 20*V^2 + V^4)/20 = V*(76*I - 20*U + U^2)/20 => V := 20*II*(76*I - 20*U + U^2)^(-1);
# CMfieldnumber = 10, D = 7
# II = (931*V + 70*V^3 + V^5)/245 = V*(931*I + 70*V^2 + V^4)/245 = V*(931*I - 70*U + U^2)/245 => V := 245*II*(931*I - 70*U + U^2)^(-1);
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import * 
from time import time
from random import randint

import csv

#--------------------------
filelocation = "./"
p = sage_eval(sys.argv[1])

if p == 2:
    QA = QuaternionAlgebra(SR, -1,-1, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
else:
    S = Integers(8)
    if S(p) == 3 or S(p) == 7:
        QA = QuaternionAlgebra(SR, -1,-p, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
    elif S(p) == 5:
        QA = QuaternionAlgebra(SR, -2,-p, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
    elif S(p) == 1:
        print ("Choose a prime l = 3 mod 4 and l not a square mod p. Choose QA(-1,-l).") #QA = QuaternionAlgebra(SR, -1,-l, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)

Mlist = [sage_eval(sys.argv[2]),sage_eval(sys.argv[3]),sage_eval(sys.argv[4]), sage_eval(sys.argv[5]),sage_eval(sys.argv[6]),sage_eval(sys.argv[7]),sage_eval(sys.argv[8]), sage_eval(sys.argv[9]),sage_eval(sys.argv[10]),sage_eval(sys.argv[11]),sage_eval(sys.argv[12]), sage_eval(sys.argv[13]),sage_eval(sys.argv[14]),sage_eval(sys.argv[15]),sage_eval(sys.argv[16]), sage_eval(sys.argv[17])]
M = matrix(QQ,4,Mlist)
O = [ QA(M[0][0]+M[0][1]*i+M[0][2]*j+M[0][3]*k), QA(M[1][0]+M[1][1]*i+M[1][2]*j+M[1][3]*k), QA(M[2][0]+M[2][1]*i+M[2][2]*j+M[2][3]*k), QA(M[3][0]+M[3][1]*i+M[3][2]*j+M[3][3]*k)]
#M = matrix([O[0].coefficient_tuple(), O[1].coefficient_tuple(), O[2].coefficient_tuple(), O[3].coefficient_tuple()])

CMfieldnumber = sage_eval(sys.argv[18])
D = sage_eval(sys.argv[19])

print ("CMfieldnumber = ",CMfieldnumber,", p = ",p,", max_order = ",O,", D = ",D)
print (sys.version)
print ("----------------------------------------------------------------------")
sys.stdout.flush()


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


basis = M.transpose().inverse()


def elem_in_order(element,basis): #"""Given a basis for an order in a quaternion algebra, checks if elem is in order"""
    v = vector(QA(element).coefficient_tuple())
    s = basis*v
    for nn in range(4):
        if s[nn] not in ZZ:
            return False
    return True


def matrix_with_entries_in_R_and_ROVERn(MM,n,basis): #"""Given a matrix check if first row entries are in R and second and third row entries are in R/n
    test = False
    nMM = matrix(QA, 3, [MM[0][0], MM[0][1], MM[0][2], n*MM[1][0], n*MM[1][1], n*MM[1][2], n*MM[2][0], n*MM[2][1], n*MM[2][2]])
    if elem_in_order(nMM[0][0],basis) == True and elem_in_order(nMM[0][1],basis) == True and elem_in_order(nMM[0][2],basis) == True:
        if elem_in_order(nMM[1][0],basis) == True and elem_in_order(nMM[1][1],basis) == True and elem_in_order(nMM[1][2],basis) == True:
            if elem_in_order(nMM[2][0],basis) == True and elem_in_order(nMM[2][1],basis) == True and elem_in_order(nMM[2][2],basis) == True:
                test = True
    return test


#Read the candidate solutions: 
print ("reading Solutions file...")
sys.stdout.flush()
tread = time()


filename1 = "CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"
for m in [0..3]:
    for n in [0..3]:
        if O[m][n] ==1:
            if n == 0:
                filename1 = filename1+"1"
            if n == 1:
                filename1 = filename1+"+i"
            elif n == 2:
                filename1 = filename1+"+j"
            elif n ==3:
                filename1 = filename1+"+k"
        elif O[m][n] !=0 and O[m][n] !=1:
            if n == 0:
                if O[m][n] in ZZ:
                    filename1 = filename1+QQ(O[m][n]).str()
                else:
                    filename1 = filename1+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"
            if n == 1:
                if O[m][n] in ZZ:
                    filename1 = filename1+"+"+QQ(O[m][n]).str()+"i"
                else:
                    filename1 = filename1+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"i"
            elif n == 2:
                if O[m][n] in ZZ:
                    filename1 = filename1+"+"+QQ(O[m][n]).str()+"j"
                else:
                    filename1 = filename1+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"j"
            elif n ==3:
                if O[m][n] in ZZ:
                    filename1 = filename1+"+"+QQ(O[m][n]).str()+"k"
                else:
                    filename1 = filename1+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"k"
    filename1 = filename1+"_"
filename1 = filename1+"Solutions_for_di_Step7.csv"
reader1 = csv.DictReader(open(filelocation + filename1, 'r'))
L = []
for row in reader1: #x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3
    L += [[QQ(row["x"]), QQ(row["dOVERn"]), QQ(row["a"]), QQ(row["cOVERn"]), QQ(row["b"]), QQ(row["gamma"]), QQ(row["n"]), eval(preparse(row["d1"])), eval(preparse(row["d2"])), eval(preparse(row["d3"])), eval(preparse(row["d4"])), eval(preparse(row["d5"])), eval(preparse(row["d6"])), eval(preparse(row["d7"])), eval(preparse(row["d8"])), eval(preparse(row["d9"])), eval(preparse(row["Dd1d2"])), eval(preparse(row["Dd1d3"])), eval(preparse(row["Dd2d3"])), eval(preparse(row["Dd1d2d3"])) ]]       


LengthAllSol = len(L)
print ("Number of solutions for Lambda3 = ", LengthAllSol)


I = matrix(QQ,3,3, [1,0,0, 0,1,0, 0,0,1])


filename2 = "CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"
for m in [0..3]:
    for n in [0..3]:
        if O[m][n] ==1:
            if n == 0:
                filename2 = filename2+"1"
            if n == 1:
                filename2 = filename2+"+i"
            elif n == 2:
                filename2 = filename2+"+j"
            elif n ==3:
                filename2 = filename2+"+k"
        elif O[m][n] !=0 and O[m][n] !=1:
            if n == 0:
                if O[m][n] in ZZ:
                    filename2 = filename2+QQ(O[m][n]).str()
                else:
                    filename2 = filename2+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"
            if n == 1:
                if O[m][n] in ZZ:
                    filename2 = filename2+"+"+QQ(O[m][n]).str()+"i"
                else:
                    filename2 = filename2+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"i"
            elif n == 2:
                if O[m][n] in ZZ:
                    filename2 = filename2+"+"+QQ(O[m][n]).str()+"j"
                else:
                    filename2 = filename2+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"j"
            elif n ==3:
                if O[m][n] in ZZ:
                    filename2 = filename2+"+"+QQ(O[m][n]).str()+"k"
                else:
                    filename2 = filename2+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"k"
    filename2 = filename2+"_"
filename2 = filename2+"Solutions_for_xi.csv"

csvfile2 = open(filelocation + filename2, 'w')
csvwriter2 = csv.writer(csvfile2) #x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3
csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'x1', 'x2', 'x3', 'Dx1x2', 'Dx1x3', 'Dx2x3', 'Dx1x2x3'))


OUTPUT = []
for LL in L: #x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3
    x = QQ(LL[0])
    dOVERn = QQ(LL[1])
    a = QQ(LL[2])
    cOVERn = QQ(LL[3])
    b = QQ(LL[4])
    gamma = QQ(LL[5])
    n = QQ(LL[6])
    Lambda1 = matrix(QQ,3,3,[x,a,b, 1,0,cOVERn, 0,1,dOVERn])
    gOVERn = QQ(gamma/n)
    aOVERn = QQ(a/n)
    bOVERn = QQ(b/n)
    k1 = QQ(aOVERn*dOVERn - aOVERn*x - bOVERn)
    k2 = QQ(aOVERn*cOVERn + bOVERn*x)
    k3 = QQ(k1*cOVERn - k2*dOVERn)
    k4 = QQ(k2*cOVERn - bOVERn*b)
    k5 = QQ(aOVERn*a - k1*dOVERn - k2)
    k6 = QQ(aOVERn*b + k1*cOVERn)
    d1 = QA(LL[7])
    d2 = QA(LL[8])
    d3 = QA(LL[9])
    d4 = QA(LL[10])
    d5 = QA(LL[11])
    d6 = QA(LL[12])
    d7 = QA(LL[13])
    d8 = QA(LL[14])
    d9 = QA(LL[15])
    Lambda3 = matrix(QA, 3,3, [d1,d2,d3, d4,d5,d6, d7,d8,d9])
    if CMfieldnumber == 1:
        Lambda2 = 7*Lambda3*((I-6*Lambda1+Lambda1^2).inverse())
    elif CMfieldnumber == 2:
        Lambda2 = Lambda3*((3*I - Lambda1).inverse())
    elif CMfieldnumber == 22:
        Lambda2 = Lambda3*((9*I - 26*Lambda1 + 3*Lambda1^2).inverse())
    elif CMfieldnumber == 3:
        Lambda2 = Lambda3*((3*I - 4*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 33:
        Lambda2 = Lambda3*((5*I - 11*Lambda1 + 2*Lambda1^2).inverse())
    elif CMfieldnumber == 4:
        Lambda2 = Lambda3*((7*I - 6*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 5:
        Lambda2 = 11*Lambda3*((21*I - Lambda1).inverse())
    elif CMfieldnumber == 55:
        Lambda2 = 3*Lambda3*((63*I - 62*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 6:
        Lambda2 = 176*Lambda3*((300*I - 37*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 66:
        Lambda2 = 44*Lambda3*((70*I - 19*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 666:
        Lambda2 = 64*Lambda3*((788*I - 153*Lambda1 + 5*Lambda1^2).inverse())
    elif CMfieldnumber == 7:
        Lambda2 = 16*Lambda3*((28*I - 13*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 8:
        Lambda2 = 28*Lambda3*((21*I - Lambda1).inverse())
    elif CMfieldnumber == 88:
        Lambda2 = 96*Lambda3*((180*I - 29*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 9:
        Lambda2 = 20*Lambda3*((76*I - 20*Lambda1 + Lambda1^2).inverse())
    elif CMfieldnumber == 10:
        Lambda2 = 245*Lambda3*((931*I - 70*Lambda1 + Lambda1^2).inverse())
    x1 = QA(Lambda2[0][0])
    x2 = QA(Lambda2[0][1])
    x3 = QA(Lambda2[0][2])
    if [x, dOVERn, a, cOVERn, b, gamma, n, QA(x1),QA(x2),QA(x3)] not in OUTPUT and QA(x1).reduced_trace() == 0 and QA(x2).reduced_trace() == 0 and QA(x3).reduced_trace() == 0 and matrix_with_entries_in_R_and_ROVERn(Lambda2,n,basis) == True:
        T = False
        Lambda22 = Lambda2^2
        Lambda23 = Lambda2^3
        Lambda24 = Lambda2^4
        Lambda25 = Lambda2^5
        Lambda26 = Lambda2^6
        if CMfieldnumber == 1: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, 1/7*(Lambda25 + 6*Lambda23 + Lambda2). i = 1/7*I*(Lambda2 + 6*Lambda23 + Lambda25)
            B6 = 1/7*I*(Lambda2 + 6*Lambda23 + Lambda25)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda23,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda24,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 2: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 3*Lambda2 + Lambda23
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda23,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda24,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda25,n,basis) == True:
                T = True
        elif CMfieldnumber == 22: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 9*Lambda2 + 26*Lambda23 + 3*Lambda25
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda23,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda24,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda25,n,basis) == True:
                T = True
        elif CMfieldnumber == 3: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 3*Lambda2 + 4*Lambda23 + Lambda25
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda23,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda24,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda25,n,basis) == True:
                T = True
        elif CMfieldnumber == 33: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 5*Lambda2 + 11*Lambda23 + 2*Lambda25
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda23,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda24,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(Lambda25,n,basis) == True:
                T = True
        elif CMfieldnumber == 4: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda22 + I), B5 = 1/2*I*(Lambda24 + Lambda22 + Lambda2 + I), B6 = 1/2*I*(Lambda25 + Lambda2 + I). sqrt(-7) = 7*Lambda2 + 6*Lambda23 + Lambda25 
            B4 = 1/2*I*(Lambda23 + Lambda22 + I)
            B5 = 1/2*I*(Lambda24 + Lambda22 + Lambda2 + I)
            B6 = 1/2*I*(Lambda25 + Lambda2 + I)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 5: #basis: 1, Lambda2, Lambda22, B4 = 1/22*I*(Lambda23 + 21*Lambda2 + 11*I), B5 = 1/110*I*(Lambda24 + 65*Lambda22 + 55*Lambda2 + 66*I), B6 = 1/110*I*(Lambda25 + 55*Lambda22 + 21*Lambda2 + 55*I). sqrt(-7) = 1/11*I*(21*Lambda2 + Lambda23)
            B4 = 1/22*I*(Lambda23 + 21*Lambda2 + 11*I)
            B5 = 1/110*I*(Lambda24 + 65*Lambda22 + 55*Lambda2 + 66*I)
            B6 = 1/110*I*(Lambda25 + 55*Lambda22 + 21*Lambda2 + 55*I)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 55: #basis: 1, Lambda2, Lambda22, B4 = 1/6*I*(Lambda23 + 3*Lambda22 + 3*I), B5 = 1/30*I*(Lambda24 + 15*Lambda22 + 15*Lambda2 + 21*I), B6 = 1/30*I*(Lambda25 + 21*Lambda2 + 15*I). sqrt(-7) = 1/3*I*(63*Lambda2 + 62*Lambda23 + Lambda25)
            B4 = 1/6*I*(Lambda23 + 3*Lambda22 + 3*I)
            B5 = 1/30*I*(Lambda24 + 15*Lambda22 + 15*Lambda2 + 21*I)
            B6 = 1/30*I*(Lambda25 + 21*Lambda2 + 15*I)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 6: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/44*I*(Lambda24 + 37*Lambda22 + 36*I), B6 = 1/176*I*(Lambda25 + 37*Lambda23 + 124*Lambda2) . i = 1/176*I*(300*Lambda2 + 37*Lambda23 + Lambda25)
            B4 = 1/2*I*(Lambda23 + Lambda2)
            B5 = 1/44*I*(Lambda24 + 37*Lambda22 + 36*I)
            B6 = 1/176*I*(Lambda25 + 37*Lambda23 + 124*Lambda2)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 66: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/8*I*(Lambda24 + 5*Lambda22 + 4*I), B6 = 1/88*I*(Lambda25 + 41*Lambda23 + 48*Lambda2). i = 1/44*I*(70*Lambda2 + 19*Lambda23 + Lambda25)
            B4 = 1/2*I*(Lambda23 + Lambda2)
            B5 = 1/8*I*(Lambda24 + 5*Lambda22 + 4*I)
            B6 = 1/88*I*(Lambda25 + 41*Lambda23 + 48*Lambda2)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 666: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/4*I*(Lambda24 + Lambda22), B6 = 1/64*I*(Lambda25 + 5*Lambda23 + 4*Lambda2). i = 1/64*I*(788*Lambda2 + 153*Lambda23 + 5*Lambda25)
            B4 = 1/2*I*(Lambda23 + Lambda2)
            B5 = 1/4*I*(Lambda24 + Lambda22)
            B6 = 1/64*I*(Lambda25 + 5*Lambda23 + 4*Lambda2)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 7: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/4*I*(Lambda24 + Lambda22), B6 = 1/16*I*(Lambda25 + 5*Lambda23 + 4*Lambda2). i = 1/16*I*(28*Lambda2 + 13*Lambda23 + Lambda25)
            B4 = 1/2*I*(Lambda23 + Lambda2)
            B5 = 1/4*I*(Lambda24 + Lambda22)
            B6 = 1/16*I*(Lambda25 + 5*Lambda23 + 4*Lambda2)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 8: #basis: 1, Lambda2, Lambda22, B4 = 1/28*I*(Lambda23 + 21*Lambda2), B5 = 1/56*I*(Lambda24 + 49*Lambda22), B6 = 1/56*I*(Lambda25 + Lambda23). i = 1/28*I*(21*Lambda2 + Lambda23)
            B4 = 1/28*I*(Lambda23 + 21*Lambda2)
            B5 = 1/56*I*(Lambda24 + 49*Lambda22)
            B6 = 1/56*I*(Lambda25 + Lambda23)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 88: #basis: 1, Lambda2, Lambda22, B4 = 1/6*I*(Lambda23 + 3*Lambda2), B5 = 1/12*I*(Lambda24 + 9*Lambda22), B6 = 1/96*I*(Lambda25 + 13*Lambda23 + 36*Lambda2). i = 1/96*I*(180*Lambda2 + 29*Lambda23 + Lambda25)
            B4 = 1/6*I*(Lambda23 + 3*Lambda2)
            B5 = 1/12*I*(Lambda24 + 9*Lambda22)
            B6 = 1/96*I*(Lambda25 + 13*Lambda23 + 36*Lambda2)
            if matrix_with_entries_in_R_and_ROVERn(Lambda22,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 9: #basis: 1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/20*I*(Lambda24 + 16*I), B6 = 1/20*I*(Lambda25 + 16*Lambda2). sqrt(-2) = 1/20*I*(76*Lambda2 + 20*Lambda23 + Lambda25)
            B3 = 1/2*Lambda22
            B4 = 1/2*Lambda23
            B5 = 1/20*I*(Lambda24 + 16*I)
            B6 = 1/20*I*(Lambda25 + 16*Lambda2)
            if matrix_with_entries_in_R_and_ROVERn(B3,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        elif CMfieldnumber == 10: #basis: 1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/14*I*(Lambda23 + Lambda22 + 7*I), B5 = 1/490*I*(Lambda24 + 35*Lambda22 + 245*Lambda2 + 441*I), B6 = 1/490*I*(Lambda25 + 441*Lambda2 + 245*I). sqrt(-7) = 1/245*I*(931*Lambda2 + 70*Lambda23 + Lambda25)
            B3 = 1/7*Lambda22
            B4 = 1/14*I*(Lambda23 + Lambda22 + 7*I)
            B5 = 1/490*I*(Lambda24 + 35*Lambda22 + 245*Lambda2 + 441*I)
            B6 = 1/490*I*(Lambda25 + 441*Lambda2 + 245*I)
            if matrix_with_entries_in_R_and_ROVERn(B3,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
                T = True
        if T == True:
            x4 = QA(Lambda2[1][0])
            x5 = QA(Lambda2[1][1])
            x6 = QA(Lambda2[1][2])
            x7 = QA(Lambda2[2][0])
            x8 = QA(Lambda2[2][1])
            x9 = QA(Lambda2[2][2])
            Nx1 = QA(x1).reduced_norm()
            Nx2 = QA(x2).reduced_norm()
            Nx3 = QA(x3).reduced_norm()
            Trx1x2 = QA(x1*x2).reduced_trace()
            Trx1x3 = QA(x1*x3).reduced_trace()
            Trx2x3 = QA(x2*x3).reduced_trace()
            if x == Nx1 + gOVERn*Nx2 + aOVERn*Nx3 + bOVERn*Trx2x3 and a == k3*Nx2 + k1*Nx3 - k2*Trx2x3 - Trx1x2 and b == k4*Nx2 - k5*Nx3 - k6*Trx2x3 - Trx1x3:
                Dx1x2 = 4*Nx1*Nx2 - Trx1x2^2
                Dx1x3 = 4*Nx1*Nx3 - Trx1x3^2
                Dx2x3 = 4*Nx2*Nx3 - Trx2x3^2
                Dx1x2x3 = 4*Nx1*Nx2*Nx3 - Nx1*Trx2x3^2 - Nx2*Trx1x3^2 - Nx3*Trx1x2^2 - Trx1x2*Trx1x3*Trx2x3 # = (Trx1x2x3)^2
                if Dx1x2x3.is_square() and Dx1x2x3/p in ZZ and Dx1x2/p in ZZ and Dx1x3/p in ZZ and Dx2x3/p in ZZ:
                    if Dx1x2 !=0 or Dx1x3 !=0 or Dx2x3!=0 or (Dx1x2 == 0 and Dx1x3 == 0 and Dx2x3 == 0 and n/p in ZZ):
                        print ('[x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3] =', [x, dOVERn, a, cOVERn, b, gamma, n, QA(x1),QA(x2),QA(x3), Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
                        csvwriter2.writerow([x, dOVERn, a, cOVERn, b, gamma, n, QA(x1),QA(x2),QA(x3), Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
                        csvfile2.flush()
                        OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, QA(x1),QA(x2),QA(x3), Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
                        sys.stdout.flush()
print ("done")
print (len(OUTPUT)," solutions found")
print ("----------------------------------------------------------------------")
sys.stdout.flush()


csvfile2.close()

print ("All done in ", time()- T)



