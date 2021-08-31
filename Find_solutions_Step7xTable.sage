# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step performs a check before Step 7 of Algorithm 1 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step7xTable.sage CMfieldnumber p (q)
# e.g. sage ./Find_solutions_Step7x.sage 34 11 runs for CM field 34, p = 11, all max orders in the table BasesForPrimesTable.csv
# If p = 1 mod 8, then choose a prime q = 3 mod 4 and q not a square mod p. Ensure you choose q such that the bases appear in BasesForPrimesTable.csv.
# e.g. sage ./Find_solutions_Step7xTable.sage 34 17 3 runs for CM field 34, p = 17, q = 3, all max orders in the table BasesForPrimesTable.csv
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 7: 
# Read solutions [x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,x4,x5,x6,x7,x8,x9,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3] from Step 6
# Check if the solution lifts to a solution in the maximal order.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# For each max order, a file CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"+O+Solutions_for_di_Step7.csv" with the xi's. 
# Example: CMfield34p11MaxOrderBasis_1_+i_+frac{1}{2}i+frac{1}{2}k_frac{1}{2}+frac{1}{2}j_Solutions_for_di_Step7.csv 
################################################################################
# Explanation on how to find the basis (Magma):
# Choose A,B,C below such that sextic is of the form x^6+A*x^4+B*x^2+C.
# A := 63; B := 126; C := 63; //Modify accordingly. 
# R<x>:=PolynomialRing(Rationals());
# f:=x^6+A*x^4+B*x^2+C;
# K<Lambda2>:=NumberField(f);
# OK := MaximalOrder(K); O:= EquationOrder(K);
# Basis(OK, NumberField(OK));

################################################################################
from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv
################################################################################
filelocation = "./"

# Set-up the quaternion algebra QA
CMfieldnumber = sage_eval(sys.argv[1])
p = sage_eval(sys.argv[2])

S = Integers(8)

if p == 2:
	QA = QuaternionAlgebra(SR, -1,-1, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
	q = 1
else:
	if S(p) == 3 or S(p) == 7:
		QA = QuaternionAlgebra(SR, -1,-p, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
		q = 1
	elif S(p) == 5:
		QA = QuaternionAlgebra(SR, -2,-p, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
		q = 2
	elif S(p) == 1:
		q = sage_eval(sys.argv[3]) # Choose a prime q = 3 mod 4 and q not a square mod p. 
		QA = QuaternionAlgebra(SR, -q,-p, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)


print ("CMfieldnumber = ",CMfieldnumber,", p = ",p,", q = ",q)
print (sys.version)
print ("----------------------------------------------------------------------")
sys.stdout.flush()


# Read max order basis from Table
print ("reading max order basis...") #
sys.stdout.flush()
tread = time() 
filename4 = "BasesForPrimesTable.csv"
reader =csv.DictReader(open(filelocation+filename4,'r',encoding='utf-8-sig'))
print ("Reading time = ", time()-tread)
all_bases = []
for row in reader:
	pp = QQ(row["p"])
	if S(p) != 1 and pp == p:
		a00 = QQ(row["a00"])
		a01 = QQ(row["a01"])
		a02 = QQ(row["a02"])
		a03 = QQ(row["a03"])
		a10 = QQ(row["a10"])
		a11 = QQ(row["a11"])
		a12 = QQ(row["a12"])
		a13 = QQ(row["a13"])
		a20 = QQ(row["a20"])
		a21 = QQ(row["a21"])
		a22 = QQ(row["a22"])
		a23 = QQ(row["a23"])
		a30 = QQ(row["a30"])
		a31 = QQ(row["a31"])
		a32 = QQ(row["a32"])
		a33 = QQ(row["a33"])
		all_bases += [[a00,a01,a02,a03, a10,a11,a12,a13, a20,a21,a22,a23, a30,a31,a32,a33]]
	elif S(p) == 1 and pp == p and q == QQ(row["q"]):
		a00 = QQ(row["a00"])
		a01 = QQ(row["a01"])
		a02 = QQ(row["a02"])
		a03 = QQ(row["a03"])
		a10 = QQ(row["a10"])
		a11 = QQ(row["a11"])
		a12 = QQ(row["a12"])
		a13 = QQ(row["a13"])
		a20 = QQ(row["a20"])
		a21 = QQ(row["a21"])
		a22 = QQ(row["a22"])
		a23 = QQ(row["a23"])
		a30 = QQ(row["a30"])
		a31 = QQ(row["a31"])
		a32 = QQ(row["a32"])
		a33 = QQ(row["a33"])
		all_bases += [[a00,a01,a02,a03, a10,a11,a12,a13, a20,a21,a22,a23, a30,a31,a32,a33]]


load('CM_Fields.sage')

TT = time()

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


print ("There are ", len(all_bases), " maximal orders in ", QA)
print (" ")

for base in all_bases:
	a00 = base[0]
	a01 = base[1]
	a02 = base[2]
	a03 = base[3]
	a10 = base[4]
	a11 = base[5]
	a12 = base[6]
	a13 = base[7]
	a20 = base[8]
	a21 = base[9]
	a22 = base[10]
	a23 = base[11]
	a30 = base[12]
	a31 = base[13]
	a32 = base[14]
	a33 = base[15]
	Mlist = [a00,a01,a02,a03, a10,a11,a12,a13, a20,a21,a22,a23, a30,a31,a32,a33]
	M = matrix(QQ,4,Mlist)
	O = [ QA(M[0][0]+M[0][1]*i+M[0][2]*j+M[0][3]*k), QA(M[1][0]+M[1][1]*i+M[1][2]*j+M[1][3]*k), QA(M[2][0]+M[2][1]*i+M[2][2]*j+M[2][3]*k), QA(M[3][0]+M[3][1]*i+M[3][2]*j+M[3][3]*k)]
	basis = M.transpose().inverse()
	print (" ")
	print ("Maximal order with basis ", O, " in ", QA, ".")
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
	filename1 = filename1+"Solutions_for_xi_Step6.csv"
	reader1 = csv.DictReader(open(filelocation + filename1, 'r'))
	L = []
	for row in reader1: #x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3
		L += [[QQ(row["x"]), QQ(row["dOVERn"]), QQ(row["a"]), QQ(row["cOVERn"]), QQ(row["b"]), QQ(row["gamma"]), QQ(row["n"]), eval(preparse(row["x1"])), eval(preparse(row["x2"])), eval(preparse(row["x3"])), eval(preparse(row["Dx1x2"])), eval(preparse(row["Dx1x3"])), eval(preparse(row["Dx2x3"])), eval(preparse(row["Dx1x2x3"])) ]]	
	LengthAllSol = len(L)
	print ("Number of trivial and non-trivial solutions = ", LengthAllSol)
	I = matrix(QQ,3,3, [1,0,0, 0,1,0, 0,0,1])
	#write all remaining solutions in a csv file
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
	filename2 = filename2+"Solutions_for_xi_Step7.csv"
	csvfile2 = open(filelocation + filename2, 'w')
	csvwriter2 = csv.writer(csvfile2) #x,dOVERn,a,cOVERn,b,gamma,n,x1,..., x9,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3
	csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'Dx1x2', 'Dx1x3', 'Dx2x3', 'Dx1x2x3'))
	OUTPUT = []
	for LL in L: #x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3
		x = QQ(LL[0])
		dOVERn = QQ(LL[1])
		a = QQ(LL[2])
		cOVERn = QQ(LL[3])
		b = QQ(LL[4])
		gamma = QQ(LL[5])
		n = QQ(LL[6])
		x1 = QA(LL[7])
		x2 = QA(LL[8])
		x3 = QA(LL[9])
		Dx1x2 = QQ(LL[10])
		Dx1x3 = QQ(LL[11])
		Dx2x3 = QQ(LL[12])
		Dx1x2x3 = QQ(LL[13])
		gOVERn = QQ(gamma/n)
		aOVERn = QQ(a/n)
		bOVERn = QQ(b/n)
		k1 = QQ(aOVERn*dOVERn - aOVERn*x - bOVERn)
		k2 = QQ(aOVERn*cOVERn + bOVERn*x)
		k3 = QQ(k1*cOVERn - k2*dOVERn)
		k4 = QQ(k2*cOVERn - bOVERn*b)
		k5 = QQ(aOVERn*a - k1*dOVERn - k2)
		k6 = QQ(aOVERn*b + k1*cOVERn)
		x4 = QA(gOVERn*x2 - bOVERn*x3)
		x5 = QA(x1 + k3*x2 + k2*x3)
		x6 = QA(k4*x2 + k6*x3)
		x7 = QA(-bOVERn*x2 + aOVERn*x3)
		x8 = QA(k2*x2 + k1*x3)
		x9 = QA(x1 + k6*x2 - k5*x3)
		Lambda2 = matrix(QA, 3,3, [x1,x2,x3, x4,x5,x6, x7,x8,x9]) 
		Zero = matrix(QA, 3,3, [0,0,0, 0,0,0, 0,0,0])  # Lambda2^6 + A*Lambda2^4 + B*Lambda2^2 + C*II = O
		if matrix_with_entries_in_R_and_ROVERn(Lambda2,n,basis) == True:
			T = False
			Lambda22 = Lambda2^2
			Lambda23 = Lambda2^3
			Lambda24 = Lambda2^4
			Lambda25 = Lambda2^5
			Lambda26 = Lambda2^6
			if CMfieldnumber == 1: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, 1/7*(Lambda25 + 6*Lambda23 + Lambda2). i = 1/7*I*(Lambda2 + 6*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = 1/7*I*(Lambda2 + 6*Lambda23 + Lambda25)
			elif CMfieldnumber == 2: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 3*Lambda2 + Lambda23
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = Lambda25
			elif CMfieldnumber == 22: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 9*Lambda2 + 26*Lambda23 + 3*Lambda25
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = Lambda25
			elif CMfieldnumber == 3: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 3*Lambda2 + 4*Lambda23 + Lambda25
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = Lambda25
			elif CMfieldnumber == 33: #basis: 1, Lambda2, Lambda22, Lambda23, Lambda24, Lambda25. i = 5*Lambda2 + 11*Lambda23 + 2*Lambda25
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = Lambda25
			elif CMfieldnumber == 4: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda22 + I), B5 = 1/2*I*(Lambda24 + Lambda22 + Lambda2 + I), B6 = 1/2*I*(Lambda25 + Lambda2 + I). sqrt(-7) = 7*Lambda2 + 6*Lambda23 + Lambda25 
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*I*(Lambda23 + Lambda22 + I)
				B5 = 1/2*I*(Lambda24 + Lambda22 + Lambda2 + I)
				B6 = 1/2*I*(Lambda25 + Lambda2 + I)
			elif CMfieldnumber == 5: #basis: 1, Lambda2, Lambda22, B4 = 1/22*I*(Lambda23 + 21*Lambda2 + 11*I), B5 = 1/110*I*(Lambda24 + 65*Lambda22 + 55*Lambda2 + 66*I), B6 = 1/110*I*(Lambda25 + 55*Lambda22 + 21*Lambda2 + 55*I). sqrt(-7) = 1/11*I*(21*Lambda2 + Lambda23)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/22*I*(Lambda23 + 21*Lambda2 + 11*I)
				B5 = 1/110*I*(Lambda24 + 65*Lambda22 + 55*Lambda2 + 66*I)
				B6 = 1/110*I*(Lambda25 + 55*Lambda22 + 21*Lambda2 + 55*I)
			elif CMfieldnumber == 55 or CMfieldnumber == 5555: #basis: 1, Lambda2, Lambda22, B4 = 1/6*I*(Lambda23 + 3*Lambda22 + 3*I), B5 = 1/30*I*(Lambda24 + 15*Lambda22 + 15*Lambda2 + 21*I), B6 = 1/30*I*(Lambda25 + 21*Lambda2 + 15*I). sqrt(-7) = 1/3*I*(63*Lambda2 + 62*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/6*I*(Lambda23 + 3*Lambda22 + 3*I)
				B5 = 1/30*I*(Lambda24 + 15*Lambda22 + 15*Lambda2 + 21*I)
				B6 = 1/30*I*(Lambda25 + 21*Lambda2 + 15*I)
			elif CMfieldnumber == 6: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/44*I*(Lambda24 + 37*Lambda22 + 36*I), B6 = 1/176*I*(Lambda25 + 37*Lambda23 + 124*Lambda2) . i = 1/176*I*(300*Lambda2 + 37*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*I*(Lambda23 + Lambda2)
				B5 = 1/44*I*(Lambda24 + 37*Lambda22 + 36*I)
				B6 = 1/176*I*(Lambda25 + 37*Lambda23 + 124*Lambda2)
			elif CMfieldnumber == 66: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/8*I*(Lambda24 + 5*Lambda22 + 4*I), B6 = 1/88*I*(Lambda25 + 41*Lambda23 + 48*Lambda2). i = 1/44*I*(70*Lambda2 + 19*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*I*(Lambda23 + Lambda2)
				B5 = 1/8*I*(Lambda24 + 5*Lambda22 + 4*I)
				B6 = 1/88*I*(Lambda25 + 41*Lambda23 + 48*Lambda2)
			elif CMfieldnumber == 666 or CMfieldnumber == 6666 or CMfieldnumber == 66666: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/4*I*(Lambda24 + Lambda22), B6 = 1/64*I*(Lambda25 + 5*Lambda23 + 4*Lambda2). i = 1/64*I*(788*Lambda2 + 153*Lambda23 + 5*Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*I*(Lambda23 + Lambda2)
				B5 = 1/4*I*(Lambda24 + Lambda22)
				B6 = 1/64*I*(Lambda25 + 5*Lambda23 + 4*Lambda2)
			elif CMfieldnumber == 7 or CMfieldnumber == 77: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/4*I*(Lambda24 + Lambda22), B6 = 1/16*I*(Lambda25 + 5*Lambda23 + 4*Lambda2). i = 1/16*I*(28*Lambda2 + 13*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*I*(Lambda23 + Lambda2)
				B5 = 1/4*I*(Lambda24 + Lambda22)
				B6 = 1/16*I*(Lambda25 + 5*Lambda23 + 4*Lambda2)
			elif CMfieldnumber == 8: #basis: 1, Lambda2, Lambda22, B4 = 1/28*I*(Lambda23 + 21*Lambda2), B5 = 1/56*I*(Lambda24 + 49*Lambda22), B6 = 1/56*I*(Lambda25 + Lambda23). i = 1/28*I*(21*Lambda2 + Lambda23)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/28*I*(Lambda23 + 21*Lambda2)
				B5 = 1/56*I*(Lambda24 + 49*Lambda22)
				B6 = 1/56*I*(Lambda25 + Lambda23)
			elif CMfieldnumber == 88: #basis: 1, Lambda2, Lambda22, B4 = 1/6*I*(Lambda23 + 3*Lambda2), B5 = 1/12*I*(Lambda24 + 9*Lambda22), B6 = 1/96*I*(Lambda25 + 13*Lambda23 + 36*Lambda2). i = 1/96*I*(180*Lambda2 + 29*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/6*I*(Lambda23 + 3*Lambda2)
				B5 = 1/12*I*(Lambda24 + 9*Lambda22)
				B6 = 1/96*I*(Lambda25 + 13*Lambda23 + 36*Lambda2)
			elif CMfieldnumber == 888 or CMfieldnumber == 8888 or CMfieldnumber == 88888: #basis: 1, Lambda2, Lambda22, B4 = 1/2*I*(Lambda23 + Lambda2), B5 = 1/124*I*(Lambda24 + 89*Lambda22 + 100*I), B6 = 1/496*I*(Lambda25 + 213*Lambda23 + 348*Lambda2). i = 1/496*I*(1044*Lambda2 + 143*Lambda23 + 3*Lambda25)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*I*(Lambda23 + Lambda2)
				B5 = 1/124*I*(Lambda24 + 89*Lambda22 + 100*I)
				B6 = 1/496*I*(Lambda25 + 213*Lambda23 + 348*Lambda2)
			elif CMfieldnumber == 9: #basis: 1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/20*I*(Lambda24 + 16*I), B6 = 1/20*I*(Lambda25 + 16*Lambda2). sqrt(-2) = 1/20*I*(76*Lambda2 + 20*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/20*I*(Lambda24 + 16*I)
				B6 = 1/20*I*(Lambda25 + 16*Lambda2)
			elif CMfieldnumber == 99 or CMfieldnumber == 999: #basis: 1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/4*Lambda24, B6 = 1/20*I*(Lambda25 + 16*Lambda2). sqrt(-2) = 1/10*I*(46*Lambda2 + 15*Lambda23 + Lambda25) 
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/4*Lambda24
				B6 = 1/20*I*(Lambda25 + 16*Lambda2)
			elif CMfieldnumber == 10: #basis: 1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/14*I*(Lambda23 + Lambda22 + 7*I), B5 = 1/490*I*(Lambda24 + 35*Lambda22 + 245*Lambda2 + 441*I), B6 = 1/490*I*(Lambda25 + 441*Lambda2 + 245*I). sqrt(-7) = 1/245*I*(931*Lambda2 + 70*Lambda23 + Lambda25)
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/14*I*(Lambda23 + Lambda22 + 7*I)
				B5 = 1/490*I*(Lambda24 + 35*Lambda22 + 245*Lambda2 + 441*I)
				B6 = 1/490*I*(Lambda25 + 441*Lambda2 + 245*I)
			elif CMfieldnumber == 1010 or CMfieldnumber == 1011: #basis: 1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/14*I*(Lambda23 + 7*Lambda2 + 7*I), B5 = 1/98*I*(Lambda24 + 7*Lambda22 + 49*Lambda2*I), B6 = 1/490*I*(Lambda25 + 35*Lambda22 + 441*Lambda2 + 245*I). sqrt(-7) = 1/245*I*(1127*Lambda2 + 105*Lambda23 + 2*Lambda25)   
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/14*I*(Lambda23 + 7*Lambda2 + 7*I)
				B5 = 1/98*I*(Lambda24 + 7*Lambda22 + 49*Lambda2*I)
				B6 = 1/490*I*(Lambda25 + 35*Lambda22 + 441*Lambda2 + 245*I)
			elif CMfieldnumber == 23: #basis:  1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/42*(Lambda23 + 7*Lambda2 + 21), B5 = 1/127890*(Lambda24 + 6895*Lambda22 + 63945*Lambda2 + 117306), B6 = 1/8057070*(Lambda25 + 79975*Lambda23 + 575505*Lambda22 + 2547216*Lambda2). sqrt(-7) = (8*Lambda25 + 9485*Lambda23 + 1002393*Lambda2)/575505
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/42*(Lambda23 + 7*Lambda2 + 21*I)
				B5 = 1/127890*(Lambda24 + 6895*Lambda22 + 63945*Lambda2 + 117306*I)
				B6 = 1/8057070*(Lambda25 + 79975*Lambda23 + 575505*Lambda22 + 2547216*Lambda2)
			elif CMfieldnumber == 221 or CMfieldnumber == 2211: #basis:  1, Lambda2, 1/7*Lambda2^2, 1/14*(Lambda2^3 + 7*Lambda2 + 7), 1/98*(Lambda2^4 + 7*Lambda2^2 + 49*Lambda2), 1/98*(Lambda2^5 + 7*Lambda2^2 + 49*Lambda2 + 49)
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/14*(Lambda23 + 7*Lambda2 + 7*I)
				B5 = 1/98*(Lambda24 + 7*Lambda22 + 49*Lambda2)
				B6 = 1/98*(Lambda25 + 7*Lambda22 + 49*Lambda2 + 49*I)
			elif CMfieldnumber == 233 or CMfieldnumber == 2333: #basis:  1, Lambda2, B3 = 1/3*(Lambda22 + 1), B4 = 1/30*(Lambda23 + 5*Lambda22 + 22*Lambda2 + 5), B5 = 1/450*(Lambda24 + 47*Lambda22 + 225*Lambda2 + 325), B6 = 1/2250*(Lambda25 + 47*Lambda23 + 375*Lambda22 + 1675*Lambda2 + 1500). sqrt(-7) = (Lambda25 + 122*Lambda23 + 3325*Lambda2)/1125
				B2 = Lambda2
				B3 = 1/3*(Lambda22 + I)
				B4 = 1/30*(Lambda23 + 5*Lambda22 + 22*Lambda2 + 5*I)
				B5 = 1/450*(Lambda24 + 47*Lambda22 + 225*Lambda2 + 325*I)
				B6 = 1/2250*(Lambda25 + 47*Lambda23 + 375*Lambda22 + 1675*Lambda2 + 1500*I)
			elif CMfieldnumber == 25: #basis:  1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/42*(Lambda23 + 3*Lambda22 + 28*Lambda2 + 21), B5 = 1/16758*(Lambda24 + 1813*Lambda22 + 8379*Lambda2 + 15435), B6 = 1/1055754*(Lambda25 + 7798*Lambda23 + 1037673*Lambda2 + 527877). sqrt(-7) = (Lambda25 + 616*Lambda23 + 82467*Lambda2)/75411
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/42*(Lambda23 + 3*Lambda22 + 28*Lambda2 + 21*I)
				B5 = 1/16758*(Lambda24 + 1813*Lambda22 + 8379*Lambda2 + 15435*I)
				B6 = 1/1055754*(Lambda25 + 7798*Lambda23 + 1037673*Lambda2 + 527877*I)
			elif CMfieldnumber == 255 or CMfieldnumber == 255555: #basis:  1, Lambda2, B3 = 1/3*(Lambda22 + 1), B4 = 1/6*(Lambda23 + Lambda22 + 4*Lambda2 + 1), B5 = 1/54*(Lambda24 + 5*Lambda22 + 27*Lambda2 + 31), B6 = 1/1026*(Lambda25 + 158*Lambda23 + 157*Lambda2 + 513). sqrt(-7) = (Lambda25 + 44*Lambda23 + 385*Lambda2)/171
				B2 = Lambda2
				B3 = 1/3*(Lambda22 + I)
				B4 = 1/6*(Lambda23 + Lambda22 + 4*Lambda2 + I)
				B5 = 1/54*(Lambda24 + 5*Lambda22 + 27*Lambda2 + 31*I)
				B6 = 1/1026*(Lambda25 + 158*Lambda23 + 157*Lambda2 + 513*I)
			elif CMfieldnumber == 2555: #basis:  1, Lambda2, B3 = Lambda22, B4 = 1/6*(Lambda23 + Lambda2 + 3), B5 = 1/594*(Lambda24 + 115*Lambda22 + 297*Lambda2 + 522), B6 = 1/1782*(Lambda25 + 115*Lambda23 + 891*Lambda22 + 522*Lambda2). sqrt(-7) = (2*Lambda25 + 131*Lambda23 + 945*Lambda2)/297
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/6*(Lambda23 + Lambda2 + 3*I)
				B5 = 1/594*(Lambda24 + 115*Lambda22 + 297*Lambda2 + 522*I)
				B6 = 1/1782*(Lambda25 + 115*Lambda23 + 891*Lambda22 + 522*Lambda2)
			elif CMfieldnumber == 25555: #basis:  1, Lambda2, B3 = Lambda22, B4 = 1/18*(Lambda23 + 9*Lambda22 + 10*Lambda2 + 9), B5 = 1/2682*(Lambda24 + 2233*Lambda22 + 1341*Lambda2 + 315), B6 = 1/2682*(Lambda25 + 147*Lambda23 + 1341*Lambda22 + 911*Lambda2). sqrt(-7) = (5*Lambda25 + 586*Lambda23 + 5747*Lambda2)/1341
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/18*(Lambda23 + 9*Lambda22 + 10*Lambda2 + 9*I)
				B5 = 1/2682*(Lambda24 + 2233*Lambda22 + 1341*Lambda2 + 315*I)
				B6 = 1/2682*(Lambda25 + 147*Lambda23 + 1341*Lambda22 + 911*Lambda2)
			elif CMfieldnumber == 26: #basis:  1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/42*(Lambda23 + 7*Lambda2 + 21), B5 = 1/53802*(Lambda24 + 5509*Lambda22 + 26901*Lambda2 + 48510), B6 = 1/3389526*(Lambda25 + 66997*Lambda23 + 242109*Lambda22 + 909342*Lambda2). sqrt(-7) = (4*Lambda25 + 2821*Lambda23 + 328545*Lambda2)/242109
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/42*(Lambda23 + 7*Lambda2 + 21*I)
				B5 = 1/53802*(Lambda24 + 5509*Lambda22 + 26901*Lambda2 + 48510*I)
				B6 = 1/3389526*(Lambda25 + 66997*Lambda23 + 242109*Lambda22 + 909342*Lambda2)
			elif CMfieldnumber == 266: #basis:  1, Lambda2, B3 = Lambda22, B4 = 1/6*(Lambda23 + 3*Lambda22 + 4*Lambda2 + 3), B5 = 1/414*(Lambda24 + 325*Lambda22 + 207*Lambda2 + 225), B6 = 1/3726*(Lambda25 + 118*Lambda23 + 2709*Lambda2 + 1863). sqrt(-7) = (Lambda25 + 118*Lambda23 + 2709*Lambda2)/1863
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/6*(Lambda23 + 3*Lambda22 + 4*Lambda2 + 3)
				B5 = 1/414*(Lambda24 + 325*Lambda22 + 207*Lambda2 + 225)
				B6 = 1/3726*(Lambda25 + 118*Lambda23 + 2709*Lambda2 + 1863)
			elif CMfieldnumber == 2666 or CMfieldnumber == 26666: #basis:  1, Lambda2, B3 = Lambda22, B4 = 1/6*(Lambda23 + Lambda2 + 3), B5 = 1/6*(Lambda24 + Lambda22 + 3*Lambda2), B6 = 1/138*(Lambda25 + 6*Lambda23 + 69*Lambda22 + 29*Lambda2 + 69). sqrt(-7) = (2*Lambda25 + 127*Lambda23 + 1967*Lambda2)/69
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/6*(Lambda23 + Lambda2 + 3*I)
				B5 = 1/6*(Lambda24 + Lambda22 + 3*Lambda2)
				B6 = 1/138*(Lambda25 + 6*Lambda23 + 69*Lambda22 + 29*Lambda2 + 69*I)
			elif CMfieldnumber == 27: #basis:  1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/14*(Lambda23 + 7*Lambda2 + 7), B5 = 1/1078*(Lambda24 + 63*Lambda22 + 539*Lambda2 + 980), B6 = 1/1078*(Lambda25 + 63*Lambda23 + 77*Lambda22 + 980*Lambda2). sqrt(-7) = (2*Lambda25 + 203*Lambda23 + 2499*Lambda2)/539
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/14*(Lambda23 + 7*Lambda2 + 7*I)
				B5 = 1/1078*(Lambda24 + 63*Lambda22 + 539*Lambda2 + 980*I)
				B6 = 1/1078*(Lambda25 + 63*Lambda23 + 77*Lambda22 + 980*Lambda2)
			elif CMfieldnumber == 277 or CMfieldnumber == 2777: #basis:  1, Lambda2, B3 = 1/7*Lambda22, B4 = 1/14*(Lambda23 + Lambda22 + 7), B5 = 1/98*(Lambda24 + 7*Lambda22 + 49*Lambda2 + 49), B6 = 1/686*(Lambda25 + 42*Lambda23 + 49*Lambda2 + 343). sqrt(-7) = (Lambda25 + 42*Lambda23 + 49*Lambda2)/343
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/14*(Lambda23 + Lambda22 + 7*I)
				B5 = 1/98*(Lambda24 + 7*Lambda22 + 49*Lambda2 + 49*I)
				B6 = 1/686*(Lambda25 + 42*Lambda23 + 49*Lambda2 + 343*I)
			elif CMfieldnumber == 28 or CMfieldnumber == 288 or CMfieldnumber == 2888 or CMfieldnumber == 28888: #basis: 1, Lambda2, 1/7*Lambda2^2, 1/42*(Lambda2^3 + 3*Lambda2^2 + 28*Lambda2 + 21), 1/882*(Lambda2^4 + 91*Lambda2^2 + 441*Lambda2 + 441), 1/7938*(Lambda2^5 + 154*Lambda2^3 + 5733*Lambda2 + 3969), sqrt(-7) = 1/3969*(Lambda2^5 + 154*Lambda2^3 - 2205*Lambda2)
				B2 = Lambda2
				B3 = 1/7*Lambda22
				B4 = 1/42*(Lambda23 + 3*Lambda22 + 28*Lambda2 + 21*I)
				B5 = 1/882*(Lambda24 + 91*Lambda22 + 441*Lambda2 + 441*I)
				B6 = 1/7938*(Lambda25 + 154*Lambda23 + 5733*Lambda2 + 3969*I)
			elif CMfieldnumber == 29: #basis:  1, B2 = 1/2*Lambda2, B3 = 1/8*Lambda22, B4 = 1/16*Lambda23, B5 = 1/320*(Lambda24 + 256), B6 = 1/640*(Lambda25 + 256*Lambda2). sqrt(-8) = (Lambda25 + 80*Lambda23 + 1216*Lambda2)/320
				B2 = 1/2*Lambda2
				B3 = 1/8*Lambda22
				B4 = 1/16*Lambda23
				B5 = 1/320*(Lambda24 + 256*I)
				B6 = 1/640*(Lambda25 + 256*Lambda2)
			elif CMfieldnumber == 299: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/20*(Lambda24 + 16), B6 = 1/20*(Lambda25 + 16*Lambda2). sqrt(-8) = (Lambda25 + 20*Lambda23 + 76*Lambda2)/10
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/20*(Lambda24 + 16*I)
				B6 = 1/20*(Lambda25 + 16*Lambda2)
			elif CMfieldnumber == 2999 or CMfieldnumber == 29999: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/4*Lambda24, B6 = 1/20*(Lambda25 + 16*Lambda2). sqrt(-8) = (Lambda25 + 15*Lambda23 + 46*Lambda2)/5
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/4*Lambda24
				B6 = 1/20*(Lambda25 + 16*Lambda2)
			elif CMfieldnumber == 210: #basis:  1, B2 = 1/2*Lambda2, B3 = 1/8*Lambda22, B4 = 1/16*Lambda23, B5 = 1/64*Lambda24, B6 = 1/128*Lambda25. sqrt(-8) = (Lambda25 + 32*Lambda23 + 192*Lambda2)/64
				B2 = 1/2*Lambda2
				B3 = 1/8*Lambda22
				B4 = 1/16*Lambda23
				B5 = 1/64*Lambda24
				B6 = 1/128*Lambda25
			elif CMfieldnumber == 2102 or CMfieldnumber == 21011: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/4*Lambda24, B6 = 1/4*Lambda25. sqrt(-8) = (Lambda25 + 8*Lambda23 + 12*Lambda2)/2
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/4*Lambda24
				B6 = 1/4*Lambda25
			elif CMfieldnumber == 2103: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/4*Lambda24, B6 = 1/4*Lambda25. sqrt(-8) = Lambda25 + 11*Lambda23 + 10*V
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/4*Lambda24
				B6 = 1/4*Lambda25
			elif CMfieldnumber == 2104: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/52*(Lambda24 + 16), B6 = 1/52*(Lambda25 + 16*Lambda2). sqrt(-8) = (2*Lambda25 + 39*Lambda23 + 110*Lambda2)/13
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/52*(Lambda24 + 16*I)
				B6 = 1/52*(Lambda25 + 16*Lambda2)
			elif CMfieldnumber == 2105: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/116*(Lambda24 + 38*Lambda22 + 96), B6 = 1/116*(Lambda25 + 38*Lambda23 + 96*Lambda2). sqrt(-8) = (5*Lambda25 + 132*Lambda23 + 596*Lambda2)/58
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/116*(Lambda24 + 38*Lambda22 + 96*I)
				B6 = 1/116*(Lambda25 + 38*Lambda23 + 96*Lambda2)
			elif CMfieldnumber == 2106: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/52*(Lambda24 + 40), B6 = 1/52*(Lambda25 + 40*Lambda2). sqrt(-8) = (3*Lambda25 + 104*Lambda23 + 172*Lambda2)/26
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/52*(Lambda24 + 40*I)
				B6 = 1/52*(Lambda25 + 40*Lambda2)
			elif CMfieldnumber == 2107: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/2*Lambda23, B5 = 1/4*Lambda24, B6 = 1/52*(Lambda25 + 8*Lambda23 + 12*Lambda2). sqrt(-8) = (3*Lambda25 + 76*Lambda23 + 452*Lambda2)/26
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/2*Lambda23
				B5 = 1/4*Lambda24
				B6 = 1/52*(Lambda25 + 8*Lambda23 + 12*Lambda2)
			elif CMfieldnumber == 2108: #basis:  1, Lambda2, B3 = 1/8*(Lambda22 + 2), B4 = 1/8*(Lambda23 + 2*Lambda2), B5 = 1/64*(Lambda24 + 4*Lambda22 + 4), B6 = 1/64*(Lambda25 + 4*Lambda23 + 4*Lambda2). sqrt(-8) = (Lambda25 + 36*Lambda23 + 292*Lambda2)/16
				B2 = Lambda2
				B3 = 1/8*(Lambda22 + 2)
				B4 = 1/8*(Lambda23 + 2*Lambda2)
				B5 = 1/64*(Lambda24 + 4*Lambda22 + 4)
				B6 = 1/64*(Lambda25 + 4*Lambda23 + 4*Lambda2)
			elif CMfieldnumber == 2109: #basis:  1, Lambda2, B3 = 1/2*Lambda22, B4 = 1/14*Lambda23, B5 = 1/364*(Lambda24 + 70*Lambda22 + 168), B6 = 1/364*(Lambda25 + 18*Lambda23 + 168*Lambda2). sqrt(-8) = (Lambda25 + 44*Lambda23 + 532*Lambda2)/182
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/14*Lambda23
				B5 = 1/364*(Lambda24 + 70*Lambda22 + 168*I)
				B6 = 1/364*(Lambda25 + 18*Lambda23 + 168*Lambda2)
			elif CMfieldnumber == 21010: #basis:  1, B2 = 1/2*Lambda2, B3 = 1/8*Lambda22, B4 = 1/16*Lambda23, B5 = 1/64*Lambda24, B6 = 1/128*Lambda25. sqrt(-8) = (Lambda25 + 44*Lambda23 + 160*Lambda2)/32, sqrt(-2) = 1/64*(-Lambda25 - 26*Lambda23 - 112*Lambda2)
				B2 = 1/2*Lambda2
				B3 = 1/8*Lambda22
				B4 = 1/16*Lambda23
				B5 = 1/64*Lambda24
				B6 = 1/128*Lambda25
			elif CMfieldnumber == 2112 or CMfieldnumber == 2114: #basis: 1, Lambda2, 1/2*Lambda2^2, 1/4*(Lambda2^3 + 2*Lambda2), 1/16*(Lambda2^4 + 2*Lambda2^2), 1/64*(Lambda2^5 + 10*Lambda2^3 + 16*Lambda2)
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/4*(Lambda23 + 2*Lambda2)
				B5 = 1/16*(Lambda24 + 2*Lambda22)
				B6 = 1/64*(Lambda25 + 10*Lambda23 + 16*Lambda2)
			elif CMfieldnumber == 212 or CMfieldnumber == 2122: #basis:  1, Lambda2, B3 = 1/11*Lambda22, B4 = 1/22*(Lambda23 + Lambda22 + 11), B5 = 1/242*(Lambda24 + 11*Lambda22 + 121*Lambda2 + 121), B6 = 1/242*(Lambda25 + 121*Lambda2 + 121). sqrt(-11) = (Lambda25 + 44*Lambda23 + 363*Lambda2)/121
				B2 = Lambda2
				B3 = 1/11*Lambda22
				B4 = 1/22*(Lambda23 + Lambda22 + 11*I)
				B5 = 1/242*(Lambda24 + 11*Lambda22 + 121*Lambda2 + 121*I)
				B6 = 1/242*(Lambda25 + 121*Lambda2 + 121*I)
			elif CMfieldnumber == 213 or CMfieldnumber == 2133: #basis:  1, Lambda2, 1/22*(Lambda2^2 + 11*Lambda2), 1/22*(Lambda2^3 + 11*Lambda2), 1/10648*(Lambda2^4 + 242*Lambda2^3 + 407*Lambda2^2 + 7986*Lambda2 + 9680), 1/42592*(Lambda2^5 + 407*Lambda2^3 + 36300*Lambda2 + 21296)
				B2 = Lambda2
				B3 = 1/22*(Lambda22 + 11*Lambda2)
				B4 = 1/22*(Lambda23 + 11*Lambda2)
				B5 = 1/10648*(Lambda24 + 242*Lambda23 + 407*Lambda22 + 7986*Lambda2 + 9680*I)
				B6 = 1/42592*(Lambda25 + 407*Lambda23 + 36300*Lambda2 + 21296*I)
			elif CMfieldnumber == 214 or CMfieldnumber == 2144: #basis:  1, Lambda2, 1/11*Lambda2^2, 1/22*(Lambda2^3 + 11*Lambda2 + 11), 1/2662*(Lambda2^4 + 99*Lambda2^2 + 1331*Lambda2 + 2420), 1/2662*(Lambda2^5 + 99*Lambda2^3 + 121*Lambda2^2 + 2420*Lambda2)
				B2 = Lambda2
				B3 = 1/11*Lambda22
				B4 = 1/22*(Lambda23 + 11*Lambda2 + 11*I)
				B5 = 1/2662*(Lambda24 + 99*Lambda22 + 1331*Lambda2 + 2420*I)
				B6 = 1/2662*(Lambda25 + 99*Lambda23 + 121*Lambda22 + 2420*Lambda2)
			elif CMfieldnumber == 215 or CMfieldnumber == 2155: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2 + 1), 1/14*(Lambda2^4 + 9*Lambda2^2 + 7*Lambda2 + 11), 1/14*(Lambda2^5 + 2*Lambda2^3 + 11*Lambda2 + 7)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22 + I)
				B5 = 1/14*(Lambda24 + 9*Lambda22 + 7*Lambda2 + 11*I)
				B6 = 1/14*(Lambda25 + 2*Lambda23 + 11*Lambda2 + 7*I)
			elif CMfieldnumber == 216 or CMfieldnumber == 2166: #basis:  1, Lambda2,  1/19*Lambda2^2, 1/38*(Lambda2^3 + 19*Lambda2 + 19), 1/722*(Lambda2^4 + 19*Lambda2^2 + 361*Lambda2), 1/722*(Lambda2^5 + 19*Lambda2^2 + 361*Lambda2 + 361)
				B2 = Lambda2
				B3 = 1/19*Lambda22
				B4 = 1/38*(Lambda23 + 19*Lambda2 + 19*I)
				B5 = 1/722*(Lambda24 + 19*Lambda22 + 361*Lambda2)
				B6 = 1/722*(Lambda25 + 19*Lambda22 + 361*Lambda2 + 361*I)
			elif CMfieldnumber == 217: #basis:  1, Lambda2, B3 = 1/38*(Lambda22 + 19*Lambda2), B4 = 1/38*(Lambda23 + 19*Lambda2), B5 = 1/401432*(Lambda24 + 5282*Lambda23 + 7391*Lambda22 + 301074*Lambda2 + 375440), B6 = 1/30508832*(Lambda25 + 387695*Lambda23 + 4590476*Lambda2 + 15254416). sqrt(-19) = (9*Lambda25 + 24263*Lambda23 + 2776812*Lambda2)/802864
				B2 = Lambda2
				B3 = 1/38*(Lambda22 + 19*Lambda2)
				B4 = 1/38*(Lambda23 + 19*Lambda2)
				B5 = 1/401432*(Lambda24 + 5282*Lambda23 + 7391*Lambda22 + 301074*Lambda2 + 375440*I)
				B6 = 1/30508832*(Lambda25 + 387695*Lambda23 + 4590476*Lambda2 + 15254416*I)
			elif CMfieldnumber == 2177 or CMfieldnumber == 2177777: #basis:  1, Lambda2, B3 = 1/2*(Lambda22 + Lambda2), B4 = 1/52*(Lambda23 + 5*Lambda2 + 26), B5 = 1/416*(Lambda24 + 31*Lambda22 + 208*Lambda2 + 260), B6 = 1/832*(Lambda25 + Lambda24 + 15*Lambda23 + 239*Lambda22 + 388*Lambda2 + 676). sqrt(-19) = (Lambda23 + 57*Lambda2)/26 = V*(-U + 57*I)/26
				B2 = Lambda2
				B3 = 1/2*(Lambda22 + Lambda2)
				B4 = 1/52*(Lambda23 + 5*Lambda2 + 26*I)
				B5 = 1/416*(Lambda24 + 31*Lambda22 + 208*Lambda2 + 260*I)
				B6 = 1/832*(Lambda25 + Lambda24 + 15*Lambda23 + 239*Lambda22 + 388*Lambda2 + 676*I)
			elif CMfieldnumber == 21777: #basis:  1, Lambda2, B3 = 1/2*(Lambda22 + Lambda2), B4 = 1/12*(Lambda23 + 9*Lambda2 + 6), B5 = 1/192*(Lambda24 + 63*Lambda22 + 132), B6 = 1/768*(Lambda25 + 2*Lambda24 + 63*Lambda23 + 126*Lambda22 + 708*Lambda2 + 648). sqrt(-19) = (Lambda25 + 167*Lambda23 + 684*Lambda2)/48
				B2 = Lambda2
				B3 = 1/2*(Lambda22 + Lambda2)
				B4 = 1/12*(Lambda23 + 9*Lambda2 + 6*I)
				B5 = 1/192*(Lambda24 + 63*Lambda22 + 132*I)
				B6 = 1/768*(Lambda25 + 2*Lambda24 + 63*Lambda23 + 126*Lambda22 + 708*Lambda2 + 648*I)
			elif CMfieldnumber == 217777: #basis:  1, Lambda2, B3 = 1/2*(Lambda22 + Lambda2), B4 = 1/16*(Lambda23 + 3*Lambda2 + 8), B5 = 1/1024*(Lambda24 + 491*Lambda22 + 512*Lambda2 + 64), B6 = 1/16384*(Lambda25 + 8*Lambda24 + 1003*Lambda23 + 8024*Lambda22 + 6208*Lambda2 + 512). sqrt(-19) = (3*Lambda25 + 449*Lambda23 + 10944*Lambda2)/4096
				B2 = Lambda2
				B3 = 1/2*(Lambda22 + Lambda2)
				B4 = 1/16*(Lambda23 + 3*Lambda2 + 8*I)
				B5 = 1/1024*(Lambda24 + 491*Lambda22 + 512*Lambda2 + 64*I)
				B6 = 1/16384*(Lambda25 + 8*Lambda24 + 1003*Lambda23 + 8024*Lambda22 + 6208*Lambda2 + 512*I)
			elif CMfieldnumber == 218 or CMfieldnumber == 2188 or CMfieldnumber == 21888 or CMfieldnumber == 218888 or CMfieldnumber == 2188888: #basis:  1, Lambda2, 1/38*(Lambda2^2 + 19*Lambda2), 1/114*(Lambda2^3 + 19*Lambda2), 1/1065672*(Lambda2^4 + 4674*Lambda2^3 + 6859*Lambda2^2 + 88806*Lambda2 + 545832), 1/728919648*(Lambda2^5 + 2166247*Lambda2^3 + 47968236*Lambda2 + 364459824)
				B2 = Lambda2
				B3 = 1/38*(Lambda22 + 19*Lambda2)
				B4 = 1/114*(Lambda23 + 19*Lambda2)
				B5 = 1/1065672*(Lambda24 + 4674*Lambda23 + 6859*Lambda22 + 88806*Lambda2 + 545832*I)
				B6 = 1/728919648*(Lambda25 + 2166247*Lambda23 + 47968236*Lambda2 + 364459824*I)
			elif CMfieldnumber == 219 or CMfieldnumber == 2199: #basis:  1, Lambda2, 1/86*(Lambda2^2 + 43*Lambda2), 1/86*(Lambda2^3 + 43*Lambda2), 1/162712*(Lambda2^4 + 946*Lambda2^3 + 1591*Lambda2^2 + 122034*Lambda2 + 147920), 1/27986464*(Lambda2^5 + 266471*Lambda2^3 + 23008956*Lambda2 + 13993232)
				B2 = Lambda2
				B3 = 1/86*(Lambda22 + 43*Lambda2)
				B4 = 1/86*(Lambda23 + 43*Lambda2)
				B5 = 1/162712*(Lambda24 + 946*Lambda23 + 1591*Lambda22 + 122034*Lambda2 + 147920*I)
				B6 = 1/27986464*(Lambda25 + 266471*Lambda23 + 23008956*Lambda2 + 13993232*I)
			elif CMfieldnumber == 220 or CMfieldnumber == 2200 or CMfieldnumber == 2201: #basis:  1, Lambda2, 1/67*Lambda2^2, 1/402*(Lambda2^3 + 67*Lambda2 + 201), 1/404010*(Lambda2^4 + 3685*Lambda2^2 + 202005*Lambda2 + 242406), 1/243618030*(Lambda2^5 + 187600*Lambda2^3 + 1818045*Lambda2^2 + 24281001*Lambda2 + 121809015)
				B2 = Lambda2
				B3 = 1/67*Lambda22
				B4 = 1/402*(Lambda23 + 67*Lambda2 + 201*I)
				B5 = 1/404010*(Lambda24 + 3685*Lambda22 + 202005*Lambda2 + 242406*I)
				B6 = 1/243618030*(Lambda25 + 187600*Lambda23 + 1818045*Lambda22 + 24281001*Lambda2 + 121809015*I)
			elif CMfieldnumber == 31: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2), 1/4*(Lambda2^4 + Lambda2^2 + 2), 1/8*(Lambda2^5 + Lambda2^4 + Lambda2^3 + Lambda2^2 + 6*Lambda2 + 6)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22)
				B5 = 1/4*(Lambda24 + Lambda22 + 2*I)
				B6 = 1/8*(Lambda25 + Lambda24 + Lambda23 + Lambda22 + 6*Lambda2 + 6*I)
			elif CMfieldnumber == 311: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2), 1/2*(Lambda2^4 + Lambda2^2), 1/16*(Lambda2^5 + 4*Lambda2^4 + 5*Lambda2^3 + 12*Lambda2^2 + 14*Lambda2 + 8)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22)
				B5 = 1/2*(Lambda24 + Lambda22)
				B6 = 1/16*(Lambda25 + 4*Lambda24 + 5*Lambda23 + 12*Lambda22 + 14*Lambda2 + 8*I)
			elif CMfieldnumber == 32: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2 + 1), 1/2*(Lambda2^5 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2 + I)
				B6 = 1/2*(Lambda25 + Lambda2 + I)
			elif CMfieldnumber == 34: #basis:  1, Lambda2, Lambda2^2, 1/3*(Lambda2^3 + 2*Lambda2), 1/3*(Lambda2^4 + 2*Lambda2^2), 1/9*(Lambda2^5 + Lambda2^3 + 7*Lambda2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/3*(Lambda23 + 2*Lambda2)
				B5 = 1/3*(Lambda24 + 2*Lambda22)
				B6 = 1/9*(Lambda25 + Lambda23 + 7*Lambda2)
			elif CMfieldnumber == 35: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, Lambda2^4, Lambda2^5
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = Lambda25
			elif CMfieldnumber == 36: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2 + 1), 1/6*(Lambda2^4 + 3*Lambda2^2 + 3*Lambda2 + 2), 1/6*(Lambda2^5 + 3*Lambda2^2 + 5*Lambda2 + 3)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2 + I)
				B5 = 1/6*(Lambda24 + 3*Lambda22 + 3*Lambda2 + 2*I)
				B6 = 1/6*(Lambda25 + 3*Lambda22 + 5*Lambda2 + 3*I)
			elif CMfieldnumber == 37: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, Lambda2^4, Lambda2^5
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = Lambda25
			elif CMfieldnumber == 38: #basis:  1, Lambda2, 1/2*(Lambda2^2 + Lambda2), 1/4*(Lambda2^3 + Lambda2 + 2), 1/8*(Lambda2^4 + 3*Lambda2^2 + 4), 1/16*(Lambda2^5 + Lambda2^4 + 3*Lambda2^3 + 3*Lambda2^2 + 4*Lambda2 + 4)
				B2 = Lambda2
				B3 = 1/2*(Lambda22 + Lambda2)
				B4 = 1/4*(Lambda23 + Lambda2 + 2*I)
				B5 = 1/8*(Lambda24 + 3*Lambda22 + 4*I)
				B6 = 1/16*(Lambda25 + Lambda24 + 3*Lambda23 + 3*Lambda22 + 4*Lambda2 + 4*I)
			elif CMfieldnumber == 388: #basis:  1, Lambda2, 1/2*(Lambda2^2 + Lambda2), 1/4*(Lambda2^3 + Lambda2 + 2), 1/4*(Lambda2^4 + Lambda2^2 + 2*Lambda2), 1/8*(Lambda2^5 + Lambda2^4 + Lambda2^3 + 3*Lambda2^2 + 2*Lambda2)
				B2 = Lambda2
				B3 = 1/2*(Lambda22 + Lambda2)
				B4 = 1/4*(Lambda23 + Lambda2 + 2)
				B5 = 1/4*(Lambda24 + Lambda22 + 2*Lambda2)
				B6 = 1/8*(Lambda25 + Lambda24 + Lambda23 + 3*Lambda22 + 2*Lambda2)
			elif CMfieldnumber == 39: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, Lambda2^4, Lambda2^5
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = Lambda25
			elif CMfieldnumber == 40: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22 + Lambda2 + I)
			elif CMfieldnumber == 41: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2 + 1), 1/2*(Lambda2^5 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2 + I)
				B6 = 1/2*(Lambda25 + Lambda2 + I)
			elif CMfieldnumber == 42: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2), 1/4*(Lambda2^4 + Lambda2^2 + 2), 1/8*(Lambda2^5 + Lambda2^4 + Lambda2^3 + Lambda2^2 + 6*Lambda2 + 6)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22)
				B5 = 1/4*(Lambda24 + Lambda22 + 2)
				B6 = 1/8*(Lambda25 + Lambda24 + Lambda23 + Lambda22 + 6*Lambda2 + 6*I)
			elif CMfieldnumber == 422: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2), 1/2*(Lambda2^4 + Lambda2^2), 1/16*(Lambda2^5 + 4*Lambda2^4 + 5*Lambda2^3 + 12*Lambda2^2 + 14*Lambda2 + 8)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22)
				B5 = 1/2*(Lambda24 + Lambda22)
				B6 =  1/16*(Lambda25 + 4*Lambda24 + 5*Lambda23 + 12*Lambda22 + 14*Lambda2 + 8*I)
			elif CMfieldnumber == 43: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, Lambda2^4, 1/2*(Lambda2^5 + Lambda2^4)
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = 1/2*(Lambda25 + Lambda24)
			elif CMfieldnumber == 45: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + 1), 1/2*(Lambda2^4 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + I)
				B5 = 1/2*(Lambda24 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22)
			elif CMfieldnumber == 46: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, Lambda2^4, 1/2*(Lambda2^5 + Lambda2^4)
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = 1/2*(Lambda25 + Lambda24)
			elif CMfieldnumber == 47: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2), 1/2*(Lambda2^4 + Lambda2^2), 1/2*(Lambda2^5 + Lambda2^2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22)
				B5 = 1/2*(Lambda24 + Lambda22)
				B6 = 1/2*(Lambda25 + Lambda22)
			elif CMfieldnumber == 48: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2 + Lambda2), 1/2*(Lambda2^4 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22 + Lambda2)
				B5 = 1/2*(Lambda24 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22)
			elif CMfieldnumber == 488: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, 1/2*(Lambda2^4 + Lambda2^3 + Lambda2^2), 1/2*(Lambda2^5 + Lambda2^2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = 1/2*(Lambda24 + Lambda23 + Lambda22)
				B6 = 1/2*(Lambda25 + Lambda22)
			elif CMfieldnumber == 49: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2 + Lambda2 + 1), 1/4*(Lambda2^4 + 2*Lambda2 + 1), 1/4*(Lambda2^5 + 2*Lambda2^2 + Lambda2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22 + Lambda2 + I)
				B5 = 1/4*(Lambda24 + 2*Lambda2 + I)
				B6 = 1/4*(Lambda25 + 2*Lambda22 + Lambda2)
			elif CMfieldnumber == 50: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2^2 + Lambda2), 1/2*(Lambda2^4 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda22 + Lambda2)
				B5 = 1/2*(Lambda24 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22)
			elif CMfieldnumber == 51: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22 + Lambda2 + I)
			elif CMfieldnumber == 52: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22 + Lambda2 + I)
			elif CMfieldnumber == 53: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22 + Lambda2 + I)
			elif CMfieldnumber == 54: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2), 1/2*(Lambda2^5 + Lambda2^2 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda22 + Lambda2 + I)
			elif CMfieldnumber == 56: #basis:  1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2 + 1), 1/2*(Lambda2^4 + Lambda2^2 + Lambda2), 1/2*(Lambda2^5 + Lambda2 + 1)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2 + I)
				B5 = 1/2*(Lambda24 + Lambda22 + Lambda2)
				B6 = 1/2*(Lambda25 + Lambda2 + I)
			elif CMfieldnumber == 57: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, Lambda2^4, 1/4*(Lambda2^5 + 2*Lambda2^4 + 2*Lambda2^3 + 3*Lambda2 + 2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = 1/4*(Lambda25 + 2*Lambda24 + 2*Lambda23 + 3*Lambda2 + 2*I)
			elif CMfieldnumber == 58: #basis:  1, Lambda2, Lambda2^2, Lambda2^3, Lambda2^4, 1/2*(Lambda2^5 + Lambda2^4)
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = Lambda24
				B6 = 1/2*(Lambda25 + Lambda24)
			elif CMfieldnumber == 59: #basis: 1, Lambda2, Lambda2^2, Lambda2^3, 1/12*(Lambda2^4 + 6*Lambda2^3 + 11*Lambda2^2 + 6*Lambda2 + 10), 1/12*(Lambda2^5 + 11*Lambda2^3 + 10*Lambda2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = Lambda23
				B5 = 1/12*(Lambda24 + 6*Lambda23 + 11*Lambda22 + 6*Lambda2 + 10*I)
				B6 = 1/12*(Lambda25 + 11*Lambda23 + 10*Lambda2)
			elif CMfieldnumber == 599: #basis: 1, Lambda2, 1/2*(Lambda2^2 + 1), 1/2*(Lambda2^3 + Lambda2), 1/4*(Lambda2^4 + 3), 1/24*(Lambda2^5 + 3*Lambda2^4 + 8*Lambda2^3 + 7*Lambda2 + 21)
				B2 = Lambda2
				B3 = 1/2*(Lambda22 + I)
				B4 = 1/2*(Lambda23 + Lambda2)
				B5 = 1/4*(Lambda24 + 3*I)
				B6 = 1/24*(Lambda25 + 3*Lambda24 + 8*Lambda23 + 7*Lambda2 + 21*I)
			elif CMfieldnumber == 65: #basis: 1, Lambda2, Lambda2^2, 1/6*(Lambda2^3 + 3*Lambda2^2), 1/6*(Lambda2^4 + 3*Lambda2^2), 1/6*(Lambda2^5 + 3*Lambda2^2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/6*(Lambda23 + 3*Lambda22)
				B5 = 1/6*(Lambda24 + 3*Lambda22)
				B6 = 1/6*(Lambda25 + 3*Lambda22)
			elif CMfieldnumber == 71: #basis: 1, Lambda2, 1/2*Lambda2^2, 1/4*(Lambda2^3 + 2*Lambda2), 1/24*(Lambda2^4 + 10*Lambda2^2 + 12*Lambda2), 1/24*(Lambda2^5 + 4*Lambda2^3 + 12*Lambda2)
				B2 = Lambda2
				B3 = 1/2*Lambda22
				B4 = 1/4*(Lambda23 + 2*Lambda2)
				B5 = 1/24*(Lambda24 + 10*Lambda22 + 12*Lambda2)
				B6 = 1/24*(Lambda25 + 4*Lambda23 + 12*Lambda2)
			elif CMfieldnumber == 711: #basis: 1, Lambda2, Lambda2^2, 1/2*(Lambda2^3 + Lambda2), 1/12*(Lambda2^4 + 3*Lambda2^3 + 7*Lambda2^2 + 9*Lambda2), 1/24*(Lambda2^5 + 10*Lambda2^3 + 21*Lambda2)
				B2 = Lambda2
				B3 = Lambda22
				B4 = 1/2*(Lambda23 + Lambda2)
				B5 = 1/12*(Lambda24 + 3*Lambda23 + 7*Lambda22 + 9*Lambda2)
				B6 = 1/24*(Lambda25 + 10*Lambda23 + 21*Lambda2)
			if matrix_with_entries_in_R_and_ROVERn(B2,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B3,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B4,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B5,n,basis) == True and matrix_with_entries_in_R_and_ROVERn(B6,n,basis) == True:
				T = True
			if T == True:
				Nx1 = QA(x1).reduced_norm()
				Nx2 = QA(x2).reduced_norm()
				Nx3 = QA(x3).reduced_norm()
				Trx1x2 = QA(x1*x2).reduced_trace()
				Trx1x3 = QA(x1*x3).reduced_trace()
				Trx2x3 = QA(x2*x3).reduced_trace()
				Dx1x2 = 4*Nx1*Nx2 - Trx1x2^2
				Dx1x3 = 4*Nx1*Nx3 - Trx1x3^2
				Dx2x3 = 4*Nx2*Nx3 - Trx2x3^2
				Dx1x2x3 = 4*Nx1*Nx2*Nx3 - Nx1*Trx2x3^2 - Nx2*Trx1x3^2 - Nx3*Trx1x2^2 - Trx1x2*Trx1x3*Trx2x3 # = (Trx1x2x3)^2
				if Dx1x2x3.is_square() and Dx1x2x3/p in ZZ and Dx1x2/p in ZZ and Dx1x3/p in ZZ and Dx2x3/p in ZZ:
					if Dx1x2 != 0 or Dx1x3 != 0 or Dx2x3 != 0: #this condition shouldn't be here is we run for fields with D!!!
						print ('[x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3] =', [x, dOVERn, a, cOVERn, b, gamma, n, QA(x1),QA(x2),QA(x3), QA(x4),QA(x5),QA(x6), QA(x7),QA(x8),QA(x9), Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
						csvwriter2.writerow([x, dOVERn, a, cOVERn, b, gamma, n, QA(x1),QA(x2),QA(x3), QA(x4),QA(x5),QA(x6), QA(x7),QA(x8),QA(x9), Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
						csvfile2.flush()
						OUTPUT.append([x, dOVERn, a, cOVERn, b, gamma, n, QA(x1),QA(x2),QA(x3), QA(x4),QA(x5),QA(x6), QA(x7),QA(x8),QA(x9), Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
						sys.stdout.flush()
	print ("done")
	print (len(OUTPUT)," solutions found")
	print ("----------------------------------------------------------------------")
	sys.stdout.flush()
	csvfile2.close()
	print ("All done in ", time() - TT)

