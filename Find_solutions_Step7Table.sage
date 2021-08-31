# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step performs a check before Step 7 of Algorithm 1 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step7Table.sage CMfieldnumber p (q)
# e.g. sage ./Find_solutions_Step7.sage 2 3 runs for CM field 2, p = 3, all max orders in the table BasesForPrimesTable.csv
# If p = 1 mod 8, then choose a prime q = 3 mod 4 and q not a square mod p. Ensure you choose q such that the bases appear in BasesForPrimesTable.csv.
# e.g. sage ./Find_solutions_Step7Table.sage 26666 17 3 runs for CM field 26666, p = 17, q = 3, all max orders in the table BasesForPrimesTable.csv
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 7: 
# Read solutions [x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3] from Step 6
# Check if Lambda32 == -D*II 
# and -Lambda13 + A*Lambda12 - B*Lambda1 + C*II == Zero 
# and Lambda1*Lambda3==Lambda3*Lambda1 
# and Lambda3 has row1, n*row2, n*row3 in max order
# and if -D = 1 mod 4 check 1/2*II + 1/2*Lambda3 has row1, n*row2, n*row3 in max order
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# For each max order, a file CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"+O+Solutions_for_di_Step7.csv"
# Example: CMfield1p7MaxOrderBasis_1_+i_+frac{1}{2}i+frac{1}{2}k_frac{1}{2}+frac{1}{2}j_Solutions_for_di_Step7.csv
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv

#######################################################
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


def matrix_with_entries_in_R_and_ROVERn(MM,n,basis): #Given a matrix check if first row entries are in R and second and third row entries are in R/n
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
	print ("reading solutions file from Step 6...")
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
	filename1 = filename1+"Solutions_for_di_Step6.csv"
	reader1 = csv.DictReader(open(filelocation + filename1, 'r'))
	L = []
	for row in reader1: #x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3
		L += [[QQ(row["x"]), QQ(row["dOVERn"]), QQ(row["a"]), QQ(row["cOVERn"]), QQ(row["b"]), QQ(row["gamma"]), QQ(row["n"]), eval(preparse(row["d1"])), eval(preparse(row["d2"])), eval(preparse(row["d3"])), eval(preparse(row["Dd1d2"])), eval(preparse(row["Dd1d3"])), eval(preparse(row["Dd2d3"])), eval(preparse(row["Dd1d2d3"])) ]]		
	LengthAllSol = len(L)
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
	filename2 = filename2+"Solutions_for_di_Step7.csv"
	csvfile2 = open(filelocation + filename2, 'w')
	csvwriter2 = csv.writer(csvfile2) #x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3
	csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'Dd1d2', 'Dd1d3', 'Dd2d3', 'Dd1d2d3'))
	count = 0
	print ("Number of read solutions from Step 6 = ", LengthAllSol)
	for LL in L:
		x = QQ(LL[0])
		dOVERn = QQ(LL[1])
		a = QQ(LL[2])
		cOVERn = QQ(LL[3])
		b = QQ(LL[4])
		gamma = QQ(LL[5])
		n = QQ(LL[6])
		d1 = QA(LL[7])
		d2 = QA(LL[8])
		d3 = QA(LL[9])
		Dd1d2 = QQ(LL[10])
		Dd1d3 = QQ(LL[11])
		Dd2d3 = QQ(LL[12])
		Dd1d2d3 = QQ(LL[13])
		gOVERn = QQ(gamma/n)
		aOVERn = QQ(a/n)
		bOVERn = QQ(b/n)
		k1 = QQ(aOVERn*dOVERn - aOVERn*x - bOVERn)
		k2 = QQ(aOVERn*cOVERn + bOVERn*x)
		k3 = QQ(k1*cOVERn - k2*dOVERn)
		k4 = QQ(k2*cOVERn - bOVERn*b)
		k5 = QQ(aOVERn*a - k1*dOVERn - k2)
		k6 = QQ(aOVERn*b + k1*cOVERn)
		d4 = QA(gOVERn*d2 - bOVERn*d3)
		d5 = QA(d1 + k3*d2 + k2*d3)
		d6 = QA(k4*d2 + k6*d3)
		d7 = QA(-bOVERn*d2 + aOVERn*d3)
		d8 = QA(k2*d2 + k1*d3)
		d9 = QA(d1 + k6*d2 - k5*d3)
		Lambda1 = matrix(QQ, 3,3, [x,a,b, 1,0,cOVERn, 0,1,dOVERn])
		Lambda3 = matrix(QA, 3,3, [d1,d2,d3, d4,d5,d6, d7,d8,d9]) 
		Zero = matrix(QA, 3,3, [0,0,0, 0,0,0, 0,0,0])
		II = matrix(QA, 3,3, [1,0,0, 0,1,0, 0,0,1]) # -Lambda1^3 + A*Lambda1^2 - B*Lambda1 + C*II = O
		Lambda32 = Lambda3^2
		Lambda12 = Lambda1^2
		Lambda13 = Lambda1^3
		if Lambda32 == -D*II and -Lambda13 + A*Lambda12 - B*Lambda1 + C*II == Zero and Lambda1*Lambda3==Lambda3*Lambda1 and matrix_with_entries_in_R_and_ROVERn(Lambda3,n,basis): 
			S = Integers(4)
			if S(D) == 3:
				if matrix_with_entries_in_R_and_ROVERn(1/2*II + 1/2*Lambda3,n,basis):
					#print 'LL, Lambda3 =', LL, Lambda3
					csvwriter2.writerow([x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3])
					csvfile2.flush()
					print ([x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3])
					count = count + 1
					sys.stdout.flush()
			else:
				#print 'LL, Lambda3 =', LL, Lambda3
				csvwriter2.writerow([x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3])
				csvfile2.flush()
				print ([x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3])
				count = count + 1
				sys.stdout.flush()
	sys.stdout.flush()
	print ("Number of solutions left = ", count)
	csvfile2.close()

print (" ")
print ("All done in ", time() - TT)
print (" ")


