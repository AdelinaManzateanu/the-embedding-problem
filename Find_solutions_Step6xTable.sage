# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements Algorithm 4 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step6xTable.sage CMfieldnumber p (q)
# e.g. sage ./Find_solutions_Step6xTable.sage 34 11 runs for CM field 34, p = 11, all max orders in the table BasesForPrimesTable.csv
# If p = 1 mod 8, then choose a prime q = 3 mod 4 and q not a square mod p. Ensure you choose q such that the bases appear in BasesForPrimesTable.csv.
# e.g. sage ./Find_solutions_Step6xTable.sage 34 17 3 runs for CM field 34, p = 17, q = 3, all max orders in the table BasesForPrimesTable.csv
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 6: Check if a solution in terms of Norms and Traces (Step 4) leads to a solution in terms of elements di.
# Read solutions file from Step 4: solutions have the form [x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]
# Read di's from file from Step 5: all di's with trace 0 and norms Nx1, Nx2, Nx3.
# For each solution from Step 4, check if the xi's from Step 5, with norms Nx1, Nx2, Nx3, can give traces with values Trx1x2, Trx1x3, Trx2x3.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# For each maximal order, a file CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"+O+Solutions_for_di_Step6.csv"
# Example: CMfield34p11MaxOrderBasis_1_+i_+frac{1}{2}i+frac{1}{2}k_frac{1}{2}+frac{1}{2}j_Solutions_for_di_Step6.csv
#-------------------------------------------------------------------------------------------------------------------------------------------

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


T = time()

@parallel
def find_solutions(p, reader, csvwriter, csvfile):
	SOLS = []
	L = 0
	for row in reader:
		solutions = []
		x = int(row["x"])
		dOVERn = int(row["dOVERn"])
		a = int(row["a"])
		cOVERn = int(row["cOVERn"])
		b = int(row["b"])
		gamma = int(row["gamma"])
		n = int(row["n"])
		Nx1 = int(row["Nx1"])
		Nx2 = int(row["Nx2"])
		Nx3 = int(row["Nx3"])
		Trx1x2 = int(row["Trx1x2"])
		Trx1x3 = int(row["Trx1x3"])
		Trx2x3 = int(row["Trx2x3"])
		Dx1x2 = int(row["Dx1x2"])
		Dx1x3 = int(row["Dx1x3"])
		Dx2x3 = int(row["Dx2x3"])
		Dx1x2x3 = int(row["Dx1x2x3"])
		if Nx1 in all_xs and Nx2 in all_xs and Nx3 in all_xs:
			for x1 in all_xs[Nx1]:
				for x3 in all_xs[Nx3]:
					if Trx1x3 == QA(x1*x3).reduced_trace():
						for x2 in all_xs[Nx2]:
							if Trx1x2 == QA(x1*x2).reduced_trace() and Trx2x3 == QA(x2*x3).reduced_trace():
								solutions += [[x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]]
								csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
								csvfile.flush()
								print ("solution [x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3] = ",[x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])				
								sys.stdout.flush()
		print ("done for [x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3] =", [x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3])
		SOLS += [solutions]
		L = L + len(solutions)
	csvfile.close()
	print ("Number of solutions = ", L)
	print ("----------------------------------------------------------------------")
	sys.stdout.flush()
	return SOLS


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
	print ("Maximal order with basis ", O, " in ", QA, ".")
	print ("reading d_i file...")
	sys.stdout.flush()
	tread = time()
	filename = "CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"
	for m in [0..3]:
		for n in [0..3]:
			if O[m][n] ==1:
				if n == 0:
					filename = filename+"1"
				if n == 1:
					filename = filename+"+i"
				elif n == 2:
					filename = filename+"+j"
				elif n ==3:
					filename = filename+"+k"
			elif O[m][n] !=0 and O[m][n] !=1:
				if n == 0:
					if O[m][n] in ZZ:
						filename = filename+QQ(O[m][n]).str()
					else:
						filename = filename+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"
				if n == 1:
					if O[m][n] in ZZ:
						filename = filename+"+"+QQ(O[m][n]).str()+"i"
					else:
						filename = filename+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"i"
				elif n == 2:
					if O[m][n] in ZZ:
						filename = filename+"+"+QQ(O[m][n]).str()+"j"
					else:
						filename = filename+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"j"
				elif n ==3:
					if O[m][n] in ZZ:
						filename = filename+"+"+QQ(O[m][n]).str()+"k"
					else:
						filename = filename+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"k"
		filename = filename+"_"
	filename = filename+"x.csv"
	reader=csv.DictReader(open(filelocation+filename, 'r')) 
	all_xs = {}
	for row in reader:
		if int(row["N"]) not in all_xs: 
			all_xs[int(row["N"])] = [eval(preparse(row["x with trace 0 and norm N"]))]
		else:
			all_xs[int(row["N"])].append(eval(preparse(row["x with trace 0 and norm N"])))
	print ("done reading x_i file from Step 5, time = ",time()-tread)
	print ("----------------------------------------------------------------------")
	sys.stdout.flush()
	#write all solutions in a csv file
	filename3 = "CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"
	for m in [0..3]:
		for n in [0..3]:
			if O[m][n] ==1:
				if n == 0:
					filename3 = filename3+"1"
				if n == 1:
					filename3 = filename3+"+i"
				elif n == 2:
					filename3 = filename3+"+j"
				elif n ==3:
					filename3 = filename3+"+k"
			elif O[m][n] !=0 and O[m][n] !=1:
				if n == 0:
					if O[m][n] in ZZ:
						filename3 = filename3+QQ(O[m][n]).str()
					else:
						filename3 = filename3+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"
				if n == 1:
					if O[m][n] in ZZ:
						filename3 = filename3+"+"+QQ(O[m][n]).str()+"i"
					else:
						filename3 = filename3+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"i"
				elif n == 2:
					if O[m][n] in ZZ:
						filename3 = filename3+"+"+QQ(O[m][n]).str()+"j"
					else:
						filename3 = filename3+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"j"
				elif n ==3:
					if O[m][n] in ZZ:
						filename3 = filename3+"+"+QQ(O[m][n]).str()+"k"
					else:
						filename3 = filename3+"+"+"frac{"+QQ(O[m][n]).numerator().str()+"}{"+QQ(M[m][n]).denominator().str()+"}"+"k"
		filename3 = filename3+"_"
	filename3 = filename3+"Solutions_for_xi_Step6.csv"
	csvfile3 = open(filelocation+filename3, 'w')
	csvwriter3 = csv.writer(csvfile3)
	csvwriter3.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'x1', 'x2', 'x3', 'Dx1x2', 'Dx1x3', 'Dx2x3', 'Dx1x2x3'))
	#read potsolutions from file: x, dOVERn, a, cOVERn, b, gamma, n, Nx1, Nx2, Nx3, Trx1x2, Trx1x3, Trx2x3, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3
	print ("reading solutions file from Step 4...") #
	sys.stdout.flush()
	tread = time() 
	filename2 = "CMfield"+CMfieldnumber.str()+"Step4"+"p"+p.str()+".csv"
	reader2 =csv.DictReader(open(filelocation+filename2))
	SOLUTIONS = find_solutions(p, reader2, csvwriter3, csvfile3)
	csvfile3.close()
	print('')

print ("All done in ", time()- T)
print('')

