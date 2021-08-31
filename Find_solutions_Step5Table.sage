# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements Step 5 in Algorithm 1 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step5Table.sage CMfieldnumber p (q)
# e.g. sage ./Find_solutions_Step5Table.sage 2 3 runs for CM field 2, p = 3, all max orders in the table BasesForPrimesTable.csv
# If p = 1 mod 8, then choose a prime q = 3 mod 4 and q not a square mod p. Ensure you choose q such that the bases appear in BasesForPrimesTable.csv.
# e.g. sage ./Find_solutions_Step5Table.sage 26666 17 3 runs for CM field 26666, p = 17, q = 3, all max orders in the table BasesForPrimesTable.csv
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 5: Find elements with trace 0 and a certain norm
# Given the CM field number and the prime p, read all solutions from the file "CMfield"+CMfieldnumber+"Step4"+"p"+p.csv" created at Step 4
# Solutions have the form: [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]
# Save all norms and find all d_i in O with Tr = 0 and those norms
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# For each max order, a file "CMfield"+CMfieldnumber+"p"+p+"MaxOrderBasis_"+O+d.csv" is created.
# example: CMfield1p7MaxOrderBasis_1_+i_+frac{1}{2}i+frac{1}{2}k_frac{1}{2}+frac{1}{2}j_d.csv
#-------------------------------------------------------------------------------------------------------------------------------------------
# To find the basis for the maximal order, use Magma. For example,
# p := 3;
# Q := RationalField();
# A<i,j,k>:=QuaternionAlgebra<Q|-1,-p>; #change this accordingly
# CC:=ConjugacyClasses(MaximalOrder(A));
# #CC;
# for G in CC do
#     Generators(G);
# end for;
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv

######################################################
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
filename4 = "BasesForPrimesTable.csv"
reader =csv.DictReader(open(filelocation+filename4,'r',encoding='utf-8-sig'))
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


# Read solutions [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3] from Step 4.
print ("reading solutions from Step 4...") #
sys.stdout.flush()
tread = time() 
filename3 = "CMfield"+CMfieldnumber.str()+"Step4"+"p"+p.str()+".csv"
reader =csv.DictReader(open(filelocation+filename3,'r'))
print ("Reading time = ", time()-tread)


# Create a list all_norms of all the norms for which we need to find (if they exist) elements with trace 0 in the maximal order max_order in the quaternion algebra QA
all_norms = []

# Each row if of the form [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]
for row in reader:
	x = QQ(row["x"])
	dOVERn = QQ(row["dOVERn"])
	a = QQ(row["a"])
	cOVERn = QQ(row["cOVERn"])
	b = QQ(row["b"])
	gamma = QQ(row["gamma"])
	n = QQ(row["n"])
	Nd1 = eval(preparse(row["Nd1"]))
	Nd2 = eval(preparse(row["Nd2"]))
	Nd3 = eval(preparse(row["Nd3"]))
	Trd1d2 = eval(preparse(row["Trd1d2"]))
	Trd1d3 = eval(preparse(row["Trd1d3"]))
	Trd2d3 = eval(preparse(row["Trd2d3"]))
	Dd1d2 = eval(preparse(row["Dd1d2"]))
	Dd1d3 = eval(preparse(row["Dd1d3"]))
	Dd2d3 = eval(preparse(row["Dd2d3"]))
	Dd1d2d3 = eval(preparse(row["Dd1d2d3"]))
	if Nd1 not in all_norms:
		all_norms += [Nd1]
	if Nd2 not in all_norms:
		all_norms += [Nd2]
	if Nd3 not in all_norms:
		all_norms += [Nd3]

print ("len(all_norms) = ", len(all_norms))
all_norms.sort()
print ("max norm = ", all_norms[len(all_norms)-1])


# Load the function that finds elements with trace 0 and a certain norm in the maximal order max_order in the quaternion algebra QA
load('FindElemWithTr0AndNormNshort.sage')

sys.stdout.flush()


# This function finds all elements with trace 0 and norm in a certain list in the maximal order max_order in the quaternion algebra QA
@parallel
def find_all_x_i(p, q, a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33, norms):
	all_x_i = {}
	for Nx_i in norms:
		if Nx_i%10 == 0:
			print ('Nx_i = ', Nx_i)
			sys.stdout.flush()
		all_x_i[Nx_i] = find_elem_with_norm_and_trace_zero_short(p,q,Nx_i,a00,a01,a02,a03,a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33)
		L = len(all_x_i[Nx_i])
		if L > 0:
			for i in [0..L-1]:
				csvwriter.writerow(([Nx_i , all_x_i[Nx_i][i]]))
				csvfile.flush()
		sys.stdout.flush()
	sys.stdout.flush()
	return all_x_i

print (" ")
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
	# Name the file in which we will write all elements found
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
	filename = filename+"d.csv"
	csvfile = open(filelocation+filename, 'w')
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(('N', 'd with trace 0 and norm N'))
	# Find all elements with trace 0 and norm in the list all_norms in the maximal order max_order in the quaternion algebra QA
	print ("finding d_i...")
	tx = time() 
	all = find_all_x_i(p, q, a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33, all_norms)
	csvfile.close()
	print ("time for finding all d_i's with all norms =", time() - tx)
	L=0
	for NN in all:
		L = L + len(all[NN])
	print ("There are ", L, "elements with trace 0 and these norms.") # " in the maximal order with basis ", O, " in ", QA, ".")
	print('')










