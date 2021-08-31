# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements Algorithm 5 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step8xTable.sage CMfieldnumber p (q)
# e.g. sage ./Find_solutions_Step8xTable.sage 34 3 runs for CM field 34, p = 11, all max orders in the table BasesForPrimesTable.csv
# If p = 1 mod 8, then choose a prime q = 3 mod 4 and q not a square mod p. Ensure you choose q such that the bases appear in BasesForPrimesTable.csv.
# e.g. sage ./Find_solutions_Step93xTable.sage 34 17 3 runs for CM field 34, p = 17, q = 3, all max orders in the table BasesForPrimesTable.csv
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 8:
# Read solutions x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,x4,x5,x6,x7,x8,x9,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3 from Step 7
# Check that if they are equivalent due to Proposition 4.2
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# One file CMfield"+CMfieldnumber.str()+"p"+p.str()+"MaxOrderBasis_"+O+Solutions_for_xi_Step93.txt" with equivalence classes of solutions in the form [x, dOVERn, a, cOVERn, b, gamma, n, x1, x2, x3, x4, x5, x6, x7, x8, x9, Dx1x2, Dx1x3, Dx2x3, Dx1x2x3]
# Example: CMfield34p11MaxOrderBasis_1_+i_+frac{1}{2}i+frac{1}{2}k_frac{1}{2}+frac{1}{2}j_Solutions_for_xi_Step93.txt
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


load('CM_Fields.sage')


TT = time()


load('All_Roots.sage')

U = []
if len(all_roots[CMfieldnumber]) >1:
	for r in all_roots[CMfieldnumber]:
		for u in [0..2]:
			U.append(r.coefficient({mu:u}))
	Roots = matrix(QQ,3,U)
elif len(all_roots[CMfieldnumber]) ==1:
	Roots = matrix(QQ,1,[0,1,0])


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
	filename1 = filename1+"Solutions_for_xi_Step7.csv"
	reader1 = csv.DictReader(open(filelocation + filename1, 'r'))
	L = []
	for row in reader1: #x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,x4,x5,x6,x7,x8,x9,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3
		L += [ [ QQ(row["x"]), QQ(row["dOVERn"]), QQ(row["a"]), QQ(row["cOVERn"]), QQ(row["b"]), QQ(row["gamma"]), QQ(row["n"]), eval(preparse(row["x1"])), eval(preparse(row["x2"])), eval(preparse(row["x3"])), eval(preparse(row["x4"])), eval(preparse(row["x5"])), eval(preparse(row["x6"])), eval(preparse(row["x7"])), eval(preparse(row["x8"])), eval(preparse(row["x9"])), eval(preparse(row["Dx1x2"])), eval(preparse(row["Dx1x3"])), eval(preparse(row["Dx2x3"])), eval(preparse(row["Dx1x2x3"])) ] ]
	LengthAllSol = len(L)
	print ("Number of trivial and non-trivial solutions = ", LengthAllSol)
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
	filename2 = filename2+"Solutions_for_xi_Step93.txt"
	#csvfile2 = open(filelocation + filename2, 'w')
	#csvwriter2 = csv.writer(csvfile2) #x,dOVERn,a,cOVERn,b,gamma,n,x1,x2,x3,Dx1x2,Dx1x3,Dx2x3,Dx1x2x3
	#csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'x1', 'x2', 'x3','Dx1x2', 'Dx1x3', 'Dx2x3', 'Dx1x2x3'))
	#write all solutions in txt files
	file2 = open(filelocation+filename2, 'w')
	file2.write('sets of equivalent solutions: ' + '\n')
	OUTPUT = []
	r1 = 0
	all = copy(L)
	while r1 < len(all): #
		sol = []
		x = all[r1][0]
		dOVERn = all[r1][1]
		a = all[r1][2]
		cOVERn = all[r1][3]
		b = all[r1][4]
		gamma = all[r1][5]
		n = all[r1][6]
		for r2 in [0..(Roots.nrows()-1)]:
			a0 = Roots[r2][0]
			a1 = Roots[r2][1]
			a2 = Roots[r2][2]
			xx = a0 + a1*x + a2*x^2 + a2*a
			aa = a*(a2^2*cOVERn + (a1+a2*x)^2) + b*(a2^2*dOVERn + 2*a2*(a1+a2*x))
			dOVERnn = A - xx
			cOVERnn = dOVERnn*xx - aa - B
			bb = cOVERnn*xx + dOVERnn*aa + C
			gammaa = aa*cOVERnn + bb*dOVERnn
			nn = aa*gammaa - bb^2
			SS = [xx,dOVERnn,aa,cOVERnn,bb,gammaa,nn]
			for sols in L:
				if SS == [sols[0],sols[1],sols[2],sols[3],sols[4],sols[5],sols[6]]:
					if sols in all:
						sol.append(sols)
						all.remove(sols)
		if len(all) == 1:
			sol = [all[0]]
			all.remove(all[0])
		OUTPUT.append(sol)
		print(sol)
		if len(all) == 0:
			break
	sys.stdout.flush()
	print (len(OUTPUT)," sets of equivalent Lambda1s")
	OUTPUT2 = []
	for S in OUTPUT:
		temp = []
		for LL in S:
			x = QA(LL[0])
			dOVERn = QA(LL[1])
			a = QA(LL[2])
			cOVERn = QA(LL[3])
			b = QA(LL[4])
			gamma = QA(LL[5])
			n = QA(LL[6])
			x1 = QA(LL[7])
			x2 = QA(LL[8])
			x3 = QA(LL[9])
			x4 = QA(LL[10])
			x5 = QA(LL[11])
			x6 = QA(LL[12])
			x7 = QA(LL[13])
			x8 = QA(LL[14])
			x9 = QA(LL[15])
			Dx1x2 = QQ(LL[16])
			Dx1x3 = QQ(LL[17])
			Dx2x3 = QQ(LL[18])
			Dx1x2x3 = QQ(LL[19])
			sol = matrix(QA, 3, 3, [x1,x2,x3, x4,x5,x6, x7,x8,x9]) #Lambda2
			sol2 = matrix(QA, 3, 3, [x,a,b, 1,0,cOVERn, 0,1,dOVERn]) #Lambda1
			TQ = []
			Tlist = []
			Qlist = []
			for r in [0..2]:
				for c in [0..2]:
					t0 = sol[r][c][0]
					t1 = sol[r][c][1]
					t2 = sol[r][c][2]
					t3 = sol[r][c][3]
					q0 = sol2[r][c][0]
					q1 = sol2[r][c][1]
					q2 = sol2[r][c][2]
					q3 = sol2[r][c][3]
					Telem = matrix(QQ, 4, [t0,-t1,-t2,-t3, t1,t0,-t3,t2, t2,t3,t0,-t1, t3,-t2,t1,t0])
					Qelem = matrix(QQ, 4, [q0,-q1,-q2,-q3, q1,q0,-q3,q2, q2,q3,q0,-q1,  q3,-q2,q1,q0])
					Tlist += [Telem]
					Qlist += [Qelem]
			T = block_matrix(QQ, 3, 3, Tlist, subdivide=True)
			Q = block_matrix(QQ, 3, 3, Qlist, subdivide=True)
			TQ += [T,Q]
			temp.append(TQ)
		OUTPUT2.append(temp)
	final = []
	for TQ in OUTPUT2:
		all_indices = [0..len(TQ)-1]
		used_indices = []
		while (len(used_indices) < len(all_indices)):
			for r3 in all_indices:
				indices = []
				if r3 not in used_indices:
					equiv_class = []
					N = TQ[r3][0] #Lambda2 
					V = TQ[r3][1] #Lambda1 = [x,a,b, 1,0,cOVERn, 0,1,dOVERn]
					for pair in TQ: #(Lambda2,Lambda1)
						if (N.is_similar(pair[0])) and (V.is_similar(pair[1])):
							QApair01 = QA(pair[0][0][0] + pair[0][1][0]*i + pair[0][2][0]*j + pair[0][3][0]*k) #x1
							QApair02 = QA(pair[0][0][4] + pair[0][1][4]*i + pair[0][2][4]*j + pair[0][3][4]*k) #x2
							QApair03 = QA(pair[0][0][8] + pair[0][1][8]*i + pair[0][2][8]*j + pair[0][3][8]*k) #x3
							QApair04 = QA(pair[0][4][0] + pair[0][5][0]*i + pair[0][6][0]*j + pair[0][7][0]*k) #x4
							QApair05 = QA(pair[0][4][4] + pair[0][5][4]*i + pair[0][6][4]*j + pair[0][7][4]*k) #x5
							QApair06 = QA(pair[0][4][8] + pair[0][5][8]*i + pair[0][6][8]*j + pair[0][7][8]*k) #x6
							QApair07 = QA(pair[0][8][0] + pair[0][9][0]*i + pair[0][10][0]*j + pair[0][11][0]*k) #x7
							QApair08 = QA(pair[0][8][4] + pair[0][9][4]*i + pair[0][10][4]*j + pair[0][11][4]*k) #x8
							QApair09 = QA(pair[0][8][8] + pair[0][9][8]*i + pair[0][10][8]*j + pair[0][11][8]*k) #x9
							#QApair0 = [QApair01, QApair02, QApair03, QApair04, QApair05, QApair06, QApair07, QApair08, QApair09]
							#QApair0 = matrix(QA,3,QApair0)
							QApair11 = QA(pair[1][0][0] + pair[1][1][0]*i + pair[1][2][0]*j + pair[1][3][0]*k) #x
							QApair12 = QA(pair[1][0][4] + pair[1][1][4]*i + pair[1][2][4]*j + pair[1][3][4]*k) #a
							QApair13 = QA(pair[1][0][8] + pair[1][1][8]*i + pair[1][2][8]*j + pair[1][3][8]*k) #b
							#QApair14 = QA(pair[1][4][0] + pair[1][5][0]*i + pair[1][6][0]*j + pair[1][7][0]*k) #1
							#QApair15 = QA(pair[1][4][4] + pair[1][5][4]*i + pair[1][6][4]*j + pair[1][7][4]*k) #0
							QApair16 = QA(pair[1][4][8] + pair[1][5][8]*i + pair[1][6][8]*j + pair[1][7][8]*k) #cOVERn
							#QApair17 = QA(pair[1][8][0] + pair[1][9][0]*i + pair[1][10][0]*j + pair[1][11][0]*k) #0
							#QApair18 = QA(pair[1][8][4] + pair[1][9][4]*i + pair[1][10][4]*j + pair[1][11][4]*k) #1
							QApair19 = QA(pair[1][8][8] + pair[1][9][8]*i + pair[1][10][8]*j + pair[1][11][8]*k) #dOVERn
							#QApair1 = [QApair11, QApair12, QApair13, QApair14, QApair15, QApair16, QApair17, QApair18, QApair19]
							QApair1g = QApair12*QApair16 + QApair13*QApair19 #gamma = a*cOVERn + b*dOVERn
							QApair1n = QApair12*QApair1g-QApair13^2 #n = a*gamma - b^2
							#QApair1 = matrix(QA,3,QApair1) 
							origsol = [QApair11, QApair19, QApair12, QApair16, QApair13, QApair1g, QApair1n, QApair01, QApair02, QApair03, QApair04, QApair05, QApair06, QApair07, QApair08, QApair09] #x,dOVERn,a,cOVERn,b,gamma,n, x1, ..., x9
							#print ('solution',  origsol)
							indices += [TQ.index(pair)+1]
							used_indices += [TQ.index(pair)]
							#print ('pair', pair)
							#print ('solution',  origsol)
							equiv_class.append(origsol)
					print('equivalence class', equiv_class)
					file2.write(str(equiv_class))
					file2.write(',' + '\n')
					final.append(equiv_class)
	print('There are ', len(final), 'equivalence classes of solutions.')
	file2.close()


print (" ")
print ("All done in ", time()- TT)
print (" ")
