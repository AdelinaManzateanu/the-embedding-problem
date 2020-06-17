# run as sage ./Find_solutions_Step7.sage p max_order CMfieldnumber D
# e.g. sage ./Find_solutions_Step7.sage 3 1 0 0 0 0 1 0 0 0 1/2 0 1/2 1/2 0 1/2 0 2 1 runs for p=3, O=1 with basis [1, i, i/2 + k/2, 1/2 + j/2], CM field 2, D = 1.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Read solutions x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3 from Step6
# Check if Lambda32 == -D*II 
# and -Lambda13 + A*Lambda12 - B*Lambda1 + C*II == Zero 
# and Lambda1*Lambda3==Lambda3*Lambda1 
# and Lambda3 has row1, n*row2, n*row3 in max order O
# Save solutions x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3 in a file.
#-------------------------------------------------------------------------------------------------------------------------------------------
from sage.all_cmdline import *  
from time import time
from random import randint
import csv
#--------------------------------
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
#M = matrix(QQ,4,sage_eval([sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5], sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9], sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13], sys.argv[14],sys.argv[15],sys.argv[16],sys.argv[17]]))

O = [ QA(M[0][0]+M[0][1]*i+M[0][2]*j+M[0][3]*k), QA(M[1][0]+M[1][1]*i+M[1][2]*j+M[1][3]*k), QA(M[2][0]+M[2][1]*i+M[2][2]*j+M[2][3]*k), QA(M[3][0]+M[3][1]*i+M[3][2]*j+M[3][3]*k)]
#M = matrix([O[0].coefficient_tuple(), O[1].coefficient_tuple(), O[2].coefficient_tuple(), O[3].coefficient_tuple()])

CMfieldnumber = sage_eval(sys.argv[18])
D = sage_eval(sys.argv[19])

TT = time()

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

filename1 = filename1+"Solutions_for_di_Step6.csv"
reader1 = csv.DictReader(open(filelocation + filename1, 'r'))
L = []
for row in reader1: #x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3
	L += [[QQ(row["x"]), QQ(row["dOVERn"]), QQ(row["a"]), QQ(row["cOVERn"]), QQ(row["b"]), QQ(row["gamma"]), QQ(row["n"]), eval(preparse(row["d1"])), eval(preparse(row["d2"])), eval(preparse(row["d3"])), eval(preparse(row["Dd1d2"])), eval(preparse(row["Dd1d3"])), eval(preparse(row["Dd2d3"])), eval(preparse(row["Dd1d2d3"])) ]]		

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
filename2 = filename2+"Solutions_for_di_Step7.csv"

csvfile2 = open(filelocation + filename2, 'w')
csvwriter2 = csv.writer(csvfile2) #x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3
csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'Dd1d2', 'Dd1d3', 'Dd2d3', 'Dd1d2d3'))


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
		#print 'LL, Lambda3 =', LL, Lambda3
		csvwriter2.writerow([x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3])
		csvfile2.flush()
		print ([x,dOVERn,a,cOVERn,b,gamma,n,d1,d2,d3,d4,d5,d6,d7,d8,d9,Dd1d2,Dd1d3,Dd2d3,Dd1d2d3])
		sys.stdout.flush()

sys.stdout.flush()

csvfile2.close()


print ("All done in ", time() - TT)



