# run as sage ./Find_solutions_Step6.sage p max_order CMfieldnumber
# e.g. sage ./Find_solutions_Step6.sage 3 1 0 0 0 0 1 0 0 0 1/2 0 1/2 1/2 0 1/2 0 2 runs for p=3, O=1 with basis [1, i, i/2 + k/2, 1/2 + j/2], CM field 2.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 5 created a file containing all elements d_i with trace 0 and norms Nd1, Nd2, Nd3 that occur in the solutions in Step 4.
# Read solutions file from Step 4. Shape of solution: [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3].
# Read di's from file from Step 5.
# Check if a solution in terms of Norms and Traces (Step 4) leads to a solution in terms of elements di.
# That is, the function "find_solutions" checks for each solution from Step 4, if the corresponding di's from Step 5 with norms Nd1, Nd2, Nd3 can give traces with values Trd1d2, Trd1d3, Trd2d3.
# Save the solution in terms of di. Shape of solution: [x, dOVERn, a, cOVERn, b, gamma, n, d1, d2, d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3].
#-------------------------------------------------------------------------------------------------------------------------------------------
from sage.all_cmdline import * 
from time import time
from random import randint
import csv
#-------------------------------------------------------------------------------------------------------------------------------------------
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

print ("CMfieldnumber = ",CMfieldnumber,", p = ",p,", max_order = ",O)
print (sys.version)
print ("----------------------------------------------------------------------")
sys.stdout.flush()

T = time()

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
filename = filename+"d.csv"
reader=csv.DictReader(open(filelocation+filename, 'r')) 


all_xs = {}
for row in reader:
	if int(row["N"]) not in all_xs: 
		all_xs[int(row["N"])] = [eval(preparse(row["d with trace 0 and norm N"]))]
	else:
		all_xs[int(row["N"])].append(eval(preparse(row["d with trace 0 and norm N"])))
print ("done reading d_i file from Step 5, time = ",time()-tread)
print ("----------------------------------------------------------------------")
sys.stdout.flush()


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
		Nd1 = int(row["Nd1"])
		Nd2 = int(row["Nd2"])
		Nd3 = int(row["Nd3"])
		Trd1d2 = int(row["Trd1d2"])
		Trd1d3 = int(row["Trd1d3"])
		Trd2d3 = int(row["Trd2d3"])
		Dd1d2 = int(row["Dd1d2"])
		Dd1d3 = int(row["Dd1d3"])
		Dd2d3 = int(row["Dd2d3"])
		Dd1d2d3 = int(row["Dd1d2d3"])
		if Nd1 in all_xs and Nd2 in all_xs and Nd3 in all_xs:
			for d1 in all_xs[Nd1]:
				for d3 in all_xs[Nd3]:
					if Trd1d3 == QA(d1*d3).reduced_trace():
						for d2 in all_xs[Nd2]:
							if Trd1d2 == QA(d1*d2).reduced_trace() and Trd2d3 == QA(d2*d3).reduced_trace():
								solutions += [[x, dOVERn, a, cOVERn, b, gamma, n, d1, d2, d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]]
								csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n, d1, d2, d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
								csvfile.flush()
								print ("solution [x, dOVERn, a, cOVERn, b, gamma, n, n, d1, d2, d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3] = ",[x, dOVERn, a, cOVERn, b, gamma, n, d1, d2, d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])				
								sys.stdout.flush()
		print ("done for [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3] =", [x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
		SOLS += [solutions]
		L = L + len(solutions)
	csvfile.close()
	print ("Number of solutions = ", L)
	print ("----------------------------------------------------------------------")
	sys.stdout.flush()
	return SOLS



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
filename3 = filename3+"Solutions_for_di_Step6.csv"

csvfile3 = open(filelocation+filename3, 'w')
csvwriter3 = csv.writer(csvfile3)
csvwriter3.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'd1', 'd2', 'd3', 'Dd1d2', 'Dd1d3', 'Dd2d3', 'Dd1d2d3'))


#read potsolutions from file: x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3
print ("reading solutions file from Step 4...") #
sys.stdout.flush()
tread = time() 
filename2 = "CMfield"+CMfieldnumber.str()+"Step4"+"p"+p.str()+".csv"


reader2 =csv.DictReader(open(filelocation+filename2))
SOLUTIONS = find_solutions(p, reader2, csvwriter3, csvfile3)


csvfile3.close()


print ("All done in ", time()- T)


