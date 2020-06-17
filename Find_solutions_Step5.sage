# run as sage ./Find_solutions_Step5.sage p max_order CMfieldnumber
# e.g. sage ./Find_solutions_Step5.sage 3 1 0 0 0 0 1 0 0 0 1/2 0 1/2 1/2 0 1/2 0 2 runs p=3, O=1 with basis [1, i, i/2 + k/2, 1/2 + j/2], for CM field 2.
################################################################################
#Read solutions from Step 4
#Solutions have the form: x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3
#Save all norms and find all x_i in O with Tr = 0 and those norms to create CMfield1p7MaxOrderBasis_1_+i_+frac{1}{2}i+frac{1}{2}k_frac{1}{2}+frac{1}{2}j_d

#-------------------------------------------------------------------------------------------------------------------------------------------
# To find the basis for the maximal order, use Magma. For example,
# p := 3;
# Q := RationalField();
# A<i,j,k>:=QuaternionAlgebra<Q|-1,-p>;
# CC:=ConjugacyClasses(MaximalOrder(A));
# #CC;
# Generators(CC[1]);

################################################################################
from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv
################################################################################
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


#read solutions: x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3
print ("reading solutions from Step 4...") #
sys.stdout.flush()
tread = time() 
filename3 = "CMfield"+CMfieldnumber.str()+"Step4"+"p"+p.str()+".csv"
reader =csv.DictReader(open(filelocation+filename3,'r'))
print ("Reading time = ", time()-tread)


all_norms = []

#x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3
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

B2 = M.determinant()

C = []
for m in [0..3]:
	for n in [0..3]:
		C += [M.delete_rows([m]).delete_columns([n]).determinant()]
C = matrix(4,C) #= [[C0_wtv, C1_wtv, C2_wtv, C3_wtv], [C0_wtu, C1_wtu, C2_wtu, C3_wtu], [C0_uvt, C1_uvt, C2_uvt, C3_uvt], [C0_uvw, C1_uvw, C2_uvw, C3_uvw]]

# Ci_uvw = det of minor obtained by deleting row 3 (for t) anc column i = C[]

B4_uvw = (C[3][0]^2 + C[3][1]^2)*p + (C[3][2]^2 + C[3][3]^2) #(C0_uvw^2 + C1_uvw^2)*p + (C2_uvw^2 + C3_uvw^2)
B4_uvt = (C[2][0]^2 + C[2][1]^2)*p + (C[2][2]^2 + C[2][3]^2) #(C0_uvt^2 + C1_uvt^2)*p + (C2_uvt^2 + C3_uvt^2)
B4_wtu = (C[1][0]^2 + C[1][1]^2)*p + (C[1][2]^2 + C[1][3]^2) #(C0_wtu^2 + C1_wtu^2)*p + (C2_wtu^2 + C3_wtu^2)
B4_wtv = (C[0][0]^2 + C[0][1]^2)*p + (C[0][2]^2 + C[0][3]^2) #(C0_wtv^2 + C1_wtv^2)*p + (C2_wtv^2 + C3_wtv^2)

bdu = sqrt(abs(B4_wtv/p))/abs(B2)
bdv = sqrt(abs(B4_wtu/p))/abs(B2)
bdw = sqrt(abs(B4_uvt/p))/abs(B2)
bdt = sqrt(abs(B4_uvw/p))/abs(B2)

sys.stdout.flush()


@parallel
def find_elem_with_norm_and_trace_zero(p,O,NN):
	sqN = sqrt(NN)
	res = []
	P = p^2
	if set(O) == set([ QA(1), QA(i), QA(1/2*i + 1/2*k), QA(1/2 + 1/2*j)]): #x = QA(u + v*(1/2 + 1/2*j) + w*i + t*(1/2*i + 1/2*k)). Trx = 2*u + v. Nx = + u^2 + u*v + (p+1)/4*v^2 + w^2 + w*t + (p+1)/4*t^2. Trx = 0 => v = -2*u => Nx = p*u^2 + w^2 + w*t + (p+1)/4*t^2 
		for u in [ceil(-sqN/sqrt(p))..floor(sqN/sqrt(p))]:
			MM = NN-p*u^2
			bddt = 2*sqrt(MM/p)
			for t in [ceil(-bddt)..floor(bddt)]:
				Delta_w = -p*t^2 + 4*MM
				if Delta_w >=0:
					w1 = (-t + sqrt(Delta_w))/2
					w2 = (-t - sqrt(Delta_w))/2
					if w1 in ZZ:
						res += [QA([0, w1 + t/2, -u, t/2])] #sols.append([u, - 2*u, w, t])
					if w2 in ZZ:
						res += [QA([0, w2 + t/2, -u, t/2])]
	elif set(O) == set([ QA(1), QA(1/2 - 1/4*i + 1/4*k), QA(1/2 + 3/4*i + 1/4*k), QA(-1/2 + 1/2*i + 1/2*j)]): # Magma for p = 5. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = -t + 2*u + v + w. Nx = 1/4*p*t^2 + 1/16*p*v^2 + 1/8*p*v*w + 1/16*p*w^2 + 1/2*t^2 - t*u + u^2 - 3/4*t*v + u*v + 5/16*v^2 + 1/4*t*w + u*w + 1/8*v*w + 13/16*w^2. Trx = 0 => v = -2*u - w + t. => x = (1/4*t + 1/2*u + w)*i + 1/2*t*j + (1/4*t - 1/2*u)*k. Nx = 5/16*p*t^2 - 1/4*p*t*u + 1/4*p*u^2 + 1/16*t^2 + 1/4*t*u + 1/4*u^2 + 1/2*t*w + u*w + w^2.
		newbdu = sqrt(3/p) #bdu = sqrt(abs((p + 5)/p)).
		for u in [ceil(-newbdu*sqN)..floor(newbdu*sqN)]: # Delta_wt = (-32)*u^2*p^2 + 96*N*p >=0 => 3/p*N >= u^2
			D2 = p*u
			D1 = -D2*u
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wu  = (-32)*t^2*p^2 + 128*N*p >=0 => 4/p*N >= t^2. Same as bdt.
				Delta_w_over_4 = D1 + D2*t + (-3*p/4)*t^2 + 2*NN
				if Delta_w_over_4 >= 0: #norm_eq = + 2*w^2 + (2*u + t)*w + (p+1)/2*u^2 + (-p+1)/2*u*t + (3*p+1)/8*t^2 - N
					w1 = -(2*u+t)/4 + sqrt(Delta_w_over_4)/2
					w2 = -(2*u+t)/4 - sqrt(Delta_w_over_4)/2
					if w1 in ZZ: #x = (1/4*t + 1/2*u + w)*i + 1/2*t*j + (1/4*t - 1/2*u)*k
						res += [QA([0, t/4 + u/2 + w1, t/2, t/4 - u/2])] #sols.append([u, v, w, t]) #v = -2*u - w + t
					if w2 in ZZ: #x = (1/4*t + 1/2*u + w)*i + 1/2*t*j + (1/4*t - 1/2*u)*k
						res += [QA([0, t/4 + u/2 + w2, t/2, t/4 - u/2])] #sols.append([u, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j + k), QA(1/4*i + 1/2*j + 5/4*k), QA(j), QA(2*k)]): # 1st/Pinar's for p11 O2, Pinar's for p43 O2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Nx = 4*p*t^2 + 4*p*t*u + 5/4*p*u^2 + 5*p*t*v + 3*p*u*v + 29/16*p*v^2 + p*u*w + p*v*w + p*w^2 + 1/4*u^2 + 1/16*v^2. Trx = 0 => u = 0 => x = v/4*i + (v/2 + w)*j + (2*t + 5/4*v)*k => Nx = (29*p+1)/16*v^2 + p*v*w + 5*p*v*t + p*w^2 + 4*p*t^2.
		newbdt = 1/2*sqrt((25*p+1)/p) # bdt = 1/2*sqrt((29*p+1)/p)
		newbdw = sqrt((4*p+1)/p) #bdw = sqrt((5*p+1)/p)
		cst00 = (29*p + 1)/4
		cst0 = cst00*NN
		cst1 = (-25*P - p)/4
		cst2 = (-4*P - p)
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-29)*w^2*p^3 - w^2*p^2 + 116*N*p^3 + 33*N*p^2 + N*p >=0 => (116*p^2 + 33*p + 1)*N >= p*(29*p + 1)*w^2 => w^2 <= (4*p+1)*N/p	
			D1 = cst1*w^2
			D2 = 10*P*w
			for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: # Delta_vw = (-29)*t^2*p^3 - t^2*p^2 + 725/4*N*p^3 + 27/2*N*p^2 + 1/4*N*p >=0 => (29*p + 1)*(25*p + 1)/4*N >= p*(29*p + 1)*t^2 => t^2 <= (25*p + 1)/4*N/p
				Delta_v = cst2*t^2 + D1 + D2*t + cst0 
				if Delta_v >= 0:
					v1 = 2*((-w-5*t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-w-5*t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/4, v1/2 + w, 5*v1/4+2*t])] #sols.append([0, v, w, t]) 
					if v2 in ZZ:
						res += [QA([0, v2/4, v2/2 + w, 5*v2/4+2*t])] #sols.append([0, v, w, t]) 	
	elif set(O) == set([ QA(1/2 + 1/2*j + k), QA(1/4*i + 1/2*j + 3/4*k), QA(j), QA(2*k)]): # 2nd for p11 O2, QA(original/2nd for p43 O2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Nx = 4*p*t^2 + 4*p*t*u + 5/4*p*u^2 + 3*p*t*v + 2*p*u*v + 13/16*p*v^2 + p*u*w + p*v*w + p*w^2 + 1/4*u^2 + 1/16*v^2. Trx = 0 => u = 0 => x = v/4*i + (v/2 + w)*j + (2*t + 3/4*v)*k => Nx = (13*p+1)/16*v^2 + p*v*w + 3*p*v*t + p*w^2 + 4*p*t^2.
		newbdw = sqrt((4*p + 1)/p) # bdw = sqrt((5*p + 1)/p)
		newbdt = 1/2*sqrt((9*p + 1)/p) # bdt = 1/2*sqrt((13*p + 1)/p) 
		cst00 = (13*p+1)/4
		cst0 = cst00*NN
		cst1 = (-9*P-p)/4
		cst2 = (-4*P-p)
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-13)*w^2*p^3 - w^2*p^2 + 52*N*p^3 + 17*N*p^2 + N*p >= 0 => (52*p^2 + 17*p + 1)*N >= p*(13*p + 1)*w^2 => w^2 <= (4*p + 1)*N/p
			D1 = cst1*w^2
			D2 = 6*P*w
			for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: # Delta_vw = (-13)*t^2*p^3 - t^2*p^2 + 117/4*N*p^3 + 11/2*N*p^2 + 1/4*N*p >=0 => (117*p^2 + 22*p + 1)/4*N >= p*(13*p + 1)*t^2 => t^2 <= 1/4*(9*p + 1)*N/p
				Delta_v = cst2*t^2 + D1 + D2*t +  cst0
				if Delta_v >= 0:
					v1 = 2*((-w-3*t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-w-3*t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/4, v1/2 + w, 3*v1/4+2*t])] #sols.append([0, v, w, t]) 
					if v2 in ZZ:
						res += [QA([0, v2/4, v2/2 + w, 3*v2/4+2*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j + k), QA(1/4*i + 1/2*j + 1/4*k), QA(j), QA(2*k)]): # original/Pinar's for p19 O2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Nx = 4*p*t^2 + 4*p*t*u + 5/4*p*u^2 + p*t*v + p*u*v + 5/16*p*v^2 + p*u*w + p*v*w + p*w^2 + 1/4*u^2 + 1/16*v^2. Trx = 0 => u = 0 => x = v/4*i + (v/2 + w)*j + (2*t + v/4)*k => Nx = 4*p*t^2 + p*t*v + 5/16*p*v^2 + p*v*w + p*w^2 + 1/16*v^2.
		newbdw = sqrt((4*p + 1)/p) # bdw = sqrt((5*p + 1)/p)
		newbdt = 1/2*sqrt((p + 1)/p) # bdt = 1/2*sqrt((5*p + 1)/p) 
		cst00 = (5*p+1)/4
		cst0 = cst00*NN
		cst1 = (-4*P - p)
		cst2 = (-P -p)/4
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-5)*w^2*p^3 - w^2*p^2 + 20*N*p^3 + 9*N*p^2 + N*p >=0 => (20*p^2 + 9*p + 1)*N >= p*(5*p + 1)*w^2 => w^2 <= (4*p + 1)*N/p
			D1 = cst2*w^2
			D2 = 2*P*w
			for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: # Delta_vw = (-5)*t^2*p^3 - t^2*p^2 + 5/4*N*p^3 + 3/2*N*p^2 + 1/4*N*p >=0 =>  (5*p^2 + 6*p + 1)/4*N >= p*(5*p + 1)*t^2 => t^2 <= 1/4*(p + 1)*N/p
				Delta_v = D1 + D2*t + cst1*t^2 + cst0
				if Delta_v >= 0:
					v1 = 2*((-w-t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-w-t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/4, v1/2 + w, v1/4+2*t])] #sols.append([0, v, w, t])
					if v2 in ZZ:
						res += [QA([0, v2/4, v2/2 + w, v2/4+2*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j + k), QA(1/4*i + 1/2*j + 7/4*k), QA(j), QA(2*k)]): #2nd for p19 O2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]) => Trx = u. Nx = 4*p*t^2 + 4*p*t*u + 5/4*p*u^2 + 7*p*t*v + 4*p*u*v + 53/16*p*v^2 + p*u*w + p*v*w + p*w^2 + 1/4*u^2 + 1/16*v^2. Trx = 0 => u = 0 => x= v/4*i + (v/2 + w)*j + (2*t + 7/4*v)*k => Nx = 4*p*t^2 + 7*p*t*v + 53/16*p*v^2 + p*v*w + p*w^2 + 1/16*v^2.
		newbdw = sqrt((4*p + 1)/p) # bdw = sqrt((5*p + 1)/p)
		newbdt = 1/2*sqrt((49*p + 1)/p) # bdt = 1/2*sqrt((53*p + 1)/p) 
		cst00 = (53*p+1)/4
		cst0 = cst00*NN
		cst1 = (-4*P-p)
		cst2 = (-49*P-p)/4 #Doing loop for w first, then for t is faster.
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-53)*w^2*p^3 - w^2*p^2 + 212*N*p^3 + 57*N*p^2 + N*p >= 0 => (53*p + 1)*(4*p + 1)*N >= p*(53*p + 1)*w^2 => w^2 <= (4*p + 1)*N/p
			D1 = cst2*w^2
			D2 = 14*P*w
			for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: # Delta_vw = (-53)*t^2*p^3 - t^2*p^2 + 2597/4*N*p^3 + 51/2*N*p^2 + 1/4*N*p >= 0 => (53*p + 1)*(49*p + 1)/4*N >= p*(53*p + 1)*t^2 => t^2 <= 1/4*(49*p + 1)*N/p
				Delta_v = D1 + D2*t + cst1*t^2 + cst0
				if Delta_v >= 0:
					v1 = 2*((-w-7*t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-w-7*t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/4, v1/2 + w, 7*v1/4+2*t])] #sols.append([0, v, w, t]) 
					if v2 in ZZ:
						res += [QA([0, v2/4, v2/2 + w, 7*v2/4+2*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j), QA(1/4*i + 7/4*k), QA(j), QA(2*k)]): # original/Pinar's for p31 and p47 O2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]) => Trx = u. Nx = 4*p*t^2 + 1/4*p*u^2 + 7*p*t*v + 49/16*p*v^2 + p*u*w + p*w^2 + 1/4*u^2 + 1/16*v^2. Trx =0 => u = 0 => x = v/4*i + w*j + (2*t + 7/4*v)*k => Nx = 4*p*t^2 + 7*p*t*v + 49/16*p*v^2 + p*w^2 + 1/16*v^2
		newbdw = sqrt(1/p) # bdw = sqrt((p + 1)/p)
		cst00 = (49*p+1)/4
		cst0 = cst00*NN
		cst1 = (-49*P-p)/4
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-49)*w^2*p^3 - w^2*p^2 + 49*N*p^2 + N*p >= 0 => (49*p + 1)*N >= p*(49*p + 1)*w^2 => w^2 <= N/p
			D1 = cst1*w^2
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_vw = (-49)*t^2*p^3 - t^2*p^2 + 2401/4*N*p^3 + 49/2*N*p^2 + 1/4*N*p >= 0 => (49*p + 1)^2/4*N >= p*(49*p + 1)*t^2 => t^2 <= 1/4*(49*p + 1)*N/p = Same as bdt.
				Delta_v = D1 - t^2*p + cst0
				if Delta_v >= 0:
					v1 = 2*((-7*t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-7*t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/4, w, 7*v1/4+2*t])] #sols.append([0, v, w, t]) 
					if v2 in ZZ:
						res += [QA([0, v2/4, w, 7*v2/4+2*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j), QA(1/4*i + 1/4*k), QA(j), QA(2*k)]): # 2nd for p31, p47, p79 O2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Nx = 4*p*t^2 + 1/4*p*u^2 + p*t*v + 1/16*p*v^2 + p*u*w + p*w^2 + 1/4*u^2 + 1/16*v^2. Trx = 0 => u = 0 => x = v/4*i + w*j + (2*t + v/4)*k => Nx = 4*p*t^2 + p*t*v + 1/16*p*v^2 + p*w^2 + 1/16*v^2
		newbdw = sqrt(1/p) # bdw = sqrt((p + 1)/p)
		cst00 = (p+1)/4
		cst0 = cst00*NN
		cst1 = (-p^2-p)/4
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = -w^2*p^3 - w^2*p^2 + N*p^2 + N*p >= 0 =>  (p + 1)*N >= p*(p + 1)*w^2 => w^2 <= N/p
			D1 = cst1*w^2
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_vw =  -t^2*p^3 - t^2*p^2 + 1/4*N*p^3 + 1/2*N*p^2 + 1/4*N*p >= 0 => (p+1)^2/4*N >= p*(p + 1)*t^2 => t^2 <= 1/4*(p+1)*N/p. Same as bdt.
				Delta_v = D1 - t^2*p + cst0
				if Delta_v >= 0:
					v1 = 2*(-t*p + sqrt(Delta_v))/cst00
					v2 = 2*(-t*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/4, w, v1/4+2*t])] #sols.append([0, v, w, t]) 
					if v2 in ZZ:
						res += [QA([0, v2/4, w, v2/4+2*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j + k), QA(1/6*i + 2/3*j + 7/6*k), QA(j + 2*k), QA(3*k)]): # original for p31 O3. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Nx = 9*p*t^2 + 6*p*t*u + 5/4*p*u^2 + 7*p*t*v + 3*p*u*v + 65/36*p*v^2 + 12*p*t*w + 5*p*u*w + 6*p*v*w + 5*p*w^2 + 1/4*u^2 + 1/36*v^2. Trx =0 => u = 0 => x = v/6*i + (2/3*v + w)*j + (3*t + 7/6*v + 2*w)*k => Nx = 9*p*t^2 + 7*p*t*v + 65/36*p*v^2 + 12*p*t*w + 6*p*v*w + 5*p*w^2 + 1/36*v^2
		newbdw = sqrt((16*p+1)/p) # bdw = sqrt((17*p+1)/p)
		cst00 = (65*p+1)/9
		cst0 = cst00*NN
		cst1 = (-p^2-5*p)/9
		cst2 = (-8*p^2-4*p)/3
		cst3 = (-16*p^2-p)
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-260/9)*w^2*p^3 + (-4/9)*w^2*p^2 + 4160/9*N*p^3 + 36*N*p^2 + 4/9*N*p >= 0 => 4/9*N*(65*p + 1)*(16*p + 1)*p >= 4/9*(65*p + 1)*p^2*w^2 => w^2 <= (16*p + 1)*N/p
			D1 = cst1*w^2
			D2 = cst2*w
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_vw = (-260/9)*t^2*p^3 + (-4/9)*t^2*p^2 + 260/81*N*p^3 + 1304/81*N*p^2 + 20/81*N*p >= 0 => 4/81*N*(65*p + 1)*(p + 5)*p <= 4/9*(65*p + 1)*p^2*t^2 => t^2 <= 1/9*(p + 5)*N/p. Same as bdt.
				Delta_v = D1 + D2*t + cst3*t^2 + cst0 
				if Delta_v >= 0:
					v1 = 2*((-6*w-7*t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-6*w-7*t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/6, 2/3*v1 + w, 7*v1/6 + 2*w + 3*t])] #sols.append([0, v, w, t])
					if v2 in ZZ:
						res += [QA([0, v2/6, 2/3*v2 + w, 7*v2/6 + 2*w + 3*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j + 2*k), QA(1/6*i + 1/3*j + 13/6*k), QA(j + k), QA(3*k)]): # 2nd/Pinar's for p31 O3. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Nx = 9*p*t^2 + 12*p*t*u + 17/4*p*u^2 + 13*p*t*v + 9*p*u*v + 173/36*p*v^2 + 6*p*t*w + 5*p*u*w + 5*p*v*w + 2*p*w^2 + 1/4*u^2 + 1/36*v^2. Trx = 0 => u = 0 => x = v/6*i + (v/3 + w)*j + (3*t + 13/6*v + w)*k => Nx = 9*p*t^2 + 13*p*t*v + 173/36*p*v^2 + 6*p*t*w + 5*p*v*w + 2*p*w^2 + 1/36*v^2
		newbdw = sqrt((4*p + 1)/p) # bdw = sqrt((5*p + 1)/p)
		newbdt = 1/3*sqrt((121*p + 2)/p) # bdt = 1/3*sqrt((130*p + 2)/p) 
		cst00 = (173*p+1)/9
		cst0 = cst00*NN
		cst1 = (-121*p^2-2*p)/9
		cst2 = (44*p^2-2*p)/3
		cst3 = (-4*p^2-p)
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-692/9)*w^2*p^3 + (-4/9)*w^2*p^2 + 2768/9*N*p^3 + 236/3*N*p^2 + 4/9*N*p >= 0 =>  4/9*N*(173*p + 1)*(4*p + 1)*p <= 4/9*(173*p + 1)*p^2*w^2 => w^2 <= (4*p + 1)*N/p
			D1 = cst1*w^2
			D2 = cst2*w
			for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: # Delta_vw = (-692/9)*t^2*p^3 + (-4/9)*t^2*p^2 + 83732/81*N*p^3 + 1868/81*N*p^2 + 8/81*N*p >= 0 => 4/81*N*(173*p + 1)*(121*p + 2)*p <= 4/9*(173*p + 1)*p^2*t^2 => t^2 <= 1/9*(121*p + 2)*N/p
				Delta_v = D1 + D2*t + cst3*t^2 + cst0
				if Delta_v >= 0:
					v1 = 2*((-5*w-13*t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-5*w-13*t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/6, 1/3*v1 + w, 13*v1/6 + w + 3*t])] #sols.append([0, v, w, t]) 
					if v2 in ZZ:
						res += [QA([0, v2/6, 1/3*v2 + w, 13*v2/6 + w + 3*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1/2 + 1/2*j + 2*k), QA(1/6*i + 2/3*j + 5/6*k), QA(j + k), QA(3*k)]): # origina/2nd/Pinar's for p43 O3. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Nx = 9*p*t^2 + 12*p*t*u + 17/4*p*u^2 + 5*p*t*v + 4*p*u*v + 41/36*p*v^2 + 6*p*t*w + 5*p*u*w + 3*p*v*w + 2*p*w^2 + 1/4*u^2 + 1/36*v^2. Trx = 0 => u = 0 => x = v/6*i + (2/3*v + w)*j + (3*t + 5/6*v + w)*k => Nx = 9*p*t^2 + 5*p*t*v + 41/36*p*v^2 + 6*p*t*w + 3*p*v*w + 2*p*w^2 + 1/36*v^2
		newbdw = sqrt((16*p + 1)/p) # bdw = sqrt((17*p + 1)/p)
		newbdt = 1/3*sqrt((p + 2)/p) # bdt = 1/3*sqrt((10*p + 2)/p) 
		cst00 = (41*p+1)/9
		cst0 = cst00*NN
		cst1 = (-p^2-2*p)/9
		cst2 = (8*p^2-2*p)/3
		cst3 = (-16*p^2-p)
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vt = (-164/9)*w^2*p^3 + (-4/9)*w^2*p^2 + 2624/9*N*p^3 + 76/3*N*p^2 + 4/9*N*p >=0 => 4/9*N*(41*p + 1)*(16*p + 1)*p >= 4/9*(41*p + 1)*p^2*w^2 => w^2 <= (16*p + 1)*N/p
			D1 = cst1*w^2
			D2 = cst2*w
			for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: # Delta_vw = (-164/9)*t^2*p^3 + (-4/9)*t^2*p^2 + 164/81*N*p^3 + 332/81*N*p^2 + 8/81*N*p >=0 => 4/81*N*(41*p + 1)*(p + 2)*p >= 4/9*(41*p + 1)*p^2*t^2 => t^2 <= 1/9*(p + 2)*N/p
				Delta_v = D1 + D2*t + cst3*t^2 + cst0
				if Delta_v >= 0:
					v1 = 2*((-3*w-5*t)*p + sqrt(Delta_v))/cst00
					v2 = 2*((-3*w-5*t)*p - sqrt(Delta_v))/cst00
					if v1 in ZZ:
						res += [QA([0, v1/6, 2/3*v1 + w, 5*v1/6 + w + 3*t])] #sols.append([0, v, w, t]) 
					if v2 in ZZ:
						res += [QA([0, v2/6, 2/3*v2 + w, 5*v2/6 + w + 3*t])] #sols.append([0, v, w, t])
	elif set(O) == set([ QA(1), QA(i), QA(1/2 - 1/4*i + 1/4*k), QA(-1/2 + 1/2*i - 1/2*j)]): #p13. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = -t + 2*u + w. Nx = 1/4*p*t^2 + 1/8*p*w^2 + 3/4*t^2 - t*u + u^2 + 2*t*v + 2*v^2 - t*w + u*w - v*w + 3/8*w^2. Put t = 2*u + w. Nx = (p+2)*u^2 + 4*u*v + (p+1)*u*w + 2*v^2 + v*w + (3*p+1)/8*w^2
		newbdu = sqrt(3/p) #bdu = sqrt(abs((p + 5)/p))
		newbdw = sqrt(8/p) #bdw = 4/sqrt(abs(p))
		for u in [ceil(-newbdu*sqN)..floor(newbdu*sqN)]: #Delta_vw = (-32)*w^2*p^2 + 256*N*p >=0 => 256*N*p >= 32*w^2*p^2 => 8*N/p >= w^2
			D1 = (-4*p)*u^2 + 4*NN
			D2 = 4*p*u
			for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vu = (-32)*w^2*p^2 + 256*N*p >=0 => 256*N*p >= 32*w^2*p^2 => 8*N/p >= w^2
				Delta_v_over_2 = D1 - D2*w - 3*p/2*w^2 #norm_eq = 2*v^2 + (4*u + w)*v + u^2*p + u*w*p + 3/8*w^2*p + 2*u^2 + u*w + 1/8*w^2 - N
				if Delta_v_over_2 >= 0: 
					v1 = (-u + w/4) + sqrt(Delta_v_over_2)
					v2 = (-u + w/4) - sqrt(Delta_v_over_2)
					if v1 in ZZ: #x = (u + v + 1/4*w)*i + (-u - 1/2*w)*j + 1/4*w*k
						res += [QA([0, u + v1 + w/4, -u - w/2, w/4])] #sols.append([u, v, w, t]) #t = 2*u + w
					if v2 in ZZ: #x = (u + v + 1/4*w)*i + (-u - 1/2*w)*j + 1/4*w*k
						res += [QA([0, u + v2 + w/4, -u - w/2, w/4])] #sols.append([u, v, w, t]) #t = 2*u + w
	elif set(O) == set([QA(1/2 + 1/2*j), QA(1/4*i + 5/4*k), QA(j), QA(2*k)]): #p23 O2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (25*p+1)/16*p*v^2 + 5*p*v*t  + p*w^2 + 4*p*t^2
		newbdw = sqrt(1/p) #bdw = sqrt(abs((p + 1)/p) 
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-64)*w^2*p^3 + 64*N*p^2 >=0 => N/p >= w^2
			D1 = (25*p+1)/4*(NN-p*w^2)
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wv = = (-16)*t^2*p^3 + 100*N*p^3 + 4*N*p^2 >=0 => (25*p + 1)/(4*p)*N >= t^2. bdt = newbdt.
				Delta_v = D1 - t^2*p  #norm_eq = p*w^2 + 4*p*t^2 + 5*p*t*v + (25*p+1)/16*v^2 - N
				if Delta_v >= 0: 
					v1 = 8*(- 5*p*t + sqrt(Delta_v))/(25*p+1)
					v2 = 8*(- 5*p*t - sqrt(Delta_v))/(25*p+1)
					if v1 in ZZ: #x = 1/4*v*i + w*j + (2*t + 5/4*v)*k
						res += [QA([0, v1/4 , w, 2*t + 5*v1/4])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #1/4*v*i + w*j + (2*t + 5/4*v)*k
						res += [QA([0, v2/4 , w, 2*t + 5*v2/4])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([QA(1/2 + 3/2*j), QA(1/6*i + 7/3*j + 1/2*k), QA(3*j), QA(k)]): #p23 O3. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (205*p+1)/36*v^2 + 14*p*v*w + p*v*t + 9*p*w^2 + p*t^2 
		pp1 = (205*p + 1)/9
		pp2 = (196*p+1)*p/9
		newbdw = sqrt(pp2)/p #bdw = sqrt(abs((205*p + 1)/(9*p)))
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-16)*w^2*p^3 + 3136/9*N*p^3 + 16/9*N*p^2 >= 0 => (196*p + 1)/(9*p)*N >= w^2
			D1 = -(9*p+1)*p*w^2 + pp1*NN
			D2 = 28*p^2*w
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wv = (-144)*t^2*p^3 + 1296*N*p^3 + 144*N*p^2 >= 0 => (9*p + 1)/p*N >= t^2
				Delta_v = D1 + D2*t - pp2*t^2 #norm_eq = (205*p+1)/36*v^2 + (14*w + t)*p*v + 9*w^2*p + t^2*p - NN
				if Delta_v >= 0: 
					v1 = 2*(- (14*w + t)*p + sqrt(Delta_v))/pp1
					v2 = 2*(- (14*w + t)*p - sqrt(Delta_v))/pp1
					if v1 in ZZ: #x = 1/6*v*i + (7/3*v + 3*w)*j + (t + 1/2*v)*k
						res += [QA([0, v1/6, 7*v1/3 + 3*w, t + v1/2])] #sols.append([u, v, w, t]) #u = 0
					if v2 in ZZ: #x = 1/6*v*i + (7/3*v + 3*w)*j + (t + 1/2*v)*k
						res += [QA([0, v1/6, 7*v2/3 + 3*w, t + v2/2])] #sols.append([u, v, w, t]) #u = 0
	elif set(O) == set([QA(1/2 + 1/2*j), QA(1/6*i + 11/6*k), QA(j), QA(3*k)]): #p47 O3. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (121*p+1)/36*v^2 + 11*p*v*t + p*w^2 + 9*p*t^2
		newbdw = sqrt(1/p) #bdw = sqrt(abs((p + 1)/p)
		pp = (121*p+1)/9
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv =  (-144)*w^2*p^3 + 144*N*p^2 >=0 => 1/p*N >= w^2
			D1 = pp*(NN - p*w^2)
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wv = (-16)*t^2*p^3 + 1936/9*N*p^3 + 16/9*N*p^2 >=0 => (121*p + 1)/(9*p)*N >= t^2 bdt = newbdt.
				Delta_v = D1 - t^2*p #norm_eq = (121*p+1)/36*v^2 + 11*p*t*v + p*w^2 + 9*p*t^2 - NN
				if Delta_v >= 0: 
					v1 = 2*(-11*p*t + sqrt(Delta_v))/pp
					v2 = 2*(-11*p*t - sqrt(Delta_v))/pp
					if v1 in ZZ: #x = 1/6*v*i + w*j + (3*t + 11/6*v)*k
						res += [QA([0, v1/6 , w, 3*t + 11*v1/6])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #1/4*v*i + w*j + (2*t + 5/4*v)*k
						res += [QA([0, v2/6 , w, 3*t + 11*v2/6])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([QA(1/2 + 3/2*j), QA(1/6*i + 1/3*j + 1/2*k), QA(3*j), QA(k)]): #p47 O4. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (13*p+1)/36*v^2 + 2*p*v*w + p*v*t + 9*p*w^2 + p*t^2
		newbdw = sqrt((4*p + 1)/(9*p)) #bdw = sqrt((13*p + 1)/(9*p))
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-16)*w^2*p^3 + 64/9*N*p^3 + 16/9*N*p^2 >=0 => (4*p + 1)/(9*p)*N >= w^2
			D1 = -(9*p+1)*p*w^2 + (13*p+1)/9*NN 
			D2 = 4*w*p^2
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wv = (-144)*v^2*p^3 + 5184*N*p^3 >=0 => 36*N >= v^2. bdt = newbdt.
				Delta_v = D1 + D2*t - (4*p+1)*p/9*t^2  #norm_eq = (13*p+1)/36*v^2 + 2*p*w*v + p*v*t + 9*p*w^2 + p*t^2 - NN
				if Delta_v >= 0: 
					v1 = 18*(- 2*p*w + sqrt(Delta_v))/(13*p+1)
					v2 = 18*(- 2*p*w - sqrt(Delta_v))/(13*p+1)
					if v1 in ZZ: #x = 1/6*v*i + (1/3*v + 3*w)*j + (t + 1/2*v)*k
						res += [QA([0, v1/6 , v1/3 + 3*w, t + v1/2])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #1/6*v*i + (1/3*v + 3*w)*j + (t + 1/2*v)*k
						res += [QA([0, v2/6 , v2/3 + 3*w, t + v2/2])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([QA(1/2 + 7/2*j), QA(1/14*i + 27/7*j + 1/2*k), QA(7*j), QA(k)]): #p47 O5. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (2965*p+1)/196*p*v^2 + 54*p*v*w + p*v*t + 49*p*w^2 + p*t^2
		pp = (2916*p+1)*p/49
		newbdw = sqrt(pp)/p #bdw = sqrt((2965*p+1)/(49*p))
		pp1 = (2965*p + 1)/49
		pp2 = (49*p+1)*p
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-16)*w^2*p^3 + 46656/49*N*p^3 + 16/49*N*p^2 >=0 => (2916*p+1)/(49*p)*N >= w^2
			D1 = -pp2*w^2 + pp1*NN 
			D2 = 108*w*p^2
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wv = (-784)*t^2*p^3 + 38416*N*p^3 + 784*N*p^2 >=0 => (49*p + 1)/p*N >= t^2. bdt = newbdt.
				Delta_v = D1 + D2*t - pp*t^2  #norm_eq = (2965*p+1)/196*p*v^2 + (54*w + t)*p*v + 49*p*w^2 + p*t^2 - NN
				if Delta_v >= 0: 
					v1 = 2*p*(- (54*w + t)*p + sqrt(Delta_v))/pp1
					v2 = 2*p*(- (54*w + t)*p - sqrt(Delta_v))/pp1
					if v1 in ZZ: #x = 1/14*v*i + (27/7*v + 7*w)*j + (t + 1/2*v)*k
						res += [QA([0, v1/14 , 27/7*v2 + 7*w, t + v1/2])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #x = 1/14*v*i + (27/7*v + 7*w)*j + (t + 1/2*v)*k
						res += [QA([0, v2/14 , 27/7*v2 + 7*w, t + v2/2])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([QA(1/2 + 1/2*j + k), QA(1/6*i + 2/3*j + 13/6*k), QA(j + 2*k), QA(3*k)]): #p79 O3. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (185*p+1)/36*v^2 + 10*p*v*w + 13*p*v*t + 5*p*w^2 + 12*p*w*t + 9*p*t^2
		newbdw = sqrt((16*p + 1)/p) #bdw = sqrt((17*p + 1)/p)
		pp1 = -5*(5*p+1)*p/9
		pp2 = -(16*p+1)*p
		pp3 = 4*(10*p-1)*p/3*p
		pp4 = (185*p+1)/9*p
		for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: #Delta_wv = (-80)*t^2*p^3 + 2000/9*N*p^3 + 400/9*N*p^2 >=0 => 5*(5*p + 1)/(9*p)*N >= t^2. bdt = newbdt.
			D1 = pp2*t^2 + pp4*NN
			D2 = pp3*t
			for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_tv = (-144)*w^2*p^3 + 2304*N*p^3 + 144*N*p^2 >=0 => (16*p + 1)/p*N >= w^2
				Delta_v = pp1*w^2 + D2*w + D1  #norm_eq = (185*p+1)/36*v^2 + (10*w + 13*t)*p*v + 5*p*w^2 + 12*p*w*t + 9*p*t^2 - NN
				if Delta_v >= 0: 
					v1 = 18*(- (10*w + 13*t)*p + sqrt(Delta_v))/(185*p+1)
					v2 = 18*(- (10*w + 13*t)*p - sqrt(Delta_v))/(185*p+1)
					if v1 in ZZ: #x = 1/6*v*i + (2/3*v + w)*j + (3*t + 13/6*v + 2*w)*
						res += [QA([0, v1/6 , 2*v1/3 + w, 3*t + 13*v1/6 + 2*w])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #x = 1/6*v*i + (2/3*v + w)*j + (3*t + 13/6*v + 2*w)*
						res += [QA([0, v2/6 , 2*v2/3 + w, 3*t + 13*v2/6 + 2*w])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([QA(1/2 + 5/2*j), QA(1/10*i + 8/5*j + 1/2*k), QA(5*j), QA(k)]): #p79 O4. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (281*p+1)/100*p*v^2 + 16*p*v*w + p*v*t + 25*p*w^2 + p*t^2
		pp = (256*p+1)*p/25
		newbdw = sqrt(pp)/p #sqrt((281*p + 1)/(25*p))
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-16)*w^2*p^3 + 4096/25*N*p^3 + 16/25*N*p^2 >=0 => (256*p + 1)/(25*p)*N >= w^2
			D1 = (281*p+1)/25*NN -(25*p+1)*p*w^2 
			D2 = 32*p^2*w
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wv = (-400)*t^2*p^3 + 10000*N*p^3 + 400*N*p^2 >=0 => (25*p + 1)/p*N >= t^2. bdt = newbdt.
				Delta_v = D2*t - pp*t^2  #norm_eq = (281*p+1)/100*p*v^2 + (16*w + t)*p*v + 25*p*w^2 + p*t^2 - NN
				if Delta_v >= 0: 
					v1 = 50*(- (16*w + t)*p + sqrt(Delta_v))/(281*p+1)
					v2 = 50*(- (16*w + t)*p - sqrt(Delta_v))/(281*p+1)
					if v1 in ZZ: #x = 1/10*v*i + (8/5*v + 5*w)*j + (t + 1/2*v)*k
						res += [QA([0, v1/10 , 8*v1/5 + 5*w, t + v1/2])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #x = 1/10*v*i + (8/5*v + 5*w)*j + (t + 1/2*v)*k
						res += [QA([0, v2/10 , 8*v2/5 + 5*w, t + v2/2])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([QA(1/2 + 1/2*j), QA(1/10*i + 9/10*k), QA(j), QA(5*k)]): #p79 O5. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (81*p+1)/100*p*v^2 + 9*p*v*t + p*w^2 + 25*p*t^2
		newbdw = sqrt(1/p) #bdw = sqrt(abs((p + 1)/p))
		pp = (81*p+1)/25
		for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-400)*w^2*p^3 + 400*N*p^2 >=0 => 1/p*N >= w^2
			D1 = pp*(NN - p*w^2)
			for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: # Delta_wv = (-16)*t^2*p^3 + 1296/25*N*p^3 + 16/25*N*p^2 >=0 => (81*p + 1)/(25*p)*N >= t^2
				Delta_v = D1 - t^2*p  #norm_eq = (81*p+1)/100*p*v^2 + 9*p*t*v + p*w^2 + 25*p*t^2 - NN
				if Delta_v >= 0: 
					v1 = 2*p*(- 9*p*t + sqrt(Delta_v))/pp
					v2 = 2*p*(- 9*p*t - sqrt(Delta_v))/pp
					if v1 in ZZ: #x = 1/10*v*i + w*j + (5*t + 9/10*v)*k
						res += [QA([0, v1/10 , w, 5*t + 9*v1/10])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #x = 1/10*v*i + w*j + (5*t + 9/10*v)*k
						res += [QA([0, v2/10 , w, 5*t + 9*v2/10])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([QA(1/2 + 11/2*j), QA(1/22*i + 64/11*j + 1/2*k), QA(11*j), QA(k)]): #p79 O6. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = u. Put u = 0. Nx = (16505*p+1)/484*v^2 + 128*p*v*w + p*v*t + 121*p*w^2 + p*t^2 
		pp1 = (16384*p+1)*p/121
		pp2 = (16505*p+1)/121*NN
		pp3 = (121*p+1)*p
		newbdw = sqrt(pp1)/p #bdw = 1/11*sqrt(abs((16505*p + 1)/p))
		for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: #Delta_wv = (-1936)*t^2*p^3 + 234256*N*p^3 + 1936*N*p^2 >=0 => (121*p + 1)/p*N >= t^2. bdt = newbdt.
			D1 = - pp1*t^2 + pp2
			D2 = 256*p^2*t
			for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_tv = (-16)*w^2*p^3 + 262144/121*N*p^3 + 16/121*N*p^2 >=0 => (16384*p+1)/(121*p)*N >= w^2
				Delta_v = D1 - pp3*w^2 + D2*w   #norm_eq = (16505*p+1)/484*v^2 + (128*w + t)*p*v + 121*p*w^2 + p*t^2  - NN
				if Delta_v >= 0: 
					v1 = 2*p*(- (128*w + t)*p + sqrt(Delta_v))/pp1
					v2 = 2*p*(- (128*w + t)*p - sqrt(Delta_v))/pp1
					if v1 in ZZ: #x = 1/22*v*i + (64/11*v + 11*w)*j + (t + 1/2*v)*k
						res += [QA([0, v1/22 , 64/11*v1 + 11*w, t + v1/2 ])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #x = 1/22*v*i + (64/11*v + 11*w)*j + (t + 1/2*v)*k
						res += [QA([0, v2/22 , 64/11*v2 + 11*w, t + v2/2 ])] #sols.append([u, v, w, t]) #u=0
	elif set(O) == set([ QA(1), QA(1/2 + 1/2*i + 1/2*j + 1/2*k), QA(1/2 + 1/2*i - 1/2*j + 1/2*k), QA(1/2 - 1/2*i + 1/2*j + 1/2*k) ]): #p2. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Put t = - 2*u -v -w. Nx = 3*u^2 + 2*u*v + v^2 + 4*u*w + 2*v*w + 2*w^2.
		newbdw = sqrt(2)
		for u in [ceil(-sqN)..floor(sqN)]: #Delta_wv = (-128)*u^2 + 128*N >= 0 =>  N > u^2. #newbdu = 1
			D1 = (NN - 2*u^2)
			for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vu = (-64)*w^2 + 128*N >= 0 => 2*N > w^2
				Delta_v_over_4 = D1 - 2*u*w - w^2  #norm_eq = v^2 + 2*(u + 2*w + w)*v + 3*u^2 + 2*w^2 - NN
				if Delta_v_over_4 >= 0: 
					v1 = - (u + 2*w + w) + sqrt(Delta_v_over_4)
					v2 = - (u + 2*w + w) - sqrt(Delta_v_over_4)
					if v1 in ZZ: #x = (u + v + w)*i + (-u - w)*j + (-u)*k
						res += [QA([0, u + v1 + w , -u - w, -u])] #sols.append([u, v, w, t]) #u=0
					if v2 in ZZ: #x = (u + v + w)*i + (-u - w)*j + (-u)*k
						res += [QA([0, u + v2 + w , -u - w, -u])] #sols.append([u, v, w, t]) #u=0
	else:	
		for v in [ceil(-bdv*sqN)..floor(bdv*sqN)]:
			for w in [ceil(-bdw*sqN)..floor(bdw*sqN)]:
				for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]:
					if M[0][0] != 0:
						u = (- M[1][0]*v - M[2][0]*w - M[3][0]*t)/M[0][0]
						elem = QA(u*O[0] + v*O[1] + w*O[2] + t*O[3])
						if elem.reduced_norm() == NN:
							res += [QA(elem)]
	return list(set(res))


@parallel
def find_all_x_i(p, O, norms):
	all_x_i = {}
	for Nx_i in norms:
		if Nx_i%10 == 0:
			print ('Nx_i = ', Nx_i)
			sys.stdout.flush()
		all_x_i[Nx_i] = find_elem_with_norm_and_trace_zero(p,O,Nx_i)
		L = len(all_x_i[Nx_i])
		if L > 0:
			for i in [0..L-1]:
				csvwriter.writerow(([Nx_i , all_x_i[Nx_i][i]]))
				csvfile.flush()
		sys.stdout.flush()
	sys.stdout.flush()
	return all_x_i


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


print ("finding d_i...")
tx = time() 
all = find_all_x_i(p, O, all_norms)
csvfile.close()
print ("time for finding all x with all norms =", time() - tx)


L=0
for NN in all:
	L = L + len(all[NN])
print ("L=", L)


