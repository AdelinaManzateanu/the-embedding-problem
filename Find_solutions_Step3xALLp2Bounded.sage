# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements the arguments of Section 5.3 and Section 5.4 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step3xALLp2Bounded.sage CMfieldnumber Bound
# e.g. sage ./Find_solutions_Step3xALLp2Bounded.sage 34 1000 runs Step 3 for CMfield 34 and p < 1001.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 3: Find the bound for primes p up to which we may find solutions [x, dOVERn, a, cOVERn, b, gamma, n] and separate solutions from Step 2 by prime. 
# This step works for fields in Section 6.3.
# The bound for primes p is as described in Section 5.3.
# We use Proposition 4.1 in KLLLRS which proves that an abelian threefold with CM by the maximal order OK of a cyclic CM field has potential good reduction if p factors into 2 or 6 primes in OK. This way we discard some primes that we know are of good reduction.
#-------------------------------------------------------------------------------------------------------------------------------------------
from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv

######################################################
CMfieldnumber = sage_eval(sys.argv[1])
filelocation = "./"
T = time()

max_bdp = sage_eval(sys.argv[2])

load('CM_Fields.sage')

    
print ("CM field number ",CMfieldnumber,": ","A = ",A,", B = ",B,", C = ",C)
#print version()
print ("----------------------------------------------------------------------")
sys.stdout.flush()


print ("reading file with solutions from Step 2...")
sys.stdout.flush()
filename = "CMfield"+CMfieldnumber.str()+"Step2.csv"
reader=csv.DictReader(open(filelocation+filename))

@parallel
def find_all_bdp(CMfieldnumber):
    all_bounds = []
    solutions = []
    for row in reader:
        x = int(row["x"])
        dOVERn = int(row["dOVERn"])
        a = int(row["a"])
        cOVERn = int(row["cOVERn"])
        b = int(row["b"])
        gamma = int(row["gamma"])
        n = int(row["n"])
        bdp = 4*a*gamma*x^2
        solutions.append([x, dOVERn, a, cOVERn, b, gamma, n, bdp])
        print("[x, dOVERn, a, cOVERn, b, gamma, n, bdp] = ",[x, dOVERn, a, cOVERn, b, gamma, n, bdp])
        all_bounds.append(bdp)
        sys.stdout.flush()
    print ("----------------------------------------------------------------------")
    print (all_bounds)
    all_bounds = sorted(all_bounds)
    return all_bounds, solutions


all_bdp, sols_with_bdp = find_all_bdp(CMfieldnumber)

L = len(all_bdp)
print ("There are ", L, "solutions.")
print ("The largest bound for p is", all_bdp[L-1],".")
#max_bdp = all_bdp[L-1]


@parallel
def eliminate_potential_good_reduction_primes_in_general(CMfieldnumber,max_bdp):
    List_of_potential_bad_primes = []
    Bound = max_bdp 
    R.<x> = PolynomialRing(QQ)
    f = x^6 + A*x^4 + B*x^2 + C
    K.<eta> = NumberField(f)
    Delta = K.discriminant()
    f_disc = f.discriminant().factor()
    potential_singular_primes = [t[0] for t in f_disc]
    for p in primes(ceil(max_bdp)+1):
        if p in potential_singular_primes:
            if len(K.factor(p)) == 4 or len(K.factor(p)) == 5 or len(K.factor(p)) == 6:
                pass
            else:
                List_of_potential_bad_primes.append(p)
        else:   
            R.<t> = PolynomialRing(GF(p))
            f_over_finite = R(f)
            factors = f_over_finite.factor()
            if len(factors) == 4 or len(factors) == 5 or len(factors) == 6: 
                pass
            else:
                List_of_potential_bad_primes.append(p)
    return List_of_potential_bad_primes


@parallel
def find_solutions_for_each_p(CMfieldnumber): 
    #tm = time()
    all_bounds = all_bdp
    solutions = sols_with_bdp
    List_of_primes = eliminate_potential_good_reduction_primes_in_general(CMfieldnumber,max_bdp)
    print (List_of_primes)
    for p in List_of_primes:
        filename2 = "CMfield"+CMfieldnumber.str()+"Step3"+"p"+p.str()+".csv"
        csvfile = open(filelocation+filename2, 'w')
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n'))
        for S in solutions: #[x, dOVERn, a, cOVERn, b, gamma, n, bdp]
            if p < S[7] + 1:
                csvwriter.writerow([S[0],S[1],S[2],S[3],S[4],S[5],S[6]])
                csvfile.flush()
                sys.stdout.flush()
        csvfile.close()
        sys.stdout.flush()
    sys.stdout.flush()
    return solutions


bound_for_p = find_solutions_for_each_p(CMfieldnumber)


print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")