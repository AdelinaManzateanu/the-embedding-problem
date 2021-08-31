# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This step implements the arguments of Section 3.4.1, Section 5.3, and Section 5.4 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".


# run as sage ./Find_solutions_Step3ALLp2Bounded.sage CMfieldnumber Bound
# e.g. sage ./Find_solutions_Step3ALLp2Bounded.sage 1 100 runs Step 3 for CMfield 1 and p < 101.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 3: # Step 3: Find the bound for primes p up to which we may find solutions [x, dOVERn, a, cOVERn, b, gamma, n] and separate solutions from Step 2 by prime. 
# This step can only be performed for fields in Sections 6.1 and 6.2.
# This step cannot be performed for the cases with no imaginary quadratic subfield in Section 6.3.
# The extra automorphism condition as described in Section 3.4.1 leads to the condition that p must divide n for the solution to be valid.
# For the CM fields for which we do not have the extra automorphism condition, but there is an imaginary quadratic subfield, we use Proposition 4.1 in KLLLRS which proves that an abelian threefold with CM by the maximal order OK of a cyclic CM field has potential good reduction if p factors into 2 or 6 primes in OK. This way we discard some primes that we know are of good reduction.
#-------------------------------------------------------------------------------------------------------------------------------------------
# Files created:
# For each valid p up the specified Bound, a file "CMfield"+CMfieldnumber+"Step3"+"p"+p.csv" is created which contains solutions [x, dOVERn, a, cOVERn, b, gamma, n] corresponding to that prime, and the solutions come from Step 2.
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import *   # import sage library
from time import time
from random import randint
import csv

######################################################
CMfieldnumber = sage_eval(sys.argv[1])
filelocation = "./"
T = time()

load('CM_Fields.sage')
    
    
print ("CM field number ",CMfieldnumber,": ","A = ",A,", B = ",B,", C = ",C,", D = ",D)
#print version()
print ("----------------------------------------------------------------------")
sys.stdout.flush()


print ("reading file with solutions from Step 2...")
sys.stdout.flush()
filename = "CMfield"+CMfieldnumber.str()+"Step2.csv"
reader=csv.DictReader(open(filelocation+filename))


@parallel
def find_all_bdp(CMfieldnumber,D):
    all_bounds = []
    solutions = []
    lcm_n = 1
    for row in reader:
        x = int(row["x"])
        dOVERn = int(row["dOVERn"])
        a = int(row["a"])
        cOVERn = int(row["cOVERn"])
        b = int(row["b"])
        gamma = int(row["gamma"])
        n = int(row["n"])
        if D == 1: 
            if n > 1: # n >1 because we have p | n 
                solutions.append([x, dOVERn, a, cOVERn, b, gamma, n, n])
                print("[x, dOVERn, a, cOVERn, b, gamma, n, bdp] = ",[x, dOVERn, a, cOVERn, b, gamma, n, n])
                all_bounds.append(n)
                lcm_n = lcm(lcm_n,n)
        else: # D>1
            bdp = min(4*a*gamma*x^2,4*a*gamma*D^2)
            solutions.append([x, dOVERn, a, cOVERn, b, gamma, n, bdp])
            print("[x, dOVERn, a, cOVERn, b, gamma, n, bdp] = ",[x, dOVERn, a, cOVERn, b, gamma, n, bdp])
            all_bounds.append(bdp)
        sys.stdout.flush()
    print ("----------------------------------------------------------------------")
    print (all_bounds)
    all_bounds = sorted(all_bounds)
    return lcm_n, all_bounds, solutions


LCM_n, all_bdp, sols_with_bdp = find_all_bdp(CMfieldnumber,D)

    
L = len(all_bdp)
print ("There are ", L, "solutions.")
print ("The largest bound for p is", all_bdp[L-1],".")


max_bdp = sage_eval(sys.argv[2])
print ("Create files only for  p < ", max_bdp+1,".")


@parallel
def eliminate_potential_good_reduction_primes_cyclic_case(CMfieldnumber,D,max_bdp): 
    """
    This uses Proposition 4.1 in KLLLRS.

    - Returns True if the number of primes lying above p is not 2 or 6. 
    
    NOTE: Proposition 4.1 in KLLLRS proves that an abelian threefold with CM by the maximal order OK of 
    a cyclic CM field has potential good reduction if p factors into 2 or 6 primes in OK.
        
    
    """
    List_of_potential_bad_primes = []
    R.<x> = PolynomialRing(QQ)
    f = x^3 - A*x^2 + B*x - C
    Kplus.<mu> = NumberField(f)
    S.<t> = PolynomialRing(Kplus)
    g = t^2 + D
    Krel.<sqrtD> = NumberField(g)
    K = Krel.absolute_field('eta')
    f = K.polynomial()
    f_disc = f.discriminant().factor()
    potential_singular_primes = [t[0] for t in f_disc]
    for p in primes(ceil(max_bdp)+1):
        if p in potential_singular_primes:
            if len(K.factor(p)) == 2 or len(K.factor(p)) == 6:
                pass
            else:
                List_of_potential_bad_primes.append(p)
        else:   
            R.<t> = PolynomialRing(GF(p))
            f_over_finite = R(f)
            factors = f_over_finite.factor()
            if len(factors) == 6 or len(factors) == 2: 
                pass
            else:
                List_of_potential_bad_primes.append(p)
    return List_of_potential_bad_primes


# The function above only works for the CM fields in Section 6.1 and 6.2 (these are cyclic examples)


@parallel
def find_solutions_for_each_p(CMfieldnumber,D,max_bdp): 
    all_bounds = all_bdp
    solutions = sols_with_bdp
    lcm_n = LCM_n
    List_of_primes = eliminate_potential_good_reduction_primes_cyclic_case(CMfieldnumber, D, max_bdp)
    if D == 1:
        factors_lcm_n = []
        for f in list(factor(lcm_n)):
            factors_lcm_n.append(f[0])
        #print (lcm_n, factors_lcm_n)
        for p in factors_lcm_n:
            if p in List_of_primes:
                filename2 = "CMfield"+CMfieldnumber.str()+"Step3"+"p"+p.str()+".csv"
                csvfile = open(filelocation+filename2, 'w')
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n'))
                for S in solutions: #[x, dOVERn, a, cOVERn, b, gamma, n, bdp]
                    if S[6]/p in ZZ:
                        csvwriter.writerow([S[0],S[1],S[2],S[3],S[4],S[5],S[6]])
                        csvfile.flush()
                        sys.stdout.flush()
                csvfile.close()
                sys.stdout.flush()
        sys.stdout.flush()
    else: #if D > 1
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


bound_for_p = find_solutions_for_each_p(CMfieldnumber,D,max_bdp)


print ("total time: ", time()- T)
print ("----------------------------------------------------------------------")