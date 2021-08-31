# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This is part of Algorithm 1 in "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".

#------------------------------------------------------------------------------------------------------------------
# This document contains a function that finds elements in the maximal order OO in a certain quaternion algebra QA with a trace 0 and given norm NN.
# For optimisation reasons, several particular cases for the integral basis of the maximal order have been considered separately, but the function works for any integral basis of the maximal order. 
# The only requirement is that the basis is ordered such that the first element in the integral basis u0 + u1*i + u2*j + u3*k has u0 not equal to 0.

#------------------------------------------------------------------------------------------------------------------

@parallel
def find_elem_with_norm_and_trace_zero_short(p, q, NN, a00,a01,a02,a03,a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33):
    res = []
    if NN == 0:
        res += [QA([0, 0 , 0, 0])]
    else: #NN !=0
        MATRIX = matrix(QQ,4,[a00,a01,a02,a03,a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33])
        OO = [ QA(a00+a01*i+a02*j+a03*k), QA(a10+a11*i+a12*j+a13*k), QA(a20+a21*i+a22*j+a23*k), QA(a30+a31*i+a32*j+a33*k)]
        B2 = MATRIX.determinant()
        C = []
        for m in [0..3]:
            for n in [0..3]:
                C += [MATRIX.delete_rows([m]).delete_columns([n]).determinant()]
        C = matrix(4,C) #= [[C0_wtv, C1_wtv, C2_wtv, C3_wtv], [C0_wtu, C1_wtu, C2_wtu, C3_wtu], [C0_uvt, C1_uvt, C2_uvt, C3_uvt], [C0_uvw, C1_uvw, C2_uvw, C3_uvw]] # Ci_uvw = det of minor obtained by deleting row 3 (for t) anc column i = C[]
        B4_uvw = (q*C[3][0]^2 + C[3][1]^2)*p + (q*C[3][2]^2 + C[3][3]^2) #(q*C0_uvw^2 + C1_uvw^2)*p + (q*C2_uvw^2 + C3_uvw^2)
        B4_uvt = (q*C[2][0]^2 + C[2][1]^2)*p + (q*C[2][2]^2 + C[2][3]^2) #(q*C0_uvt^2 + C1_uvt^2)*p + (q*C2_uvt^2 + C3_uvt^2)
        B4_wtu = (q*C[1][0]^2 + C[1][1]^2)*p + (q*C[1][2]^2 + C[1][3]^2) #(q*C0_wtu^2 + C1_wtu^2)*p + (q*C2_wtu^2 + C3_wtu^2)
        B4_wtv = (q*C[0][0]^2 + C[0][1]^2)*p + (q*C[0][2]^2 + C[0][3]^2) #(q*C0_wtv^2 + C1_wtv^2)*p + (q*C2_wtv^2 + C3_wtv^2)
        bdu = sqrt(abs(B4_wtv)/(q*p))/abs(B2)
        bdv = sqrt(abs(B4_wtu)/(q*p))/abs(B2)
        bdw = sqrt(abs(B4_uvt)/(q*p))/abs(B2)
        bdt = sqrt(abs(B4_uvw)/(q*p))/abs(B2)
        sqN = sqrt(NN)
        P = p^2
        if a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a10 == 1/2 and a11 == 1/2 and a12 == 1/2 and a13 == 1/2 and a20 == 1/2 and a21 == 1/2 and a22 == -1/2 and a23 == 1/2 and a30 == 1/2 and a31 == -1/2 and a32 == 1/2 and a33 == 1/2 and p == 2: #p2. x = QA(u*OO[0]+ v*OO[1]+w*OO[2]+ t*OO[3]). Put t = - 2*u -v -w. Nx = 3*u^2 + 2*u*v + v^2 + 4*u*w + 2*v*w + 2*w^2.
            print ('Case 1')
            newbdw = sqrt(2)
            for u in [ceil(-sqN)..floor(sqN)]: #Delta_wv = (-128)*u^2 + 128*N >= 0 =>  N > u^2. #newbdu = 1
                D0 = - 2*u^2 + NN
                D1 = 2*u
                for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: # Delta_vu = (-64)*w^2 + 128*N >= 0 => 2*N > w^2
                    Delta_v_over_4 = - w^2 - D1*w + D0 #norm_eq = v^2 + 2*(u + w)*v + 4*u*w + 3*u^2 + 2*w^2 - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = - (u + w) + sqrt(Delta_v_over_4)
                        v2 = - (u + w) - sqrt(Delta_v_over_4)
                        if v1 in ZZ: #x = (u + v + w)*i + (-u - w)*j + (-u)*k
                            res += [QA([0, u + v1 + w , -u - w, -u])] #sols.append([u, v, w, t]) #u=0
                        if v2 in ZZ: #x = (u + v + w)*i + (-u - w)*j + (-u)*k
                            res += [QA([0, u + v2 + w , -u - w, -u])] #sols.append([u, v, w, t]) #u=0
        elif a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a10 == 0 and a11 == 1 and a12 == 0 and a13 == 0 and a20 == 0 and a21 != 0 and a22 == 0 and a23 != 0 and a30 != 0 and a31 == 0 and a32 != 0 and a33 == 0: #O = [ QA(1), QA(i), QA(a21*i + a23*k), QA(a30 + a32*j) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a30*t + 2*u. Put t = -u/a30. Nx = a23^2*p*q*w^2 + a21^2*q*w^2 + 2*a21*q*v*w + a32^2*p*u^2/a30^2 + q*v^2
            print ('Case 2')
            newbdu = sqrt(a30^2/(a32^2*p)) #< sqrt((p + (a30/a32)^2)/p) = bdu
            pp1 = a23^2*p*q^2
            D0 = NN*q
            pp2 = a32^2*p*q/a30^2
            pp3 = 1/(a23^2*p*q)
            pp4 = a32^2/(a30^2*a23^2*q*NN)
            for u in [ceil(-newbdu*sqN)..floor(newbdu*sqN)]: #Delta_wv = -64*(a23^2*p + a21^2)*a23^2*a32^2*p^2*q^3*u^2/a30^2 + 64*(a23^2*p + a21^2)*NN*a23^2*p*q^3 >= 0 => a30^2/(a32^2*p)*NN >= u^2
                newbdw = sqrt(pp3 - pp4*u^2) # Delta_v = -4*a23^2*p*q^2*w^2 - 4*a32^2*p*q*u^2/a30^2 + 4*NN*q >= 0
                D1 = D0 - pp2*u^2
                for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_uv = -64*a23^2*a32^4*p^3*q^2*w^2/a30^4 + 64*NN*a32^4*p^2*q/a30^4 >= 0 => (1/a23^2*p*q)*NN >= w^2 => newbdw = sqrt(1/(a23^2*p*q)) = bdw
                    Delta_v_over_4 = - pp1*w^2 + D1 #norm_eq = q*v^2 + 2*a21*q*w*v + a23^2*p*q*w^2 + a21^2*q*w^2 + a32^2*p*u^2/a30^2 - NN
                    if Delta_v_over_4 >=0:
                        v1 = ( - a21*q*w + sqrt(Delta_v_over_4))/q
                        v2 = ( - a21*q*w - sqrt(Delta_v_over_4))/q
                        if v1 in ZZ: #x = (a21*w + v)*i + (-a32*u/a30)*j + a23*w*k
                            res += [QA([0, (a21*w + v1),  -a32*u/a30, a23*w])] #sols.append([u, v1, w, t]) #t = -u/a30
                        if v2 in ZZ: #x = 
                            res += [QA([0, (a21*w + v2),  -a32*u/a30, a23*w])] #sols.append([u, v2, w, t]) #t = -u/a30
        elif a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a10 != 0 and a11 != 0 and a12 == 0 and a13 != 0 and a20 != 0 and a21 != 0 and a22 == 0 and a23 != 0 and a30 != 0 and a31 != 0 and a32 != 0 and a33 == 0 and (a13*a21 - a11*a23) !=0 : #O = [ QA(1), QA(a10 + a11*i + a13*k), QA(a20 + a21*i + a23*k), QA(a30 + a31*i + a32*j) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a30*t + 2*a10*v + 2*a20*w + 2*u. Put v = - u/a10 - a20/a10*w - a30/a10*t. Nx = a13^2*a30^2*p*q*t^2/a10^2 + 2*a13^2*a20*a30*p*q*t*w/a10^2 - 2*a13*a23*a30*p*q*t*w/a10 + a13^2*a20^2*p*q*w^2/a10^2 - 2*a13*a20*a23*p*q*w^2/a10 + a23^2*p*q*w^2 + a32^2*p*t^2 + a11^2*a30^2*q*t^2/a10^2 - 2*a11*a30*a31*q*t^2/a10 + a31^2*q*t^2 + 2*a13^2*a30*p*q*t*u/a10^2 + 2*a11^2*a20*a30*q*t*w/a10^2 - 2*a11*a21*a30*q*t*w/a10 - 2*a11*a20*a31*q*t*w/a10 + 2*a21*a31*q*t*w + 2*a13^2*a20*p*q*u*w/a10^2 - 2*a13*a23*p*q*u*w/a10 + a11^2*a20^2*q*w^2/a10^2 - 2*a11*a20*a21*q*w^2/a10 + a21^2*q*w^2 + 2*a11^2*a30*q*t*u/a10^2 - 2*a11*a31*q*t*u/a10 + a13^2*p*q*u^2/a10^2 + 2*a11^2*a20*q*u*w/a10^2 - 2*a11*a21*q*u*w/a10 + a11^2*q*u^2/a10^2
            print ('Case 3')
            pp0 = (a13*a21 - a11*a23)^2
            pp1 = ((a13*a20 - a10*a23)*a31 - (a13*a21 - a11*a23)*a30)^2*q + (a11*a20 - a10*a21)^2*a32^2 + (a13*a20 - a10*a23)^2*a32^2*p
            pp2 = p*q*pp0*a32^2
            newbdu = sqrt( pp1/pp2 ) 
            pp3 = ((a13*a20 - a10*a23)^2*a32^2*p + (a13*a21 - a11*a23)^2*a30^2*q -2*(a13*a20 - a10*a23)*(a13*a21 - a11*a23)*a30*a31*q + (a13*a20 - a10*a23)^2*a31^2*q + (a11*a20 - a10*a21)^2*a32^2)*p*q/a10^2
            pp4 = 2*(a13*a21*a30 - a11*a23*a30 - a13*a20*a31 + a10*a23*a31)*(a13*a21 - a11*a23)*p*q^2/a10^2
            pp5 = pp0*p*q^2/a10^2
            D0 = (a13^2*a20^2*p - 2*a10*a13*a20*a23*p + a10^2*a23^2*p + a11^2*a20^2 - 2*a10*a11*a20*a21 + a10^2*a21^2)*NN*q/a10^2
            pp6 = (a13^2*a20*a30*p - a10*a13*a23*a30*p + a11^2*a20*a30 - a10*a11*a21*a30 - a10*a11*a20*a31 + a10^2*a21*a31)*q/a10^2
            denom = ((a13*a20 - a10*a23)^2*p + (a11*a20 - a10*a21)^2)*q/a10^2
            pp7 = (a13^2*a20*p - a10*a13*a23*p + a11^2*a20 - a10*a11*a21)*q/a10^2
            for u in [ceil(-newbdu*sqN)..floor(newbdu*sqN)]: #Delta_tw = -64*(a13^2*a30^2*p*q + a10^2*a32^2*p + a11^2*a30^2*q - 2*a10*a11*a30*a31*q + a10^2*a31^2*q)*(a13*a21 - a11*a23)^2*a32^2*p^2*q^2*u^2/a10^4
                D1 = - pp5*u^2 + D0
                D2 = pp4*u
                for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: #Delta_uw = -64*(a13^2*p + a11^2)*(a13*a21 - a11*a23)^2*a32^2*p^2*q^3*t^2/a10^4 + 64*(a13^2*p + a11^2)*(a13*a21 - a11*a23)^2*NN*p*q^3/a10^4 >= 0
                    Delta_w_over_4 = - pp3*t^2 - D2*t + D1
                    #norm_eq = + ((a13*a20 - a10*a23)^2*p + (a11*a20 - a10*a21)^2)*q*w^2/a10^2  
                    # + 2*(a13^2*a20*a30*p*t - a10*a13*a23*a30*p*t + a11^2*a20*a30*t - a10*a11*a21*a30*t - a10*a11*a20*a31*t + a10^2*a21*a31*t + a13^2*a20*p*u - a10*a13*a23*p*u + a11^2*a20*u - a10*a11*a21*u)*q*w/a10^2 
                    #+ 2*(a13^2*a30*p + a11^2*a30 - a10*a11*a31)*q*t*u/a10^2 + (a13^2*a30^2*p*q + a10^2*a32^2*p + a11^2*a30^2*q - 2*a10*a11*a30*a31*q + a10^2*a31^2*q)*t^2/a10^2 + (a13^2*p + a11^2)*q/a10^2*u^2 - NN
                    if Delta_w_over_4 >= 0: 
                        w1 = (- (pp6*t + pp7*u) + sqrt(Delta_w_over_4))/denom
                        w2 = (- (pp6*t + pp7*u) - sqrt(Delta_w_over_4))/denom
                        if w1 in ZZ: #x = (a31*t + a21*w - (a30*t + a20*w + u)*a11/a10)*i + a32*t*j + (a23*w - (a30*t + a20*w + u)*a13/a10)*k
                            res += [QA([0, a31*t + a21*w1 - (a30*t + a20*w1 + u)*a11/a10,  a32*t, a23*w1 - (a30*t + a20*w1 + u)*a13/a10])] #sols.append([u, v, w1, t]) #v = - u/a10 - a20/a10*w - a30/a10*t
                        if w2 in ZZ: #x = (a31*t + a21*w - (a30*t + a20*w + u)*a11/a10)*i + a32*t*j + (a23*w - (a30*t + a20*w + u)*a13/a10)*k
                            res += [QA([0, a31*t + a21*w1 - (a30*t + a20*w2 + u)*a11/a10,  a32*t, a23*w2 - (a30*t + a20*w2 + u)*a13/a10])] #sols.append([u, v, w2, t]) #v = - u/a10 - a20/a10*w - a30/a10*t
        elif a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a10 != 0 and a11 != 0 and a12 == 0 and a13 != 0 and a20 != 0 and a21 != 0 and a22 != 0 and a23 != 0 and a30 != 0 and a31 != 0 and a32 != 0 and a33 != 0 and (-a13*a22*a31 + a13*a21*a32 - a11*a23*a32 + a11*a22*a33) !=0 : #O = [ QA(1), QA(a10 + a11*i + a13*k), QA(a20 + a21*i + a22*j + a23*k), QA(a30 + a31*i + a32*j + a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a30*t + 2*a10*v + 2*a20*w + 2*u. Put v = - u/a10 - a20/a10*w - a30/a10*t. Nx = a13^2*a30^2*p*q*t^2/a10^2 - 2*a13*a30*a33*p*q*t^2/a10 + a33^2*p*q*t^2 + 2*a13^2*a20*a30*p*q*t*w/a10^2 - 2*a13*a23*a30*p*q*t*w/a10 - 2*a13*a20*a33*p*q*t*w/a10 + 2*a23*a33*p*q*t*w + a13^2*a20^2*p*q*w^2/a10^2 - 2*a13*a20*a23*p*q*w^2/a10 + a23^2*p*q*w^2 + a32^2*p*t^2 + a11^2*a30^2*q*t^2/a10^2 - 2*a11*a30*a31*q*t^2/a10 + a31^2*q*t^2 + 2*a13^2*a30*p*q*t*u/a10^2 - 2*a13*a33*p*q*t*u/a10 + 2*a22*a32*p*t*w + 2*a11^2*a20*a30*q*t*w/a10^2 - 2*a11*a21*a30*q*t*w/a10 - 2*a11*a20*a31*q*t*w/a10 + 2*a21*a31*q*t*w + 2*a13^2*a20*p*q*u*w/a10^2 - 2*a13*a23*p*q*u*w/a10 + a22^2*p*w^2 + a11^2*a20^2*q*w^2/a10^2 - 2*a11*a20*a21*q*w^2/a10 + a21^2*q*w^2 + 2*a11^2*a30*q*t*u/a10^2 - 2*a11*a31*q*t*u/a10 + a13^2*p*q*u^2/a10^2 + 2*a11^2*a20*q*u*w/a10^2 - 2*a11*a21*q*u*w/a10 + a11^2*q*u^2/a10^2
            print ('Case 4')
            pp0 = (a13*a21 - a11*a23)*a30 - (a13*a31 - a11*a33)*a20 + (a23*a31 - a21*a33)*a10 
            pp1 = a11*a22*a30 - a11*a20*a32 - (a22*a31 - a21*a32)*a10
            pp2 = a13*a22*a30 - a13*a20*a32 + (a23*a32 - a22*a33)*a10
            pp3 = a13*a22*a31 - a13*a21*a32 + a11*a23*a32 - a11*a22*a33
            pp4 = pp0^2*q + pp1^2 + pp2^2*p
            newbdu = sqrt(pp4/(p*q))/abs(pp3)
            denom = 2*((a13*a20 - a10*a23)^2*p*q + a10^2*a22^2*p + (a11*a20 - a10*a21)^2*q)
            D0 = 2*demon*NN/a10^2 
            D1 = 4*(a13^2*a22^2*p + (a13*a21 - a11*a23)^2*q + a11^2*a22^2)*p*q/a10^2
            D2 = 8*(pp2*a13*a22*p + pp0*(a13*a21 - a11*a23)*q + pp1*a11*a22)*p*q/a10^2 
            D3 = 4*pp4*p*q/a10^2
            pp5 = ((a13*a20 - a10*a23)*(a13*a30 - a10*a33)*p*q + a10^2*a22*a32*p + (a11*a20 - a10*a21)*(a11*a30 - a10*a31)*q)/a10^2
            pp6 = ((a13*a20 - a10*a23)*a13*p + (a11*a20 - a10*a21)*a11)*q/a10^2
            for u in [ceil(-newbdu*sqN)..floor(newbdu*sqN)]: #Delta_wt = 64*(((a23*a31 - a21*a33)*a10 - (a13*a31 - a11*a33)*a20 + (a13*a21 - a11*a23)*a30)^2*q + (a11*a22*a30 - a11*a20*a32 - (a22*a31 - a21*a32)*a10)^2 + (a13*a22*a30 - a13*a20*a32 + (a23*a32 - a22*a33)*a10)^2*p)*(a13^2*a20^2*p*q - 2*a10*a13*a20*a23*p*q + a10^2*a23^2*p*q + a10^2*a22^2*p + a11^2*a20^2*q - 2*a10*a11*a20*a21*q + a10^2*a21^2*q)*NN*p*q/a10^4 -64*(a13*a22*a31 - a13*a21*a32 + a11*a23*a32 - a11*a22*a33)^2*(a13^2*a20^2*p*q - 2*a10*a13*a20*a23*p*q + a10^2*a23^2*p*q + a10^2*a22^2*p + a11^2*a20^2*q - 2*a10*a11*a20*a21*q + a10^2*a21^2*q)*p^2*q^2/a10^4*u^2
                DD1 = D1*u^2
                DD2 = D2*u
                for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: #Delta_uw = 64*(a13^2*a22^2*p + a11^2*a22^2 + (a13*a21 - a11*a23)^2*q)*(a13^2*p + a11^2)*NN*p*q^2/a10^4 - 64*(a13*a22*a31 - a13*a21*a32 + a11*a23*a32 - a11*a22*a33)^2*(a13^2*p + a11^2)*p^2*q^3/a10^4*t^2 >= 0
                    Delta_w = - D3*t^2 - DD2*t - DD1 + D0
                    #norm_eq = + ((a13*a20 - a10*a23)^2*p*q + a10^2*a22^2*p + (a11*a20 - a10*a21)^2*q)/a10^2*w^2 
                    # + 2*((a13*a20 - a10*a23)*(a13*a30 - a10*a33)*p*q + a10^2*a22*a32*p + (a11*a20 - a10*a21)*(a11*a30 - a10*a31)*q)/a10^2*t*w 
                    # + 2*((a13*a20 - a10*a23)*a13*p + (a11*a20 - a10*a21)*a11)*q/a10^2*u*w 
                    # + ((a13*a30 - a10*a33)^2*p*q + a10^2*a32^2*p + (a11*a30 - a10*a31)^2*q)*t^2/a10^2 + (a13^2*p + a11^2)*q/a10^2*u^2 + 2*((a13*a30 - a10*a33)*a13*p + (a11*a30 - a10*a31)*a11)*q/a10^2*u*t - NN
                    if Delta_w >= 0: 
                        w1 = a10^2*(- 2*(pp5*t + pp6*u) + sqrt(Delta_w))/denom
                        w2 = a10^2*(- 2*(pp5*t + pp6*u) - sqrt(Delta_w))/denom
                        if w1 in ZZ: #x = (a31*t + a21*w - (a30*t + a20*w + u)*a11/a10)*i + (a32*t + a22*w)*j + (a33*t + a23*w - (a30*t + a20*w + u)*a13/a10)*k
                            res += [QA([0, a31*t + a21*w1 - (a30*t + a20*w1 + u)*a11/a10,  a32*t + a22*w1, a33*t + a23*w1 - (a30*t + a20*w1 + u)*a13/a10])] #sols.append([u, v, w1, t]) #v = - u/a10 - a20/a10*w - a30/a10*t
                        if w2 in ZZ: #(a31*t + a21*w - (a30*t + a20*w + u)*a11/a10)*i + (a32*t + a22*w)*j + (a33*t + a23*w - (a30*t + a20*w + u)*a13/a10)*k
                            res += [QA([0, a31*t + a21*w1 - (a30*t + a20*w2 + u)*a11/a10,  a32*t + a22*w2, a33*t + a23*w2 - (a30*t + a20*w2 + u)*a13/a10])] #sols.append([u, v, w2, t]) #v = - u/a10 - a20/a10*w - a30/a10*t
        elif a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a10 == 0 and a11 == 1 and a12 == 0 and a13 == 0 and a21 != 0 and a22 == 0 and a23 != 0 and a31 != 0 and a32 != 0 and a33 == 0: #O[ QA(1), QA(i), QA(a20 + a21*i + a23*k), QA(a30 + a31*i + a32*j) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a30*t + 2*a20*w + 2*u. Put u = - a30*t - a20*w. Nx = a23^2*p*q*w^2 + a32^2*p*t^2 + a31^2*q*t^2 + 2*a21*a31*q*t*w + a21^2*q*w^2 + 2*a31*q*t*v + 2*a21*q*v*w + q*v^2. a20 and a30 may be 0.
            print ('Case 5')
            D0 = q*NN
            pp0 = a23^2*p*q
            pp1 = a32^2*p
            pp2 = pp0*q
            pp3 = pp1*q
            newbdw = sqrt(1/pp0)
            newbdt = sqrt(1/pp1)
            for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_vt = (-64)*w^2*a23^2*a32^2*p^2*q^3 + 64*a32^2*p*q^2*NN >= 0 => 1/(a23^2*p*q)*NN >= w^2 => newbdw = bdw
                D1 = - pp2*w^2 + D0
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = (-64)*t^2*a23^2*a32^2*p^2*q^3 + 64*a23^2*p*q^3*NN >= 0 => 1/(a32^2*p)*NN >= t^2 => newbdt = bdt
                    Delta_v_over_4 = - pp3*t^2 + D1 #norm_eq = q*v^2 + 2*q*(w*a21 + t*a31)*v + q*(a23^2*p + a21^2)*w^2 + 2*a21*a31*q*w*t + (a32^2*p + a31^2*q)*t^2 - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- q*(w*a21 + t*a31) + sqrt(Delta_v_over_4))/q
                        v2 = (- q*(w*a21 + t*a31) - sqrt(Delta_v_over_4))/q
                        if v1 in ZZ: #x = (a31*t + a21*w + v)*i + a32*t*j + a23*w*k
                            res += [QA([0, a31*t + a21*w + v1, a32*t, a23*w])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = (a31*t + a21*w + v)*i + a32*t*j + a23*w*k
                            res += [QA([0, a31*t + a21*w + v2, a32*t, a23*w])] #sols.append([u, v2, w, t]) #u=0
        elif a00 != 0 and a01 == 0 and a03 == 0 and a10 == 0 and a11 != 0 and a12 == 0 and a13 != 0 and a20 == 0 and a21 == 0 and a22 != 0 and a23 == 0 and a30 == 0 and a31 == 0 and a32 == 0 and a33 != 0 : #O = [ QA(a00 + a02*j), QA(a11*i + a13*k), QA(a22*j), QA(a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a13*a33*p*q*t*v + a13^2*p*q*v^2 + a11^2*q*v^2 + a22^2*p*w^2. a02 may be 0
            print ('Case 6')
            pp6 = a13^2*p + a11^2
            pp0 = pp6*q
            pp1 = a13*a33*p*q
            pp2 = a22^2*p
            pp3 = a11^2*a33^2*p*q
            newbdw = sqrt(1/pp2) #< bdw = sqrt((a02^2*p + a00^2)/(a00^2*a22^2*p))
            newbdt = sqrt(pp6/pp3) #= bdt = (a13^2*p + a11^2)/(p*q*a33^2**a11^2)
            D0 = pp0*NN
            pp4 = pp0*pp2
            pp5 = pp3*q
            for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = -64*w^2*a11^2*a22^2*a33^4*p^3*q^3 + 64*a11^2*a33^4*p^2*q^3*NN >= 0 => 1/(a22^2*p)*NN >= w^2
                D1 = - pp4*w^2 + D0
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = -64*p^2*a33^2*a22^2*a11^2*q^3*(a13^2*p + a11^2)*t^2 + 64*NN*p*q^2*a22^2*(a13^2*p + a11^2)^2 >= 0 => newbdt = bdt = sqrt((a13^2*p + a11^2)/(a11^2*a33^2*p*q))
                    Delta_v_over_4 = - pp5*t^2 + D1 #norm_eq = q*(a13^2*p + a11^2)*v^2 + 2*a13*a33*p*q*t*v + t^2*a33^2*p*q + w^2*a22^2*p - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- pp1*t + sqrt(Delta_v_over_4))/pp0
                        v2 = (- pp1*t - sqrt(Delta_v_over_4))/pp0
                        if v1 in ZZ: #x = a11*v*i + a22*w*j + (a33*t + a13*v)*k
                            res += [QA([0, a11*v1, a22*w, a33*t + a13*v1])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = a11*v*i + a22*w*j + (a33*t + a13*v)*k
                            res += [QA([0, a11*v2, a22*w, a33*t + a13*v2])] #sols.append([u, v2, w, t]) #u=0
        elif a00 != 0 and a10 == 0 and a11 != 0 and a12 == 0 and a13 != 0 and a20 == 0 and a21 == 0 and a22 != 0 and a23 != 0 and a30 == 0 and a31 == 0 and a32 == 0 and a33 != 0 : #O = [ QA(a00 + a01*i + a02*j + a03*k), QA(a11*i + a13*k), QA(a22*j + a23*k), QA(a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a13*a33*p*q*t*v + a13^2*p*q*v^2 + 2*a23*a33*p*q*t*w + 2*a13*a23*p*q*v*w + a23^2*p*q*w^2 + a11^2*q*v^2 + a22^2*p*w^2. a01, a02, a03 may be 0
            print ('Case 7')
            pp0 = a22^2*p
            pp1 = a13^2*p + a11^2
            pp2 = pp1*a22^2 + a11^2*a23^2*q
            pp3 = pp0*q*a33^2*a11^2
            newbdt = sqrt(pp2/pp3) #< bdw = sqrt((a02^2*p + a00^2)/(a00^2*a22^2*p))
            newbdw = sqrt(1/pp0) #< bdt = sqrt((a02^2*p + a00^2)/(a00^2*a22^2*p))
            D0 =  pp1*NN*q
            pp4 = pp2*p*q
            pp5 = 2*a11^2*a23*a33*p*q^2
            pp6 =  a11^2*a33^2*p*q^2
            pp7 = a33*a13*p*q
            pp8 = a23*a13*p*q
            for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = -64*w^2*a11^2*a22^2*a33^4*p^3*q^3 + 64*a11^2*a33^4*p^2*q^3*NN >= 0 => 1/(a22^2*p)*NN >= w^2
                D1 = - pp2*p*q*w^2 + D0
                D2 = pp5*w
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = -64*p^2*a33^2*a22^2*a11^2*t^2*q^3*(a13^2*p + a11^2) + 64*NN*p*q^2*(a13^2*p + a11^2)*(a13^2*a22^2*p + a11^2*a23^2*q + a11^2*a22^2) >= 0 =>  (a13^2*a22^2*p + a11^2*a23^2*q + a11^2*a22^2)/(p*q*a33^2*a22^2*a11^2)*NN > = t^2
                    Delta_v_over_4 = - pp6*t^2 - D2*t + D1 #norm_eq = (a13^2*p + a11^2)*q*v^2 + 2*(a33*t + a23*w)*a13*p*q*v + a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a22^2*p*w^2 - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- (pp7*t + pp8*w) + sqrt(Delta_v_over_4))/(pp1*q)
                        v2 = (- (pp7*t + pp8*w) - sqrt(Delta_v_over_4))/(pp1*q)
                        if v1 in ZZ: #x = a11*v*i + a22*w*j + (a33*t + a13*v + a23*w)*k
                            res += [QA([0, a11*v1, a22*w, a33*t + a13*v1 + a23*w])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = a11*v*i + a22*w*j + (a33*t + a13*v + a23*w)*k
                            res += [QA([0, a11*v2, a22*w, a33*t + a13*v2 + a23*w])] #sols.append([u, v2, w, t]) #u=0
        elif a00 != 0 and a10 == 0 and a11 != 0 and a12 != 0 and a13 != 0 and a20 == 0 and a21 == 0 and a22 != 0 and a23 == 0 and a30 == 0 and a31 == 0 and a32 == 0 and a33 != 0 : #O = [ QA(a00 + a01*i + a02*j + a03*k), QA(a11*i + a12*j + a13*k), QA(a22*j), QA(a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a13*a33*p*q*t*v + a13^2*p*q*v^2 + a12^2*p*v^2 + a11^2*q*v^2 + 2*a12*a22*p*v*w + a22^2*p*w^2. a01, a02, a03 may be 0.
            print ('Case 8')
            pp0 = a12^2*p + a11^2*q
            pp1 = a13^2*p + a11^2
            pp2 = a11^2*a22^2*p*q
            pp6 = a11^2*a33^2*p*q
            newbdw = sqrt(pp0/pp2) # < bdw = sqrt(abs((a00^2*a11^2*q + (a00^2*a12^2 + (a02*a11 - a01*a12)^2*q)*p)/(p*q*a00^2*a11^2*a22^2)))
            newbdt = sqrt(pp1/pp6) # < bdt = sqrt((a00^2*a11^2 + (a00^2*a13^2 + (a03*a11 - a01*a13)^2*q)*p)/(p*q*a00^2*a11^2*a33^2))
            pp3 = a13^2*p*q + pp0
            D0 = pp3*NN
            pp4 = a13*a33*p*q
            pp5 = a12*a22*p
            pp7 = pp1*a22^2*p*q
            pp8 = 2*pp5*pp4
            pp9 = pp0*a33^2*p*q
            for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-64)*w^2*a11^2*a22^2*a33^4*p^3*q^3 + 64*(a12^2*p + a11^2*q)*NN*a33^4*p^2*q^2 >= 0
                D1 = - pp7*w^2 + D0
                D2 = pp8*w
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = -64*(a13^2*p*q + a12^2*p + a11^2*q)*a11^2*a22^2*a33^2*p^2*q^2*t^2 + 64*(a13^2*p*q + a12^2*p + a11^2*q)*(a13^2*p + a11^2)*NN*a22^2*p*q >= 0
                    Delta_v_over_4 = - pp9*t^2 + D2*t + D1 #norm_eq = (a13^2*p*q + a12^2*p + a11^2*q)*v^2 + 2*(a13*a33*q*t + a12*a22*w)*p*v + t^2*a33^2*p*q + w^2*a22^2*p - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- (pp4*t + pp5*w) + sqrt(Delta_v_over_4))/pp3
                        v2 = (- (pp4*t + pp5*w) - sqrt(Delta_v_over_4))/pp3
                        if v1 in ZZ: #x = a11*v*i + (a12*v + a22*w)*j + (a33*t + a13*v)*k
                            res += [QA([0, a11*v1, a12*v1 + a22*w, a33*t + a13*v1])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = a11*v*i + (a12*v + a22*w)*j + (a33*t + a13*v)*k
                            res += [QA([0, a11*v2, a12*v2 + a22*w, a33*t + a13*v2])] #sols.append([u, v2, w, t]) #u=0
        elif a00 != 0 and a10 == 0 and a11 != 0 and a12 != 0 and a13 != 0 and a20 == 0 and a21 == 0 and a22 != 0 and a23 != 0 and a30 == 0 and a31 == 0 and a32 == 0 and a33 != 0 : #O = [QA(a00 + a01*i + a02*j + a03*k), QA(a11*i + a12*j + a13*k), QA(a22*j + a23*k), QA(a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a13*a33*p*q*t*v + a13^2*p*q*v^2 + 2*a23*a33*p*q*t*w + 2*a13*a23*p*q*v*w + a23^2*p*q*w^2 + a12^2*p*v^2 + a11^2*q*v^2 + 2*a12*a22*p*v*w + a22^2*p*w^2. a01, a02, a03 may be 0
            print ('Case 9')
            pp0 = a11^2*a22^2*p*q
            pp1 = a12^2*p + a11^2*q
            pp2 = (a13*a22 - a12*a23)^2*p + a11^2*a23^2*q + a11^2*a22^2
            newbdw = sqrt(pp1/pp0) # < bdw = sqrt((a00^2*a11^2*q + (a00^2*a12^2 + (a02*a11 - a01*a12)^2*q)*p)/(p*q*a00^2*a11^2*a22^2))
            newbdt = sqrt(pp2/(pp0*a33^2)) # < bdt = sqrt(abs((a00^2*a11^2*a23^2*q + a00^2*a11^2*a22^2 + ((a13*a22 - a12*a23)^2*a00^2 + ((a13*a22 - a12*a23)*a01 - (a03*a22 - a02*a23)*a11)^2*q)*p)/(p*q)))/abs(a00*a11*a22*a33)
            pp3 = a13^2*p*q + pp1
            D0 = pp3*NN
            pp4 = pp2*p*q
            pp5 = 2*(a12*a13*a22*p - a12^2*a23*p - a11^2*a23*q)*a33*p*q
            pp6 = (a12^2*p + a11^2*q)*a33^2*p*q
            pp7 = a13*a33*p*q
            pp8 = a13*a23*q + a12*a22
            for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = (-64)*w^2*a11^2*a22^2*a33^4*p^3*q^3 + 64*(a12^2*p + a11^2*q)*NN*a33^4*p^2*q^2 >= 0
                D1 = - pp4*w^2 + D0
                D2 = pp5*w
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = -64*(a13^2*p*q + a12^2*p + a11^2*q)*a11^2*a22^2*a33^2*p^2*q^2*t^2 + 64*(a13^2*a22^2*p - 2*a12*a13*a22*a23*p + a12^2*a23^2*p + a11^2*a23^2*q + a11^2*a22^2)*(a13^2*p*q + a12^2*p + a11^2*q)*NN*p*q >= 0
                    Delta_v_over_4 = - pp6*t^2 + D2*t + D1  #norm_eq = (a13^2*p*q + a12^2*p + a11^2*q)*v^2 + 2*(a13*a33*p*q*t + a13*a23*p*q*w + a12*a22*p*w)*v + a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a22^2*p*w^2 - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- (pp7*t + pp8*p*w) + sqrt(Delta_v_over_4))/pp3
                        v2 = (- (pp7*t + pp8*p*w) - sqrt(Delta_v_over_4))/pp3
                        if v1 in ZZ: #x = a11*v*i + (a12*v + a22*w)*j + (a33*t + a13*v + a23*w)*k
                            res += [QA([0, a11*v1, a12*v1 + a22*w, a33*t + a13*v1 + a23*w])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = a11*v*i + (a12*v + a22*w)*j + (a33*t + a13*v + a23*w)*k
                            res += [QA([0, a11*v2, a12*v2 + a22*w, a33*t + a13*v2 + a23*w])] #sols.append([u, v2, w, t]) #u=0
        elif a00 != 0 and a10 == 0 and a11 != 0 and a12 != 0 and a13 != 0 and a20 == 0 and a21 == 0 and a22 != 0 and a23 != 0 and a30 == 0 and a31 == 0 and a32 != 0 and a33 != 0 and (a23*a32 - a22*a33) != 0 : #O = [ QA(a00 + a01*i + a02*j + a03*k), QA(a11*i + a12*j + a13*k), QA(a22*j + a23*k), QA(a32*j + a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a13*a33*p*q*t*v + a13^2*p*q*v^2 + 2*a23*a33*p*q*t*w + 2*a13*a23*p*q*v*w + a23^2*p*q*w^2 + a32^2*p*t^2 + 2*a12*a32*p*t*v + a12^2*p*v^2 + a11^2*q*v^2 + 2*a22*a32*p*t*w + 2*a12*a22*p*v*w + a22^2*p*w^2. a01, a02, a03 may be 0.
            print ('Case 10')
            pp0 = (a23*a32 - a22*a33)^2*a11^2*p*q
            pp1 = (a13*a32 - a12*a33)^2*p + a11^2*(a33^2*q + a32^2)
            pp2 = (a13*a22 - a12*a23)^2*p + a11^2*(a23^2*q + a22^2)
            newbdw = sqrt(pp1/pp0) #< bdw = sqrt(abs((a00^2*a11^2*a33^2*q + a00^2*a11^2*a32^2 + ((a13*a32 - a12*a33)^2*a00^2 + ((a13*a32 - a12*a33)*a01 - (a03*a32 - a02*a33)*a11)^2*q)*p)/(p*q)))/abs(-a00*a11*a23*a32 + a00*a11*a22*a33)
            newbdt = sqrt(pp2/pp0) #< bdt = sqrt(abs((a00^2*a11^2*a23^2*q + a00^2*a11^2*a22^2 + ((a13*a22 - a12*a23)^2*a00^2 + ((a13*a22 - a12*a23)*a01 - (a03*a22 - a02*a23)*a11)^2*q)*p)/(p*q)))/abs(-a00*a11*a23*a32 + a00*a11*a22*a33)
            pp3 = a13^2*p*q + a12^2*p + a11^2*q
            D0 = pp3*NN
            pp4 = pp2*p*q
            pp5 = 2*(a13^2*a22*a32*p - a12*a13*a23*a32*p - a12*a13*a22*a33*p + a12^2*a23*a33*p + a11^2*a23*a33*q + a11^2*a22*a32)*p*q
            pp6 = pp1*p*q
            pp7 = a13*a33*q + a12*a32
            pp8 = a13*a23*q + a12*a22
            for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_tv = -64*(a33^2*q + a32^2)*(a23*a32 - a22*a33)^2*a11^2*p^3*q^2*w^2 + 64*(a13^2*a32^2*p - 2*a12*a13*a32*a33*p + a12^2*a33^2*p + a11^2*a33^2*q + a11^2*a32^2)*(a33^2*q + a32^2)*NN*p^2*q >= 0
                D1 = - pp4*w^2 + D0
                D2 = pp5*w
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = -64*(a13^2*p*q + a12^2*p + a11^2*q)*(a23*a32 - a22*a33)^2*a11^2*p^2*q^2*t^2 + 64*(a13^2*a22^2*p - 2*a12*a13*a22*a23*p + a12^2*a23^2*p + a11^2*a23^2*q + a11^2*a22^2)*(a13^2*p*q + a12^2*p + a11^2*q)*NN*p*q >= 0
                    Delta_v_over_4 = - pp6*t^2 - D2*t + D1 #norm_eq = (a13^2*p*q + a12^2*p + a11^2*q)*v^2 + 2*((a13*a33*q + a12*a32)*t + (a13*a23*q + a12*a22)*w)*p*v + (a33^2*q + a32^2)*p*t^2 + (a23^2*q + a22^2)*p*w^2 + 2*(a23*a33*q + a22*a32)*p*t*w - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- (pp7*t + pp8*w)*p + sqrt(Delta_v_over_4))/pp3
                        v2 = (- (pp7*t + pp8*w)*p - sqrt(Delta_v_over_4))/pp3
                        if v1 in ZZ: #x = a11*v*i + (a32*t + a12*v + a22*w)*j + (a33*t + a13*v + a23*w)*k
                            res += [QA([0, a11*v1, a12*v1 + a22*w + a32*t, a33*t + a13*v1 + a23*w])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = a11*v*i + (a32*t + a12*v + a22*w)*j + (a33*t + a13*v + a23*w)*k
                            res += [QA([0, a11*v2, a12*v2 + a22*w + a32*t, a33*t + a13*v2 + a23*w])] #sols.append([u, v2, w, t]) #u=0                          
        elif a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a10 != 0 and a11 != 0 and a12 == 0 and a13 == 0 and a20 == 0 and a21 == 0 and a22 != 0 and a23 != 0 and a30 == 0 and a31 == 0 and a32 != 0 and a33 != 0 and (a23*a32 - a22*a33) != 0 : #O = [ QA(1), QA(a10 + a11*i), QA(a22*j + a23*k), QA(a32*j + a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a10*v + 2*u. Put v = - u/a10. Nx = a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a32^2*p*t^2 + 2*a22*a32*p*t*w + a22^2*p*w^2 + a11^2*q*u^2/a10^2
            print ('Case 11')
            pp0 = a11^2*q
            pp3 = a23^2*q + a22^2
            pp2 = pp0*pp3*p/a10^2
            pp4 = - (a23*a32 - a22*a33)^2*p^2*q
            pp1 = pp3*p
            newbdu = sqrt(a10^2/pp0)
            D0 = NN*a10^2*pp3*p/a10^2
            pp5 = (a23*a33*q + a22*a32)*p
            for u in [ceil(-newbdu*sqN)..floor(newbdu*sqN)]: #Delta_wt = 64*(a23^2*q + a22^2)*(a23*a32 - a22*a33)^2*NN*p^3*q - 64*(a23^2*q + a22^2)*(a23*a32 - a22*a33)^2*a11^2*p^3*q^2/a10^2*u^2 >= 0 
                D1 = - pp2*u^2 + D0
                for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: #Delta_uw = 64*(a23^2*q + a22^2)*NN*a11^4*p*q^2/a10^4 - 64*(a23*a32 - a22*a33)^2*a11^4*p^2*q^3/a10^4*t^2 => 0
                    Delta_w_over_4 = pp4*t^2 + D1 
                    #norm_eq = (a23^2*q + a22^2)*p*w^2 + 2*(a23*a33*p*q*t + a22*a32*p*t)*w + a33^2*p*q*t^2 + a32^2*p*t^2 + a11^2*q*u^2/a10^2 - NN
                    if Delta_w_over_4 >= 0: 
                        w1 = ( pp5*t + sqrt(Delta_w_over_4))/pp1
                        w2 = ( pp5*t - sqrt(Delta_w_over_4))/pp1
                        if w1 in ZZ: #x = (-a11*u/a10)*i + (a32*t + a22*w)*j + (a33*t + a23*w)*k
                            res += [QA([0, -a11*u/a10,  a32*t + a22*w1, a33*t + a23*w1])] #sols.append([u, v, w1, t]) #v = - u/a10 - a20/a10*w - a30/a10*t
                        if w2 in ZZ: #(-a11*u/a10)*i + (a32*t + a22*w)*j + (a33*t + a23*w)*k
                            res += [QA([0, -a11*u/a10,  a32*t + a22*w2, a33*t + a23*w2])] #sols.append([u, v, w2, t]) #v = - u/a10 - a20/a10*w - a30/a10*t
        elif a00 != 0 and a10 == 0 and a11 != 0 and a12 == 0 and a13 == 0 and a20 == 0 and a21 == 0 and a22 != 0 and a23 != 0 and a30 == 0 and a31 == 0 and a32 == 0 and a33 != 0 : #O = [ QA(a00 + a01*i + a02*j + a03*k), QA(a11*i), QA(a22*j + a23*k), QA(a33*k)  ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a11^2*q*v^2 + a22^2*p*w^2. a01, a02, a03 may be 0
            print ('Case 12')
            pp4 = a23^2*q + a22^2
            pp0 = pp4*p
            pp1 = a11^2*q
            pp5 = a22^2*a33^2*p*q
            newbdv = sqrt(1/pp1) #< bdv = sqrt((a01^2*q + a00^2)/(a00^2*a11^2q))
            newbdt = sqrt(pp4/pp5) #= bdt = sqrt((a23^2*q + a22^2)/(p*q*a22^2*a33^2))
            D0 = pp0*NN
            pp2 = a23*a33*p*q
            pp3 = pp5*p
            for v in [ceil(-newbdv*sqN)..floor(newbdv*sqN)]: #Delta_tw = (-64)*v^2*a11^2*a22^2*a33^4*p^3*q^3 + 64*a22^2*a33^4*p^3*q^2*NN >= 0 => 1/(a11^2*q)*NN >= v^2
                D1 = - pp0*pp1*v^2 + D0
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = (-64)*t^2*a11^4*a22^2*a33^2*p^2*q^3 + 64*(a23^2*q + a22^2)*NN*a11^4*p*q^2 >= 0 => (a23^2*q + a22^2)/(a22^2*a33^2*p*q)*NN >= t^2
                    Delta_w_over_4 =  - pp3*t^2 + D1 #norm_eq = (a23^2*q + a22^2)*p*w^2 + 2*a23*a33*p*q*t*w + t^2*a33^2*p*q + v^2*a11^2*q - NN
                    if Delta_w_over_4 >= 0: 
                        w1 = (- pp2*t + sqrt(Delta_w_over_4))/pp0
                        w2 = (- pp2*t - sqrt(Delta_w_over_4))/pp0
                        if w1 in ZZ: #x = a11*v*i + a22*w*j + (a33*t + a23*w)*k
                            res += [QA([0, a11*v, a22*w1, a33*t + a23*w1])] #sols.append([u, v1, w, t]) #u=0
                        if w2 in ZZ: #x = a11*v*i + a22*w*j + (a33*t + a23*w)*k
                            res += [QA([0, a11*v, a22*w2, a33*t + a23*w2])] #sols.append([u, v2, w, t]) #u=0
        elif a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a11 != 0 and a12 == 0 and a13 == 0 and a21 != 0 and a22 != 0 and a23 != 0 and a31 != 0 and a32 == 0 and a33 != 0 : #O = [QA(1), QA(a10 + a11*i), QA(a20 + a21*i + a22*j + a23*k), QA(a30 + a31*i + a33*k) ] x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a31^2*q*t^2 + 2*a11*a31*q*t*v + a11^2*q*v^2 + 2*a21*a31*q*t*w + 2*a11*a21*q*v*w + a22^2*p*w^2 + a21^2*q*w^2. a10, a20, a30 may be 0.
            print ('Case 13')
            pp0 = a11^2*q
            pp1 = a11*a31*q
            pp2 = a11*a21*q
            D0 = pp0*NN
            pp3 = (a23^2*q + a22^2)*p*pp0
            pp4 = 2*a23*a33*p*pp0*q
            for w in [ceil(-bdw*sqN)..floor(bdw*sqN)]: #Delta_vt = (-64)*w^2*a11^4*a22^2*a33^2*p^2*q^3 + 64*a11^4*a33^2*p*q^3*NN >= 0 => 1/(a22^2*p)*NN >= w^2 
                D1 = - pp3*w^2 + D0
                D2 = pp4*w
                for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]: #Delta_vw = (-64)*t^2*a11^4*a22^2*a33^2*p^2*q^3 + 64*(a23^2*q + a22^2)*NN*a11^4*p*q^2 >= 0 => (a23^2*q + a22^2)/(a22^2*a33^2*p*q)*NN >= t^2
                    Delta_v_over_4 = - a11^2*a33^2*p*q^2*t^2 - D2*t + D1 #norm_eq = a11^2*q*v^2 +  2*(a11*a31*q*t + a11*a21*q*w)*v + a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a31^2*q*t^2 + 2*a21*a31*q*t*w + a22^2*p*w^2 + a21^2*q*w^2 - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- (pp1*t + pp2*w) + sqrt(Delta_v_over_4))/pp0
                        v2 = (- (pp1*t + pp2*w) - sqrt(Delta_v_over_4))/pp0
                        if v1 in ZZ: #x = (a11*v + a21*w + a31*t)*i + a22*w*j + (a33*t + a23*w)*k
                            res += [QA([0, a11*v1 + a21*w + a31*t, a22*w, a33*t + a23*w])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = (a11*v + a21*w + a31*t)*i + a22*w*j + (a33*t + a23*w)*k
                            res += [QA([0, a11*v2 + a21*w + a31*t, a22*w, a33*t + a23*w])] #sols.append([u, v2, w, t]) #u=0
        elif a00 == 1 and a01 == 0 and a02 == 0 and a03 == 0 and a11 != 0 and a12 == 0 and a13 == 0 and a21 != 0 and a22 == 0 and a23 != 0 and a31 != 0 and a32 != 0 and a33 != 0 : #O = [QA(1), QA(a10 + a11*i), QA(a20 + a21*i + a23*k), QA(a30 + a31*i + a32*j + a33*k) ]. x = QA(u*O[0]+ v*O[1]+w*O[2]+ t*O[3]). Trx = 2*a00*u. Put u = 0. Nx = a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a32^2*p*t^2 + a31^2*q*t^2 + 2*a11*a31*q*t*v + a11^2*q*v^2 + 2*a21*a31*q*t*w + 2*a11*a21*q*v*w + a21^2*q*w^2. a10, a20, a30 may be 0.
            print ('Case 14')
            pp0 = a11^2*q
            pp1 = a33^2*q + a32^2
            pp2 = a32^2*p
            pp3 = a23^2*pp2*q
            newbdw = sqrt(pp1/pp3) #= bdw = sqrt((a33^2*q + a32^2)/(p*q*a23^2*a32^2))
            newbdt = sqrt(1/pp2) # = bdt =  sqrt(1/(p*a32^2))
            D0 = pp0*NN
            pp4 = pp0*a23^2*p*q
            pp5 = 2*a23*a33*p*q*pp0
            pp6 = pp1*p*pp0
            pp7 = a11*a31*q
            pp8 = a11*a21*q
            for w in [ceil(-newbdw*sqN)..floor(newbdw*sqN)]: #Delta_vt = (-64)*w^2*a11^4*a23^2*a32^2*p^2*q^3 + 64*(a33^2*q + a32^2)*NN*a11^4*p*q^2 >= 0 => (a33^2*q + a32^2)/(a23^2*a32^2*p*q)*NN >= w^2 
                D1 = - pp4*w^2 + D0
                D2 = pp5*w
                for t in [ceil(-newbdt*sqN)..floor(newbdt*sqN)]: #Delta_vw = (-64)*t^2*a11^4*a23^2*a32^2*p^2*q^3 + 64*a11^4*a23^2*p*q^3*NN >= 0 => 1/(a32^2*p)*NN >= t^2
                    Delta_v_over_4 =  - pp6*t^2 - D2*t + D1 #norm_eq = a11^2*q*v^2 + 2*(a11*a31*q*t + a11*a21*q*w)*v + a33^2*p*q*t^2 + 2*a23*a33*p*q*t*w + a23^2*p*q*w^2 + a32^2*p*t^2 + a31^2*q*t^2 + 2*a21*a31*q*t*w + a21^2*q*w^2 - NN
                    if Delta_v_over_4 >= 0: 
                        v1 = (- (pp7*t + pp8*w) + sqrt(Delta_v_over_4))/pp0
                        v2 = (- (pp7*t + pp8*w) - sqrt(Delta_v_over_4))/pp0
                        if v1 in ZZ: #x = (a11*v + a21*w + a31*t )*i + a32*t*j + (a33*t + a23*w)*k
                            res += [QA([0, a11*v1 + a21*w + a31*t, a32*t, a33*t + a23*w])] #sols.append([u, v1, w, t]) #u=0
                        if v2 in ZZ: #x = (a11*v + a21*w + a31*t )*i + a32*t*j + (a33*t + a23*w)*k
                            res += [QA([0, a11*v2 + a21*w + a31*t, a32*t, a33*t + a23*w])] #sols.append([u, v2, w, t]) #u=0
        elif a00 != 0:
            print ('Case 15')
            for v in [ceil(-bdv*sqN)..floor(bdv*sqN)]:
                for w in [ceil(-bdw*sqN)..floor(bdw*sqN)]:
                    for t in [ceil(-bdt*sqN)..floor(bdt*sqN)]:
                        u = (- a10*v - a20*w - a30*t)/a00
                        elem = QA(u*OO[0] + v*OO[1] + w*OO[2] + t*OO[3])
                        if elem.reduced_norm() == NN:
                            res += [QA(elem)]
        else:
            print ('Reorder basis such that first element in the basis u0 + u1*i + u2*j + u3*k is such that u0 is not equal to 0.')
    return list(set(res))


