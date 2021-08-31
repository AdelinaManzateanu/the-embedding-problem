# Authors: Sorina Ionica, Pınar Kılıçer, Kristin Lauter, Elisa Lorenzo Garcia, Adelina Mânzățeanu and Christelle Vincent.

# This document contains a list of CM fields corresponding to the curves studied as part of the project "DETERMINING THE PRIMES OF BAD REDUCTION OF CM CURVES OF GENUS 3".

#--------------------------------------------------------------------------------------------------------------------------------------------------

# CM fields 1, 2 (= 22), 3 (= 33), 4, 5 (= 55, 555, 5555), 6 (= 66, 666, 6666, 66666), 7 (= 77), 8 (= 88, 888, 8888, 88888) correspond to the hyperelliptic genus 3 curves H1 - H8 in Section 7.1 of [IKLLMV].

# CM fields 9/29 (= 99/2999, 999/29999, 29), 10 (= 1010, 1011), 23 (= 233, 2333), 25 (= 255, 2555, 25555, 255555), 26 (= 266, 2666, 26666), 27 (= 277, 2777), 28 (= 288, 2888), 210 (=2102, 2103, 2104, 2105, 2106, 2107, 2108, 2109, 21010, 21011), 211 (= 2112, 2113, 2114), 212 (= 2122), 217 (= 2177, 21777, 217777, 2177777) correspond to the CM non-hyperelliptic curves of genus 3 defined over the rationals X9, X1, X3, X5, X6, X8, X7, X10, X11, X12, X17, respectively in Section 7.2 of [IKLLMV]. See also Table 4 in https://arxiv.org/pdf/1803.05816.pdf and page 21 of https://arxiv.org/pdf/1701.06489.pdf

# CM fields 31 (= 311), 32, 34, 35, 36, 37, 38 (= 388), 39, 40, 41, 42 (= 422), 43, 45, 46, 47, 48 (= 488), 49, 50, 51, 52, 53, 54, 56, 57 are examples of sextic CM fields with no imaginary quadraric subfield. Some of these are presented in Section 7.3 of [IKLLMV]: G1 (32), G2 (41), G3 (43), G4 (50), G5 (54).

#--------------------------------------------------------------------------------------------------------------------------------------------------


if CMfieldnumber == 1:          # curve H1. indices: 1 in K+, 7 in K (best A, best indices)
    A = 13; B = 50; C = 49; D = 1
if CMfieldnumber == 2:          # curve H2. indices: 1 in K+, 1 in K (best A, best indices)
    A = 6; B = 9; C = 1; D = 1
if CMfieldnumber == 22:         # curve H2. indices: 1 in K+, 1 in K (best indices)
    A = 9; B = 6; C = 1; D = 1
if CMfieldnumber == 3:          # curve H3. indices: 1 in K+, 1 in K (best A, best indices)
    A = 5; B = 6; C = 1; D = 1 
if CMfieldnumber == 33:         # curve H3. indices: 1 in K+, 1 in K (best A, best indices)
    A = 6; B = 5; C = 1; D = 1 
if CMfieldnumber == 4:          # curve H4. indices: 1 in K+, 8 in K (best A, best indices)
    A = 7; B = 14; C = 7; D = 7
if CMfieldnumber == 5:          # curve H5. indices: 55 in K+, 266200 in K (best A)
    A = 42; B = 441; C = 847; D = 7
if CMfieldnumber == 55:         # curve H5. indices: 15 in K+, 5400 in K (best indices ?)
    A = 63; B = 126; C = 63; D = 7
if CMfieldnumber == 555:        # curve H5. indices: 4983 in K+, 4966057800 in K (index in K+ not divisible by 5)
    A = 210; B = 2541; C = 4375; D = 7
if CMfieldnumber == 5555:       # curve H5. cubic is isom, but sextic not. index: 1 in K+
    A = 12; B = 27; C = 15; D = 7
if CMfieldnumber == 6:          # curve H6. indices: 44 in K+, 15488 in K (best A)
    A = 29; B = 180; C = 64; D = 1
if CMfieldnumber == 66:         # curve H6. indices: 8 in K+, 1408 in K
    A = 30; B = 257; C = 484; D = 1
if CMfieldnumber == 666:        # curve H6. indices: 4 in K+, 512 in K (best indices ?)
    A = 37; B = 356; C = 1024; D = 1
if CMfieldnumber == 6666:       # curve H6. cubic is isom, but sextic not. index: 4 in K+
    A = 19; B = 20; C = 4; D = 1
if CMfieldnumber == 66666:      # curve H6. cubic is isom, but sextic not. index: 2 in K+
    A = 10; B = 19; C = 2; D = 1
if CMfieldnumber == 7:          # curve H7. indices: 4 in K+, 128 in K (best indices ?)
    A = 21; B = 116; C = 64; D = 1
if CMfieldnumber == 77:         # curve H7. indices: 4 in K+, 128 in K (smallest A, best indices ?)
    A = 11; B = 30; C = 16; D = 1
if CMfieldnumber == 8:          # curve H8. indices: 56 in K+, 87808 in K (best A)
    A = 42; B = 441; C = 784; D = 1
if CMfieldnumber == 88:         # curve H8. indices: 12 in K+, 6912 in K (best indices ?)
    A = 45; B = 612; C = 2304; D = 1
if CMfieldnumber == 888:        # curve H8. indices: 124 in K+, 1230008 in K
    A = 45; B = 276; C = 64; D = 1
if CMfieldnumber == 8888:       # curve H8. cubic is isom, but sextic not. index: 2 in K+
    A = 33; B = 90; C = 64; D = 1
if CMfieldnumber == 88888:      # curve H8. cubic is isom, but sextic not. index: 2 in K+
    A = 12; B = 27; C = 8; D = 1    
if CMfieldnumber == 9:          # curve X9. indices: 40 in K+, 1600 in K
    A = 18; B = 56; C = 8; D = 2 
if CMfieldnumber == 99:         # curve X9. indices: 8 in K+, 320 in K
    A = 20; B = 116; C = 200; D = 2 
if CMfieldnumber == 999:        # curve X9. cubic is isom, but sextic not. index: 1 in K+
    A = 7; B = 12; C = 5; D = 2 
if CMfieldnumber == 10:         # curve X1. indices: 1715 in K+, 23529800 in K
    A = 63; B = 686; C = 343; D = 7 
if CMfieldnumber == 1010:       # curve X1. indices: 343 in K+, 4705960 in K
    A = 70; B = 1421; C = 8575; D = 7
if CMfieldnumber == 1011:       # curve X1. cubic is isom, but sextic not. index: 1 in K+
    A = 7; B = 12; C = 5; D = 7 
if CMfieldnumber == 221:         # curve X2. indices: 343 in K+, 941192 in K
    A = 42; B = 441; C = 343; D = 7
if CMfieldnumber == 2211:         # curve X2. icubic is isom, but sextic not. index: 1 in K+
    A = 6; B = 9; C = 1; D = 7
if CMfieldnumber == 23:         # curve X3. indices: 447615 in K+, 302943092596200 in K
    A = 1162; B = 106281; C = 250047; D = 7
if CMfieldnumber == 233:        # curve X3. indices: 675 in K+, 91125000 in K
    A = 147; B = 5250; C = 4375; D = 7
if CMfieldnumber == 2333:       # curve X3. cubic is isom, but sextic not. index: 3 in K+
    A = 20; B = 61; C = 15; D = 7
if CMfieldnumber == 25:         # curve X5. indices: 58653 in K+, 5201543706408 in K
    A = 427; B = 41454; C = 250047; D = 7
if CMfieldnumber == 255:        # curve X5. indices: 81 in K+, 997272 in K
    A = 63; B = 1050; C = 2527; D = 7
if CMfieldnumber == 2555:       # curve X5. indices: 297 in K+, 6351048 in K
    A = 70; B = 693; C = 567; D = 7
if CMfieldnumber == 25555:      # curve X5. indices: 1341 in K+, 129476232 in K
    A = 119; B = 1414; C = 567; D = 7
if CMfieldnumber == 255555:     # curve X5. cubic is isom, but sextic not. index: 3 in K+
    A = 16; B = 55; C = 27; D = 7
if CMfieldnumber == 26:         # curve X6. indices: 188307 in K+, 53614803688488 in K
    A = 658; B = 63945; C = 250047; D = 7
if CMfieldnumber == 266:        # curve X6. indices: 207 in K+, 9255384 in K
    A = 91; B = 1386; C = 5103; D = 7
if CMfieldnumber == 2666:       # curve X6. indices: 3 in K+, 4968 in K
    A = 98; B = 3157; C = 33327; D = 7
if CMfieldnumber == 26666:      # curve X6. cubic is isom, but sextic not. index: 3 in K+
    A = 19; B = 76; C = 57; D = 7
if CMfieldnumber == 27:         # curve X8. indices: 3773 in K+, 113884232 in K
    A = 98; B = 1029; C = 343; D = 7
if CMfieldnumber == 277:        # curve X8. indices: 343 in K+, 6588344 in K
    A = 91; B = 2450; C = 16807; D = 7
if CMfieldnumber == 2777:       # curve X8. indices: cubic is isom, but sextic not. index: 1 in K+
    A = 8; B = 15; C = 7; D = 7
if CMfieldnumber == 28:         # curve X7. indices: 3087 in K+, 2058386904 in K
    A = 343; B = 30870; C = 250047; D = 7
if CMfieldnumber == 288:        # curve X7. cubic is isom, but sextic not. index: 3 in K+
    A = 14; B = 41; C = 7; D = 7
if CMfieldnumber == 2888:       # curve X7. cubic is isom, but sextic not. index: 3 in K+
    A = 17; B = 72; C = 63; D = 7
if CMfieldnumber == 28888:       # curve X7. cubic is isom, but sextic not. index: 3 in K+
    A = 16; B = 61; C = 3; D = 7
if CMfieldnumber == 29:         # curve X9. indices: 2560 in K, 52428800 in K
    A = 72; B = 896; C = 512; D = 2
if CMfieldnumber == 299:        # curve X9. indices: 40 in K+, 1600 in K (Same as CM 9)
    A = 18; B = 56; C = 8; D = 2
if CMfieldnumber == 2999:       # curve X9. indices: 8 in K+, 320 in K (Same as CM 99)
    A = 20; B = 116; C = 200; D = 2
if CMfieldnumber == 29999:      # curve X9. cubic is isom, but sextic not. index: 1 in K+ (Same as CM 999)
    A = 7; B = 12; C = 5; D = 2
if CMfieldnumber == 210:        # curve X10. indices: 512 in K+, 2097152 in K
    A = 40; B = 384; C = 512; D = 2
if CMfieldnumber == 2102:       # curve X10. indices: 8 in K+, 64 in K
    A = 10; B = 24; C = 8; D = 2
if CMfieldnumber == 2103:       # curve X10. indices: 8 in K+, 64 in K
    A = 12; B = 20; C = 8; D = 2
if CMfieldnumber == 2104:       # curve X10. indices: 104 in K+, 10816 in K
    A = 20; B = 68; C = 8; D = 2
if CMfieldnumber == 2105:       # curve X10. indices: 232 in K+, 53824 in K
    A = 26; B = 104; C = 8; D = 2
if CMfieldnumber == 2106:       # curve X10. indices: 104 in K+, 10816 in K
    A = 34; B = 40; C = 8; D = 2
if CMfieldnumber == 2107:       # curve X10. indices: 8 in K+, 832 in K
    A = 34; B = 376; C = 1352; D = 2
if CMfieldnumber == 2108:       # curve X10. indices: 512 in K+, 262144 in K
    A = 38; B = 332; C = 8; D = 2
if CMfieldnumber == 2109:       # curve X10. indices: 728 in K+, 3709888 in K
    A = 42; B = 392; C = 392; D = 2
if CMfieldnumber == 21010:      # curve X10. indices: 512 in K+, 2097152 in K
    A = 48; B = 320; C = 512; D = 2
if CMfieldnumber == 21011:      # curve X10. cubic is isom, but sextic not. index: 1 in K+
    A = 5; B = 6; C = 1; D = 2
if CMfieldnumber == 211:        # curve X11. indices: 2048 in K+, 268435456 in K
    A = 168; B = 7424; C = 32768; D = 2
if CMfieldnumber == 2112:       # curve X11. indices: 32 in K+, 8192 in K (smallest A with sextic isom)
    A = 42; B = 464; C = 512; D = 2
if CMfieldnumber == 2113:       # curve X11. cubic is isom, but sextic not. index: 4 in K+
    A = 12; B = 17; C = 2; D = 2
if CMfieldnumber == 2114:       # curve X11. cubic is isom, but sextic not. index: 2 in K+
    A = 13; B = 46; C = 32; D = 2
if CMfieldnumber == 212:        # curve X12. indices: 1331 in K+, 14172488 in K
    A = 55; B = 726; C = 1331; D = 11
if CMfieldnumber == 2122:       # curve X12. cubic is isom, but sextic not. index: 1 in K+
    A = 5; B = 6; C = 1; D = 11
if CMfieldnumber == 215:        # curve X15. indices: 7 in K+, 392 in K
    A = 19; B = 38; C = 19; D = 19
if CMfieldnumber == 2155:       # curve X15. cubic is isom, but sextic not. index: 1 in K+
    A = 8; B = 15; C = 7; D = 19
if CMfieldnumber == 216:        # curve X16. indices: 6859 in K+, 376367048 in K
    A = 114; B = 3249; C = 6859; D = 19
if CMfieldnumber == 2166:       # curve X16. cubic is isom, but sextic not. index: 1 in K+
    A = 6; B = 9; C = 1; D = 19
if CMfieldnumber == 217:        # curve X17. indices: 3813604 in K+, 17684987770080256 in K
    A = 2679; B = 272916; C = 438976; D = 19
if CMfieldnumber == 2177:       # curve X17. indices: 416 in K+, 35995648 in K
    A = 114; B = 3249; C = 12844; D = 19
if CMfieldnumber == 21777:      # curve X17. indices: 192 in K+, 3538944 in K
    A = 171; B = 1368; C = 2736; D = 19
if CMfieldnumber == 217777:     # curve X17. indices: 1024 in K+, 536870912 in K
    A = 171; B = 7296; C = 77824; D = 19
if CMfieldnumber == 2177777:    # curve X17. cubic is isom, but sextic not. index: 2 in K+
    A = 18; B = 51; C = 26; D = 19
if CMfieldnumber == 31:         # no imaginary quadratic subfield. indices: 4, 64
    A = 12; B = 17; C = 2
if CMfieldnumber == 311:        # no imaginary quadratic subfield. indices: 2, 64
    A = 13; B = 46; C = 32
if CMfieldnumber == 32:         # curve G1. no imaginary quadratic subfield. indices: 1, 8
    A = 15; B = 14; C = 3
if CMfieldnumber == 34:         # no imaginary quadratic subfield. indices: 3, 81
    A = 22; B = 139; C = 243
if CMfieldnumber == 35:         # no imaginary quadratic subfield. indices: 1, 1
    A = 8; B = 15; C = 7
if CMfieldnumber == 36:         # no imaginary quadratic subfield. indices: 3, 72
    A = 14; B = 41; C = 7
if CMfieldnumber == 37:         # no imaginary quadratic subfield. indices: 1, 1
    A = 10; B = 27; C = 11
if CMfieldnumber == 38:         # no imaginary quadratic subfield. indices: 8, 1024
    A = 18; B = 65; C = 44
if CMfieldnumber == 388:        # no imaginary quadratic subfield. indices: 4, 256
    A = 22; B = 61; C = 44
if CMfieldnumber == 39:         # no imaginary quadratic subfield. indices: 1, 1
    A = 9; B = 24; C = 19
if CMfieldnumber == 40:         # no imaginary quadratic subfield. indices: 1, 8
    A = 10; B = 21; C = 11
if CMfieldnumber == 41:         # curve G2. no imaginary quadratic subfield. indices: 1, 8
    A = 11; B = 34; C = 31
if CMfieldnumber == 42:         # no imaginary quadratic subfield. indices: 4, 64
    A = 12; B = 17; C = 2
if CMfieldnumber == 422:        # no imaginary quadratic subfield. indices: 2, 64
    A = 13; B = 46; C = 32
if CMfieldnumber == 43:         # curve G3. no imaginary quadratic subfield. indices: 1, 2
    A = 7; B = 10; C = 2
if CMfieldnumber == 45:         # no imaginary quadratic subfield. indices: 1, 8
    A = 8; B = 12; C = 3
if CMfieldnumber == 46:         # no imaginary quadratic subfield. indices: 1, 2
    A = 9; B = 16; C = 2
if CMfieldnumber == 47:         # no imaginary quadratic subfield. indices: 2, 8
    A = 9; B = 14; C = 4
if CMfieldnumber == 48:         # no imaginary quadratic subfield. indices: 1, 8
    A = 9; B = 17; C = 8
if CMfieldnumber == 488:        # no imaginary quadratic subfield. indices: 1, 4
    A = 15; B = 11; C = 2
if CMfieldnumber == 49:         # no imaginary quadratic subfield. indices: 2, 32
    A = 9; B = 19; C = 7
if CMfieldnumber == 50:         # curve G4. no imaginary quadratic subfield. indices: 1, 8
    A = 9; B = 21; C = 8
if CMfieldnumber == 51:         # no imaginary quadratic subfield. indices: 1, 8
    A = 10; B = 17; C = 3
if CMfieldnumber == 52:         # no imaginary quadratic subfield. indices: 1, 8
    A = 10; B = 17; C = 7
if CMfieldnumber == 53:         # no imaginary quadratic subfield. indices: 1, 8
    A = 10; B = 21; C = 3
if CMfieldnumber == 54:         # curve G5. no imaginary quadratic subfield. indices: 1, 8
    A = 10; B = 29; C = 23
if CMfieldnumber == 56:         # no imaginary quadratic subfield. indices: 1, 8
    A = 11; B = 14; C = 3
if CMfieldnumber == 57:         # no imaginary quadratic subfield. indices: 2, 16
    A = 10; B = 21; C = 4
if CMfieldnumber == 577:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+ 
    A = 9; B = 16; C = 6
if CMfieldnumber == 58:         # no imaginary quadratic subfield. indices: 1, 12
    A = 14; B = 43; C = 36
if CMfieldnumber == 588:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+ 
    A = 11; B = 18; C = 6
if CMfieldnumber == 59:         # no imaginary quadratic subfield. indices: 6, 144
    A = 21; B = 60; C = 4
if CMfieldnumber == 599:        # no imaginary quadratic subfield. indices: 8, 384
    A = 23; B = 79; C = 9
if CMfieldnumber == 5999:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 14 in K+ 
    A = 32; B = 182; C = 196
if CMfieldnumber == 59999:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+ 
    A = 13; B = 32; C = 14
if CMfieldnumber == 60:         # no imaginary quadratic subfield. indices: 6, 432
    A = 23; B = 112; C = 36
if CMfieldnumber == 600:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 6 in K+
    A = 26; B = 162; C = 36
if CMfieldnumber == 6000:       # no imaginary quadratic subfield. indices: 3, 384 
    A = 26; B = 187; C = 324
if CMfieldnumber == 60000:      # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 3 in K+
    A = 19; B = 82; C = 6   
if CMfieldnumber == 600000:     # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 3 in K+
    A = 20; B = 95; C = 46    
if CMfieldnumber == 61:         # no imaginary quadratic subfield. indices: 2, 128 
    A = 26; B = 177; C = 128
if CMfieldnumber == 611:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 19; B = 42; C = 22
if CMfieldnumber == 62:         # no imaginary quadratic subfield. indices: 2, 256 
    A = 29; B = 246; C = 512
if CMfieldnumber == 622:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 2 in K+
    A = 15; B = 34; C = 8
if CMfieldnumber == 63:         # no imaginary quadratic subfield. indices: 6, 216 
    A = 29; B = 226; C = 252
if CMfieldnumber == 633:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 3 in K+
    A = 18; B = 84; C = 63
if CMfieldnumber == 6333:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 4 in K+
    A = 23; B = 140; C = 112
if CMfieldnumber == 63333:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 11; B = 19; C = 6
if CMfieldnumber == 633333:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 14; B = 44; C = 37
if CMfieldnumber == 64:         # no imaginary quadratic subfield. indices: 8, 2560 
    A = 30; B = 169; C = 200
if CMfieldnumber == 644:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 10 in K+
    A = 44; B = 530; C = 1800
if CMfieldnumber == 6444:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 14; B = 43; C = 20
if CMfieldnumber == 64444:        # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 16; B = 63; C = 10
if CMfieldnumber == 65:       # no imaginary quadratic subfield. indices: 6, 216 
    A = 33; B = 342; C = 1116
if CMfieldnumber == 655:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 12; B = 36; C = 31
if CMfieldnumber == 6555:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 9; B = 15; C = 6
if CMfieldnumber == 67:       # no imaginary quadratic subfield. indices: 1, 8 
    A = 20; B = 16; C = 3
if CMfieldnumber == 677:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 11; B = 17; C = 4
if CMfieldnumber == 68:       # no imaginary quadratic subfield. indices: 65840, 8669811200 
    A = 1027; B = 4795; C = 25
if CMfieldnumber == 688:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 8 in K+
    A = 30; B = 208; C = 400
if CMfieldnumber == 6888:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 2 in K+
    A = 12; B = 25; C = 12
if CMfieldnumber == 69:       # no imaginary quadratic subfield. indices: 210768, 177692599296 
    A = 1768; B = 19893; C = 18
if CMfieldnumber == 699:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 59 in K+
    A = 48; B = 295; C = 18
if CMfieldnumber == 6999:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 13; B = 34; C = 10
if CMfieldnumber == 70:       # no imaginary quadratic subfield. indices: 3167731, 1846351622658424 
    A = 4332; B = 164512; C = 12167
if CMfieldnumber == 700:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 27 in K+
    A = 45; B = 213; C = 92
if CMfieldnumber == 7000:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 16; B = 34; C = 13
if CMfieldnumber == 71:       # no imaginary quadratic subfield. indices: 24, 4608 
    A = 28; B = 180; C = 144
if CMfieldnumber == 711:       # no imaginary quadratic subfield. indices: 6, 576 
    A = 34; B = 261; C = 576
if CMfieldnumber == 712:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 56 in K+
    A = 36; B = 196; C = 144
if CMfieldnumber == 713:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 16 in K+
    A = 39; B = 319; C = 729
if CMfieldnumber == 714:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 2 in K+
    A = 13; B = 41; C = 21
if CMfieldnumber == 715:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 3 in K+
    A = 14; B = 45; C = 18
if CMfieldnumber == 72:       # no imaginary quadratic subfield. indices: 29, 40368 
    A = 53; B = 761; C = 1692
if CMfieldnumber == 722:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 27 in K+
    A = 54; B = 846; C = 3807
if CMfieldnumber == 7222:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 12; B = 34; C = 17
if CMfieldnumber == 73:       # no imaginary quadratic subfield. indices: 192, 294912 
    A = 389; B = 179; C = 7
if CMfieldnumber == 733:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 12 in K+
    A = 53; B = 732; C = 1008
if CMfieldnumber == 7333:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 8 in K+
    A = 58; B = 980; C = 4032
if CMfieldnumber == 73333:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 16; B = 50; C = 21
if CMfieldnumber == 74:       # no imaginary quadratic subfield. indices: 121, 117128 
    A = 52; B = 704; C = 847
if CMfieldnumber == 744:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 214 in K+
    A = 59; B = 970; C = 3388
if CMfieldnumber == 7444:       # no imaginary quadratic subfield. cubic is isom, but sextic not. index: 1 in K+
    A = 7; B = 9; C = 2
    