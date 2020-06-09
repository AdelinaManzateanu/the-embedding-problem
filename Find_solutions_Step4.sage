# run as sage ./Find_solutions_Step4.sage CMfieldnumber p D
# e.g. sage ./Find_solutions_Step4.sage 1 7 1 runs Step 4 for CMfield 1 and prime p = 7 D = 1
#-------------------------------------------------------------------------------------------------------------------------------------------
# Step 4: Read solutions [x, dOVERn, a, cOVERn, b, gamma, n] from Step 3. For each, 
# find solutions [x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3] that satisfy the following:
# 1) the 3 equations
# 2) Nd1, Nd2, Nd3, N(d1+d2), N(d1+d3), N(d2+d3), N(d1+d2+d3)  are achievable norms of elements with trace 0 in the maximal order (i.e. satisfy certain congruence conditions mod p)
# 3) |Tr(di*dj)| \leq floor(2*sqrt(Ndi*Ndj)) for all i,j=1,2,3
# 4) Ddidj = 4*Ndi*Ndj - Trdidj^2 are integers divisible by p for all i,j=1,2,3
# 5) Dd1d2d3 = 4*Nd1*Nd2*Nd3 - Nd1*Trd2d3^2 - Nd2*Trd1d2^2 - Nd3*Trd1d3^2 - Trd1d2*Trd1d3*Trd2d3 is an integer divisible by p and is a square
# 6) if Ndi = 0 then Trdidj = Trdjdi = 0 for all i,j=1,2,3.
#-------------------------------------------------------------------------------------------------------------------------------------------
#Know: d1,d2,d3 in max order O in quaternion algebra QA with Tr(d1)=Tr(d2)=Tr(d3)=0, N(d1)=Nd1, N(d2)=Nd2, N(d3)=Nd3, Tr(d1*d2)=Trd1d2, Tr(d1*d3)=Trd1d3, Tr(d2*d3)=Trd2d3.
#Know: Tr(d1*d2)= N(d1)+N(d2)-N(d1+d2) and Tr(d1*d2) + Tr(d1*d3) + Tr(d2*d3) = N(d1) + N(d2) + N(d3) - N(d1+d2+d3)
#Know: Tr(d1*d2*d3) = Tr(d2*d3*d1) = Tr(d3*d1*d2) = - Tr(d1*d3*d2) = - Tr(d3*d2*d1) = - Tr(d2*d1*d3).
#Check: N(d1)+N(d2)-Tr(d1*d2), N(d1)+N(d3)-Tr(d1*d3), N(d2)+N(d3)-Tr(d2*d3), N(d1) + N(d2) + N(d3) - Tr(d1*d2) - Tr(d1*d3) - Tr(d2*d3) in GoodNorms.
#-------------------------------------------------------------------------------------------------------------------------------------------

from sage.all_cmdline import *  
from time import time
from random import randint

import csv

#-------------------------------------------------
T = time()

filelocation = "./"
CMfieldnumber = sage_eval(sys.argv[1])
p = sage_eval(sys.argv[2])
D = sage_eval(sys.argv[3])


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


print ("CMfieldnumber = ",CMfieldnumber,", p = ",p,", D = ",D)
print (sys.version)
print ("----------------------------------------------------------------------")
sys.stdout.flush()


R = Integers(p)


def good_norms(p): #Lists congruence conditions for norms of elems with trace 0. See Explanation for good norms below.
    goodnorms = []
    if p == 3: #x = (1/2*t + w)*i + (-u)*j + 1/2*t*k => Nx = 3*u^2 + w^2 + w*t + t^2 => Nx = w^2 + w*t + t^2 = 0, 1 (mod 3) 
        goodnorms = [R(0), R(1)]
    elif p == 5: #x = (t/4 + u/2 + w)*i + t/2*j + (t/4 - u/2)*k. Nx = 2*t^2 - 2*t*u + 3*u^2 + t*w + 2*u*w + 2*w^2 = 0, 2, 3 (mod 5) 
        goodnorms = [R(0), R(2), R(3)]
    elif p == 7: #x = (1/2*t + w)*i + (-u)*j + 1/2*t*k => Nx = 7*u^2 + w^2 + w*t + 2*t^2 => Nx = w^2 + w*t + 2*t^2 = 0, 1, 2, 4 (mod 7)
        goodnorms = [R(0), R(1), R(2), R(4)]
    elif p == 11: #x = (1/2*t + w)*i + (-u)*j + 1/2*t*k => Nx = 11*u^2 + w^2 + w*t + 3*t^2 => Nx =  w^2 + w*t + 3*t^2 = 0, 1, 3, 4, 5, 9 (mod 11)
        goodnorms = [R(0), R(1), R(3), R(4), R(5), R(9)]
    elif p == 13:
        goodnorms = [R(0), R(2), R(5), R(6), R(7), R(8), R(11)]
    elif p == 19: #x = (1/2*t + w)*i + (-u)*j + 1/2*t*k => Nx = 19*u^2 + w^2 + w*t + 5*t^2 => Nx = w^2 + w*t + 5*t^2 = 0, 1, 4, 5, 6, 7, 9, 11, 16, 17 (mod 19)
        goodnorms = [R(0), R(1), R(4), R(5), R(6), R(7), R(9), R(11), R(16), R(17)]
    elif p == 23:
        goodnorms = [R(0), R(1), R(2), R(3), R(4), R(6), R(8), R(9), R(12), R(13), R(16), R(18)]
    elif p == 31: #x = (1/2*t + w)*i + (-u)*j + 1/2*t*k => Nx = 31*u^2 + w^2 + w*t + 8*t^2 => Nx = w^2 + w*t + 8*t^2 = 0, 1, 2, 4, 5, 7, 8, 9, 10, 14, 16, 18, 19, 20, 25, 28 (mod 31)
        goodnorms = [R(0), R(1), R(2), R(4), R(5), R(7), R(8), R(9), R(10), R(14), R(16), R(18), R(19), R(20), R(25), R(28)]
    elif p == 43: #x = (1/2*t + w)*i + (-u)*j + 1/2*t*k => Nx = 43*u^2 + w^2 + w*t + 11*t^2 => Nx = w^2 + w*t + 11*t^2 = 0, 1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17, 21, 23, 24, 25, 31, 35, 36, 38, 40, 41 (mod 43)
        goodnorms = [R(0), R(1), R(4), R(6), R(9), R(10), R(11), R(13), R(14), R(15), R(16), R(17), R(21), R(23), R(24), R(25), R(31), R(35), R(36), R(38), R(40), R(41)]
    elif p == 47:
        goodnorms = [R(0), R(1), R(2), R(3), R(4), R(6), R(7), R(8), R(9), R(12), R(14), R(16), R(17), R(18), R(21), R(24), R(25), R(27), R(28), R(32), R(34), R(36), R(37), R(42)]
    elif p == 53:
        goodnorms = [R(0), R(2), R(3), R(5), R(8), R(12), R(14), R(18), R(19), R(20), R(21), R(22), R(23), R(26), R(27), R(30), R(31), R(32), R(33), R(34), R(35), R(39), R(41), R(45), R(48), R(50), R(51)]
    elif p == 71:
        goodnorms = [R(0), R(1), R(2), R(3), R(4), R(5), R(6), R(8), R(9), R(10), R(12), R(15), R(16), R(18), R(19), R(20), R(24), R(25), R(27), R(29), R(30), R(32), R(36), R(37), R(38), R(40), R(43), R(45), R(48), R(49), R(50), R(54), R(57), R(58), R(60), R(64)]
    elif p == 79:
        goodnorms = [R(0), R(1), R(2), R(4), R(5), R(8), R(9), R(10), R(11), R(13), R(16), R(18), R(19), R(20), R(21), R(22), R(23), R(25), R(26), R(31), R(32), R(36), R(38), R(40), R(42), R(44), R(45), R(46), R(49), R(50), R(51), R(52), R(55), R(62), R(64), R(65), R(67), R(72), R(73), R(76)]
    elif p == 127:
        goodnorms = [R(0), R(1), R(2), R(4), R(8), R(9), R(11), R(13), R(15), R(16), R(17), R(18), R(19), R(21), R(22), R(25), R(26), R(30), R(31), R(32), R(34), R(35), R(36), R(37), R(38), R(41), R(42), R(44), R(47), R(49), R(50), R(52), R(60), R(61), R(62), R(64), R(68), R(69), R(70), R(71), R(72), R(73), R(74), R(76), R(79), R(81), R(82), R(84), R(87), R(88), R(94), R(98), R(99), R(100), R(103), R(104), R(107), R(113), R(115), R(117), R(120), R(121), R(122), R(124)]
    elif p == 863:
        goodnorms = [R(0), R(1), R(2), R(3), R(4), R(6), R(8), R(9), R(12), R(16), R(17), R(18), R(19), R(24), R(25), R(27), R(29), R(31), R(32), R(34), R(35), R(36), R(37), R(38), R(41), R(43), R(48), R(49), R(50), R(51), R(53), R(54), R(55), R(57), R(58), R(59), R(61), R(62), R(64), R(65), R(68), R(70), R(71), R(72), R(74), R(75), R(76), R(77), R(81), R(82), R(86), R(87), R(91), R(93), R(96), R(98), R(100), R(102), R(103), R(105), R(106), R(107), R(108), R(109), R(110), R(111), R(113), R(114), R(115), R(116), R(118), R(121), R(122), R(123), R(124), R(127), R(128), R(129), R(130), R(136), R(140), R(142), R(143), R(144), R(147), R(148), R(149), R(150), R(151), R(152), R(153), R(154), R(159), R(161), R(162), R(163), R(164), R(165), R(169), R(171), R(172), R(174), R(177), R(181), R(182), R(183), R(186), R(191), R(192), R(195), R(196), R(199), R(200), R(204), R(206), R(210), R(212), R(213), R(214), R(216), R(218), R(220), R(222), R(223), R(225), R(226), R(228), R(229), R(230), R(231), R(232), R(235), R(236), R(239), R(242), R(243), R(244), R(246), R(248), R(253), R(254), R(256), R(257), R(258), R(260), R(261), R(269), R(272), R(273), R(279), R(280), R(281), R(283), R(284), R(286), R(288), R(289), R(293), R(294), R(296), R(298), R(299), R(300), R(302), R(304), R(306), R(307), R(308), R(309), R(311), R(313), R(315), R(318), R(321), R(322), R(323), R(324), R(326), R(327), R(328), R(329), R(330), R(331), R(333), R(335), R(337), R(338), R(339), R(342), R(344), R(345), R(348), R(353), R(354), R(359), R(361), R(362), R(363), R(364), R(365), R(366), R(369), R(372), R(373), R(381), R(382), R(383), R(384), R(387), R(389), R(390), R(392), R(395), R(397), R(398), R(400), R(408), R(409), R(412), R(415), R(420), R(421), R(424), R(425), R(426), R(428), R(429), R(432), R(433), R(436), R(440), R(441), R(444), R(445), R(446), R(447), R(449), R(450), R(452), R(453), R(456), R(457), R(458), R(459), R(460), R(461), R(462), R(464), R(467), R(469), R(470), R(472), R(475), R(477), R(478), R(483), R(484), R(485), R(486), R(487), R(488), R(489), R(492), R(493), R(495), R(496), R(503), R(505), R(506), R(507), R(508), R(511), R(512), R(513), R(514), R(516), R(517), R(520), R(522), R(523), R(527), R(529), R(531), R(538), R(543), R(544), R(546), R(547), R(549), R(551), R(553), R(558), R(560), R(562), R(566), R(568), R(571), R(572), R(573), R(576), R(578), R(581), R(585), R(586), R(587), R(588), R(589), R(592), R(593), R(595), R(596), R(597), R(598), R(599), R(600), R(601), R(604), R(608), R(611), R(612), R(613), R(614), R(616), R(618), R(622), R(623), R(625), R(626), R(629), R(630), R(636), R(639), R(642), R(644), R(646), R(648), R(652), R(654), R(655), R(656), R(658), R(660), R(661), R(662), R(665), R(666), R(669), R(670), R(673), R(674), R(675), R(676), R(678), R(679), R(683), R(684), R(685), R(687), R(688), R(690), R(693), R(695), R(696), R(697), R(703), R(705), R(706), R(707), R(708), R(717), R(718), R(722), R(724), R(725), R(726), R(728), R(729), R(730), R(731), R(732), R(737), R(738), R(743), R(744), R(746), R(751), R(759), R(762), R(764), R(766), R(768), R(769), R(771), R(773), R(774), R(775), R(778), R(779), R(780), R(783), R(784), R(785), R(790), R(794), R(796), R(797), R(800), R(803), R(807), R(811), R(816), R(817), R(818), R(819), R(821), R(823), R(824), R(830), R(833), R(835), R(837), R(840), R(841), R(842), R(843), R(848), R(849), R(850), R(852), R(853), R(856), R(858)]
    elif p == 1013:
        goodnorms = [R(0), R(2), R(3), R(5), R(7), R(8), R(12), R(17), R(18), R(20), R(22), R(26), R(27), R(28), R(29), R(30), R(31), R(32), R(33), R(37), R(38), R(39), R(41), R(42), R(45), R(46), R(47), R(48), R(50), R(55), R(57), R(59), R(61), R(63), R(65), R(67), R(68), R(69), R(70), R(72), R(75), R(77), R(80), R(86), R(88), R(91), R(95), R(98), R(101), R(102), R(103), R(104), R(105), R(106), R(107), R(108), R(109), R(112), R(115), R(116), R(120), R(124), R(125), R(128), R(129), R(131), R(132), R(133), R(137), R(139), R(142), R(146), R(147), R(148), R(151), R(152), R(153), R(156), R(158), R(159), R(161), R(162), R(164), R(166), R(168), R(170), R(174), R(175), R(178), R(179), R(180), R(184), R(186), R(187), R(188), R(191), R(192), R(194), R(198), R(200), R(213), R(215), R(219), R(220), R(221), R(222), R(226), R(227), R(228), R(234), R(236), R(237), R(238), R(239), R(242), R(243), R(244), R(245), R(246), R(249), R(252), R(254), R(255), R(260), R(261), R(263), R(265), R(267), R(268), R(269), R(270), R(272), R(276), R(277), R(279), R(280), R(282), R(286), R(288), R(290), R(291), R(293), R(297), R(298), R(300), R(301), R(307), R(308), R(310), R(313), R(314), R(317), R(319), R(320), R(323), R(326), R(330), R(333), R(334), R(338), R(339), R(341), R(342), R(343), R(344), R(346), R(349), R(351), R(352), R(354), R(355), R(357), R(359), R(362), R(363), R(364), R(365), R(366), R(369), R(370), R(371), R(377), R(378), R(380), R(381), R(383), R(386), R(389), R(390), R(391), R(392), R(394), R(395), R(398), R(401), R(402), R(403), R(404), R(405), R(406), R(407), R(408), R(409), R(410), R(412), R(414), R(415), R(416), R(418), R(420), R(421), R(422), R(423), R(424), R(425), R(428), R(429), R(432), R(434), R(435), R(436), R(439), R(443), R(445), R(446), R(447), R(448), R(450), R(451), R(457), R(458), R(460), R(462), R(463), R(464), R(465), R(466), R(467), R(470), R(471), R(480), R(481), R(482), R(485), R(489), R(494), R(495), R(496), R(497), R(499), R(500), R(501), R(502), R(506), R(507), R(511), R(512), R(513), R(514), R(516), R(517), R(518), R(519), R(524), R(528), R(531), R(532), R(533), R(542), R(543), R(546), R(547), R(548), R(549), R(550), R(551), R(553), R(555), R(556), R(562), R(563), R(565), R(566), R(567), R(568), R(570), R(574), R(577), R(578), R(579), R(581), R(584), R(585), R(588), R(589), R(590), R(591), R(592), R(593), R(595), R(597), R(598), R(599), R(601), R(603), R(604), R(605), R(606), R(607), R(608), R(609), R(610), R(611), R(612), R(615), R(618), R(619), R(621), R(622), R(623), R(624), R(627), R(630), R(632), R(633), R(635), R(636), R(642), R(643), R(644), R(647), R(648), R(649), R(650), R(651), R(654), R(656), R(658), R(659), R(661), R(662), R(664), R(667), R(669), R(670), R(671), R(672), R(674), R(675), R(679), R(680), R(683), R(687), R(690), R(693), R(694), R(696), R(699), R(700), R(703), R(705), R(706), R(712), R(713), R(715), R(716), R(720), R(722), R(723), R(725), R(727), R(731), R(733), R(734), R(736), R(737), R(741), R(743), R(744), R(745), R(746), R(748), R(750), R(752), R(753), R(758), R(759), R(761), R(764), R(767), R(768), R(769), R(770), R(771), R(774), R(775), R(776), R(777), R(779), R(785), R(786), R(787), R(791), R(792), R(793), R(794), R(798), R(800), R(813), R(815), R(819), R(821), R(822), R(825), R(826), R(827), R(829), R(833), R(834), R(835), R(838), R(839), R(843), R(845), R(847), R(849), R(851), R(852), R(854), R(855), R(857), R(860), R(861), R(862), R(865), R(866), R(867), R(871), R(874), R(876), R(880), R(881), R(882), R(884), R(885), R(888), R(889), R(893), R(897), R(898), R(901), R(904), R(905), R(906), R(907), R(908), R(909), R(910), R(911), R(912), R(915), R(918), R(922), R(925), R(927), R(933), R(936), R(938), R(941), R(943), R(944), R(945), R(946), R(948), R(950), R(952), R(954), R(956), R(958), R(963), R(965), R(966), R(967), R(968), R(971), R(972), R(974), R(975), R(976), R(980), R(981), R(982), R(983), R(984), R(985), R(986), R(987), R(991), R(993), R(995), R(996), R(1001), R(1005), R(1006), R(1008), R(1010), R(1011)]
    elif p == 1087:
        [R(0), R(1), R(2), R(4), R(8), R(9), R(15), R(16), R(17), R(18), R(21), R(25), R(30), R(32), R(33), R(34), R(35), R(36), R(39), R(41), R(42), R(43), R(49), R(50), R(55), R(57), R(60), R(64), R(65), R(66), R(68), R(69), R(70), R(71), R(72), R(73), R(77), R(78), R(79), R(81), R(82), R(83), R(84), R(86), R(87), R(91), R(93), R(95), R(98), R(100), R(101), R(103), R(107), R(109), R(110), R(111), R(114), R(115), R(120), R(121), R(128), R(130), R(132), R(133), R(135), R(136), R(137), R(138), R(139), R(140), R(141), R(142), R(143), R(144), R(145), R(146), R(151), R(153), R(154), R(155), R(156), R(157), R(158), R(159), R(161), R(162), R(163), R(164), R(166), R(168), R(169), R(172), R(173), R(174), R(177), R(181), R(182), R(183), R(185), R(186), R(189), R(190), R(191), R(196), R(200), R(201), R(202), R(203), R(206), R(209), R(211), R(214), R(217), R(218), R(220), R(222), R(223), R(225), R(227), R(228), R(229), R(230), R(235), R(239), R(240), R(241), R(242), R(247), R(253), R(255), R(256), R(257), R(259), R(260), R(264), R(265), R(266), R(267), R(269), R(270), R(271), R(272), R(274), R(276), R(277), R(278), R(280), R(282), R(284), R(286), R(288), R(289), R(290), R(291), R(292), R(295), R(297), R(299), R(302), R(305), R(306), R(307), R(308), R(310), R(311), R(312), R(314), R(315), R(316), R(317), R(318), R(319), R(322), R(324), R(326), R(328), R(329), R(332), R(335), R(336), R(338), R(339), R(341), R(344), R(346), R(348), R(351), R(354), R(357), R(361), R(362), R(364), R(366), R(367), R(369), R(370), R(371), R(372), R(375), R(377), R(378), R(380), R(381), R(382), R(383), R(387), R(389), R(392), R(393), R(397), R(400), R(401), R(402), R(403), R(404), R(406), R(407), R(412), R(413), R(418), R(419), R(421), R(422), R(425), R(427), R(428), R(433), R(434), R(436), R(437), R(440), R(441), R(444), R(445), R(446), R(447), R(450), R(454), R(456), R(458), R(460), R(461), R(469), R(470), R(478), R(479), R(480), R(481), R(482), R(484), R(485), R(487), R(491), R(494), R(495), R(499), R(501), R(506), R(510), R(512), R(513), R(514), R(517), R(518), R(520), R(521), R(525), R(528), R(529), R(530), R(532), R(534), R(537), R(538), R(540), R(541), R(542), R(544), R(548), R(551), R(552), R(554), R(556), R(560), R(561), R(563), R(564), R(565), R(568), R(571), R(572), R(576), R(578), R(579), R(580), R(582), R(583), R(584), R(585), R(587), R(589), R(590), R(591), R(594), R(595), R(597), R(598), R(599), R(601), R(604), R(610), R(611), R(612), R(613), R(614), R(615), R(616), R(619), R(620), R(621), R(622), R(623), R(624), R(625), R(628), R(630), R(632), R(634), R(635), R(636), R(638), R(639), R(644), R(645), R(648), R(649), R(652), R(655), R(656), R(657), R(658), R(661), R(663), R(664), R(667), R(670), R(671), R(672), R(673), R(676), R(677), R(678), R(679), R(682), R(688), R(689), R(691), R(692), R(693), R(696), R(697), R(699), R(701), R(702), R(703), R(708), R(711), R(713), R(714), R(719), R(722), R(724), R(727), R(728), R(729), R(731), R(732), R(734), R(735), R(737), R(738), R(740), R(742), R(744), R(745), R(747), R(750), R(753), R(754), R(756), R(757), R(760), R(762), R(764), R(766), R(767), R(774), R(778), R(783), R(784), R(786), R(787), R(789), R(791), R(793), R(794), R(800), R(802), R(804), R(806), R(808), R(812), R(814), R(819), R(824), R(825), R(826), R(829), R(833), R(835), R(836), R(837), R(838), R(839), R(841), R(842), R(843), R(844), R(849), R(850), R(851), R(853), R(854), R(855), R(856), R(861), R(863), R(866), R(868), R(871), R(872), R(874), R(875), R(877), R(879), R(880), R(882), R(883), R(888), R(889), R(890), R(892), R(893), R(894), R(895), R(899), R(900), R(903), R(907), R(908), R(909), R(911), R(912), R(916), R(917), R(920), R(922), R(927), R(935), R(937), R(938), R(939), R(940), R(953), R(956), R(958), R(960), R(961), R(962), R(963), R(964), R(965), R(968), R(969), R(970), R(971), R(974), R(975), R(979), R(981), R(982), R(983), R(985), R(988), R(990), R(991), R(993), R(995), R(997), R(998), R(999), R(1002), R(1007), R(1011), R(1012), R(1013), R(1020), R(1024), R(1025), R(1026), R(1028), R(1029), R(1031), R(1033), R(1034), R(1035), R(1036), R(1039), R(1040), R(1041), R(1042), R(1043), R(1047), R(1049), R(1050), R(1056), R(1058), R(1059), R(1060), R(1061), R(1063), R(1064), R(1065), R(1067), R(1068), R(1073), R(1074), R(1075), R(1076), R(1077), R(1080), R(1081), R(1082), R(1084)]
    else:
        goodnorms = list(R)
    return goodnorms


GoodNorms = good_norms(p)

@parallel
def find_potential_solutions_new_and_elim_bad(reader, csvwriter, csvfile):
    tm = time()
    OUTPUT = []
    for row in reader: #x, dOVERn, a, cOVERn, b, gamma, n
        x = QQ(row["x"])
        dOVERn = QQ(row["dOVERn"])
        a = QQ(row["a"])
        cOVERn = QQ(row["cOVERn"])
        b = QQ(row["b"])
        gamma = QQ(row["gamma"])
        n = QQ(row["n"])
        gOVERn = QQ(gamma/n)
        aOVERn = QQ(a/n)
        bOVERn = QQ(b/n)
        k1 = QQ(aOVERn*dOVERn - aOVERn*x - bOVERn)
        k2 = QQ(aOVERn*cOVERn + bOVERn*x)
        k3 = QQ(k1*cOVERn - k2*dOVERn)
        k4 = QQ(k2*cOVERn - bOVERn*b)
        k5 = QQ(aOVERn*a - k1*dOVERn - k2)
        k6 = QQ(aOVERn*b + k1*cOVERn)
        l = dOVERn*x - x^2 + cOVERn
        l1 = a*cOVERn + b*x
        l2 = a*dOVERn - a*x - b
        l3 = l - a
        for Nd1 in [0..D]: # 
            if Nd1 == 0: # see Explanation for case Nd1 = 0 below
                Trd1d2 = 0
                Trd1d3 = 0
                Nd2 = -(l-a)*D #Nd2 = -(dOVERn*x - x^2 + cOVERn - a)*D
                Nd3 = (cOVERn*l + b*x)*D #Nd3 = (cOVERn*(dOVERn*x - x^2 + cOVERn) + b*x)*D
                Trd2d3 = (dOVERn*l - a*x - b)*D #Trd2d3 = (dOVERn*(dOVERn*x - x^2 + cOVERn) - a*x - b)*D
                if Nd2 >=0 and Nd2 in ZZ and Nd3>=0 and Nd3 in ZZ and R(Nd3) in GoodNorms and R(Nd2) in GoodNorms:
                    bd23 = floor(2*sqrt(Nd2*Nd3))
                    Dd2d3 = 4*Nd2*Nd3 - Trd2d3^2 
                    Dd1d2 = 0
                    Dd1d3 = 0
                    Dd1d2d3 = 0
                    if abs(Trd2d3) <= bd23 and Dd2d3/p in ZZ and R(Nd2+Nd3-Trd2d3) in GoodNorms:
                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                        csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                        csvfile.flush()
                        OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                        sys.stdout.flush()
            else:  # Nd1 != 0
                if R(Nd1) in GoodNorms:
                    upNd2 = a*(D - Nd1)
                    l12 = l1*(D - Nd1)
                    l22 = l2*(D - Nd1)
                    l32 = l3*(D - Nd1) 
                    for Nd2 in [0..upNd2]:
                        if Nd2 == 0: # see Explanation for case Nd1 != 0 and Nd2 = 0:
                            Trd1d2 = 0
                            Trd2d3 = 0
                            if k1 != 0 and Nd1 == D:
                                Nd3 = 0
                                Trd1d3 = 0
                                Dd1d2 = 0
                                Dd2d3 = 0
                                Dd1d3 = 0
                                Dd1d2d3 = 0
                                print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                                csvfile.flush()
                                OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                sys.stdout.flush()
                            elif k1 == 0:
                                Nd3 = (D-Nd1)*n/a
                                if R(Nd3) in GoodNorms:
                                    Trd1d3 = - k5*Nd3
                                    bd13 = floor(2*sqrt(Nd1*Nd3))
                                    Dd1d2 = 0
                                    Dd2d3 = 0
                                    Dd1d3 = 4*Nd1*Nd3 - Trd1d3^2 
                                    Dd1d2d3 = 0
                                    if abs(Trd1d3) <= bd13 and Dd1d3/p in ZZ and R(Nd1+Nd3-Trd1d3) in GoodNorms:
                                        print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                        csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                                        csvfile.flush()
                                        OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                        sys.stdout.flush()
                        else: # Nd1 !=0 and Nd2 !=0
                            if R(Nd2) in GoodNorms:
                                l123 = l12 - cOVERn*Nd2
                                l223 = l22 - dOVERn*Nd2
                                l323 = l32 + Nd2
                                bd12 = floor(2*sqrt(Nd1*Nd2))
                                for Trd1d2 in [-bd12..bd12]: # see Explanation for Nd1 !=0 and Nd2 !=0
                                    Dd1d2 = 4*Nd1*Nd2 - Trd1d2^2
                                    Nd3 = l123 + b*Trd1d2 # Nd3 = (a*cOVERn + b*x)*(D-Nd1) + Trd1d2*b - cOVERn*Nd2
                                    if Nd3 >=0 and Nd3 in ZZ and R(Nd3) in GoodNorms  and Dd1d2/p in ZZ:   
                                        upNd3 = min(floor(2/a*(n*(D-Nd1)+(b^2-n)*Nd2/a)),gamma*(D-Nd1))
                                        if Nd3 <= upNd3:
                                            Trd2d3 = l223 - a*Trd1d2 # Trd2d3 = (a*dOVERn - a*x - b)*(D-Nd1) - dOVERn*Nd2 - Trd1d2*a
                                            if (Nd3 == 0 and Trd2d3 == 0) or (Nd3 != 0):
                                                Dd2d3 = 4*Nd2*Nd3 - Trd2d3^2
                                                bd23 = floor(2*sqrt(Nd2*Nd3))
                                                upNd2 = floor(2/gamma*(n*(x-Nd1)+(b^2-n)*Nd3/gamma))
                                                if Trd2d3 in ZZ and abs(Trd2d3) <= bd23 and Dd2d3/p in ZZ and R(Nd2+Nd3-Trd2d3) in GoodNorms and Nd2 <= upNd2:
                                                    Trd1d3 = l323 + (dOVERn - x)*Trd1d2 #Trd1d3 = (dOVERn*x - x^2 - a + cOVERn)*(D - Nd1) + Nd2 + (dOVERn - x)*Trd1d2
                                                    if (Nd3 == 0 and Trd1d3 == 0) or (Nd3 != 0):
                                                        Dd1d3 = 4*Nd1*Nd3 - Trd1d3^2
                                                        bd13 = floor(2*sqrt(Nd1*Nd3))        
                                                        if Trd1d3 in ZZ and abs(Trd1d3) <= bd13 and Dd1d3/p in ZZ and R(Nd1+Nd3-Trd1d3) in GoodNorms and R(Nd1+Nd2+Nd3-Trd1d2-Trd2d3-Trd1d3) in GoodNorms:
                                                            Dd1d2d3 = 4*Nd1*Nd2*Nd3 - Nd1*Trd2d3^2 - Nd2*Trd1d3^2 - Nd3*Trd1d2^2 - Trd1d2*Trd1d3*Trd2d3 # = (Trd1d2d3)^2
                                                            if Dd1d2d3.is_square() and Dd1d2d3/p in ZZ:
                                                                print ('[[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]] = ', [[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd2, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                                                csvwriter.writerow([x, dOVERn, a, cOVERn, b, gamma, n, Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3, Dd1d2, Dd1d3, Dd2d3, Dd1d2d3])
                                                                csvfile.flush()
                                                                OUTPUT.append([[x, dOVERn, a, cOVERn, b, gamma, n], [Nd1, Nd2, Nd3, Trd1d2, Trd1d3, Trd2d3], [Dd1d2, Dd1d3, Dd2d3, Dd1d2d3]])
                                                                sys.stdout.flush()
    print ("done")
    print (len(OUTPUT)," solutions found")
    print ("total time: ", time() - tm)
    print ("----------------------------------------------------------------------")
    sys.stdout.flush()
    return OUTPUT


print ("reading solutions from Step 3...") #x, dOVERn, a, cOVERn, b, gamma, n
sys.stdout.flush()
tread = time()
filename3 = "CMfield"+CMfieldnumber.str()+"Step3"+"p"+p.str()+".csv"
reader3 = csv.DictReader(open(filelocation+filename3, 'r'))
print ("Reading time = ", time()-tread)


#write all solutions in csv files
filename2 = "CMfield"+CMfieldnumber.str()+"Step4"+"p"+p.str()+".csv"
csvfile2 = open(filelocation+filename2, 'w')
csvwriter2 = csv.writer(csvfile2)
csvwriter2.writerow(('x', 'dOVERn', 'a', 'cOVERn', 'b', 'gamma', 'n', 'Nd1', 'Nd2', 'Nd3', 'Trd1d2', 'Trd1d3', 'Trd2d3', 'Dd1d2', 'Dd1d3', 'Dd2d3', 'Dd1d2d3'))

potsolutions = find_potential_solutions_new_and_elim_bad(reader3, csvwriter2, csvfile2)

csvfile2.close()

print ("All done in ", time()- T)


#------------------------------------------------------------------------------------------------------------------
# Explanation for good norms:
# Choose a prime p.
# if p == 2:
#     QA = QuaternionAlgebra(SR, -1,-1, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
# else:
#     S = Integers(8)
#     if S(p) == 3 or S(p) == 7:
#         QA = QuaternionAlgebra(SR, -1,-p, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
#     elif S(p) == 5:
#         QA = QuaternionAlgebra(SR, -2,-p, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
#     elif S(p) == 1:
#         print ("Choose a prime l = 3 mod 4 and l not a square mod p. Choose QA(-1,-l).") #QA = QuaternionAlgebra(SR, -1,-l, names=('i', 'j', 'k',) ); (i, j, k,) = QA._first_ngens(3)
# #Use Magma to find basis for the maximal order: 
# p := 11; //modify accordingly
# Q := RationalField();
# A<i,j,k>:=QuaternionAlgebra<Q|-1,-p>; //modify accordingly
# CC:=ConjugacyClasses(MaximalOrder(A));
# #CC;
# for g in [1..#CC] do
# Generators(CC[g]);
# end for;
# #Back to sage: for example, the basis is [ 1, i, 1/2*i + 1/2*k, 1/2 + 1/2*j ].
# u,v,w,t = var('u,v,w,t')
# x = QA(u*1 + v*i + w*(1/2*i + 1/2*k) + t*(1/2 + 1/2*j))
# x.reduced_trace()
# #But x has Tr(x)=0 => t + 2*u = 0 => t = -2*u 
# t = -2*u
# x = QA(u*1 + v*i + w*(1/2*i + 1/2*k) + t*(1/2 + 1/2*j))
# x.reduced_norm().expand()
# #This gives Nx = p*u^2 + (p+1)/4*p*w^2 + v^2 + v*w. 
# #For example, if p = 11, then Nx = 11*u^2 + 3*w^2 + v^2 + v*w.
# R = Integers(p) 
# good_norms = []
# for u in R:
#     for v in R:
#         for w in R:
#             N = R(11*u^2 + 3*w^2 + v^2 + v*w) #change this accordingly
#             if N not in good_norms:
#                 good_norms += [N]
#                 
# good_norms.sort()

#------------------------------------------------------------------------------------------------------------------
# Explanation for case Nd1 = 0:
# x, dOVERn,a,cOVERn,b,gamma,n,D = var('x, dOVERn,a,cOVERn,b,gamma,n,D')
# aOVERn = a/n
# bOVERn = b/n
# gammaOVERn = gamma/n
# gamma = a*cOVERn + b*dOVERn
# k1 = aOVERn*dOVERn - aOVERn*x - bOVERn
# k2 = aOVERn*cOVERn + bOVERn*x
# A = matrix(SR,3,3,[-gammaOVERn,-aOVERn,-bOVERn, k1*cOVERn-k2*dOVERn,k1,-k2, k2*cOVERn-b*bOVERn, -a*aOVERn+k1*dOVERn+k2, -a*bOVERn-k1*cOVERn])
# B = matrix(SR,3,1,[-D,0,0])
# A.determinant().expand().factor()
# = (a*cOVERn*dOVERn*x + b*dOVERn^2*x - a*cOVERn*x^2 - b*dOVERn*x^2 + a*cOVERn^2 + b*cOVERn*dOVERn - dOVERn*gamma*x + gamma*x^2 - b^2 + a*gamma - cOVERn*gamma)*(a^2*cOVERn + a*b*dOVERn - b^2)/n^3
# (a*cOVERn*dOVERn*x + b*dOVERn^2*x - a*cOVERn*x^2 - b*dOVERn*x^2 + a*cOVERn^2 + b*cOVERn*dOVERn - dOVERn*gamma*x + gamma*x^2 - b^2 + a*gamma - cOVERn*gamma).expand()
# = a^2*cOVERn + a*b*dOVERn - b^2
# Then det(A) = (a^2*cOVERn + a*b*dOVERn - b^2)^2/n^3.
# Use b^2= a*gamma - n = a^2*cOVERn + a*b*dOVERn - n. Thus, n = a^2*cOVERn + a*b*dOVERn - b^2.
# Then det(A)=1/n. So A is always invertible.
# S=(A.inverse()*B)
# Nd2=S[0]
# Nd3=S[1]
# Trd2d3=S[2] 
# S[0]
# = (-D*(((cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*n/gamma + (dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*(dOVERn*(a*dOVERn/n - a*x/n - b/n) - (cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*a/gamma - a^2/n + a*cOVERn/n + b*x/n)*n/(((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n)*gamma))*(b/gamma - a*((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*b/gamma - a*cOVERn/n - b*x/n)/(((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n)*gamma))/(cOVERn*(a*dOVERn/n - a*x/n - b/n) + (dOVERn*(a*dOVERn/n - a*x/n - b/n) - (cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*a/gamma - a^2/n + a*cOVERn/n + b*x/n)*((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*b/gamma - a*cOVERn/n - b*x/n)/((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n) + (cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*b/gamma + a*b/n) + (dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a*n/(((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*a/gamma + a*dOVERn/n - a*x/n - b/n)*gamma^2) - n/gamma))
# (copy what S[0] displays).expand().factor()
# = -(dOVERn*x - x^2 - a + cOVERn)*D*n/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Nd2 = -(dOVERn*x - x^2 + cOVERn - a)*D
# Similarly, 
# S[1] gives (cOVERn*dOVERn*x - cOVERn*x^2 + cOVERn^2 + b*x)*D*n/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Nd3 = (cOVERn*dOVERn*x - cOVERn*x^2 + cOVERn^2 + b*x)*D
# Thus, Nd3 = (cOVERn*(dOVERn*x - x^2 + cOVERn) + b*x)*D = cOVERn*(dOVERn*x - x^2 + cOVERn)*D + b*x*D = cOVERn*(a*D-Nd2) + b*x*D = -cOVERn*Nd2 + (cOVERn*a + b*x)*D
# and
# S[2] gives (dOVERn^2*x - dOVERn*x^2 + cOVERn*dOVERn - a*x - b)*D*n/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Trd2d3 = (dOVERn^2*x - dOVERn*x^2 + cOVERn*dOVERn - a*x - b)*D
# Thus, Trd2d3 = (dOVERn*(dOVERn*x - x^2 + cOVERn) - a*x - b)*D = dOVERn*(a*D-Nd2) - (a*x + b)*D = -dOVERn*Nd2 + (dOVERn*a - a*x - b)*D

#------------------------------------------------------------------------------------------------------------------
# Explanation for case Nd1 != 0 and Nd2 = 0:
# Then Trd1d2 = Trd2d3 = 0
# The 3 equations become: D = Nd1 + (a/n)*Nd3, 0 = k1*Nd3, 0 = - (a^2/n - k1*d/n - k2)*Nd3 - Trd1d3 = -k5*Nd3 - Trd1d3.
# If k1 = 0 then Nd3 = (D-Nd1)*n/a and Trd1d3 = - k5*Nd3.
# Otherwise, if k1 !=0 (which is actually always the case), then Nd3 = 0, D = Nd1, Trd1d3 = 0.

#------------------------------------------------------------------------------------------------------------------
# Explanation for Nd1 !=0 and Nd2 !=0
# x, dOVERn,a,cOVERn,b,gamma,n,D = var('x, dOVERn,a,cOVERn,b,gamma,n,D')
# aOVERn = a/n
# bOVERn = b/n
# gammaOVERn = gamma/n
# gamma = a*cOVERn + b*dOVERn
# k1 = aOVERn*dOVERn - aOVERn*x - bOVERn
# k2 = aOVERn*cOVERn + bOVERn*x
# A = matrix(SR,3,3,[aOVERn,bOVERn,0, k1,-k2,0, -a*aOVERn+k1*dOVERn+k2, -a*bOVERn-k1*cOVERn,-1])
# Nd1,Nd2,Trd1d2 = var('Nd1,Nd2,Trd1d2')
# B = matrix(SR,3,1,[D-Nd1-gammaOVERn*Nd2, -(k1*cOVERn-k2*dOVERn)*Nd2 + Trd1d2, -(k2*cOVERn-b*bOVERn)*Nd2])
# A.determinant().expand().factor()
# = (a^2*cOVERn + a*b*dOVERn - b^2)/n^2
# Use b^2= a*gamma - n = a^2*cOVERn + a*b*dOVERn - n. Thus, n = a^2*cOVERn + a*b*dOVERn - b^2.
# Then det(A)= 1/n. So A is always invertible.
# S = (A.inverse()*B)
# Nd3 = S[0]
# Trd2d3 = S[1]
# Trd1d3 = S[2] = k4*Nd2 - k5*Nd3 - k6*Trd2d3 - b

# S[0]
# = ((D - Nd1 - Nd2*gamma/n)*(n/a - b*(a*dOVERn/n - a*x/n - b/n)*n/(a^2*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))) + ((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*Nd2 + Trd1d2)*b/(a*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n)))
# (copy what S[0] displays).expand().factor()
# = -(Nd2*a^2*cOVERn^2 + Nd2*a*b*cOVERn*dOVERn - Nd2*b^2*cOVERn - D*a*cOVERn*n + Nd1*a*cOVERn*n - D*b*n*x + Nd1*b*n*x - Trd1d2*b*n)/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Nd3 = (D*a*cOVERn*n + D*b*n*x + Trd1d2*b*n - (a*cOVERn*n + b*n*x)*Nd1 - (a^2*cOVERn^2 + a*b*cOVERn*dOVERn - b^2*cOVERn)*Nd2)/n
# = ((a*cOVERn + b*x)*D*n + Trd1d2*b*n - (a*cOVERn + b*x)*n*Nd1 - (a^2*cOVERn + a*b*dOVERn - b^2)*cOVERn*Nd2)/n
# = ((a*cOVERn + b*x)*D*n + Trd1d2*b*n - (a*cOVERn + b*x)*n*Nd1 - n*cOVERn*Nd2)/n
# = (a*cOVERn + b*x)*D + Trd1d2*b - (a*cOVERn + b*x)*Nd1 - cOVERn*Nd2
# Thus, Nd3 = (a*cOVERn + b*x)*(D - Nd1) - cOVERn*Nd2 + Trd1d2*b
# Put l1 = a*cOVERn + b*x
# Thus, Nd3 = l1*(D-Nd1) - cOVERn*Nd2 + Trd1d2*b
# Put l12 = l1*(D-Nd1)
# Thus, Nd3 = l12 - cOVERn*Nd2 + Trd1d2*b
# Put l123 = l12 - cOVERn*Nd2
# Thus, Nd3 = l123 + b*Trd1d2

# S[1]
# = ((D - Nd1 - Nd2*gamma/n)*(a*dOVERn/n - a*x/n - b/n)*n/(a*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n)) - ((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*Nd2 + Trd1d2)/(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))
# (copy what S[1] displays).expand().factor()
# = -(Nd2*a^2*cOVERn*dOVERn + Nd2*a*b*dOVERn^2 - Nd2*b^2*dOVERn - D*a*dOVERn*n + Nd1*a*dOVERn*n + D*a*n*x - Nd1*a*n*x + Trd1d2*a*n + D*b*n - Nd1*b*n)/(a^2*cOVERn + a*b*dOVERn - b^2)
# Thus, Trd2d3 = -(Nd2*a^2*cOVERn*dOVERn + Nd2*a*b*dOVERn^2 - Nd2*b^2*dOVERn - D*a*dOVERn*n + Nd1*a*dOVERn*n + D*a*n*x - Nd1*a*n*x + Trd1d2*a*n + D*b*n - Nd1*b*n)/n
# = ((a*dOVERn - a*x - b)*D*n - (a*dOVERn - a*x - b)*n*Nd1 - (a^2*cOVERn + a*b*dOVERn - b^2)*dOVERn*Nd2 - Trd1d2*a*n)/n
# = ((a*dOVERn - a*x - b)*D*n - (a*dOVERn - a*x - b)*n*Nd1 - n*dOVERn*Nd2 - Trd1d2*a*n)/n
# = (a*dOVERn - a*x - b)*D - (a*dOVERn - a*x - b)*Nd1 - dOVERn*Nd2 - Trd1d2*a
# Thus, Trd2d3 = (a*dOVERn - a*x - b)*(D - Nd1) - dOVERn*Nd2 - Trd1d2*a
# Put l2 = a*dOVERn - a*x - b
# Thus, Trd2d3 = l2*(D-Nd1) - dOVERn*Nd2 - Trd1d2*a
# Put l22 = l2*(D - Nd1)
# Thus, Trd2d3 = l22 - dOVERn*Nd2 - a*Trd1d2 
# Put l223 = l22 - dOVERn*Nd2
# Thus, Trd2d3 = l223 - a*Trd1d2 

# S[2]
# = ((cOVERn*(a*cOVERn/n + b*x/n) - b^2/n)*Nd2 + (D - Nd1 - Nd2*gamma/n)*((dOVERn*(a*dOVERn/n - a*x/n - b/n) - a^2/n + a*cOVERn/n + b*x/n)*n/a - (cOVERn*(a*dOVERn/n - a*x/n - b/n) + (dOVERn*(a*dOVERn/n - a*x/n - b/n) - a^2/n + a*cOVERn/n + b*x/n)*b/a + a*b/n)*(a*dOVERn/n - a*x/n - b/n)*n/(a*(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))) + ((dOVERn*(a*cOVERn/n + b*x/n) - cOVERn*(a*dOVERn/n - a*x/n - b/n))*Nd2 + Trd1d2)*(cOVERn*(a*dOVERn/n - a*x/n - b/n) + (dOVERn*(a*dOVERn/n - a*x/n - b/n) - a^2/n + a*cOVERn/n + b*x/n)*b/a + a*b/n)/(b*(a*dOVERn/n - a*x/n - b/n)/a + a*cOVERn/n + b*x/n))
# (copy what S[2] displays).expand().factor()
# = (Nd2*a^2*cOVERn + Nd2*a*b*dOVERn + D*dOVERn*n*x - Nd1*dOVERn*n*x - D*n*x^2 + Nd1*n*x^2 - Nd2*b^2 - D*a*n + Nd1*a*n + D*cOVERn*n - Nd1*cOVERn*n + Trd1d2*dOVERn*n - Trd1d2*n*x)/n
# Thus, Trd1d3 = (Nd2*a^2*cOVERn + Nd2*a*b*dOVERn + D*dOVERn*n*x - Nd1*dOVERn*n*x - D*n*x^2 + Nd1*n*x^2 - Nd2*b^2 - D*a*n + Nd1*a*n + D*cOVERn*n - Nd1*cOVERn*n + Trd1d2*dOVERn*n - Trd1d2*n*x)/n
# = ((dOVERn*x - x^2 - a + cOVERn)*D*n - (dOVERn*x - x^2 - a + cOVERn)*n*Nd1 + (a^2*cOVERn + a*b*dOVERn - b^2)*Nd2 + (dOVERn*n - n*x)*Trd1d2)/n
# = ((dOVERn*x - x^2 - a + cOVERn)*D*n - (dOVERn*x - x^2 - a + cOVERn)*n*Nd1 + n*Nd2 + (dOVERn - x)*n*Trd1d2)/n
# Thus, Trd1d3 = (dOVERn*x - x^2 - a + cOVERn)*(D - Nd1) + Nd2 + (dOVERn - x)*Trd1d2
# Put l3 = dOVERn*x - x^2 - a + cOVERn = l - a
# Thus, Trd1d3 = l3*(D - Nd1) + Nd2 + (dOVERn - x)*Trd1d2
# Put l32 = l3*(D - Nd1) 
# Thus, Trd1d3 = l32 + Nd2 + (dOVERn - x)*Trd1d2
# Put l323 = l32 + Nd2
# Thus, Trd1d3 = l323 + (dOVERn - x)*Trd1d2
#------------------------------------------------------------------------------------------------------------------
