    function [AFCData_qi] = AFC_dataqi_v1()
     AFCData_qi(1) = 0.848; %CH3
     AFCData_qi(2)  = 0.54 ;%CH2
     AFCData_qi(3)  = 0.228 ; %CH
     AFCData_qi(4)  = 0 ;%C
     AFCData_qi(5)  = 1.176 ;%CH2=CH
     AFCData_qi(6)  = 0.867 ;%CH=CH
     AFCData_qi(7)  = 0.988 ;%CH2=C
     AFCData_qi(8)  = 0.676 ;%CH=C
     AFCData_qi(9)  = 0.485 ;%C=C
     AFCData_qi(10)  = 0.4 ;%ACH
     AFCData_qi(11)  = 0.12 ;%AC
     AFCData_qi(12)  = 0.968 ;%ACCH3
     AFCData_qi(13)  = 0.66 ;%ACCH2
     AFCData_qi(14)  = 0.348 ;%ACCH
     AFCData_qi(15)  = 1.2 ;%OH
     AFCData_qi(16)  = 1.432 ;%CH3OH
     AFCData_qi(17)  = 1.4 ;%H2O
     AFCData_qi(18)  = 0.68 ;%ACOH
     AFCData_qi(19)  = 1.488 ;%CH3CO
     AFCData_qi(20)  = 1.18 ;%CH2CO
     AFCData_qi(21)  = 0.948 ;%CHO
     AFCData_qi(22)  = 1.728 ;%CH3COO
     AFCData_qi(23)  = 1.42 ;%CH2COO
     AFCData_qi(24)  = 1.188 ;%HCOO
     AFCData_qi(25)  = 1.088 ;%CH3O
     AFCData_qi(26)  = 0.78 ;%CH2O
     AFCData_qi(27)  = 0.468 ;%CHO
     AFCData_qi(28)  = 1.1 ;%THF
     AFCData_qi(29)  = 1.544 ;%CH3NH2
     AFCData_qi(30)  = 1.236 ;%CH2NH2
     AFCData_qi(31)  = 0.924 ;%CHNH2
     AFCData_qi(32)  = 1.244 ; %CH3NH
     AFCData_qi(33)  = 0.936 ;%CH2NH
     AFCData_qi(34)  = 0.624 ;%CHNH
     AFCData_qi(35)  = 0.94 ;%CH3N
     AFCData_qi(36)  = 0.632 ;%CH2N
     AFCData_qi(37)  = 0.816 ;%ACNH2
     AFCData_qi(38)  = 2.113 ;%C5H5N
     AFCData_qi(39)  = 1.833 ;%C5H4N
     AFCData_qi(40)  = 1.553 ; %C5H3N
     AFCData_qi(41)  = 1.724 ;%CH3CN
     AFCData_qi(42)  = 1.416 ;%CH2CN
     AFCData_qi(43)  = 1.224 ;%COOH
     AFCData_qi(44)  = 1.532 ;%HCOOH
     AFCData_qi(45)  = 1.264 ;%CH2CL
     AFCData_qi(46)  = 0.952 ;%CHCL
     AFCData_qi(47)  = 0.724 ;%CCL
     AFCData_qi(48)  = 1.988 ;%CH2CL2
     AFCData_qi(49)  = 1.684 ;%CHCL2
     AFCData_qi(50)  = 1.448 ;%CCL2
     AFCData_qi(51)  = 2.41 ;%CHCL3
     AFCData_qi(52)  = 2.184 ;%CCL3
     AFCData_qi(53)  = 2.91 ;%CCL4
     AFCData_qi(54)  = 0.844 ;%ACCL
     AFCData_qi(55)  = 1.868 ;%CH3NO2
     AFCData_qi(56)  = 1.56 ;%CH2NO2
     AFCData_qi(57)  = 1.248 ;%CHNO2
     AFCData_qi(58)  = 1.104 ;%ACNO2
     AFCData_qi(59)  = 1.65 ;%CS2
     AFCData_qi(60)  = 1.676 ;%CH3SH
     AFCData_qi(61)  = 1.368 ;%CH2SH
     AFCData_qi(62)  = 2.481 ;%FURFURAL
     AFCData_qi(63)  = 2.248 ;%DOH
     AFCData_qi(64)  = 0.992 ;%I
     AFCData_qi(65)  = 0.832 ;%Br
     AFCData_qi(66)  = 1.088 ;%CH=-C
     AFCData_qi(67)  = 0.784 ;%C=-C
     AFCData_qi(68)  = 2.472 ;%DMSO
     AFCData_qi(69)  = 2.052 ;%Acrylonitrate
     AFCData_qi(70)  = 0.724 ;%Cl-(C=C)
     AFCData_qi(71)  = 0.524 ;%ACF
     AFCData_qi(72)  = 2.736 ;%DMF
     AFCData_qi(73)  = 2.12 ;%HCON(CH2)2
     AFCData_qi(74)  = 1.38 ;%CF3
     AFCData_qi(75)  = 0.92 ;%CF2
     AFCData_qi(76)  = 0.46 ;%CF
     AFCData_qi(77)  = 1.2 ;%COO
     AFCData_qi(78)  = 1.2632 ;%SiH3
     AFCData_qi(79)  = 1.0063 ;%SiH2
     AFCData_qi(80)  = 0.7494 ;%SiH
     AFCData_qi(81)  = 0.4099 ;%Si
     AFCData_qi(82)  = 1.0621 ;%SiH2O
     AFCData_qi(83)  = 0.7639 ;%SiHO
     AFCData_qi(84)  = 0.4657 ; %SiO
     AFCData_qi(85)  = 0.64 ; %CO
     %extension to new functional groups not covered in previous unifac
     %versions
     AFCData_qi(86)  = 3.2  ; %NMP
     AFCData_qi(87)  = 2.644  ; %CCl3F
     AFCData_qi(88)  = 1.916  ; %CCl2F
     AFCData_qi(89)  = 2.116  ; %HCCl2F
     AFCData_qi(90)  = 1.416  ; %HCClF
     AFCData_qi(91)  = 1.648  ; %CClF2
     AFCData_qi(92)  = 1.828  ; %HCClF2
     AFCData_qi(93)  = 2.1  ; %CClF3
     AFCData_qi(94)  = 2.376  ; %CCl2F2
     AFCData_qi(95)  = 1.248  ; %CONH2
     AFCData_qi(96)  = 1.796  ; %CONHCH3
     AFCData_qi(97)  = 1.488  ; %CONHCH2
     AFCData_qi(98)  = 2.428  ; %CON(CH3)2
     AFCData_qi(99)  = 2.12  ; %CONCH3CH2
     AFCData_qi(100)  = 1.812  ; %CON(CH2)2
     AFCData_qi(101)  = 1.904  ; %C2H5O2
     AFCData_qi(102)  = 1.592  ; %C2H4O2
     AFCData_qi(103)  = 1.368  ; %CH3S
     AFCData_qi(104)  = 1.060  ; %CH2S
     AFCData_qi(105)  = 0.748  ; %CHS
     AFCData_qi(106)  = 2.796  ; %MORPH
     AFCData_qi(107)  = 2.140  ; %C4H4S
     AFCData_qi(108)  = 1.860  ; %C4H3S
     AFCData_qi(109)  = 1.580 ; %C4H2S     
     %--- Electrolytes ----
     AFCData_qi(110)  = 2.7  ; %H +%Zeund et al 2008 - dynamic hydration
     AFCData_qi(111)  = 0.78  ; %NH4 +%Zeund et al 2008 - dynamic hydration
     AFCData_qi(112)  = 0.62  ; %Na +%Zeund et al 2008 - dynamic hydration
     AFCData_qi(113)  = 0.58  ; %K +%Zeund et al 2008 - dynamic hydration
     AFCData_qi(114)  = 8.35  ; %Mg 2+%Zeund et al 2008 - dynamic hydration
     AFCData_qi(115)  = 3.4  ; %Ca 2+%Zeund et al 2008 - dynamic hydration
     AFCData_qi(116)  = 0.99  ; %Li + %Zeund et al 2008 - dynamic hydration
     AFCData_qi(117)  = 3.0  ; %Ba 2+
     AFCData_qi(118)  = 1.0  ; %Sr 2+
     AFCData_qi(119)  = 1.0  ; %Co 2+
     AFCData_qi(120)  = 1.0  ; %Ni 2+
     AFCData_qi(121)  = 1.0  ; %Cu 2+
     AFCData_qi(122)  = 1.0  ; %Zn 2+
     AFCData_qi(123)  = 3.0  ; %Hg 2+
     AFCData_qi(124)  = 3.96  ; %SO4 2-%Zeund et al 2008 - dynamic hydration
     AFCData_qi(125)  = 0.97  ; %NO3 - %Zeund et al 2008 - dynamic hydration
     AFCData_qi(126)  = 0.99  ; %Cl - %Zeund et al 2008 - dynamic hydration
     AFCData_qi(127)  = 0.5597  ; %F -
     AFCData_qi(128)  = 1.16  ; %Br - %Zeund et al 2008 - dynamic hydration
     AFCData_qi(129)  = 1.4118  ; %I -
     AFCData_qi(130)  = 1.9  ; %CH3COO -
     AFCData_qi(131)  = 1.1506  ; %SCN -
     AFCData_qi(132)  = 1.4  ; %HSO4-%Zeund et al 2008 - dynamic hydration
     