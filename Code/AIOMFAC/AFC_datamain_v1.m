    function [AFCData_Main,charge] = AFC_datamain_v1()
    %--- Non-eletrolyte functional groups ---
     AFCData_Main(1)  = 1 ; %CH3
     AFCData_Main(2)  = 1 ; %CH2
     AFCData_Main(3)  = 1 ; %CH
     AFCData_Main(4) = 1 ; %C
     AFCData_Main(5)  = 2 ; %CH2=CH
     AFCData_Main(6)  = 2 ; %CH=CH
     AFCData_Main(7)  = 2 ; %CH2=C
     AFCData_Main(8)  = 2 ; %CH=C
     AFCData_Main(9)  = 2 ; %C=C
     AFCData_Main(10)  = 3 ; %ACH
     AFCData_Main(11)  = 3 ; %AC
     AFCData_Main(12)  = 4 ; %ACCH3
     AFCData_Main(13)  = 4 ; %ACCH2
     AFCData_Main(14)  = 4 ; %ACCH
     AFCData_Main(15)  = 5 ; %OH
     AFCData_Main(16)  = 6 ; %CH3OH
     AFCData_Main(17)  = 7 ; %H2O
     AFCData_Main(18)  = 8 ; %ACOH
     AFCData_Main(19)  = 9 ; %CH3CO
     AFCData_Main(20)  = 9 ; %CH2CO
     AFCData_Main(21)  = 10 ; %CHO
     AFCData_Main(22)  = 11 ; %CH3COO
     AFCData_Main(23)  = 11 ; %CH2COO
     AFCData_Main(24)  = 12 ; %HCOO
     AFCData_Main(25)  = 13 ; %CH3O
     AFCData_Main(26)  = 13 ; %CH2O
     AFCData_Main(27)  = 13 ; %CHO
     AFCData_Main(28)  = 13 ; %THF
     AFCData_Main(29)  = 14 ; %CH3NH2
     AFCData_Main(30)  = 14 ; %CH2NH2
     AFCData_Main(31)  = 14 ; %CHNH2
     AFCData_Main(32)  = 15 ; %CH3NH
     AFCData_Main(33)  = 15 ; %CH2NH
     AFCData_Main(34)  = 15 ; %CHNH
     AFCData_Main(35)  = 16 ; %CH3N
     AFCData_Main(36)  = 16 ; %CH2N
     AFCData_Main(37)  = 17 ; %ACNH2
     AFCData_Main(38)  = 18 ; %C5H5N
     AFCData_Main(39)  = 18 ; %C5H4N
     AFCData_Main(40)  = 18 ; %C5H3N
     AFCData_Main(41)  = 19 ; %CH3CN
     AFCData_Main(42)  = 19 ; %CH2CN
     AFCData_Main(43)  = 20 ; %COOH
     AFCData_Main(44)  = 20 ; %HCOOH
     AFCData_Main(45)  = 21 ; %CH2CL
     AFCData_Main(46)  = 21 ; %CHCL
     AFCData_Main(47)  = 21 ; %CCL
     AFCData_Main(48)  = 22 ; %CH2CL2
     AFCData_Main(49)  = 22 ; %CHCL2
     AFCData_Main(50)  = 22 ; %CCL2
     AFCData_Main(51)  = 23 ; %CHCL3
     AFCData_Main(52)  = 23 ; %CCL3
     AFCData_Main(53)  = 24 ; %CCL4
     AFCData_Main(54)  = 25 ; %ACCL
     AFCData_Main(55)  = 26 ; %CH3NO2
     AFCData_Main(56)  = 26 ; %CH2NO2
     AFCData_Main(57)  = 26 ; %CHNO2
     AFCData_Main(58)  = 27 ; %ACNO2
     AFCData_Main(59)  = 28 ; %CS2
     AFCData_Main(60)  = 29 ; %CH3SH
     AFCData_Main(61)  = 29 ; %CH2SH
     AFCData_Main(62)  = 30 ; %FURFURAL
     AFCData_Main(63)  = 31 ; %DOH
     AFCData_Main(64)  = 32 ; %I
     AFCData_Main(65)  = 33 ; %Br
     AFCData_Main(66)  = 34 ; %CH=-C
     AFCData_Main(67)  = 34 ; %C=-C
     AFCData_Main(68)  = 35 ; %DMSO
     AFCData_Main(69)  = 36 ; %Acrylonitrate
     AFCData_Main(70)  = 37 ; %Cl-(C=C)
     AFCData_Main(71)  = 38 ; %ACF
     AFCData_Main(72)  = 39 ; %DMF
     AFCData_Main(73)  = 39 ; %HCON(CH2)2
     AFCData_Main(74)  = 40 ; %CF3
     AFCData_Main(75)  = 40 ; %CF2
     AFCData_Main(76)  = 40 ; %CF
     AFCData_Main(77)  = 41 ; %COO
     AFCData_Main(78)  = 42 ; %SiH3
     AFCData_Main(79)  = 42 ; %SiH2
     AFCData_Main(80)  = 42 ; %SiH
     AFCData_Main(81)  = 42 ; %Si
     AFCData_Main(82)  = 43 ; %SiH2O
     AFCData_Main(83)  = 43 ; %SiHO
     AFCData_Main(84)  = 43 ; %SiO
     AFCData_Main(85)  = 9  ; %CO
     %extension to new functional groups not covered in previous unifac
     %versions
     AFCData_Main(86)  = 44  ; %NMP
     AFCData_Main(87)  = 45  ; %CCl3F
     AFCData_Main(88)  = 45  ; %CCl2F
     AFCData_Main(89)  = 45  ; %HCCl2F
     AFCData_Main(90)  = 45  ; %HCClF
     AFCData_Main(91)  = 45  ; %CClF2
     AFCData_Main(92)  = 45  ; %HCClF2
     AFCData_Main(93)  = 45  ; %CClF3
     AFCData_Main(94)  = 45  ; %CCl2F2
     AFCData_Main(95)  = 46  ; %CONH2
     AFCData_Main(96)  = 46  ; %CONHCH3
     AFCData_Main(97)  = 46  ; %CONHCH2
     AFCData_Main(98)  = 46  ; %CON(CH3)2
     AFCData_Main(99)  = 46  ; %CONCH3CH2
     AFCData_Main(100)  = 46  ; %CON(CH2)2
     AFCData_Main(101)  = 47  ; %C2H5O2
     AFCData_Main(102)  = 47  ; %C2H4O2
     AFCData_Main(103)  = 48  ; %CH3S
     AFCData_Main(104)  = 48  ; %CH2S
     AFCData_Main(105)  = 48  ; %CHS
     AFCData_Main(106)  = 49  ; %MORPH
     AFCData_Main(107)  = 50  ; %C4H4S
     AFCData_Main(108)  = 50  ; %C4H3S
     AFCData_Main(109)  = 50  ; %C4H2S     
     %--- Electrolytes ----
     AFCData_Main(110)  = 51  ; %H +
     AFCData_Main(111)  = 52  ; %NH4 +
     AFCData_Main(112)  = 53  ; %Na +
     AFCData_Main(113)  = 54  ; %K +
     AFCData_Main(114)  = 55  ; %Mg 2+
     AFCData_Main(115)  = 56  ; %Ca 2+
     AFCData_Main(116)  = 57  ; %Li +
     AFCData_Main(117)  = 58  ; %Ba 2+
     AFCData_Main(118)  = 59  ; %Sr 2+
     AFCData_Main(119)  = 60  ; %Co 2+
     AFCData_Main(120)  = 61  ; %Ni 2+
     AFCData_Main(121)  = 62  ; %Cu 2+
     AFCData_Main(122)  = 63  ; %Zn 2+
     AFCData_Main(123)  = 64  ; %Hg 2+
     AFCData_Main(124)  = 65  ; %SO4 2-
     AFCData_Main(125)  = 66  ; %NO3 -
     AFCData_Main(126)  = 67  ; %Cl -
     AFCData_Main(127)  = 68  ; %F -
     AFCData_Main(128)  = 69  ; %Br -
     AFCData_Main(129)  = 70  ; %I -
     AFCData_Main(130)  = 71  ; %CH3COO -
     AFCData_Main(131)  = 72  ; %SCN -
     AFCData_Main(132)  = 73  ; %HSO4 -
     
     charge=identify_ions();
     %sub function to retrieve ionic identifiers for ion-ion interactions
        function [charge] = identify_ions()
            %0-non ionic component
            %1-cation
            %2-annion
            charge=sparse(1,132);
            charge(110)  = 1  ; %H +
            charge(111)  = 1  ; %NH4 +
            charge(112)  = 1  ; %Na +
            charge(113)  = 1  ; %K +
            charge(114)  = 2  ; %Mg 2+
            charge(115)  = 2  ; %Ca 2+
            charge(116)  = 1  ; %Li +
            charge(117)  = 2  ; %Ba 2+
            charge(118)  = 2  ; %Sr 2+
            charge(119)  = 2  ; %Co 2+
            charge(120)  = 2  ; %Ni 2+
            charge(121)  = 2  ; %Cu 2+
            charge(122)  = 2  ; %Zn 2+
            charge(123)  = 2 ; %Hg 2+
            charge(124)  = -2  ; %SO4 2-
            charge(125)  = -1  ; %NO3 -
            charge(126)  = -1  ; %Cl -
            charge(127)  = -1  ; %F -
            charge(128)  = -1  ; %Br -
            charge(129)  = -1  ; %I -
            charge(130)  = -1  ; %CH3COO -
            charge(131)  = -1  ; %SCN -
            charge(132)  = -1  ; %HSO4 -
     
     
     
     
     