      REAL FUNCTION O2ABS2(TEMP,PRES,VAPDEN,FREQ,FAC,FAC2,O2C)
!C  Copyright (c) 2009 Massachusetts Institute of Technology
!C
!C     RETURNS POWER ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
!C     IN NEPERS/KM.  MULTIPLY O2ABS BY 4.343 TO CONVERT TO DB/KM.
!C
!C      5/1/95  P. Rosenkranz 
!C      11/5/97  P. Rosenkranz - 1- line modification.
!c      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
!c      8/21/02  pwr - revised width at 425
!c      3/20/03  pwr - 1- line mixing and width revised
!c      9/29/04  pwr - new widths and mixing, using HITRAN intensities
!c                     for all lines
!c      6/12/06  pwr - chg. T dependence of 1- line to 0.8
!C      10/14/08  pwr - moved isotope abundance back into intensities, 
!c                     added selected O16O18 lines.
!c      6/5/09  pwr - remove common block, add weak lines and second-order mixing
!C
      IMPLICIT NONE
!C
!C     ARGUMENTS:
      REAL TEMP,PRES,VAPDEN,FREQ,FAC,FAC2,O2C
!C
!C     NAME    UNITS    DESCRIPTION        VALID RANGE
!C
!C     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
!c                                          valid for atmosphere
!C     PRES   MILLIBARS PRESSURE           3 TO 1000
!C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
!C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
!C     FREQ    GHZ      FREQUENCY          0 TO 900
!C     FAC2    none     FACTOR MULTIPLYING SECOND-ORDER COEFFICIENTS
!C	   O2C		none    concentration of O2 (0.2095 for standard air, 1 for pure O2)
!C
!C     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
!C     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
!C      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
!c     G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, 
!c      pp.282-287 (2003).
!C     M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
!C     M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
!c     B.J. Drouin, JQSRT v.105, pp.450-458 (2007).
!C     Line intensities from HITRAN2004.
!c     Non-resonant intensity from JPL catalog.
!C     Second-order coefficients from 
!c      E.W. Smith, J.Chem.Phys. v.74, pp.6658-6673 (1981).
!c
!c     note:
!c     1. The mm line-width and mixing coefficients are from Tretyakov et al;
!c        submm line-widths from Golubiatnikov & Krupnov (except 
!c        234 GHz from Drouin)
!c     2. The same temperature dependence (X) is used for submillimeter 
!c        line widths as in the 60 GHz band: (1/T)**X
!c     3. The sign of DNU in the shape factor is corrected.
!C
!c     Local variables:
      INTEGER K,NL
      PARAMETER (NL=49)
      REAL TH,TH1,B,PRESWV,PRESDA,DEN,DFNR,SUM,STR,Y,SF1,SF2
      REAL DEL1,DEL2,GFAC,PE2,D1,D2,DF
      REAL X,WB300,W300(NL),F(NL),Y300(NL),S300(NL),V(NL),BE(NL), &
      G(NL),DNU(NL),GT(NL),DNUT(NL),YO300(NL),GO(NL),DNUO(NL),WN(NL),WO(NL)	  
!C      LINES ARE ARRANGED 1-,1+,...33-,33+ IN SPIN-ROTATION SPECTRUM;
!c      BY FREQUENCY IN SUBMM SPECTRUM.
      DATA F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,	&
       59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,&
       56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685, &
       55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241, &
       53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, &
       52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310,  &
       50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630, &
       487.2493, 566.8956, 715.3929, 731.1866,	&
       773.8395, 834.1455, 895.0710/
      DATA S300/ &
      0.2906E-14,0.7957E-15,0.2444E-14,0.2194E-14, &
      0.3301E-14,0.3243E-14,0.3664E-14,0.3834E-14, &
      0.3588E-14,0.3947E-14,0.3179E-14,0.3661E-14, &
      0.2590E-14,0.3111E-14,0.1954E-14,0.2443E-14, &
      0.1373E-14,0.1784E-14,0.9013E-15,0.1217E-14, &
      0.5545E-15,0.7766E-15,0.3201E-15,0.4651E-15, &
      0.1738E-15,0.2619E-15,0.8880E-16,0.1387E-15, &
      0.4272E-16,0.6923E-16,0.1939E-16,0.3255E-16, &
      0.8301E-17,0.1445E-16,0.3356E-17,0.6049E-17, &
      0.1280E-17,0.2394E-17,					   &
      0.3287E-16,0.6463E-15,0.1334E-16,0.7049E-14, &
      0.3011E-14,0.1797E-16,0.1826E-14,0.2193E-16, &
      0.1153E-13,0.3974E-14,0.2512E-16/
      DATA BE /&
	  .010, .014, 2*.083,& 
	  2*.207, 2*.387, &
	  2*.621, 2*.910, &
	  2*1.255, 2*1.654, &
	  2*2.109, 2*2.618, &
      2*3.182, 2*3.800, &
	  2*4.474, 2*5.201, &
	  2*5.983, 2*6.819, &
      2*7.709, 2*8.653, &
	  2*9.651,	&
      .019, .048, .045, .044, .049, .084, .145, .136, .141, .145, .201/
!C      WIDTHS IN MHZ/MB
      DATA WB300/.56/, X/.754/
      DATA W300/ 1.664, 1.703, 1.513, 1.491, 1.415, 1.408, &
      1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217,&
      1.189, 1.174, 1.134, 1.134, 1.089, 1.088, 1.037,1.038,&
      2*0.996, 2*0.955, 2*0.906, 2*0.858, 2*0.811, 2*0.764, &
      2*0.717, 2*0.669,&
      1.65, 3*1.64, 4*1.60, 1.62, 2*1.47/
	  DATA WN  / 1.684, 1.733, 1.524, 1.506, 1.428, 1.425, &
	  1.365, 1.352, 1.305, 1.301, 1.267, 1.268, 1.228, 1.218, &
	  1.188, 1.169, 1.127, 1.128, 1.089, 1.066, 1.065, 1.065, &
	  1.037, 1.013, 0.986, 0.983, 0.930, 0.930, &
	  0.828, 0.828, 0.765, 0.765, 0.720, 0.720, &
	  0.668, 0.668, 0.623, 0.623, &
	  11*0./
	  DATA WO  / 1.673, 1.679, 1.55, 1.508, 1.44, 1.418, &
	  1.376, 1.362, 1.322, 1.325, 1.307, 1.31 , 1.27 , 1.276,& 
      1.255, 1.251, 1.218, 1.214, 1.19 , 1.183, 1.156, 1.157,&
      1.14 , 1.113, 1.098, 1.095, 1.05, 1.049, &
      0.947, 0.947, 0.885, 0.885, 0.84, 0.84 , &
      0.795, 0.795, 0.75 , 0.75, &
      11*0./

      DATA Y300/ -0.0360, 0.2547, -0.3655,  0.5495,&
      -0.5696,  0.6181, -0.4252,  0.3517, -0.1496,  0.0430,&
       0.0640, -0.1605,  0.2906, -0.3730,  0.4169, -0.4819,&
       0.4963, -0.5481,  0.5512, -0.5931,  0.6212, -0.6558,&
       0.6920, -0.7208,  0.7312, -0.7550,  0.7555, -0.7751,&
       0.7914, -0.8073,  0.8307, -0.8431,  0.8676, -0.8761,&
       0.9046, -0.9092,  0.9416, -0.9423,  11*0./
	  DATA YO300/-0.035, 0.0915, -0.3273,   0.5059,&
      -0.4691,  0.5054, -0.3215,  0.2548, -0.0977,  0.006621,&
       0.1166, -0.2026,  0.2509, -0.3224,  0.3305, -0.3882, &
       0.3490, -0.3962,  0.3877, -0.4276,  0.4393, -0.4740, &
       0.5575, -0.5871,  0.6025, -0.6276,  0.6183, -0.6402, &
       0.6252, -0.6443,  0.6703, -0.6861,  0.7313, -0.7433, &
       0.7803, -0.7893,  0.8161, -0.8218, 11*0./

      DATA V/  0.0079, -0.0978,  0.0844, -0.1273, &
       0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,&
       0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,&
       0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,&
       0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,&
       0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545,&
       0.680,  -0.660,   0.685,  -0.665,   11*0./
      DATA G/ &
      -0.0001, -0.0501, -0.0295, -0.0765, -0.0396, -0.0281, -0.0230, &
       0.0222,  0.0421,  0.0451,  0.0318,  0.0320,  0.0206,  0.0065, &
       0.0017, -0.0046, -0.0111, -0.0182, -0.0230, -0.0244, -0.0308, &
      -0.0352, -0.0338, -0.0347, -0.0517, -0.0509, -0.0438, -0.0428, &
      -0.0574, -0.0555, -0.0491, -0.0469, -0.0613, -0.0585, -0.0543, &
      -0.0515, -0.0479, -0.0452, 11*0./
      DATA DNU/  &
      -0.0004,  0.0118, -0.0152,  0.0267, -0.0217,  0.0247, -0.0183, &
       0.0172, -0.0122,  0.0105, -0.0064,  0.0046, -0.0029,  0.0021, &
      -0.0014,  0.0009, -0.0005,  0.0002,  0.0000, -0.0003,  0.0003, &
      -0.0005,  0.0005, -0.0006,  0.0006, -0.0007,  0.0007, -0.0007, &
       0.0007, -0.0007,  0.0007, -0.0007,  0.0007, -0.0007,  0.0006, &
      -0.0006,  0.0007, -0.0007, 11*0./
	  DATA GT /49*0./
	  DATA DNUT /49*0./
	  DATA GO/49*0./
	  DATA DNUO/49*0./
!C
1     FORMAT(2(e20.9,'  '))
2	  FORMAT(2(f20.9))


	  OPEN(6,FILE='y_1new.txt',status='OLD')
	  OPEN(7,FILE='g_1new.txt',status='OLD')
	  OPEN(8,FILE='dn_1new.txt',status='OLD')
	  OPEN(5,FILE='log.txt',status='unknown')
	  OPEN(13, FILE='ygdno.txt',status='OLD') !1st and 2nd order parameters for pure O2 only at room T, might be not accurate

	  DO 16 K=1,38 
		READ(6,*)Y300(K),V(K)
		READ(7,*)G(K),GT(K)
		READ(8,*)DNU(K),DNUT(K)
		READ(13,*)YO300(K),GO(K),DNUO(K)
16	  CONTINUE

	  CLOSE(6)
	  CLOSE(7)
	  CLOSE(8)
	  CLOSE(13)
	  TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      PRESWV = VAPDEN*TEMP/217.
      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DFNR = WB300*DEN
      PE2 = DEN*DEN*FAC2
!c  1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
      SUM = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      DO 32 K=1,NL
      DF = DEN*(O2C*WO(K)+(1-O2C)*WN(K))*0.98931903 
!0.98931903 is (296/300)**0.8 as widths are given at 296K, not 300

  	  IF (O2C.GE.0.95) THEN !for pure or almost pure O2 instead of air
		DNU(K)=DNUO(K)
		G(K)=GO(K)
		DNUT(K)=0.
		GT(K)=0.
		Y300(K)=O2C*YO300(K)+(1-O2C)*(Y300(K)-0.2096*YO300(K))/(1-0.2096)
		V(K)=0.
		ENDIF
      Y = DEN*(Y300(K)+V(K)*TH1) 
      STR = S300(K)*EXP(-BE(K)*TH1)
      DEL1 = FREQ -F(K) -PE2*(DNU(K)+DNUT(K)*TH1) 
      DEL2 = FREQ +F(K) +PE2*(DNU(K)+DNUT(K)*TH1) 
      GFAC = 1. + PE2*(G(K)+GT(K)*TH1) 
      D1 = DEL1*DEL1 + DF*DF
      D2 = DEL2*DEL2 + DF*DF
      SF1 = (DF*GFAC + DEL1*Y*FAC)/D1
      SF2 = (DF*GFAC - DEL2*Y*FAC)/D2
32    SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))**2
!c   .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
      O2ABS2 = 1.6097E11*SUM*PRESDA*TH**3
      O2ABS2 = AMAX1(O2ABS2,0.)
	  CLOSE(5)
      RETURN
      END function O2ABS2