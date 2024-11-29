********************************************************************************
** CONSTANT-NVT MONTE CARLO FOR LENNARD JONES ATOMS                           **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************
C    *******************************************************************
C    ** MONTE CARLO SIMULATION PROGRAM IN THE CONSTANT-NVT ENSEMBLE.  **
C    **                                                               **
C    ** THIS PROGRAM TAKES A CONFIGURATION OF LENNARD JONES ATOMS     **
C    ** AND PERFORMS A CONVENTIONAL NVT MC SIMULATION. THE BOX IS OF  **
C    ** UNIT LENGTH, -0.5 TO +0.5 AND THERE ARE NO LOOKUP TABLES.     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF MOLECULES               **
C    ** INTEGER NSTEP               MAXIMUM NUMBER OF CYCLES          **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    DENS                REDUCED DENSITY                   **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    SIGMA               REDUCED LJ DIAMETER               **
C    ** REAL    RMIN                MINIMUM REDUCED PAIR SEPARATION   **
C    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
C    ** REAL    DRMAX               REDUCED MAXIMUM DISPLACEMENT      **
C    ** REAL    V                   THE POTENTIAL ENERGY              **
C    ** REAL    W                   THE VIRIAL                        **
C    ** REAL    PRES                THE PRESSURE                      **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE PROGRAM TAKES IN A CONFIGURATION OF ATOMS                 **
C    ** AND RUNS A MONTE CARLO SIMULATION AT THE GIVEN TEMPERATURE    **
C    ** FOR THE SPECIFIED NUMBER OF CYCLES.                           **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THE PROGRAM USES LENNARD-JONES UNITS FOR USER INPUT AND       **
C    ** OUTPUT BUT CONDUCTS THE SIMULATION IN A BOX OF UNIT LENGTH.   **
C    ** FOR EXAMPLE, FOR A BOXLENGTH L, AND LENNARD-JONES PARAMETERS  **
C    ** EPSILON AND SIGMA, THE UNITS ARE:                             **
C    **                                                               **
C    **     PROPERTY       LJ  UNITS            PROGRAM UNITS         **
C    **                                                               **
C    **     TEMP           EPSILON/K            EPSILON/K             **
C    **     PRES           EPSILON/SIGMA**3     EPSILON/L**3          **
C    **     V              EPSILON              EPSILON               **
C    **     DENS           1/SIGMA**3           1/L**3                **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE SUMUP ( RCUT, RMIN, SIGMA, OVRLAP, V, W )          **
C    **    CALCULATES THE TOTAL POTENTIAL ENERGY FOR A CONFIGURATION  **
C    ** SUBROUTINE ENERGY ( RXI, RYI, RZI, I, RCUT, SIGMA, V, W )     **
C    **    CALCULATES THE POTENTIAL ENERGY OF ATOM I WITH ALL THE     **
C    **    OTHER ATOMS IN THE LIQUID                                  **
C    ** SUBROUTINE READCN (CNFILE )                                   **
C    **    READS IN A CONFIGURATION                                   **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT A CONFIGURATION                                 **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM NUMBER BETWEEN ZERO AND ONE       **
C    *******************************************************************
        PROGRAM MCNVT
        IMPLICIT REAL*8 (A-H,O-Z)                                  
        PARAMETER ( NP = 10976 , IL=1029, IC=221591, IM=1048576)
        PARAMETER ( N0 = 1000 )
        PARAMETER ( PI = 3.1415927 )
        COMMON / BLOCK1 / RX(NP), RY(NP), RZ(NP), N
        COMMON / BLOCK3 / GRANG(N0) , DRING
        COMMON / BLOCK4 / SIGMA, SIGSQ, RCUT, RCUTSQ
        COMMON / BLOCK5 / DRMAX, DOTMIN, DENSLJ
c        COMMON / GAUSSP / EG,AG,RG,BAG,BRG,WAG
        REAL ETIME
        REAL ELAPSED(2)

        LOGICAL     OVRLAP
        CHARACTER   TITLE*80, CNFILE*70
        SAVE        ISEED
        DATA        ISEED / 0 /

C** READ INPUT DATA **
        WRITE(*,'('' =========================================== '')')
        WRITE(*,'(''      CONSTANT-NVT MONTE CARLO PROGRAM       '')')
        WRITE(*,'(''         FOR CORE-SOFTENED PARTICLES         '')')
        WRITE(*,'('' =========================================== '')')
        WRITE(*,'('' ENTER THE RUN TITLE                         '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER NUMBER OF CYCLES                      '')')
        READ (*,*) NSTEP
        WRITE(*,*) NSTEP
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN OUTPUT LINES  '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN DATA SAVES    '')')
        READ (*,*) ISAVE
        WRITE(*,'('' ENTER INTERVAL FOR UPDATE OF MAX. DISPL.    '')')
        READ (*,*) IRATIO
        WRITE(*,'('' ENTER THE CONFIGURATION FILE NAME           '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER THE FOLLOWING IN LENNARD-JONES UNITS'',/)')
        WRITE(*,'('' ENTER THE DENSITY                           '')')
        READ (*,*) DENS
        WRITE (*,'(F12.5)') DENS
        WRITE(*,'('' ENTER THE TEMPERATURE                       '')')
        READ (*,*) TEMP
        WRITE (*,'(F12.5)') TEMP


*angular correlations
        write(*,*)'IF CORRELATION FUNCTIONS  ENTER MODULO (NO=0)'
        read(*,*)MO
        IANG=0
        DRANG=0.02

C** READ INITIAL CONFIGURATION **
        CALL READCN ( CNFILE )

C** CONVERT INPUT DATA TO PROGRAM UNITS **
        BETA   = 1.0 / TEMP
        SIGMA  = ( DENS / REAL ( N ) ) ** ( 1.0 / 3.0 )
        SIGSQ  = SIGMA * SIGMA
        RMIN   = 0.70 * SIGMA
        RCUT   = 0.5
        RCUTSQ = RCUT * RCUT
        DRMAX  = 0.15 * SIGMA
        DRING = DRANG * SIGMA
        DENSLJ = DENS
        DENS   = DENS / ( SIGMA ** 3 )


C** WRITE INPUT DATA **
        WRITE(*,'(       //1X                    ,A     )') TITLE
        WRITE(*,'('' NUMBER OF ATOMS           '',I10   )') N
        WRITE(*,'('' NUMBER OF CYCLES          '',I10   )') NSTEP
        WRITE(*,'('' OUTPUT FREQUENCY          '',I10   )') IPRINT
        WRITE(*,'('' SAVE FREQUENCY            '',I10   )') ISAVE
        WRITE(*,'('' RATIO UPDATE FREQUENCY    '',I10   )') IRATIO
        WRITE(*,'('' CONFIGURATION FILE  NAME  '',A     )') CNFILE
        WRITE(*,'('' TEMPERATURE               '',F10.4 )') TEMP
        WRITE(*,'('' DENSITY                   '',F10.4 )') DENSLJ
        WRITE(*,'('' POTENTIAL CUTOFF          '',F10.4 )') RCUT

C** ZERO ACCUMULATORS **
        ACV    = 0.0
        ACVSQ  = 0.0
        ACP    = 0.0
        ACPSQ  = 0.0
        FLV    = 0.0
        FLP    = 0.0
        ACM    = 0.0
        ACATMA = 0.0

C** CALCULATE LONG RANGE CORRECTIONS    **
C** SPECIFIC TO THE LENNARD JONES FLUID **
        SR3 = ( SIGMA / RCUT ) ** 3
        SR9 = SR3 ** 3

        VLRC12 =   8.0 * PI * DENSLJ * REAL ( N ) * SR9 / 9.0
        VLRC6  = - 8.0 * PI * DENSLJ * REAL ( N ) * SR3 / 3.0
        VLRC   =   VLRC12 + VLRC6
        WLRC12 =   4.0  * VLRC12
        WLRC6  =   2.0  * VLRC6
        WLRC   =   WLRC12 + WLRC6

C** WRITE OUT SOME USEFUL INFORMATION **
        WRITE(*,'('' SIGMA/BOX              =  '',F10.4)')  SIGMA
        WRITE(*,'('' RMIN/BOX               =  '',F10.4)')  RMIN
        WRITE(*,'('' RCUT/BOX               =  '',F10.4)')  RCUT
        WRITE(*,'('' LRC FOR <V>            =  '',F10.4)')  VLRC
        WRITE(*,'('' LRC FOR <W>            =  '',F10.4)')  WLRC

C** CALCULATE INITIAL ENERGY AND CHECK FOR OVERLAPS **
        CALL SUMUP ( RMIN, OVRLAP, V, W )

        IF ( OVRLAP ) STOP ' OVERLAP IN INITIAL CONFIGURATION '

        VS = ( V + VLRC ) / REAL ( N )
        WS = ( W + WLRC ) / REAL ( N )
        PS = DENS * TEMP + W + WLRC
        PS = PS * SIGMA ** 3

        WRITE(*,'('' INITIAL V              =  '', F10.4 )' ) VS
        WRITE(*,'('' INITIAL W              =  '', F10.4 )' ) WS
        WRITE(*,'('' INITIAL P              =  '', F10.4 )' ) PS

        WRITE(*,'(//'' START OF MARKOV CHAIN               ''//)')
        WRITE(*,'(''  NMOVE     RATIO       V/N            P''/)')

C*******************************************************************
C** LOOPS OVER ALL CYCLES AND ALL MOLECULES                       **
C*******************************************************************

        DO 100 ISTEP = 1, NSTEP

           DO 99 I = 1, N

              RXIOLD = RX(I)
              RYIOLD = RY(I)
              RZIOLD = RZ(I)

C** CALCULATE THE ENERGY OF I IN THE OLD CONFIGURATION **
        CALL ENERGY ( RXIOLD, RYIOLD, RZIOLD, I, VOLD, WOLD )

C** MOVE I AND PICKUP THE CENTRAL IMAGE **
        ISEED = MOD ( ISEED * IL + IC, IM )
        RANFD = REAL ( ISEED ) / IM
              RXINEW = RXIOLD + ( 2.0 * RANFD - 1.0 ) * DRMAX
        ISEED = MOD ( ISEED * IL + IC, IM )
        RANFD = REAL ( ISEED ) / IM
              RYINEW = RYIOLD + ( 2.0 * RANFD - 1.0 ) * DRMAX
        ISEED = MOD ( ISEED * IL + IC, IM )
        RANFD = REAL ( ISEED ) / IM
              RZINEW = RZIOLD + ( 2.0 * RANFD - 1.0 ) * DRMAX

              RXINEW = RXINEW - ANINT ( RXINEW )
              RYINEW = RYINEW - ANINT ( RYINEW )
              RZINEW = RZINEW - ANINT ( RZINEW )

C** CALCULATE THE ENERGY OF I IN THE NEW CONFIGURATION **
         CALL ENERGY ( RXINEW, RYINEW, RZINEW, I,  VNEW, WNEW )

C** CHECK FOR ACCEPTANCE **
              DELTV  = VNEW - VOLD
              DELTW  = WNEW - WOLD
              DELTVB = BETA * DELTV

              IF ( DELTVB .LT. 75.0 ) THEN

        ISEED = MOD ( ISEED * IL + IC, IM )
        RANFD = REAL ( ISEED ) / IM

                 IF ( DELTV .LE. 0.0 ) THEN

                    V      = V + DELTV
                    W      = W + DELTW
                    RX(I)  = RXINEW
                    RY(I)  = RYINEW
                    RZ(I)  = RZINEW
                    ACATMA = ACATMA + 1.0

                 ELSEIF ( EXP ( - DELTVB ) .GT. RANFD ) THEN

                    V      = V + DELTV
                    W      = W + DELTW
                    RX(I)  = RXINEW
                    RY(I)  = RYINEW
                    RZ(I)  = RZINEW
                    ACATMA = ACATMA + 1.0

                 ENDIF

              ENDIF

              ACM = ACM + 1.0

C** CALCULATE INSTANTANEOUS VALUES **
              VN     = ( V + VLRC ) / REAL ( N )
              PRES   = DENS * TEMP + W + WLRC

C** CONVERT PRESSURE TO LJ UNITS **
              PRES   = PRES * SIGMA ** 3

C** ACCUMULATE AVERAGES **
              ACV    = ACV   + VN
              ACP    = ACP   + PRES
              ACVSQ  = ACVSQ + VN * VN
              ACPSQ  = ACPSQ + PRES * PRES

C*************************************************************
C* ENDS LOOP OVER ATOMS                                    **
C*************************************************************

99         CONTINUE

C** PERFORM PERIODIC OPERATIONS  **
           IF ( MOD ( ISTEP, IRATIO ) .EQ. 0 ) THEN

C** ADJUST MAXIMUM DISPLACEMENT **
              RATIO = ACATMA / REAL ( N * IRATIO )
              IF ( RATIO .GT. 0.5 ) THEN
                 DRMAX  = DRMAX  * 1.05
              ELSE
                 DRMAX  = DRMAX  * 0.95
              ENDIF
              ACATMA = 0.0
           ENDIF

C** WRITE OUT RUNTIME INFORMATION **
           IF ( MOD ( ISTEP, IPRINT ) .EQ. 0 ) THEN
              WRITE(*,'(I8,3F12.6)') INT(ACM), RATIO, VN, PRES
           ENDIF

C** WRITE OUT THE CONFIGURATION AT INTERVALS **
           IF ( MOD ( ISTEP, ISAVE ) .EQ. 0 ) THEN
              CALL WRITCN ( CNFILE )
           ENDIF

C** ANGULAR FUNCTIONS
           IF(MO.NE.0.AND.MOD(ISTEP,MO).EQ.0) THEN
              IANG=IANG+1
              CALL ANGCOR
            ENDIF

100     CONTINUE

C*******************************************************************
C** ENDS THE LOOP OVER CYCLES                                     **
C*******************************************************************

        WRITE(*,'(//'' END OF MARKOV CHAIN          ''//)')

C** CHECKS FINAL VALUE OF THE POTENTIAL ENERGY IS CONSISTENT **
        CALL SUMUP ( RMIN, OVRLAP, VEND, WEND )

        IF ( ABS ( VEND - V ) .GT. 1.0E-03 ) THEN
           WRITE(*,'('' PROBLEM WITH ENERGY,'')')
           WRITE(*,'('' VEND              = '', E20.6)') VEND
           WRITE(*,'('' V                 = '', E20.6)') V
        ENDIF

C** WRITE OUT THE FINAL CONFIGURATION FROM THE RUN **
        CALL WRITCN ( CNFILE )

C** CALCULATE AND WRITE OUT RUNNING AVERAGES **
        AVV   = ACV / ACM
        ACVSQ = ( ACVSQ / ACM ) - AVV ** 2
        AVP   = ACP / ACM
        ACPSQ = ( ACPSQ / ACM ) - AVP ** 2

C** CALCULATE FLUCTUATIONS **
        IF ( ACVSQ .GT. 0.0 ) FLV = SQRT ( ACVSQ )
        IF ( ACPSQ .GT. 0.0 ) FLP = SQRT ( ACPSQ )

        WRITE(*,'(/'' AVERAGES ''/ )')
        WRITE(*,'('' <V/N>   = '',2F10.6)') AVV,AVV/TEMP
        WRITE(*,'('' <P>     = '',2F10.6)') AVP,AVP/TEMP/DENSLJ

        WRITE(*,'(/'' FLUCTUATIONS ''/)')

        WRITE(*,'('' FLUCTUATION IN <V/N> = '',2F10.6)') FLV
     $, FLV/TEMP
        WRITE(*,'('' FLUCTUATION IN <P>   = '',2F10.6)') FLP
     $, FLP/TEMP/DENSLJ
        WRITE(*,'(/'' END OF SIMULATION '')')

*===================================================================
*OUTPUT ANGULAR CORRELATIONS
        IF(MO.NE.0) THEN
          
***** AVERAGE FOR ANGULAR PROJECTIONS
      COCO=0.5/DRING
c     ZIZ=8.0/(4*PI /COCO**3  )/(IANG*N*(N-1.)/2.0)
*below replace the cursed N(N-1) by N**2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ZIZ=8.0/(4*PI /COCO**3  )/(FLOAT(IANG)/2.0)
      ZIZA=ZIZ/N**2

      NZAZ=0.70/DRANG-1
      NZIZ=RCUT/DRING
      IF(NZIZ.GT.N0) NZIZ=N0
      write(*,*)'   '
c     write(*,*)ziz,nzaz,nziz,iang,coco,dring,drang
*COMPUTE RADIAL G000(R)
      DO  LL=1,NZIZ
      GRANG(LL)=GRANG(LL)*ZIZA/LL/LL
      ENDDO
1133  FORMAT(3X,F8.4, F12.7,I8)
c      DO 1122 LL=1,NZIZ
c1122  WRITE(*,1133)DRANG*LL,GRANG(LL),LL 
*save to file
      open(1,file='gr.data')
      write(1,'("# gr for LJ  N=",i6," rho*=",f8.4)')N,dens
      do ir=1,nziz
       write(1,'(f12.7,f15.9)')drang*ir,grang(ir)
      enddo
      close(1)
      write(*,*)'G(R) saved to file "gr.data" '

      ENDIF

      ATIME=ETIME(ELAPSED)
      WRITE(*,*)'Elaspsed time user-system ',elapsed
      WRITE(*,'(/'' END OF SIMULATION '')')


      STOP
      END
C*******************************************************************
C** CALCULATES THE TOTAL POTENTIAL ENERGY FOR A CONFIGURATION.    **
C**                                                               **
C*******************************************************************
        SUBROUTINE SUMUP ( RMIN, OVRLAP, V, W )
        IMPLICIT REAL*8 (A-H,O-Z)                                  
        PARAMETER ( NP = 10976 )
        COMMON / BLOCK1 / RX(NP), RY(NP), RZ(NP), N
        COMMON / BLOCK4 / SIGMA, SIGSQ, RCUT, RCUTSQ
c        COMMON / GAUSSP / EG,AG,RG,BAG,BRG,WAG
        LOGICAL     OVRLAP

        OVRLAP = .FALSE.
        RCUTSQ = RCUT * RCUT
        RMINSQ = RMIN * RMIN
        SIGSQ  = SIGMA * SIGMA

        V      = 0.0
        W      = 0.0
c        VG     = 0.0
c        WG     = 0.0


C** LOOP OVER ALL THE PAIRS IN THE LIQUID **
        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)

           DO 99 J = I + 1, N

              RXIJ  = RXI - RX(J)
              RYIJ  = RYI - RY(J)
              RZIJ  = RZI - RZ(J)

              RXIJ  = RXIJ - ANINT ( RXIJ )
              RYIJ  = RYIJ - ANINT ( RYIJ )
              RZIJ  = RZIJ - ANINT ( RZIJ )
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

              IF ( RIJSQ .LT. RMINSQ ) THEN
                 OVRLAP = .TRUE.
                 RETURN

              ELSEIF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2 = SIGSQ / RIJSQ
                 SR6 = SR2 * SR2 * SR2
                 VIJ = SR6 * ( SR6 - 1.0 )
                 WIJ = SR6 * ( SR6 - 0.5 )
                 V   = V + VIJ
                 W   = W + WIJ

              ENDIF

99         CONTINUE

100     CONTINUE
 
        V = 4.0 * V 
        W = 48.0 * W / 3.0  

        RETURN
        END
C*******************************************************************
C** RETURNS THE POTENTIAL ENERGY OF ATOM I WITH ALL OTHER ATOMS.  **
C*******************************************************************
        SUBROUTINE ENERGY ( RXI, RYI, RZI, I, V, W )
        IMPLICIT REAL*8 (A-H,O-Z)                                  
        PARAMETER ( NP = 10976 )
        COMMON / BLOCK1 / RX(NP), RY(NP), RZ(NP), N
        COMMON / BLOCK4 / SIGMA, SIGSQ, RCUT, RCUTSQ
c        COMMON / GAUSSP / EG,AG,RG,BAG,BRG,WAG

        RCUTSQ = RCUT * RCUT
        SIGSQ  = SIGMA * SIGMA

        V      = 0.0
        W      = 0.0
c        VG     = 0.0
c        WG     = 0.0

C** LOOP OVER ALL MOLECULES EXCEPT I  **
        DO 100 J = 1, N

           IF ( I .NE. J ) THEN

              RXIJ  = RXI - RX(J)
              RYIJ  = RYI - RY(J)
              RZIJ  = RZI - RZ(J)

              RXIJ  = RXIJ - ANINT ( RXIJ )
              RYIJ  = RYIJ - ANINT ( RYIJ )
              RZIJ  = RZIJ - ANINT ( RZIJ )

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

              IF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2 = SIGSQ / RIJSQ
                 SR6 = SR2 * SR2 * SR2
                 VIJ = SR6 * ( SR6 - 1.0 )
                 WIJ = SR6 * ( SR6 - 0.5 )
                 V   = V + VIJ
                 W   = W + WIJ

              ENDIF

           ENDIF

100     CONTINUE
     
        V = 4.0 * V
        W = 48.0 * W / 3.0

        RETURN
        END
C*******************************************************************
C** SUBROUTINE TO READ IN THE CONFIGURATION FROM UNIT 1           **
C*******************************************************************
        SUBROUTINE READCN ( CNFILE )
        IMPLICIT REAL*8 (A-H,O-Z)                                  
        PARAMETER ( NP = 10976 )
        COMMON / BLOCK1 / RX(NP), RY(NP), RZ(NP), N
        COMMON / BLOCK5 / DRMAX, DOTMIN, DENS
        CHARACTER*70   CNFILE

        OPEN (1, FILE=CNFILE)
        READ (1,*)NN, DRMAX
        DO I=1,NN
         READ ( 1,*) RX(I), RY(I), RZ(I)
        ENDDO
        N = NN
        CLOSE ( 1 )
        RETURN
        END
C*******************************************************************
C** SUBROUTINE TO WRITE OUT THE CONFIGURATION TO UNIT 2           **
C*******************************************************************
        SUBROUTINE WRITCN ( CNFILE )
        IMPLICIT REAL*8 (A-H,O-Z)                                  
        PARAMETER ( NP = 10976 )
        COMMON / BLOCK1 / RX(NP), RY(NP), RZ(NP), N
        COMMON / BLOCK5 / DRMAX, DOTMIN, DENS
        CHARACTER*70   CNFILE

        OPEN (2, FILE=CNFILE)
        WRITE(2, '(I6,3F12.5)' ) N, DRMAX, DENS
        DO  I=1,N
          WRITE ( 2,'(3F12.7)') RX(I), RY(I), RZ(I)
        ENDDO
        CLOSE ( 2 )
        RETURN
        END
********************************************************************
* ROUTINE TO COMPUTE G(R) 
********************************************************************
       SUBROUTINE ANGCOR              
       IMPLICIT REAL*8 (A-H,O-Z)                                  
       PARAMETER ( NP = 10976 , N0=1000)    

        COMMON / BLOCK1 / RX(NP), RY(NP), RZ(NP), N
        COMMON / BLOCK3 / GRANG(N0) , DR
        COMMON / BLOCK4 / SIGMA, SIGSQ, RCUT, RCUTSQ
        COMMON / LOCAL0 / LOLO
      DIMENSION GR(N0)         
      DIMENSION NPANG(N0) 
                       
      IF(LOLO.NE.-189) THEN                                  
       DO 1 J=1,N0                                           
 1     GRANG(J)=0.0                                         
          PIVAL=4.*ATAN(1.0)                    
          DK=PIVAL/1.                                  
          DDK=DK/SQRT(3.0)
          LOLO=-189                            
      ENDIF                                                
            
      NSP = N                                                       
      NSEP=INT(0.5/DR+0.5)                                     
      IF(NSEP.GT.N0) NSEP=N0                                 
      NSPM=NSP-1                                     


      DO 5 I=1,N0                                    
       NPANG(I)=0                                                 
 5    CONTINUE                                                
      DO 6 L=1,N0                                       
       GR(L)=0                                        
 6    CONTINUE  
*LOOP ON PARTICLE PAIRS                                       
      DO 100 I=1,NSPM                       
      IP=I+1                                                     
      IM=I-1                                         
                            
      DO 110 J=IP,NSP                                              
      JM=J-1                                                    

        XDIF=RX(I)-RX(J) 
        YDIF=RY(I)-RY(J)  
        ZDIF=RZ(I)-RZ(J)                            
                          
        X12=XDIF-ANINT(XDIF) 
        Y12=YDIF-ANINT(YDIF)                     
        Z12=ZDIF-ANINT(ZDIF)                  
                                                             
        R212=X12*X12+Y12*Y12+Z12*Z12  
                                  
        IF(R212.LE.RCUTSQ) THEN                                      
*WE ACCUMULATE ONLY WITHIN CUTOFF RADIUS                    

          RR=SQRT(R212)                              
          L=INT(RR/DR+0.5)  
                                          
          NPANG(L)=NPANG(L)+1                             

C(000)                         
        GR(L)=FLOAT(NPANG(L))                               
       ENDIF 
                                                                
 110   CONTINUE 
 100   CONTINUE 
*COMPUTE G000(R) 
       DO 60 L=1,NSEP 
 60    GRANG(L)=GRANG(L) +GR(L) 
    
       RETURN  
       END    
