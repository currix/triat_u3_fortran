      DOUBLE PRECISION FUNCTION CHISQRE(P11,P21,P22,P23,P31,P32,P33,P41,
     *     P42,P43,P44,P45,P46,P47)
C
C  $Id: chisqre_u3_min_gen.f,v 1.8 2011/05/04 14:36:00 curro Exp $
C
      IMPLICIT NONE
C
C     DEFINITION OF VARIABLES
C
C     
C     HAMILTONIAN PARAMETERS
      DOUBLE PRECISION P11,P21,P22,P23,P31,P32,P33,P41,
     *     P42,P43,P44,P45,P46,P47
C
C     LOCAL VARIABLES
C
      INTEGER N2MAX, DMAX, DATMAX, NPMAX, DIMSCR
C     
C     MAXIMUM NUMBER OF U(3) BOSONS (N_2)  
      PARAMETER (N2MAX = 4000) 
C     DIMENSION OF SCRATCH VECTOR
      PARAMETER (DIMSCR = 12000) 
C     MAXIMUM DIMENSIONS OF MATRICES
      PARAMETER (DMAX = 2001)
C     MAXIMUM NUMBER OF EXPERIMENTAL DATA TO CONSIDER
      PARAMETER (DATMAX = 200)
C     MAXIMUM NUMBER OF PARAMETERS IN THE HAMILTONIAN
      PARAMETER (NPMAX = 15)
C     
C     IF .T. BENT ELSE LINEAR CASE
      LOGICAL BENT
C
C     IF .T. THEN ENERGIES REFERRED THE MINIMUM OF EACH L BLOCK, 
C     ELSE REFERRED TO THE GROUND STATE (MINIMUM L = 0)
      LOGICAL EMINL
C
C     NAME OF EXPERIMENTAL DATA FILE
      CHARACTER*30 DTFL
C
C     REPRESENTATION OF U(3)
      INTEGER N2
C
C     MAXIMUM VIBRATIONAL ANGULAR MOMENTUM AND EXPERIMENTAL NUMBER OF QUANTA 
C     TO BE CONSIDERED
      INTEGER LMAX, VMAX
C
C     EXPERIMENTAL ENERGIES AND ERRORS MATRICES
      DOUBLE PRECISION EXPEN(DATMAX), EXPERR(DATMAX)
C     
C     EXPERIMENTAL ASSIGNMENTS AND POINTER
C     EXPAS(n,1) v_n EXPAS(n,2) l_n 
      INTEGER EXPAS(DATMAX,2), LXP(0:DATMAX)
C
C     HAMILTONIAN MATRIX
      DOUBLE PRECISION HAM(DMAX,DMAX)
C
C     POINTER TO COMPARE WITH EXPERIMENTAL RESULTS
      INTEGER SPNTEX(DATMAX)
C
C     STATISTICAL VARIABLES
      DOUBLE PRECISION CHSQP, CHSQT
C     COMPUTED ENERGIES MATRIX
      DOUBLE PRECISION EIGEN(DMAX)
C     COMPUTED W2 AND SQUARED W2 (SO(3) CASIMIR) MATRICES
      DOUBLE PRECISION W2MAT(DMAX,DMAX), W4MAT(DMAX,DMAX) 
C     COMPUTED W2·n+n·W2 (SO(3) CASIMIR x U(2) CASIMIR) MATRIX
CC      DOUBLE PRECISION W2NMAT(DMAX,DMAX)
C     COMPUTED W2·n2+n2·W2 (SO(3) CASIMIR x U(2) QUAD CASIMIR) MATRIX
CC      DOUBLE PRECISION W2N2MAT(DMAX,DMAX)
C     COMPUTED W2·WBAR2+WBAR2·W2 (SO(3) CASIMIR x barSO(3) CASIMIR) MATRIX
      DOUBLE PRECISION W2W2BARMAT(DMAX,DMAX)
C
C     SCRATCH VECTORS AND VARIABLES FOR DIAGONALIZATION
C     VECTORS IN COMMON
      DOUBLE PRECISION FV1
      COMMON/SCRATCH/ FV1(DIMSCR)
      INTEGER IERR
C
C
C     COMMON BLOCK FOR HAMILTONIAN PARAMETERS
C
      DOUBLE PRECISION HPAR
C
C
      COMMON/HAMPAR/ HPAR(NPMAX)
C     
C     CONTROL OF DISPLAYED OUTPUT
      INTEGER IPRINT
      COMMON/GRAF/ IPRINT
C
C     TOTAL NUMBER OF EXPERIMENTAL DATA
      INTEGER TOTDAT
      COMMON/NUMBER_DATA/TOTDAT
C     IFLAG
      INTEGER MIFLAG
      COMMON/MINUIT_IFLAG/MIFLAG
C
C     VECTOR WITH n's(v's) ASSIGNED TO COMPUTED EIGENVECTORS 
      INTEGER BLAS(DMAX)
C     LOCAL VARIABLES                                                      
      INTEGER I,JMIN,L,DIM,NDAT
      DOUBLE PRECISION AVEZERO(DMAX), EMIN
      INTEGER VEXPAS(DATMAX,2)
      DOUBLE PRECISION VEXPEN(DATMAX), VEXPERR(DATMAX)
C     LOGICAL VARIABLE THAT CONTROL IF ALL EIGENVECTORS ARE ASSIGNED
C     HAS TO BE MIFLAG = 3 AND IPRINT > 0
      LOGICAL ASSIGNALL 
cc      DOUBLE PRECISION TCPU(5)
C     LEVELS IN THE EXPERIMENTAL ENERGY LIST BUT NOT INCLUDED IN CHI^2
      INTEGER NOT_FIT_LEVELS
C
C
C     NAMELIST DEFINITIONS
C
      NAMELIST/INP0/ BENT, DTFL
      NAMELIST/INP1/ N2, LMAX, VMAX, EMINL
      NAMELIST/INP2/ IPRINT   
C
C     OPEN INPUT FILE
      OPEN(UNIT=10,FILE='minuit_u3_fit.inp',STATUS='OLD')
C
C     READING INPUT
C
      READ(10,INP0)
      READ(10,INP1)
      READ(10,INP2)
C
      CLOSE(10)
C
      IF (IPRINT.GT.2) THEN
         WRITE(*,*) 'tSTARTING CHISQRE CALCULATION'
         WRITE(*,*) "CHISQRE :: IFLAG = ", MIFLAG
      ENDIF
C
C
C     TESTS 
C    
      IF (LMAX.GT.N2) STOP 'LMX > N2, SAYONARA '
C
C ONE BODY
C     n
      HPAR(1) = P11
C
C TWO BODY
C     n^2
      HPAR(2) = P21
C     L^2
      HPAR(3) = P22
C     W^2
      HPAR(4) = P23
C
C THREE BODY
C     n^3
      HPAR(5) = P31
C     n·l^2
      HPAR(6) = P32
C     n·W^2 + W^2·n
      HPAR(7) = P33
C
C FOUR BODY
C     n^4
      HPAR(8) = P41
C     n^2·l^2
      HPAR(9) = P42
C     l^4
      HPAR(10) = P43
C     W^2·l^2
      HPAR(11) = P44
C     n^2·W^2 + W^2·n^2
      HPAR(12) = P45
C     W^4 
      HPAR(13) = P46
C     (W^2·Wb^2 + Wb^2·W^2)/2
      HPAR(14) = P47
C     
C     
C     READING EXPERIMENTAL DATA FILE 
C     
      CALL RDXENRU3
     *     (BENT,DTFL,LMAX,VMAX,EXPEN,EXPERR,DATMAX,EXPAS,LXP)
C
C     
C     INITIALIZE CHI SQUARE
      CHSQT = 0.0D0
C     
C     INITIALIZE TOTAL NUMBER OF DATA
      TOTDAT = 0
C     
C     INITIALIZE EMIN
      EMIN = 0.0D0
C     
C     MAIN LOOP FOR CHISQRE CALCULATION
C     
      DO L = 0, LMAX
C
cc      IF (IPRINT.GT.0) THEN
cc         CALL CPU_TIME(TCPU(1))
ccC         WRITE(*,*) TCPU(1)
cc      ENDIF

C     
         IF (IPRINT.GT.3) THEN
            WRITE(*,*)
            WRITE(*,*) 'N = ',N2,' L = ',L
            WRITE(*,*)
         ENDIF
C     SELECT EXPERIMENTAL DATA FOR VIBRATIONAL ANGULAR MOMENTUM L
         CALL SLCTLU3(L,DATMAX,EXPAS,VEXPAS,EXPEN,EXPERR,
     *        LXP,VEXPEN,VEXPERR,NDAT)
C
cc         IF (IPRINT.GT.0) THEN
cc            CALL CPU_TIME(TCPU(2))
cc            WRITE(*,*) 'DELTA T_0', TCPU(2)-TCPU(1)
cc         ENDIF
C
C     BUILD HAMILTONIAN MATRIX
C     ALGEBRAIC HAMILTONIAN
         CALL HBLDU3GEN(N2,L,DMAX,HAM,W2MAT,W4MAT,W2W2BARMAT)  
C     
C     DIAGONALIZE HAMILTONIAN MATRIX
         DIM = (N2-MOD(N2-L,2)-L)/2+1
         CALL DSYEV('V', 'U', DIM, HAM, DMAX, EIGEN,  FV1,  DIMSCR,
     *        IERR)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (iprint.ge.3) then
            write(*,*) 'Hamiltonian eigenvalues'
            do i = 1, dim
               write(*,*) i, eigen(i)
            enddo
         endif
cccc         stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc         IF (IPRINT.GT.0) THEN
cc            CALL CPU_TIME(TCPU(3))
cc            WRITE(*,*) 'DELTA T_1', TCPU(3)-TCPU(2)
cc         ENDIF
C
C     ASSIGN  COMPUTED DATA TO LOCAL BASIS STATES
C
         IF (MIFLAG.EQ.3.AND.IPRINT.GT.0) THEN 
            ASSIGNALL = .TRUE.
         ELSE
            ASSIGNALL = .FALSE.
         ENDIF
         CALL ASSGNU3(BENT,N2,L,DMAX,W2MAT,W4MAT,W2W2BARMAT,
     *        HAM,EIGEN,DATMAX,VEXPAS,NDAT,ASSIGNALL,BLAS)
C
C
cc         IF (IPRINT.GT.0) THEN
cc            CALL CPU_TIME(TCPU(4))
cc            WRITE(*,*) 'DELTA T_2', TCPU(4)-TCPU(3)
cc         ENDIF
C
C     REFER ENERGIES TO GROUND L=0 STATE
C
         IF (EMINL) THEN
            JMIN = 0
            DO I = 1, DIM
               IF (BENT) THEN 
                  IF (BLAS(I).EQ.0) JMIN = I
               ELSE
                  IF (BLAS(I).EQ.L) JMIN = I
               ENDIF
            ENDDO
            EMIN = EIGEN(JMIN)
            DO I = 1, DIM
               AVEZERO(I) = HAM(I,JMIN)
            ENDDO
CC
            IF (IPRINT.GE.2) WRITE(*,*) 'J(E_MIN) = ', JMIN
CC
         ELSE IF (L.EQ.0) THEN
            JMIN = 0
            DO I = 1, DIM
               IF (BLAS(I).EQ.0) JMIN = I
            ENDDO
            EMIN = EIGEN(JMIN)
            DO I = 1, DIM
               AVEZERO(I) = HAM(I,JMIN)
            ENDDO
C     
            IF (IPRINT.GE.2) WRITE(*,*) 'J(E_MIN) = ', JMIN
C     
         ENDIF
C
         DO I = 1, DIM
            EIGEN(I) = EIGEN(I) - EMIN
         ENDDO
C
C     
C     COMPARE WITH EXPERIMENT AND COMPUTE STATISTICAL PARAMETERS
         CALL CPXPDTU3(N2,L,DATMAX,VEXPAS,EIGEN,BLAS,
     *        VEXPEN,VEXPERR,NDAT,SPNTEX,CHSQP)
C     
C     DISPLAY RESULTS
         IF (IPRINT.GT.0.OR.MIFLAG.EQ.3) CALL DSDTU3(BENT,N2,L,EIGEN,
     *        VEXPEN,VEXPERR,NDAT,DATMAX,VEXPAS,DMAX,HAM,BLAS,SPNTEX,
     *        CHSQP)
C     
C     EXCLUDE NOT FIT LEVELS
         NOT_FIT_LEVELS = 0
         DO I = 1, NDAT
            IF (VEXPERR(I).EQ.0) NOT_FIT_LEVELS = NOT_FIT_LEVELS + 1
         ENDDO
C
C     
         TOTDAT = TOTDAT + NDAT - NOT_FIT_LEVELS
C     
         CHSQT = CHSQT + CHSQP
C     
cc         IF (IPRINT.GT.0) THEN
cc            CALL CPU_TIME(TCPU(5))
cc            WRITE(*,*) 'DELTA T_3', TCPU(5)-TCPU(4)
cc         ENDIF
C
      ENDDO
C
      CHISQRE = CHSQT
C
      IF (IPRINT.GE.0) THEN
         WRITE(*,*) "NUMBER OF DATA = ", TOTDAT, "CHI2 = ", CHISQRE, 
     *        " SDEV = ", SQRT(CHISQRE/TOTDAT)
      ENDIF
C
      IF (IPRINT.GT.0.OR.MIFLAG.EQ.3)  THEN
C     HPAR(1)   P11
C     HPAR(2)   P21
C     HPAR(3)   P22
C     HPAR(4)   P23
C     HPAR(5)   P31
C     HPAR(6)   P32
C     HPAR(7)   P33
C     HPAR(8)   P41
C     HPAR(9)   P42
C     HPAR(10)  P43
C     HPAR(11)  P44
C     HPAR(12)  P45
C     HPAR(13)  P46
C     HPAR(14)  P47
C
         WRITE(*,10) HPAR(1)
         WRITE(*,20) HPAR(2),HPAR(3), HPAR(4)
         WRITE(*,30) HPAR(5),HPAR(6),HPAR(7)
         WRITE(*,40) HPAR(8),HPAR(9),HPAR(10)
         WRITE(*,41) HPAR(11),HPAR(12),HPAR(13)
         WRITE(*,42) HPAR(14)
      ENDIF
C
 10   FORMAT (1X,"ONE BODY :: ", " P11 = ", G15.5)
 20   FORMAT (1X,"TWO BODY :: ", " P21 = ", G15.6, 
     *     " P22 = ", G15.6, " P23 = ", G15.6)
 30   FORMAT (1X,"THREE BODY :: ", " P31 = ", G15.6, 
     *     " P32 = ", G15.6, " P33 = ", G15.6)
 40   FORMAT (1X,"FOUR BODY (i):: ", " P41 = ", G15.6, 
     *     " P42 = ", G15.6, " P43 = ", G15.6)
 41   FORMAT (1X,"FOUR BODY (ii):: ", " P44 = ", G15.6, 
     *     " P45 = ", G15.6, " P46 = ", G15.6)
 42   FORMAT (1X,"FOUR BODY (iii):: ", " P47 = ", G15.6)
C
      RETURN 
C
      END
      
      



