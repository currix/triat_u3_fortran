      SUBROUTINE HBLDU3GEN(N2,L,LD,HAM,W2MAT,W4MAT,W2W2BARMAT)
C     
C     SUBROUTINE THAT BUILDS THE GENERAL 1-,2-,3-, AND 4-BODY 
C     HAMILTONIAN MATRIX IN THE U(3) MODEL FOR BENDING VIBRATIONS
C
C     INPUT
C     N2   : U(3) IRREP LABEL (BENDING)
C     L    : VIBRATIONAL ANGULAR MOMENTUM LABEL
C     LD   : LEADING DIMENSION OF HAM AND BL
C
C     OUTPUT
C     HAM     : HAMILTONIAN MATRIX
C     W2MAT   : ARRAY WITH THE SO(3) CASIMIR (W^2) BLOCK)
C     W4MAT   : ARRAY WITH THE SQUARED SO(3) CASIMIR (W^4 BLOCK) 
C     W2W2BARMAT : ARRAY WITH THE OPERATOR W^2·WBAR^2+WBAR^2·W^2 
C     
C     BENDING HAMILTONIAN
C
C     H = P11 n + 
C         P21 n^2 + P22 l^2 + P23 W^2 +  
C         P31 n^3 + P32 n·l^2 + P33 (n·W^2 + W^2·n) +
C         P41 n^4 + P42 n^2·l^2 + P43 l^4 + P44 l^2·W^2 + 
C         P45 (n^2·W^2 + W^2·n^2) + P46 W^4 + P47 (W^2·Wbar^2 + Wbar^2·W^2)/2
C
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
C     
C     by Currix TM
C     
      IMPLICIT NONE
C
      INTEGER NPMAX
C     MAXIMUM NUMBER OF PARAMETERS IN THE HAMILTONIAN
      PARAMETER (NPMAX = 15)
C     
C     DEFINITION OF VARIABLES
C     
      INTEGER N2, L, LD
C
      DOUBLE PRECISION HAM(LD,*)
      DOUBLE PRECISION W2MAT(LD,*), W4MAT(LD,*), W2W2BARMAT(LD,*)
CCC     * W2NMAT(LD,*), W2N2MAT(LD,*) 
C
C     COMMON BLOCK FOR HAMILTONIAN PARAMETERS
C
      DOUBLE PRECISION HPAR
C
C
      COMMON/HAMPAR/ HPAR(NPMAX)
C     
      INTEGER IPRINT
C
C     CONTROL OUTPUT DISPLAYED
      COMMON/GRAF/ IPRINT
C
C     LOCAL VARIABLES
      INTEGER I, J, I1, DIM
C     TEMPORAL VALUES TO STORE INTEGERS AS DOUBLE PREC. REALS
      DOUBLE PRECISION V2, VL
C     
C     
      IF (IPRINT.GT.2) 
     *     WRITE(*,*) 'HAMILTONIAN BUILDING SUBROUTINE STARTS HERE'
C
C     MATRIX DIMENSIONS
      DIM = (N2-MOD(N2-L,2)-L)/2+1
C     INITIALIZATION
      DO I = 1, DIM
         DO J = 1, DIM
            HAM(J,I) = 0.0D0
         ENDDO
      ENDDO
C     
      VL = DFLOAT(L)
C
C     NON-DIAGONAL INTERACTIONS
C
C     BUILDING W2 BLOCK
C
      CALL SO3CASBUILD(LD,W2MAT,N2,L)
C
C     HAM MULTIPLYING W2 TIMES P23 PARAMETER
      DO I1 = 1, DIM
         HAM(I1,I1) = HPAR(4)*W2MAT(I1,I1)
      ENDDO
      DO I1 = 1, DIM-1
         HAM(I1+1,I1) = HPAR(4)*W2MAT(I1+1,I1)
         HAM(I1,I1+1) = HPAR(4)*W2MAT(I1,I1+1)
      ENDDO
C
C     BUILDING W^2·n + n·W^2 BLOCK
C
      IF (HPAR(7).NE.0) THEN
C
C     HAM MULTIPLYING (W^2·n + n·W^2) TIMES P33 PARAMETER
         DO I1 = 1, DIM
            V2 = DFLOAT(N2 - (2*I1-2+MOD(N2-L,2)))
            HAM(I1,I1) = HAM(I1,I1) + HPAR(7)*2.0D0*V2*W2MAT(I1,I1)
         ENDDO
         DO I1 = 1, DIM-1
            V2 = DFLOAT(N2 - (2*I1-2+MOD(N2-L,2)))
            HAM(I1+1,I1) = HAM(I1+1,I1) + 
     *           HPAR(7)*2.0D0*(V2-1.0D0)*W2MAT(I1+1,I1)
            HAM(I1,I1+1) = HAM(I1,I1+1) + 
     *           HPAR(7)*2.0D0*(V2-1.0D0)*W2MAT(I1,I1+1)
         ENDDO
      ENDIF
C
C     HAM MULTIPLYING W^2·l^2 TIMES P44 PARAMETER
      IF (HPAR(11).NE.0) THEN
         DO I1 = 1, DIM
            HAM(I1,I1) = HAM(I1,I1) + HPAR(11)*VL*VL*W2MAT(I1,I1)
         ENDDO
         DO I1 = 1, DIM-1
            HAM(I1+1,I1)=HAM(I1+1,I1) + HPAR(11)*VL*VL*W2MAT(I1+1,I1)
            HAM(I1,I1+1)=HAM(I1,I1+1) + HPAR(11)*VL*VL*W2MAT(I1,I1+1)
         ENDDO
      ENDIF
C
C     BUILDING W^2·n^2 + n^2·W^2 BLOCK
C
      IF (HPAR(12).NE.0) THEN
C
C     HAM MULTIPLYING (W^2·n^2 + n^2·W^2) TIMES P45 PARAMETER
         DO I1 = 1, DIM
            V2 = DFLOAT(N2 - (2*I1-2+MOD(N2-L,2)))
            HAM(I1,I1) = HAM(I1,I1) + HPAR(12)*2.0D0*V2*V2*W2MAT(I1,I1)
         ENDDO
         DO I1 = 1, DIM-1
            V2 = DFLOAT(N2 - (2*I1-2+MOD(N2-L,2)))
            HAM(I1+1,I1) = HAM(I1+1,I1) + 
     *           HPAR(12)*2.0D0*(V2*V2-2.0D0*V2+2.0D0)*W2MAT(I1+1,I1)
            HAM(I1,I1+1) = HAM(I1,I1+1) + 
     *           HPAR(12)*2.0D0*(V2*V2-2.0D0*V2+2.0D0)*W2MAT(I1,I1+1)
         ENDDO
      ENDIF
C
C     BUILDING W2·WBAR2 + WBAR2·W2 BLOCK (HAS TO BE EVALUATED BEFORE W4 BLOCK)
C
      IF (HPAR(14).NE.0) THEN
         CALL SO3SO3BARBUILD(LD,W2MAT,W4MAT,W2W2BARMAT,N2,L)
C
C     HAM MULTIPLYING W2·WBAR2 + WBAR2·W2 TIMES P47 PARAMETER
         DO I1 = 1, DIM
            HAM(I1,I1) = HAM(I1,I1) + HPAR(14)*W2W2BARMAT(I1,I1)
         ENDDO
         DO I1 = 1, DIM-1
            HAM(I1+1,I1) = HAM(I1+1,I1) + HPAR(14)*W2W2BARMAT(I1+1,I1)
            HAM(I1,I1+1) = HAM(I1,I1+1) + HPAR(14)*W2W2BARMAT(I1,I1+1)
         ENDDO
         DO I1 = 1, DIM-2
            HAM(I1+2,I1) = HAM(I1+2,I1) + HPAR(14)*W2W2BARMAT(I1+2,I1)
            HAM(I1,I1+2) = HAM(I1,I1+2) + HPAR(14)*W2W2BARMAT(I1,I1+2)
         ENDDO
      ENDIF
C
C     BUILDING W4 BLOCK
C
      IF (HPAR(13).NE.0) THEN
         CALL SO32CASBUILD(LD,W2MAT,W4MAT,N2,L)
C     
C     HAM MULTIPLYING P^2 TIMES P46 PARAMETER
         DO I1 = 1, DIM
            HAM(I1,I1) = HAM(I1,I1) + HPAR(13)*W4MAT(I1,I1)
         ENDDO
         DO I1 = 1, DIM-1
            HAM(I1+1,I1) = HAM(I1+1,I1) + HPAR(13)*W4MAT(I1+1,I1)
            HAM(I1,I1+1) = HAM(I1,I1+1) + HPAR(13)*W4MAT(I1,I1+1)
         ENDDO
         DO I1 = 1, DIM-2
            HAM(I1+2,I1) = HAM(I1+2,I1) + HPAR(13)*W4MAT(I1+2,I1)
            HAM(I1,I1+2) = HAM(I1,I1+2) + HPAR(13)*W4MAT(I1,I1+2)
         ENDDO
      ENDIF
C
C     DIAGONAL INTERACTIONS
C     
      DO I = 1, DIM
C
         V2 = DFLOAT(N2 - (2*I-2+MOD(N2-L,2)))
C     
         HAM(I,I) = HAM(I,I) +
C     
     *        HPAR(1)*V2 +
C     
     *        HPAR(2)*V2*V2 +
C     
     *        HPAR(3)*VL*VL +
C     
     *        HPAR(5)*V2*V2*V2 + 
C     
     *        HPAR(6)*V2*VL*VL +
C     
     *        HPAR(8)*V2*V2*V2*V2 +
C     
     *        HPAR(9)*V2*V2*VL*VL +
C     
     *        HPAR(10)*VL*VL*VL*VL 
C     
      ENDDO
C     
CCCCCCCCCCCCCC
      if (iprint.gt.4) then
         WRITE(*,*)
         WRITE(*,*) 'Hamiltonian matrix'
         DO I1 = 1, DIM
            WRITE(*,10) (ham(I,I1),I=1,DIM)
         ENDDO
         WRITE(*,*) 
      endif
 10   format(1x,10f12.5)
CCCCCCCCCCCCCC
C
      IF (IPRINT.GT.2) 
     *     WRITE(*,*)'HAMILTONIAN BUILDING SUBROUTINE ENDS HERE'
C
      RETURN
      END
      
