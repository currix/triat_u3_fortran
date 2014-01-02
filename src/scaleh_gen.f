      SUBROUTINE SCLHAM(FAC,BENT,N2,L,LD,HAM,W2MAT,W4MAT,W2W2BARMAT)
C
C     SUBROUTINE THAT SCALES DOWN (x1/2) THE MIXING PARAMETERS AND
C     RECOMPUTES HAMILTONIAN IN THE U(3) MODEL.
C     
C     INPUT
C     FAC  : FACTOR THAT RESCALES THE MIXING PARAMETERS
C     BENT : .T. BENT MOLECULE, .F. LINEAR MOLECULE
C     N2   : U(3) IRREP LABEL (BENDING)
C     L    : VIBRATIONAL ANGULAR MOMENTUM LABEL
C     LD   : LEADING DIMENSION OF MATRICES
C     W2MAT : ARRAY WITH THE SO(3) CASIMIR P BLOCK 
C     W4MAT : ARRAY WITH THE SO(3) CASIMIR P^2 BLOCK 
C     W2W2BARMAT : ARRAY WITH THE OPERATOR W2·WBAR2 + WBAR2·W2 BLOCK 
C     
C
C     OUTPUT
C     HAM      : RECOMPUTED HAMILTONIAN MATRIX
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
      DOUBLE PRECISION FAC
      LOGICAL BENT
      INTEGER N2
      INTEGER L, LD
      DOUBLE PRECISION HAM(LD,*)
      DOUBLE PRECISION W2MAT(LD,*), W4MAT(LD,*), W2W2BARMAT(LD,*)
C
C
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
C
      COMMON/GRAF/ IPRINT
C
cC     LOCAL VARIABLE
c      INTEGER I
C
      IF (IPRINT.GT.2) WRITE(*,*) 'SUBROUTINE SCALE MIXING STARTS HERE'
C
C     RESCALE PARAMETERS
C
      IF (.NOT.BENT) THEN
         HPAR(4) = HPAR(4)*FAC         
         HPAR(11) = HPAR(11)*FAC
         HPAR(13) = HPAR(13)*FAC
      ELSE
         HPAR(1) = HPAR(1)*FAC
         HPAR(2) = HPAR(2)*FAC
         HPAR(5) = HPAR(5)*FAC
         HPAR(6) = HPAR(6)*FAC
         HPAR(8) = HPAR(8)*FAC
         HPAR(9) = HPAR(9)*FAC
      ENDIF
         HPAR(7) = HPAR(7)*FAC ! Both cases
         HPAR(12) = HPAR(12)*FAC ! Both cases
         HPAR(14) = HPAR(14)*FAC ! Both cases
         CALL HBLDU3GEN(N2,L,LD,HAM,W2MAT,W4MAT,W2W2BARMAT)  
C
      RETURN
      END
