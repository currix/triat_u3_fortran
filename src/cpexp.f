      SUBROUTINE CPXPDTU3
     *     (N2,L,LD2,EXPAS,EIGEN,BLAS,EXPEN,EXPERR,NDAT,
     *     SPNTEX,CHSQP)
C
C $Id: cpexp.f,v 1.4 2011/05/04 14:19:20 curro Exp $
C
C     SUBROUTINE THAT BASED IN EXPERIMENTAL ASSIGNMENTS
C     COMPARE THE DIAGONALIZATION RESULTS WITH THE
C     EXPERIMENTAL ENERGIES IN THE BENDING SU(3) MODEL.
C     IN ADDITION IT COMPUTES STATISTICAL PARAMETERS.
C
C     
C     INPUT
C     N2      : SU(3) IRREP LABEL (BENDING)
C     L       : VIBRATIONAL ANGULAR MOMENTUM
C     LD2     : LEADING DIMENSION OF MATRICES
C     EIGEN   : EIGENVALUES VECTOR
C     BLAS    : LOCAL BASIS ASSIGNMENTS
C     EXPEN   : VECTOR WITH EXPERIMENTAL ENERGIES FOR POLYAD V
C     EXPERR  : VECTOR WITH EXPERIMENTAL ENERGY ERRORS FOR POLYAD V
C     EXPASS  : MATRIX WITH EXPERIMENTAL ASSIGNMENTS FOR POLYAD V
C     NDAT    : NUMBER OF EXPERIMENTAL DATA IN POLYAD V
C
C     OUTPUT
C     SPNTX   : POINTER WITH ASSIGNMENTS TO EXPERIMENTAL DATA
C     CHSQP   : PARTIAL CHI-SQUARE
C
C     OLD :: 15-02-01                     by Currix TM
C     
      IMPLICIT NONE
C
C     DEFINITION OF VARIABLES
C     
      INTEGER N2, L, LD2, EXPAS(LD2,*),  BLAS(*), NDAT
      DOUBLE PRECISION EIGEN(*)
      DOUBLE PRECISION EXPEN(*), EXPERR(*)
C
      INTEGER SPNTEX(*)
      DOUBLE PRECISION CHSQP
C
      INTEGER IPRINT
C
C     CONTROL OUTPUT DISPLAYED
      COMMON/GRAF/ IPRINT
C     
C     TEMPORAL VARIABLES                                     
      INTEGER DIM
      INTEGER I, J
      INTEGER NOT_FIT_LEV
C
C
      IF (IPRINT.GT.2) WRITE(*,*) 'STARTING SUBROUTINE CPXPDTU3'
C
      CHSQP = 0.0D0
      NOT_FIT_LEV = 0
C
C     POLYAD DIMENSION
      DIM = (N2-MOD(N2-L,2)-L)/2+1
C
      DO I = 1, NDAT
         DO J = 1, DIM
C     
C     COMPARE EXPERIMENTAL ASSIGNMENTS AND LOCAL BASIS ASSIGN,
            IF (BLAS(J).EQ.EXPAS(I,1)) THEN
C     POINTER
               SPNTEX(I) = J
C     CHI SQUARE               
               IF (EXPERR(I).NE.0) THEN
                  CHSQP = CHSQP + 
     *                 (EIGEN(J)-EXPEN(I))*(EIGEN(J)-EXPEN(I))/
     *                 (EXPERR(I)*EXPERR(I))
               ELSE
                  NOT_FIT_LEV = NOT_FIT_LEV + 1
               ENDIF
C
            ENDIF
C     
         ENDDO
      ENDDO
C
      IF (IPRINT.GE.2) WRITE(*,*) 'EXPER. LEVELS', NDAT, 
     *     "NOT FITTED LEVELS", NOT_FIT_LEV
      IF (IPRINT.GT.2) WRITE(*,*) 'LEAVING SUBROUTINE CPXPDTU3'
C     
      RETURN 
C     
      END
      




