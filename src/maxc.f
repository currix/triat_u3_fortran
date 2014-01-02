      SUBROUTINE MAXCU3(BENT,N2,L,LD,EVEC,DIM,LMC,FLAG)
C
C     SUBROUTINE THAT LOOKS FOR MAXIMAL TYPE I BASIS (LINEAR)
C     COMPONENTS IN THE EIGENVECTORS OF THE BENDING SU(3)
C     MODEL HAMILTONIAN.
C
C     INPUT
C     BENT  : .T. BENT MOLECULE, .F. LINEAR MOLECULE
C     N2    : SU(3) IRREP LABEL (BENDING)
C     L     : VIBRATIONAL ANGULAR MOMENTUM
C     LD    : LEADING DIMENSION OF EXPAS
C     EVEC  : EIGENVECTORS MATRIX                   
C     DIM   : POLYAD DIMENSION               
C
C     OUTPUT
C     FLAG  :  .TRUE. IF THERE IS AMBIGUITY IN THE ASSIGNMENT
C     LMC   :  LOCAL BASIS ASSIGNMENTS 
C
C     16-02-01                     by Currix TM
C     
      IMPLICIT NONE
C     
C     DEFINITION OF VARIABLES
C     
      LOGICAL BENT
      INTEGER N2, L, LD, DIM
      DOUBLE PRECISION EVEC(LD,*)
C     
      LOGICAL FLAG
      INTEGER LMC(*)
C     
      INTEGER IPRINT
C     
      COMMON/GRAF/ IPRINT
C     
      INTEGER I, I2, J, JMAX
      DOUBLE PRECISION VMAX,VT
C     
C     COMPUTATION OF MAXIMAL COMPONENTS
C     
      IF (IPRINT.GT.2) WRITE(*,*) 'SUBROUTINE MAXC'
C
      FLAG = .FALSE.
C
      IF (BENT) THEN
         DO I=1, DIM
            VMAX = 0.0D0
            LMC(I) = -1
            DO J=1, DIM
               VT = EVEC(J,I) * EVEC(J,I)
               IF (VT.GT.VMAX) THEN
                  JMAX = J
                  VMAX = VT
               END IF
            END DO
            DO I2=1, I-1
               IF (LMC(I2).EQ.(JMAX - 1)) THEN
                  FLAG = .TRUE.
               END IF
            END DO
            IF (VMAX.LT.0.5D0) THEN
               FLAG = .TRUE.
            END IF
            LMC(I) = JMAX - 1
         END DO
      ELSE
         DO I=1, DIM
            VMAX = 0.0D0
            LMC(I) = -1
            DO J=1, DIM
               VT = EVEC(J,I) * EVEC(J,I)
               IF (VT.GT.VMAX) THEN
                  JMAX = J
                  VMAX = VT
               END IF
            END DO
            DO I2=1, I-1
               IF (LMC(I2).EQ.(N2 - 2*JMAX+2-MOD(N2-L,2))) THEN
                  FLAG = .TRUE.
               END IF
            END DO
            IF (VMAX.LT.0.5D0) THEN
               FLAG = .TRUE.
            END IF
            LMC(I) =  N2 - 2*JMAX+2-MOD(N2-L,2)
         END DO
      ENDIF
C     
      RETURN
C     
      END




