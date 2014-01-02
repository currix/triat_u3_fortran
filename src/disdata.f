      SUBROUTINE DSDTU3(BENT,N2,L,EIGEN,EXPEN,EXPERR,NDAT,LD,EXPAS,
     *     LD2,AVEC,BLAS,SPNTX,CHSQP)
C
C     SUBROUTINE THAT DISPLAYS THE RESULTS OF THE HAMILTONIAN 
C     DIAGONALIZATION IN THE U(3) ALGEBRAIC MODEL.
C
C     
C     INPUT
C     BENT   : .T. BENT MOLECULE, .F. LINEAR MOLECULE
C     L      : VIBRATIONAL ANGULAR MOMENTUM
C     EIGEN  : EIGENVALUES VECTOR
C     EXPEN  : VECTOR WITH EXPERIMENTAL ENERGIES
C     EXPERR : VECTOR WITH EXPERIMENTAL ERRORS
C     NDAT   : NUMBER OF EXPERIMENTAL DATA TO CONSIDER                
C     LD     : LEADING DIMENSION OF MATRICES
C     EXPAS  : MATRIX WITH EXPERIMENTAL ASSIGNMENTS
C     AVEC   : EIGENVECTORS
C     BLAS   : POINTER TO LOCAL BASIS
C     SPNTX  : POINTER WITH ASSIGNMENTS TO EXPERIMENTAL DATA
C     CHSQP  : PARTIAL CHI-SQUARE
C
C     by Currix TM
C     
      IMPLICIT NONE
C
C     DEFINITION OF VARIABLES
C     
      LOGICAL BENT
      INTEGER N2, L, LD2, LD, EXPAS(LD,*)
      INTEGER BLAS(*), SPNTX(*), NDAT
      DOUBLE PRECISION EIGEN(*),  AVEC(LD2,*)
      DOUBLE PRECISION CHSQP
      DOUBLE PRECISION EXPEN(*), EXPERR(*)
C
      INTEGER IPRINT
C
C     CONTROL OUTPUT DISPLAYED
      COMMON/GRAF/ IPRINT
C     
C     TEMPORAL VARIABLES                                     
      INTEGER I, J, DIM
C
C
      IF (IPRINT.GT.2) WRITE(*,*) 'SUBROUTINE DSDT'
C
C     POLYAD DIMENSION
      DIM = (N2-MOD(N2-L,2)-L)/2+1
C
C     
      IF (BENT) THEN 
         WRITE(*,1061) N2, L
      ELSE
         WRITE(*,1060) N2, L 
      ENDIF
      WRITE(*,*)
      WRITE(*,1070) NDAT, CHSQP
      WRITE(*,*)
C     
      WRITE(*,1080)
C     
      IF (BENT) THEN
C     
C     BENT CASE
C     
         IF (IPRINT.GE.0.AND.NDAT.NE.0) THEN
            WRITE(*,1086)
            DO I = 1, NDAT
               WRITE(*,1090) (EXPAS(I,J),J=1,2), EIGEN(SPNTX(I)),
     *              EXPEN(I),EXPERR(I), EIGEN(SPNTX(I))-EXPEN(I)
            ENDDO
         ENDIF
C
         IF (IPRINT.GE.1) THEN
C
            WRITE(*,1080)
            DO I = 1, DIM
               WRITE(*,1110) EIGEN(I), BLAS(I), L
            ENDDO
         ENDIF     
C         
         IF (IPRINT.GE.2) THEN
C
            WRITE(*,1080)
            DO I = 1, DIM
               WRITE(*,1110) EIGEN(I), BLAS(I), L
               WRITE(*,*)
               DO J = 1, DIM
                  WRITE(*,1120) AVEC(J,I), N2-(2*J-2+MOD(N2-L,2)), L
               ENDDO
               WRITE(*,*)
            ENDDO
         ENDIF
C     
      ELSE
C     
C     LINEAR CASE
C     
         IF (IPRINT.GE.0.AND.NDAT.NE.0) THEN
            WRITE(*,1085)
            DO I = 1, NDAT
               WRITE(*,1090) (EXPAS(I,J),J=1,2), EIGEN(SPNTX(I)),
     *              EXPEN(I),EXPERR(I),EIGEN(SPNTX(I))-EXPEN(I)
            ENDDO
         ENDIF
C
         IF (IPRINT.GE.1) THEN
C
            WRITE(*,1080)
            DO I = 1, DIM
               WRITE(*,1110) EIGEN(I), BLAS(I), L
            ENDDO
         ENDIF     
C         
         IF (IPRINT.GE.2) THEN
C
            WRITE(*,1080)
            DO I = 1, DIM
               WRITE(*,1110) EIGEN(I), BLAS(I), L
               WRITE(*,*)
               DO J = 1, DIM
                  WRITE(*,1120) AVEC(J,I), N2-(2*J-2+MOD(N2-L,2)), L
               ENDDO
               WRITE(*,*)
            ENDDO
         ENDIF
C     
      ENDIF
C     
      WRITE(*,1080)
      WRITE(*,*)
      WRITE(*,*)
C
      RETURN
C      
 1060 FORMAT
     *     (4X,'LINEAR CASE',6X,'N2 = ',I4,
     *     4X,'VIBR. ANG. MOM. = ',I4)
 1061 FORMAT
     *     (4X,'BENT CASE',6X,'N2 = ',I4,
     *     4X,'VIBR. ANG. MOM. = ',I4)
 1070 FORMAT(4X,'NUMBER OF EXP. DATA = ', I4, 6X, 'PARTIAL CHI SQ.',
     *     D15.5)
 1080 FORMAT(4X,
     *'***********************************************************
     *',/,
     *     4X,
     *'***********************************************************')
 1085 FORMAT(5X,'EXP.(n  L)',4X,'Calc.Energy',6X,'Exp.Energy',
     *     3X,'Calc.-Exp.')
 1086 FORMAT(5X,'EXP.(v  K)',4X,'Calc.Energy',6X,'Exp.Energy',
     *     3X,'Calc.-Exp.')
 1090 FORMAT(7X,'(',I3,I3,')',2X,F12.4,2X,F10.2,'(',F4.1,')',2X,F10.5)
 1110 FORMAT(7X,D15.8,2X,'(',I3,I3,')')
 1120 FORMAT(7X,D15.8,7X,'|',I3,I3,' >')
C
      END






