      SUBROUTINE MAXCU3EXP(BENT,LDEXP,VEXPAS,NDAT,N2,L,LD,EVEC,DIM,LMC,
     *     FLAG)
C
C $Id: maxc_exp.f,v 1.3 2011/05/04 17:30:54 curro Exp $
C
C
C     INPUT
C     BENT  : .T. BENT MOLECULE, .F. LINEAR MOLECULE
C     LDEXP : LEADING DIMENSION OF MATRIX VEXPAS
C     VEXPAS: EXPERIMENTAL ASSIGNMENTS FOR POLYAD (L) FORMAT:(V l) OR (n l)
C     NDAT  : NUMBER OF EXPERIMENTAL DATA FOR POLYAD (L)
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
C     OLD :: 16-02-01                     by Currix TM
C     
      IMPLICIT NONE
C     
C     DEFINITION OF VARIABLES
C     
      LOGICAL BENT
      INTEGER LDEXP, VEXPAS(LDEXP,*), NDAT, N2, L, LD, DIM
      DOUBLE PRECISION EVEC(LD,*)
C     
      LOGICAL FLAG
      INTEGER LMC(*)
C     
      INTEGER IPRINT
C     
      COMMON/GRAF/ IPRINT
C     
      INTEGER I, I2, I3, J, JMAX, NASSIGNED
      DOUBLE PRECISION VMAX,VT
C     
C     COMPUTATION OF MAXIMAL COMPONENTS
C     
      IF (IPRINT.GT.2) WRITE(*,*) 'SUBROUTINE MAXC EXP'
C
      FLAG = .FALSE.
C
      NASSIGNED = 0
C
      DO I = 1, DIM
         LMC(I) = -1
      ENDDO
C
      IF (BENT) THEN
         DO I=1, DIM
            VMAX = 0.0D0
            DO J=1, DIM
               VT = EVEC(J,I) * EVEC(J,I)
               IF (VT.GT.VMAX) THEN
                  JMAX = J
                  VMAX = VT
               END IF
            END DO
C     ASSIGN GROUND STATE OF EACH POLYAD
            IF ((JMAX - 1).EQ.0) THEN
               LMC(I) = 0
C
               IF (VMAX.LT.0.5D0) FLAG = .TRUE.
C
            ENDIF
C
C     COMPARE WITH EXPERIMENTAL DATA
            DO I2=1, NDAT
               IF (VEXPAS(I2,1).EQ.(JMAX - 1)) THEN
                  NASSIGNED = NASSIGNED + 1
                  LMC(I) = VEXPAS(I2,1)
C
                  DO I3 = 1, I-1
                     IF (LMC(I3).EQ.(JMAX - 1)) FLAG = .TRUE.
                  END DO
C     
                  IF (VMAX.LT.0.5D0) FLAG = .TRUE.
C
                  EXIT
C
               ENDIF
            ENDDO

C
            IF (NASSIGNED.EQ.NDAT) EXIT
C
         END DO
C
      ELSE
C
         DO I=1, DIM
            VMAX = 0.0D0
            DO J=1, DIM
               VT = EVEC(J,I) * EVEC(J,I)
               IF (VT.GT.VMAX) THEN
                  JMAX = J
                  VMAX = VT
               END IF
            END DO
C     ASSIGN GROUND STATE OF EACH POLYAD
            IF ((N2 - 2*JMAX+2-MOD(N2-L,2)).EQ.L) THEN
               LMC(I) = L
C
               IF (VMAX.LT.0.5D0) FLAG = .TRUE.
C     
            ENDIF
C     
C     COMPARE WITH EXPERIMENTAL DATA
            DO I2=1, NDAT
               IF (VEXPAS(I2,1).EQ.(N2 - 2*JMAX+2-MOD(N2-L,2))) THEN
                  NASSIGNED = NASSIGNED + 1
                  LMC(I) = VEXPAS(I2,1)
C
                  DO I3=1, I-1
                     IF (LMC(I3).EQ.(N2 - 2*JMAX+2-MOD(N2-L,2))) 
     *                    FLAG = .TRUE.
                  END DO
C     
                  IF (VMAX.LT.0.5D0) FLAG = .TRUE.
C
                  EXIT
               ENDIF
            ENDDO
C     
            IF (NASSIGNED.EQ.NDAT) EXIT
C
         END DO
C
      ENDIF
C FLAG = .TRUE. IF NOT ALL EXPERIMENTAL ENERGY LEVELS ARE ASSIGNED
      IF (NASSIGNED.NE.NDAT) FLAG = .TRUE.
C     
      RETURN
C     
      END




