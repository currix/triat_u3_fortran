      SUBROUTINE 
     *     RDXENRU3(BENT,DTFL,LMAX,VMAX,EXPEN,EXPERR,LD,EXPAS,LXPOIN)
C
C     SUBROUTINE THAT READS EXPERIMENTAL DATA (ENERGIES, ERRORS 
C     ASSIGNMENTS IN THE U(3) MODEL FOR BENDING VIBRATIONS.
C
C     FORMAT OF THE INPUT FILE:
C     -------------------------
C     N   # NUMBER OF DATA
C     E1   ERR1   V1  L1
C     E2   ERR2   V2  L2
C             .................
C
C     EN   ERRN   VN2 LN
C
C     THE VALUES HAVE NOT TO BE IN ANY DETERMINED ORDER 
C     AND IN CASE ERRORS ARE UNKNOWN ALL MUST BE GIVEN 
C     A 1.0 VALUE.
C
C     INPUT
C     BENT : .T. BENT MOLECULE, .F. LINEAR MOLECULE
C     DTFL : FILE WITH EXPERIMENTAL DATA SET
C     LMAX : MAXIMUM VIBRATIONAL ANGULAR MOMENTUM LABEL
C     VMAX : MAXIMUM NUMBER OF QUANTA OF VIBRATION
C     LD   : LEADING DIMENSION OF EXPAS
C
C     OUTPUT
C     EXPEN  : VECTOR WITH EXPERIMENTAL ENERGIES
C     EXPERR : VECTOR WITH EXPERIMENTAL ENERGY ERRORS
C     EXPAS  : MATRIX WITH EXPERIMENTAL ASSIGNMENTS
C     LXPOIN : POINTER LOCATING VALUES OF L
C     
C
C     by Currix TM
C     
      IMPLICIT NONE
C     
C     DEFINITION OF VARIABLES
C     
      LOGICAL BENT
      CHARACTER*30 DTFL
      INTEGER LMAX, VMAX, LD, EXPAS(LD,*), LXPOIN(0:*)
      DOUBLE PRECISION EXPEN(*), EXPERR(*)
C
      INTEGER IPRINT
C
C     CONTROL OUTPUT DISPLAYED
      COMMON/GRAF/ IPRINT
C
C     
C     TEMPORAL VARIABLES 
      INTEGER I, J, K, M, INDX, JT, NT, COUNT, LTEMP
      DOUBLE PRECISION TMP1, TMP2
      INTEGER TMP(3)
C
      IF (IPRINT.GT.2) WRITE(*,*) 'SUBROUTINE READENERG'
C
C     INITIALIZE POINTER
      DO I = 0, LD
         LXPOIN(I) = 0
      ENDDO
C
C     READ ALL EXPERIMENTAL DATA SELECTING RIGHT V'S AND L
C
      OPEN(UNIT=30,FILE=DTFL,STATUS='OLD')
C
      READ(30,*) NT
C
C
      COUNT = 0
      DO I = 1, NT
         READ(30,*) TMP1, TMP2, (TMP(J), J = 1, 2)
C  
C     TEST LINEAR LABELS
         IF (.NOT.BENT.AND.MOD(TMP(1),2).NE.MOD(TMP(2),2)) THEN
            WRITE(*,*) "ERROR: LINEAR CASE WITH n = ",TMP(1),
     *           " l = ", TMP(2)
            STOP 'ERROR READING EXPERIMENTAL ENERGIES'
         ENDIF
C
         IF (TMP(2).LE.LMAX.AND.TMP(1).LE.VMAX) THEN
            COUNT = COUNT+1
            EXPEN(COUNT) = TMP1
            EXPERR(COUNT) = TMP2
            DO J = 1, 2
               EXPAS(COUNT,J) = TMP(J)
            ENDDO
            LXPOIN(TMP(2)+1) = LXPOIN(TMP(2)+1) + 1
         ENDIF
      ENDDO
      DO I = 1, LMAX+1
         LXPOIN(I) = LXPOIN(I)+LXPOIN(I-1)
      ENDDO
C
C     ORDERING DATA ACCORDING TO VIBRATIONAL ANGULAR MOMENTUM L
C     
      DO 10 J = 2, COUNT
         LTEMP = EXPAS(J,2)
C     
         TMP1 = EXPEN(J)
         TMP2 = EXPERR(J)
         DO K = 1, 2
            TMP(K) = EXPAS(J,K)
         ENDDO
C     
         JT = J
C     
         DO I = J-1, 1, -1
            IF (EXPAS(I,2).GT.LTEMP) THEN
               EXPEN(JT) = EXPEN(I)
               EXPERR(JT) = EXPERR(I)
               DO K = 1, 2
                  EXPAS(JT,K) = EXPAS(I,K)
               ENDDO
               JT = I
            ELSE
               GOTO 5
            ENDIF
         ENDDO
 5       EXPEN(JT) = TMP1
         EXPERR(JT) = TMP2
         DO K = 1, 2
            EXPAS(JT,K) = TMP(K)
         ENDDO
 10   CONTINUE
C     
C     ORDERING IN INCREASING ENERGY FOR EACH L 
C     
      DO I = 0, LMAX
         DO 20 J = 2, LXPOIN(I+1)-LXPOIN(I)
C     INDEX-1 OF THE FIRST ELEMENT ELEMENT OF THE I-TH POLYAD
            INDX = LXPOIN(I)
C     
            TMP1 = EXPEN(INDX+J)
            TMP2 = EXPERR(INDX+J)
            DO K = 1, 2
               TMP(K) = EXPAS(INDX+J,K)
            ENDDO
C     
            JT = J
C     
            DO M = J-1, 1, -1
               IF (EXPEN(INDX+M).GT.TMP1) THEN
                  EXPEN(INDX+JT) = EXPEN(INDX+M)
                  EXPERR(INDX+JT) = EXPERR(INDX+M)
                  DO K = 1, 2
                     EXPAS(INDX+JT,K) = EXPAS(INDX+M,K)
                  ENDDO
                  JT = M
               ELSE
                  GOTO 15
               ENDIF
            ENDDO
 15         EXPEN(INDX+JT) = TMP1
            EXPERR(INDX+JT) = TMP2
            DO K = 1, 2
               EXPAS(INDX+JT,K) = TMP(K)
            ENDDO
 20      CONTINUE
      ENDDO
C
      IF (IPRINT.GT.2) THEN
         WRITE(*,*)
         WRITE(*,*) 'EXPERIMENTAL ENERGIES'
         WRITE(*,*)
         DO I = 0, LMAX
            WRITE(*,*) 'VIBRATIONAL ANGULAR MOMENTUM', I
            WRITE(*,*) 'LXPOIN', LXPOIN(I+1), LXPOIN(I)
            DO J = 1, LXPOIN(I+1) - LXPOIN(I)
               WRITE(*,*)(EXPAS(LXPOIN(I)+J,K),K=1,2),
     *              EXPEN(LXPOIN(I)+J), EXPERR(LXPOIN(I)+J)
            ENDDO
         ENDDO
         WRITE(*,*)
      ENDIF
C
      CLOSE(30)
C
      RETURN
C
      END
      





