      SUBROUTINE FCN(NPAR,GRAD,FVAL,XVAL,IFLAG,CHISQRE)
C     
      IMPLICIT NONE
      DOUBLE PRECISION GRAD(*),XVAL(*), FVAL
      INTEGER IFLAG, NPAR
      DOUBLE PRECISION CHISQRE
      EXTERNAL CHISQRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CONTROL OF DISPLAYED OUTPUT
      INTEGER IPRINT
      COMMON/GRAF/ IPRINT
C     TOTAL NUMBER OF EXPERIMENTAL DATA
      INTEGER TOTDAT
      COMMON/NUMBER_DATA/TOTDAT
C     IFLAG
      INTEGER MIFLAG
      COMMON/MINUIT_IFLAG/MIFLAG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      MIFLAG = IFLAG
      IF (IFLAG .EQ. 1)  THEN
         IF (IPRINT.GT.1) WRITE(*,*) " FCN iflag = ", IFLAG
C     READ INPUT DATA
c     add here the namelist and experimental energy reading ...
C     calculate any necessary constants, etc.
      ENDIF
C     Always calculate the value of the function, FVAL,
C     which is usually a chisquare or log likelihood.
C     Optionally, calculation of FVAL may involve
ccccccccccccccccccccccc
C      if (iflag.eq.3.and.iprint.lt.0) iprint = 0
ccccccccccccccccccccccc
      FVAL = CHISQRE(XVAL(1),XVAL(2),XVAL(3),XVAL(4),XVAL(5),XVAL(6),
     *     XVAL(7),XVAL(8),XVAL(9),XVAL(10),XVAL(11),
     *     XVAL(12),XVAL(13),XVAL(14))
C     It is responsability of user to pass
C     any parameter values needed by FUTIL,
C     either through arguments, or in a COMMON block
C     XVAL(1)   P11
C     XVAL(2)   P21
C     XVAL(3)   P22
C     XVAL(4)   P23
C     XVAL(5)   P31
C     XVAL(6)   P32
C     XVAL(7)   P33
C     XVAL(8)   P41
C     XVAL(9)   P42
C     XVAL(10)  P43
C     XVAL(11)  P44
C     XVAL(12)  P45
C     XVAL(13)  P46
C     XVAL(14)  P47
C
      IF (IFLAG .EQ. 3)  THEN
C
         WRITE(*,10) XVAL(1)
         WRITE(*,20) XVAL(2),XVAL(3),XVAL(4)
         WRITE(*,30) XVAL(5),XVAL(6),XVAL(7)
         WRITE(*,40) XVAL(8),XVAL(9),XVAL(10)
         WRITE(*,41) XVAL(11),XVAL(12),XVAL(13)
         WRITE(*,42) XVAL(14)
         WRITE(*,60) FVAL, SQRT(FVAL/TOTDAT)
C
      ENDIF
C
      RETURN
C
C
 10   FORMAT (3X,"ONE BODY :: ", " P11 = ", G15.5)
 20   FORMAT (3X,"TWO BODY :: ", " P21 = ", G15.6, 
     *     " P22 = ", G15.6, " P23 = ", G15.6)
 30   FORMAT (3X,"THREE BODY :: ", " P31 = ", G15.6, 
     *     " P32 = ", G15.6, " P33 = ", G15.6)
 40   FORMAT (3X,"FOUR BODY (i):: ", " P41 = ", G15.6, 
     *     " P42 = ", G15.6, " P43 = ", G15.6)
 41   FORMAT (3X,"FOUR BODY (ii):: ", " P44 = ", G15.6, 
     *     " P45 = ", G15.6, " P46 = ", G15.6)
 42   FORMAT (3X,"FOUR BODY (iii):: ", " P47 = ", G15.6)
C
 60   FORMAT (1X, "CHISQRE = ", F15.6, 2X, "SDEV = ", F15.6)
      END

