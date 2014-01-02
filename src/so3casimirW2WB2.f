      SUBROUTINE SO3SO3BARBUILD(DIMBMAX,W2BLOCK,WB2BLOCK,W2WB2BLOCK,
     *     NN,L)
C
C     SUBROUTINE THAT BUILD THE SO(3) CASIMIR (W^2·Wbar^2 + Wbar^2·W2)/2 MATRIX IN
C     THE CYLINDRICAL OSCILLATOR BASIS. 
C     SO(3) = {L,D+,D-}
C     SObar(3) = {L,R+,R-}
C
C     PROGRAM OF THE U(3) ALGEBRAIC MODEL FOR BENDING DYNAMICS
C     
C     INPUT 
C     DIMBMAX   : LEADING DIMENSION OF W2BLOCK
C     W2BLOCK   ; LOCAL BASIS MATRIX OF W2
C     NN        : U(3) REPRESENTATION (TOTAL NUMBER OF U(3) BOSONS)
C     L         : VIBRATIONAL ANGULAR MOMENTUM
C     WB2BLOCK  : ARRAY WHERE TO BUILD THE SObar(3) CASIMIR
C
C     OUTPUT
C     W2WB2BLOCK : LOCAL BASIS MATRIX OF (W^2·Wbar^2 + Wbar^2·W2)/2
C     
C
C     
C     by Currix TM
C
C     
      IMPLICIT NONE
C     
C     DEFINITION OF VARIABLES
C     
      INTEGER DIMBMAX, NN, L
      DOUBLE PRECISION W2BLOCK(DIMBMAX,*), WB2BLOCK(DIMBMAX,*), 
     *     W2WB2BLOCK(DIMBMAX,*)
C
      INTEGER IPRINT
C
C     CONTROL OUTPUT DISPLAYED
      COMMON/GRAF/ IPRINT
C
C     TEMPORAL VARIABLES                                     
      INTEGER I1, I2, DIMB
      DOUBLE PRECISION VN2, VNN, VL, VTMP
C
      IF (IPRINT.GT.2) 
     *     WRITE(*,*) 'SUBROUTINE SO3SO3BARBUILD STARTS HERE'
C
ccccccccccccccccccccc
      if (iprint.gt.4) write(*,*) 'L = ',L
ccccccccccccccccccccc
C
      DIMB = (NN - MOD(NN-L,2) - L)/2 + 1
C
C     INITIALIZING
C
      DO I1 = 1, DIMB
         DO I2 = 1, DIMB 
            W2BLOCK(I2,I1) = 0.0D0
            WB2BLOCK(I2,I1) = 0.0D0
            W2WB2BLOCK(I2,I1) = 0.0D0
         ENDDO
      ENDDO
C
      VNN = DFLOAT(NN)
      VL = DFLOAT(L)
C
C     NON-DIAGONAL PART
      DO I1 = 1, DIMB-1
         VN2 = DFLOAT(NN - MOD(NN-L,2) - 2*(I1-1))
         VTMP =
     *       DSQRT((VNN-VN2+2.0D0)*(VNN-VN2+1.0D0)*(VN2+VL)*(VN2-VL))
         W2BLOCK(I1+1,I1) = -VTMP
         W2BLOCK(I1,I1+1) = -VTMP
         WB2BLOCK(I1+1,I1) = VTMP
         WB2BLOCK(I1,I1+1) = VTMP
      ENDDO
C
C     DIAGONAL PART
      DO I1 = 1, DIMB                     
         VN2 = DFLOAT(NN - MOD(NN-L,2) - 2*(I1-1))
         W2BLOCK(I1,I1)=(VNN-VN2)*(VN2+2.0D0)+(VNN-VN2+1.0D0)*VN2+VL*VL
         WB2BLOCK(I1,I1) = W2BLOCK(I1,I1)
      ENDDO
C
C
C     MATRIX PRODUCT  (W2·Wbar2)/2
      CALL DGEMM('N','N', DIMB, DIMB, DIMB, 0.5D0,
     $     W2BLOCK, DIMBMAX, WB2BLOCK, DIMBMAX,
     $     0.0D0, W2WB2BLOCK, DIMBMAX)      
C
cccccccccccccccccc
      if (iprint.gt.4) then
         WRITE(*,*)
         WRITE(*,*) 'operator matrix (1)'
         DO I1 = 1, DIMb
            WRITE(*,10) (w2wb2block(I1,I2),I2=1,DIMb)
         ENDDO
         WRITE(*,*)
      endif
cccccccccccccccccc
C     MATRIX PRODUCT  +(Wbar2·W2)/2
      CALL DGEMM('N','N', DIMB, DIMB, DIMB, 0.5D0,
     $     WB2BLOCK, DIMBMAX, W2BLOCK, DIMBMAX,
     $     1.0D0, W2WB2BLOCK, DIMBMAX)      
C
cccccccccccccccccc
      if (iprint.gt.4) then
         WRITE(*,*)
         WRITE(*,*) 'operator matrix'
         DO I1 = 1, DIMb
            WRITE(*,10) (w2wb2block(I2,I1),I2=1,DIMb)
         ENDDO
         WRITE(*,*)
 10   format(1x,10f12.5)
      endif
cccccccccccccccccc
C
      IF (IPRINT.GT.2) 
     *     WRITE(*,*) 'SUBROUTINE SO3SO3BARBUILD ENDS HERE'
C
      RETURN
C
      END







