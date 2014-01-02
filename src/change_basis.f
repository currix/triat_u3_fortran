      SUBROUTINE CHANGE_BASIS(N2, L, LD, AVECU2, W2MAT, AVECSO3, IPRINT)
C     
C     SUBROUTINE THAT PERFORMS A CHANGE OF BASIS FROM A SET OF VECTORS 
C     EXPRESSED IN THE U(2) BASIS (CYLINDRICAL OSCILLATOR) TO THE SO(3)
C     DISPLACED OSCILLATOR BASIS IN THE U(3) MODEL FOR BENDING VIBRATIONS
C
C     INPUT
C     N2     : U(3) IRREP LABEL (BENDING)
C     L      : VIBRATIONAL ANGULAR MOMENTUM LABEL
C     LD     : LEADING DIMENSION OF AVECU2, AVECSO3, AND W2MAT
C     AVECU2 : VECTORS IN THE U(2) BASIS
C     IPRINT : CONTROL VERBOSITY. DEBUGGING PURPOSES
C
C     OUTPUT
C     W2MAT   : ARRAY WITH THE SO(3) - U(2) TRANSF. BRACKETS
C     AVECSO3 : VECTORS TRANSFORMED TO THE SO(3) BASIS 
C     
C
C     by Currix TM
C
      IMPLICIT NONE
C     DEFINITION OF PARAMETERS
C     DIMENSION OF SCRATCH VECTOR
      INTEGER DIMSCR
      PARAMETER (DIMSCR = 18000) 
C
      INTEGER N2, L, LD, IPRINT
C
      DOUBLE PRECISION AVECU2(LD,*), AVECSO3(LD,*), W2MAT(LD,*)
C
C     SCRATCH VECTORS AND VARIABLES FOR DIAGONALIZATION
C     VECTORS IN COMMON
      DOUBLE PRECISION FV1
      COMMON/SCRATCH/ FV1(DIMSCR)
      INTEGER IERR
C
C     LOCAL VARIABLES
      INTEGER I, DIM, DIMSCALC, OMEGA
      DOUBLE PRECISION AVAL
ccc      integer j
C
C
      IF (IPRINT.GT.1) WRITE(*,*) 'SUBROUTINE CHANGE_BASIS STARTS HERE'
C
C     BUILD W^2 OPERATOR 
      CALL SO3CASBUILD(LD,W2MAT,N2,L)
C
C     BLOCK DIMENSION
      DIM = (N2-MOD(N2-L,2)-L)/2 + 1
C
C     OPTIMUM SCRATCH DIMENSION
      CALL DSYEV('V','U',DIM,W2MAT,LD,FV1(DIMSCR-DIM),FV1,-1,IERR)
      IF (FV1(1).LT.DIMSCR) THEN
         DIMSCALC = FV1(1)
      ELSE
         DIMSCALC = DIMSCR - DIM - 1
      ENDIF
C         
      IF (IPRINT.GT.2) WRITE(*,*) 'DIMSCALC', DIMSCALC
C
C     DIAGONALIZE W^2 OPERATOR 
      CALL DSYEV('V','U',DIM,W2MAT,LD,FV1(DIMSCR-DIM),FV1,DIMSCALC,IERR)
C
      IF (IPRINT.GT.1) THEN
         WRITE(*,*)
         WRITE(*,*) "CHECKING W^2 EIGENVALUES"
         WRITE(*,*)
         WRITE(*,*) "  EIGENVALUE                        W           V"
         DO I = 1, DIM
            AVAL = FV1(DIMSCR-DIM+I-1)
            OMEGA = INT(0.5D0*(SQRT(4*AVAL+1)-1)+0.5D0)
            WRITE(*,*) AVAL, OMEGA, (N2-OMEGA)/2
         ENDDO
      ENDIF
C     
cccc      do i = 1, dim
cccc         write(*,*) (w2mat(i,j),j=1,dim)
cccc      enddo
cccc      write(*,*)
C
C     CHECK PHASE
      CALL CHECK_PHASE(N2, L, LD, W2MAT, IPRINT)
c      do j = 1, dim
c         write(*,*) (avecu2(i,j),i=1,dim)
c      enddo
C     TRANSFORMING TO THE DISPLACED OSCILLATOR BASIS
cccccccccccccccccccccccccccccccccccccccccccccccccc
c      CALL DGEMM('T','N',DIM,DIM,DIM,
c     *     1.0D0,W2MAT,LD,W2MAT,LD,0.0D0,AVECSO3,LD)
cc      do i = 1, dim
cc         write(*,*) (w2mat(i,j),j=1,dim)
cc      enddo
cc      write(*,*)
cc      do i = 1, dim
cc         write(*,*) (w2mat(j,i),j=1,dim)
cc      enddo
cc      write(*,*)
cc      do i = 1, dim
cc         write(*,*) (avecu2(i,j),j=1,dim)
cc      enddo
C
C     CHANGE OF BASIS
      CALL DGEMM('T','N',DIM,DIM,DIM,
     *     1.0D0,W2MAT,LD,AVECU2,LD,0.0D0,AVECSO3,LD)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do j = 1, dim
c         write(*,*) (w2mat(i,j),i=1,dim)
c      enddo
c      do j = 1, dim
c         write(*,*) (avecSO3(i,j),i=1,dim)
c      enddo
C     
      IF (IPRINT.GT.1) WRITE(*,*) 'SUBROUTINE CHANGE_BASIS ENDS HERE'
C
      RETURN  
C     
      END
