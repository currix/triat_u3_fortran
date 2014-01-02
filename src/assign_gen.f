      SUBROUTINE ASSGNU3(BENT,N2,L,LD,W2,W4,W2W2B,HAM,EIGEN,LDEXP,VEXPAS
     *     ,NDAT,ASSIGNALL,BLAS)
C
C   $Id: assign_gen.f,v 1.6 2011/05/04 17:30:13 curro Exp $
C
C     SUBROUTINE THAT LOOKS FOR THE BEST ASSIGNATIONS  
C     FOR THE CALCULATED LEVELS TO
C     LOCAL BASIS STATES IN THE U(3) MODEL.
C
C     IF LINEAR IT ASSIGNS TO A n,l BASIS IF BENT TO A v,K BASIS
C
C
C     IN CASE OF HIGH MIXING THE ASSIGNMENT COULD BE AMBIGUOUS AND 
C     THEN IT IS FOLLOWED BY SCALING DOWN THE MIXING PARAMETERS UNTIL 
C     AN UNAMBIGUOUS ASSIGNMENT CAN BE FOUND AND THEN SCALE THEM UP 
C     AGAIN TO THE ORIGINAL VALUE, PROJECTING IN EACH STEP ONTO THE 
C     PREVIOUS ONES (PROCEDURE USED TO RESOLVE  STARK MIXED ROTATIONAL
C     LEVELS, PROGRAMMED BY THOMAS MULLER). 
C     
C     INPUT
C     BENT     : .T. BENT MOLECULE, .F. LINEAR MOLECULE
C     N2       : U(3) IRREP LABEL (BENDING)
C     L        : VIBRATIONAL ANGULAR MOMENTUM LABEL
C     LD       : LEADING DIMENSION OF MATRICES HAM,EIGEN,BLAS AND W2 
C     W2       : ARRAY WITH THE SO(3) CASIMIR BLOCK 
C     W4       : ARRAY WITH THE SO(3) CASIMIR BLOCK 
C     W2W2B    : ARRAY WITH THE SO(3) W2·W2B + W2B·W2 CASIMIR BLOCK 
C     HAM      : HAMILTONIAN MATRIX (EIGENVECTOR MATRIX IN INPUT)
C     EIGEN    : EIGENVALUES VECTOR
C     LDEXP    : LEADING DIMENSION OF MATRIX VEXPAS
C     VEXPAS   : EXPERIMENTAL ASSIGNMENTS FOR POLYAD (L) FORMAT:(V l) OR (n l)
C     NDAT     : NUMBER OF EXPERIMENTAL DATA FOR POLYAD (L)
C     ASSIGNALL: IF .F. ASSIGN UNAMBIGUOUSLY ONLY DATA WITH EXPERIMENTAL INFO  
C
C     OUTPUT
C     BLAS  :  LOCAL BASIS (BENT OR LINEAR) VECTOR WITH ASSIGNMENTS
C
C
C     by Currix TM
C     
      IMPLICIT NONE
C     
C     
C     DEFINITION OF PARAMETERS
C
      INTEGER NPMAX
C     
C     MAXIMUM NUMBER OF PARAMETERS IN THE HAMILTONIAN
      PARAMETER (NPMAX = 15)
C     
      INTEGER DMAX
C     
C     MAXIMUM DIMENSIONS OF MATRICES (LINEAR CASE) 
      PARAMETER (DMAX = 2001)
C     
      INTEGER DIMSCR
C     
C     DIMENSIONS OF SCRATCH VECTOR
      PARAMETER (DIMSCR = 12000)
C
C     DEFINITION OF VARIABLES
C     
      LOGICAL BENT, ASSIGNALL
      INTEGER N2
      INTEGER L, LD, LDEXP, VEXPAS(LDEXP,*), NDAT, BLAS(LD)
      DOUBLE PRECISION HAM(LD,*), EIGEN(*), 
     *     W2(LD,*), W4(LD,*), W2W2B(LD,*)
C
C
C     COMMON BLOCK FOR HAMILTONIAN PARAMETERS
C
      DOUBLE PRECISION HPAR
C
C
      COMMON/HAMPAR/ HPAR(NPMAX)
C
C     SCRATCH VECTORS AND VARIABLES FOR DIAGONALIZATION
C     VECTORS IN COMMON(**)
      DOUBLE PRECISION FV1
      COMMON/SCRATCH/ FV1(DIMSCR)
      INTEGER IERR
C
      INTEGER IPRINT
C
C     CONTROL OUTPUT DISPLAYED
      COMMON/GRAF/ IPRINT
C
C     BLAS EXTERNAL FUNCTION
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C
C     TEMPORAL VARIABLES                                     
      INTEGER DIM
      LOGICAL FLAG
      DOUBLE PRECISION HPART(NPMAX)
      INTEGER MIX, PTEMP(DMAX)
      DOUBLE PRECISION FAC, MIXSTEP
      DOUBLE PRECISION HAM2(DMAX,DMAX), VMAX, VT
      DOUBLE PRECISION TBRACK(DMAX,DMAX)
      INTEGER I, J, J2, JMAX, ICOUNT
C
C
      IF (IPRINT.GT.2) WRITE(*,*) 'SUBROUTINE ASSIGN STARTS HERE'
C
C     POLYAD DIMENSION
      DIM = (N2-MOD(N2-L,2)-L)/2 + 1
C
C      INITIALIZE COUNTER
C
      ICOUNT = 0
C
      IF (BENT) THEN
C
C     BENT CASE
C
C     PROJECTION OF THE EIGENVECTORS TO THE BENT BASIS
C     (I) DIAGONALIZATION OF THE PAIRING OPERATOR TO GET TRANSF. BRACK.
C     
         DO I = 1, DIM 
            HAM2(I,1) = 0.0D0
         ENDDO
C
         CALL SO3CASBUILD(LD,W2,N2,L)
C
C     CHANGE SIGN TO W2 TO GET THE CORRECT GROUND STATE WF
         DO I = 1, DIM 
            DO J = 1, DIM 
               W2(J,I) = -W2(J,I)
            ENDDO 
         ENDDO
C
C     
         CALL DSYEV('V','U',DIM,W2,LD,HAM2(1,1),FV1,DIMSCR,IERR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c         do i = 1, dim 
c            write(*,*) i, HAM2(i,1),( 
c     *           2.d0*dfloat(N2)+1.0d0-
c     *           dsqrt((2.0d0*dfloat(N2)+1.0d0)**2-
c     *           4.0d0*HAM2(i,1)))/4.0d0,i-1
c         enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         CALL DLACPY('F', DIM, DIM, W2, DMAX, TBRACK, DMAX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c         do i = 1, dim 
c            write(*,*) i, (TBRACK(j,i),j=1,dim)
c         enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     (II) MATRIX PRODUCT
         CALL DGEMM('T','N',DIM,DIM,DIM,
     *        1.0D0,TBRACK,LD,HAM,LD,0.0D0,W2,DMAX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c        do i = 1, dim 
c           write(*,*) i, (W2(j,i),j=1,dim)
c        enddo
c        WRITE(*,*) 'HAM'
c        do i = 1, dim 
c           write(*,*) i, (HAM(j,i),j=1,dim)
c        enddo
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF (ASSIGNALL) THEN
            CALL MAXCU3(BENT,N2,L,LD,W2,DIM,BLAS,FLAG)
         ELSE
            CALL MAXCU3EXP(BENT,LDEXP,VEXPAS,NDAT,N2,L,LD,W2,DIM,BLAS
     *           ,FLAG)
         ENDIF
C
         IF (IPRINT.GE.3) THEN
            WRITE(*,*) 'MAXC', FLAG
            WRITE(*,*) (BLAS(I),I=1,DIM)
         ENDIF
C     
C     IF FLAG = .T. PROBLEMS WITH ASSIGNMENT
C     
         IF (FLAG) THEN
C     SAVE COPY OF INITIAL PARAMETERS
            DO I = 1, NPMAX
               HPART(I) = HPAR(I)
            ENDDO
C     
C     MIX: NUMBER OF ITERATIONS UP AND DOWN
C     
            MIX = 0
C     
C     SCALE DOWN
C     
            FAC = 1.0D0
C     
 11         CONTINUE
C     
            FAC = FAC/2.0D0
C     
            MIX = MIX + 1
C     
C     
C     RECOVER COPY OF INITIAL PARAMETERS
            DO I = 1, NPMAX
               HPAR(I) = HPART(I)
            ENDDO
C     
            IF (IPRINT.GE.2) WRITE(*,*)'SCALING DOWN MIXING BY ', FAC
C     
            ICOUNT = ICOUNT + 1
C     
            CALL SCLHAM(FAC,BENT,N2,L,LD,HAM,W2,W4,W2W2B)
C     
C     DIAGONALIZE HAMILTONIAN
C     
            DIM = (N2-MOD(N2-L,2)-L)/2 + 1            
C     
            DO I = 1, DIM 
               EIGEN(I) = 0.0D0
            ENDDO
C     
            CALL DSYEV('V','U',DIM,HAM,DMAX,EIGEN,FV1,DIMSCR,IERR)
C     
C     LOOK FOR MAXIMUM COMPONENT
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     (II) MATRIX PRODUCT
            CALL DGEMM('T','N',DIM,DIM,DIM,
     *           1.0D0,TBRACK,LD,HAM,LD,0.0D0,W2,DMAX)
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     do i = 1, dim 
c     write(*,*) i, (W2(j,i),j=1,dim)
c     enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            IF (ASSIGNALL) THEN
               CALL MAXCU3(BENT,N2,L,LD,W2,DIM,BLAS,FLAG)
            ELSE
               CALL MAXCU3EXP(BENT,LDEXP,VEXPAS,NDAT,N2,L,LD,W2,DIM,BLAS
     *              ,FLAG)
            ENDIF
C     
            IF (IPRINT.GE.4) THEN
               write(*,*) 'maxc', flag
               WRITE(*,*) (BLAS(I),I=1,DIM)
            ENDIF
C     
            IF (FLAG) THEN 
               GOTO 11
            ELSE
               GOTO 21
            ENDIF
C
         ELSE
C
            IF (IPRINT.GE.1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
C            
            RETURN
C
         ENDIF
C     
 21      CONTINUE
C     
C     SCALE UP
C     
         MIXSTEP = 1.0D0/(1.5D0*2.0D0**MIX)
         FAC = 1.0D0/(2.0D0**MIX)      
C     
C     
         CALL DLACPY('F', DIM, DIM, W2, DMAX, HAM, DMAX)         
C     
C     
 26      CONTINUE
C     
C     
C     RECOVER COPY OF INITIAL PARAMETERS
         DO I = 1, NPMAX
            HPAR(I) = HPART(I)
         ENDDO
C     
C     
         MIXSTEP = 1.5D0*MIXSTEP
         FAC = FAC + MIXSTEP
C     
CVCVCVCVCVCVCV         CALL DLACPY('F', DIM, DIM, W2, DMAX, HAM, DMAX)         
C     
         IF (FAC.GT.1.0D0) FAC = 1.0D0
C     
         IF (IPRINT.GE.2) WRITE(*,*) 'SCALING MIXING UP BY ',FAC 
C     
         ICOUNT = ICOUNT + 1
C     
         CALL SCLHAM(FAC,BENT,N2,L,LD,HAM2,W2,W4,W2W2B)
C     
C     DIAGONALIZE HAMILTONIAN
C     
         DIM = (N2-MOD(N2-L,2)-L)/2 + 1            
C     
         DO I = 1, DIM 
            EIGEN(I) = 0.0D0
         ENDDO
C     
         CALL DSYEV('V','U',DIM,HAM2,DMAX,EIGEN,FV1,250+DMAX,IERR)
C
C     (II) MATRIX PRODUCT
         CALL DGEMM('T','N',DIM,DIM,DIM,
     *        1.0D0,TBRACK,LD,HAM2,LD,0.0D0,W2,DMAX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c         do i = 1, dim 
c            write(*,*) i, (W2(j,i),j=1,dim)
c         enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     PROJECTING ONTO FORMER EIGENVECTORS
C     THE EIGENVECTORS ARE PROJECTED ONTO THE PREVIOUS ONES AND 
C     IF EITHER SOME OF THE PROJECTIONS SQUARED IS LESS THAN 0.5
C     OR TWO HAVE THE SAME MAXIMAL PROJECTIONS THE PARAMETERS ARE 
C     RESCALED DOWN AGAIN.
C     
         IF (ASSIGNALL) THEN
C     SCALAR PRODUCT
            DO I = 1, DIM 
               VMAX = 0.0D0
               PTEMP(I) = 0
               DO J = 1, DIM
                  VT = DDOT(DIM,HAM(1,J),1,W2(1,I),1)
C     
C     MAXIMAL COMPONENT
                  VT = VT * VT
                  IF (VT.GT.VMAX) THEN
                     VMAX = VT
                     JMAX = BLAS(J)
                  ENDIF
               ENDDO
               PTEMP(I) = JMAX 
C     
C     CHECK FOR AMBIGUITIES 
               IF (VMAX.LT.0.5D0) THEN
                  FAC = FAC - MIXSTEP
                  MIXSTEP = MIXSTEP/3.0D0
                  GOTO 26
               ELSE
                  DO J2 = 1, I-1
                     IF (PTEMP(J2).EQ.JMAX) THEN
                        FAC = FAC - MIXSTEP
                        MIXSTEP = MIXSTEP/3.0D0
                        GOTO 26
                     ENDIF
                  ENDDO
               ENDIF
C
            ENDDO
C
         ELSE
C
C     COMPARING ONLY TO STATES WITH BLAS(I).NE.-1
C     SCALAR PRODUCT
            DO I = 1, DIM
               PTEMP(I) = -1
            ENDDO
            DO I = 1, DIM
               IF (BLAS(I).NE.-1) THEN
                  VMAX = 0.0D0
C     
                  DO J = 1, DIM
C     
                     VT = DDOT(DIM,HAM(1,I),1,W2(1,J),1)
C     
C     MAXIMAL COMPONENT
                     VT = VT * VT
                     IF (VT.GT.VMAX) THEN
                        VMAX = VT
                        JMAX = J
                     ENDIF
                  ENDDO
C     
                  PTEMP(JMAX) = BLAS(I)
C
                  IF (IPRINT.GT.2) WRITE(*,*) I, BLAS(I), JMAX, VMAX
C     
C     CHECK FOR AMBIGUITIES 
                  IF (VMAX.LT.0.5D0) THEN
                     FAC = FAC - MIXSTEP
                     MIXSTEP = MIXSTEP/3.0D0
                     GOTO 26
                  ENDIF
                  DO J2 = 1, DIM
                     IF (PTEMP(J2).EQ.BLAS(I).AND.J2.NE.JMAX) THEN
                        FAC = FAC - MIXSTEP
                        MIXSTEP = MIXSTEP/3.0D0
                        GOTO 26
                     ENDIF
                  ENDDO
C
               ENDIF
C
            ENDDO
C     
         ENDIF
C
         CALL DLACPY('F', DIM, DIM, W2, DMAX, HAM, DMAX)         
C     
         DO I = 1, DIM
            BLAS(I) = PTEMP(I)
         ENDDO
C
         IF (FAC.EQ.1.0D0) THEN 
            CALL SCLHAM(FAC,BENT,N2,L,LD,HAM,W2,W4,W2W2B)
C     
            DIM = (N2-MOD(N2-L,2)-L)/2 + 1            
C     
            DO I = 1, DIM 
               EIGEN(I) = 0.0D0
            ENDDO
C     
            CALL DSYEV('V','U',DIM,HAM,DMAX,EIGEN,FV1,250+DMAX,IERR)
C     
            IF (IPRINT.GE.1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
C     
            RETURN
C     
         ELSE         
C
            GOTO 26
C
         ENDIF
C
      ELSE
C     
C     LINEAR CASE
C     
C     LOOK FOR MAXIMUM COMPONENT
C     
         IF (ASSIGNALL) THEN
            CALL MAXCU3(BENT,N2,L,LD,HAM,DIM,BLAS,FLAG)
         ELSE
            CALL MAXCU3EXP(BENT,LDEXP,VEXPAS,NDAT,N2,L,LD,HAM,DIM,BLAS
     *           ,FLAG)
         ENDIF
C     
         IF (IPRINT.GE.4) THEN
            WRITE(*,*) 'MAXC', FLAG
            WRITE(*,*) (BLAS(I),I=1,DIM)
         ENDIF
C     
C     IF FLAG = .T. PROBLEMS WITH ASSIGNMENT
C     
         IF (FLAG) THEN
C     SAVE COPY OF INITIAL PARAMETERS
            DO I = 1, NPMAX
               HPART(I) = HPAR(I)
            ENDDO
C     
C     MIX: NUMBER OF ITERATIONS UP AND DOWN
C     
            MIX = 0
C     
C     SCALE DOWN
C     
            FAC = 1.0D0
C     
 10         CONTINUE
C     
            FAC = FAC/2.0D0
C     
            MIX = MIX + 1
C     
C     
C     RECOVER COPY OF INITIAL PARAMETERS
            DO I = 1, NPMAX
               HPAR(I) = HPART(I)
            ENDDO
C     
            IF (IPRINT.GE.2) WRITE(*,*)'SCALING DOWN MIXING BY ', FAC
C     
            ICOUNT = ICOUNT + 1
C     
            CALL SCLHAM(FAC,BENT,N2,L,LD,HAM,W2,W4,W2W2B)
C     
C     DIAGONALIZE HAMILTONIAN
C     
            DIM = (N2-MOD(N2-L,2)-L)/2 + 1            
C     
            DO I = 1, DIM 
               EIGEN(I) = 0.0D0
            ENDDO
C     
            CALL DSYEV('V','U',DIM,HAM,DMAX,EIGEN,FV1,DIMSCR,IERR)
C     
C     LOOK FOR MAXIMUM COMPONENT
C     
            IF (ASSIGNALL) THEN
               CALL MAXCU3(BENT,N2,L,LD,HAM,DIM,BLAS,FLAG)
            ELSE
               CALL MAXCU3EXP(BENT,LDEXP,VEXPAS,NDAT,N2,L,LD,HAM,DIM,
     *              BLAS,FLAG)
            ENDIF
C     
            IF (FLAG) THEN 
               GOTO 10
            ELSE
               GOTO 20
            ENDIF
         ELSE
C
            IF (IPRINT.GE.1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
C            
            RETURN
C
         ENDIF
C     
 20      CONTINUE
C     
C     SCALE UP
C     
         MIXSTEP = 1.0D0/(1.5D0*2.0D0**MIX)
         FAC = 1.0D0/(2.0D0**MIX)      
C     
C     
 25      CONTINUE
C     
C     
C     RECOVER COPY OF INITIAL PARAMETERS
         DO I = 1, NPMAX
            HPAR(I) = HPART(I)
         ENDDO
C     
C     
         MIXSTEP = 1.5D0*MIXSTEP
         FAC = FAC + MIXSTEP
C     
         
C
         IF (FAC.GT.1.0D0) FAC = 1.0D0
C     
         IF (IPRINT.GE.2) WRITE(*,*) 'SCALING MIXING UP BY ',FAC 
C     
         ICOUNT = ICOUNT + 1
C          
         CALL SCLHAM(FAC,BENT,N2,L,LD,HAM2,W2,W4,W2W2B)
C     
C     DIAGONALIZE HAMILTONIAN
C     
         DIM = (N2-MOD(N2-L,2)-L)/2 + 1            
C     
         DO I = 1, DIM 
            EIGEN(I) = 0.0D0
         ENDDO
C     
         CALL DSYEV('V','U',DIM,HAM2,DMAX,EIGEN,FV1,DIMSCR,IERR)
C     
C     PROJECTING ONTO FORMER EIGENVECTORS
C     THE EIGENVECTORS ARE PROJECTED ONTO THE PREVIOUS ONES AND 
C     IF EITHER SOME OF THE PROJECTIONS SQUARED IS LESS THAN 0.5
C     OR TWO HAVE THE SAME MAXIMAL PROJECTIONS THE PARAMETERS ARE 
C     RESCALED DOWN AGAIN.
C     
         IF (ASSIGNALL) THEN
C     SCALAR PRODUCT
            DO I = 1, DIM 
               VMAX = 0.0D0
               PTEMP(I) = 0
               DO J = 1, DIM
C     
                  VT = DDOT(DIM,HAM(1,J),1,HAM2(1,I),1)
C     
C     MAXIMAL COMPONENT
C     
                  VT = VT * VT
                  IF (VT.GT.VMAX) THEN
                     VMAX = VT
                     JMAX = BLAS(J)
                  ENDIF
               ENDDO
               PTEMP(I) = JMAX
C     
C     CHECK FOR AMBIGUITIES
               IF (VMAX.LT.0.5D0) THEN
                  FAC = FAC - MIXSTEP
                  MIXSTEP = MIXSTEP/3.0D0
                  GOTO 25
               ELSE
                  DO J2 = 1, I - 1
                     IF (PTEMP(J2).EQ.JMAX) THEN
                        FAC = FAC - MIXSTEP
                        MIXSTEP = MIXSTEP/3.0D0
                        GOTO 25
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
C
         ELSE
C
C     COMPARING ONLY TO STATES WITH BLAS(I).NE.-1
C     SCALAR PRODUCT
            DO I = 1, DIM
               PTEMP(I) = -1
            ENDDO
C
            DO I = 1, DIM
               IF (BLAS(I).NE.-1) THEN
                  VMAX = 0.0D0
C     
                  DO J = 1, DIM
C     
                     VT = DDOT(DIM,HAM(1,I),1,HAM2(1,J),1)
C     
C     MAXIMAL COMPONENT
                     VT = VT * VT
                     IF (VT.GT.VMAX) THEN
                        VMAX = VT
                        JMAX = J
                     ENDIF
                  ENDDO
C     
                  PTEMP(JMAX) = BLAS(I)
C
                  IF (IPRINT.GT.2) WRITE(*,*) I, BLAS(I), JMAX, VMAX
C     
C     CHECK FOR AMBIGUITIES 
                  IF (VMAX.LT.0.5D0) THEN
                     FAC = FAC - MIXSTEP
                     MIXSTEP = MIXSTEP/3.0D0
                     GOTO 25
                  ENDIF
                  DO J2 = 1, DIM
                     IF (PTEMP(J2).EQ.BLAS(I).AND.J2.NE.JMAX) THEN
                        FAC = FAC - MIXSTEP
                        MIXSTEP = MIXSTEP/3.0D0
                        GOTO 25
                     ENDIF
                  ENDDO
C
               ENDIF
C
            ENDDO
C
         ENDIF
C
         CALL DLACPY('F', DIM, DIM, HAM2, DMAX, HAM, DMAX)
C
CCCC         DO I = 1, DIM
CCC            DO J = 1, DIM
CCC               HAM(J,I) = HAM2(J,I)
CCC            ENDDO
CCCC         ENDDO
         DO I = 1, DIM
            BLAS(I) = PTEMP(I)
         ENDDO
C
         IF (FAC.EQ.1.0D0) THEN 
C     
            IF (IPRINT.GE.1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
C
            IF (IPRINT.GT.2) WRITE(*,*) 'SUBROUTINE ASSIGN ENDS HERE'
C
C     
            RETURN
         ELSE         
C
            GOTO 25
C
         ENDIF
C     
      ENDIF
C
      END
      
