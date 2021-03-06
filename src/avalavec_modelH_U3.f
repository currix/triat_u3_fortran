      PROGRAM AVALAVEC_U3
C
C     PROGRAM TO COMPUTE THE ENERGIES AND EIGENVALUES FOR THE U(3) MODEL HAMILTONIAN
C
C                    by Currix TM
C     
      IMPLICIT NONE
C
C     DEFINITION OF PARAMETERS
C
      INTEGER N2MAX, DMAX, NPMAX, DIMSCR
C     
C     MAXIMUM NUMBER OF U(3) BOSONS (N_2)  
      PARAMETER (N2MAX = 6000) 
C     DIMENSION OF SCRATCH VECTOR
      PARAMETER (DIMSCR = 18000) 
C     MAXIMUM DIMENSIONS OF MATRICES
      PARAMETER (DMAX = 3001)
C     MAXIMUM NUMBER OF PARAMETERS IN THE HAMILTONIAN
      PARAMETER (NPMAX = 2)
C     
C     DEFINITION OF VARIABLES
C
C     REPRESENTATION OF U(3)
      INTEGER N2
C
C     VIBRATIONAL ANGULAR MOMENTUM
      INTEGER LVAL
C     OPTIONS:
C     IOPTS = 0 OUTPUTS ONLY GROUND STATE ENERGY
C     IOPTS = 1 OUTPUTS ALL ENERGIES
C     IOPTS >= 2 OUTPUTS ALL ENERGIES WITH G.S. ENERGY = 0
      INTEGER IOPTS
CC     LOGICAL EMINL
C
C     HAMILTONIAN MATRIX
      DOUBLE PRECISION HAM(DMAX,DMAX)
C
C     COMPUTED ENERGIES MATRIX
      DOUBLE PRECISION EIGEN(DMAX)
C     COMPUTED W2 (O(3) CASIMIR) BLOCKS
      DOUBLE PRECISION W2MAT(DMAX,DMAX)
C
C     SCRATCH VECTORS AND VARIABLES FOR DIAGONALIZATION
C     VECTORS IN COMMON
      DOUBLE PRECISION FV1
      COMMON/SCRATCH/ FV1(DIMSCR)
      INTEGER IERR
C
C     HAMILTONIAN PARAMETERS
C
C     ALGEBRAIC MODEL
C     PARAMETERS 
C
      DOUBLE PRECISION SCALE, XI
      DOUBLE PRECISION EPSILON, A
C
C     
C     CONTROL OF OUTPUT SAVED
C     EIGENVALUES AND EIGENSTATES SAVED IF TRUE
      LOGICAL LSAVE 
C
C     COMMON BLOCK FOR HAMILTONIAN PARAMETERS
C
      DOUBLE PRECISION HPAR
C
C    
      COMMON/HAMPAR/ HPAR(NPMAX)
C     
C     CONTROL OF OUTPUT DISPLAYED 
      INTEGER IPRINT
C
C     CONTROL DISPLAYED OUTPUT 
      COMMON/GRAF/ IPRINT
C
C     LOCAL VARIABLES                                                      
      INTEGER I,L,DIM
      INTEGER J
      DOUBLE PRECISION  EMIN
C
C     READ VARIABLES
      READ(5,*) N2 
      READ(5,*) LVAL
      READ(5,*) IOPTS
      READ(5,*) SCALE
      READ(5,*) XI
      READ(5,*) LSAVE
cc
c     Debugging variable
c     Minimum output
      IPRINT = 0
c     Saves Hamiltonian matrix
cc      iprint = -1
cc      write(*,*) iprint
cc
C
      IF (IPRINT.GT.2) WRITE(*,*) 'MAIN PROGRAM'
C
C
C     TESTS 
C    
      IF (LVAL.GT.N2) STOP 'LMX > N2, SAYONARA '
      IF (N2.GT.N2MAX) STOP 'N2 > N2MAX, SAYONARA '
C
C
      EPSILON = SCALE*(1.0D0-XI)
      A = SCALE*XI/DFLOAT(N2-1)
C     
      HPAR(1) = EPSILON
      HPAR(2) = A
C     
C     INITIALIZE EMIN
      EMIN = 0.0D0
C     
C     MAIN LOOP
C     
      L = LVAL
C     
      IF (IPRINT.GT.3) THEN
         WRITE(*,*)
         WRITE(*,*) 'L',L
         WRITE(*,*)
      ENDIF
C     BUILD HAMILTONIAN MATRIX
C     ALGEBRAIC HAMILTONIAN
      CALL HBLDU3MOD(N2,L,DMAX,HAM,W2MAT)  
C     
C
C     DIAGONALIZE HAMILTONIAN MATRIX
C
      DIM = (N2-MOD(N2-L,2)-L)/2+1
C
      CALL DSYEV('V', 'U', DIM, HAM, DMAX, EIGEN,  FV1,  3*DIM,
     *     IERR)
C     
C     REFER ENERGIES TO G.S.
C     
      EMIN = 0.0D0
C     
      WRITE(*,*) DIM
      IF (IOPTS.EQ.0) THEN
         WRITE(*,*) EIGEN(1)
         DO J = 1, DIM
            WRITE(*,*) HAM(J,1), N2 - (2*J-2+MOD(N2-L,2)), L
         ENDDO
      ELSE         
         DO I = 1, DIM
            WRITE(*,*) EIGEN(I) - EMIN
            DO J = 1, DIM
               WRITE(*,*) HAM(J,I), N2 - (2*J-2+MOD(N2-L,2)), L
            ENDDO
         ENDDO
      ENDIF
C     
C     SAVE EIGENVALUES AND EIGENSTATES
      IF (LSAVE) THEN
         OPEN(UNIT=100, 
     *        FILE='u3_model_Ham_results.dat', STATUS='UNKNOWN')
         WRITE(100,*) DIM
         IF (IOPTS.EQ.0) THEN
            WRITE(100,*) EIGEN(1)
            DO J = 1, DIM
               WRITE(100,*) HAM(J,1), N2 - (2*J-2+MOD(N2-L,2)), L
            ENDDO
         ELSE         
            DO I = 1, DIM
               WRITE(100,*) EIGEN(I) - EMIN
               DO J = 1, DIM
                  WRITE(100,*) HAM(J,I), N2 - (2*J-2+MOD(N2-L,2)), L
               ENDDO
            ENDDO
         ENDIF
      ENDIF
C     
      END
CC     input template
CC  10 	# N
CC  1  	# l
CC  1  	# iopts
CC  10.0d0 	# scale
CC  0.0d0	# xi
CC  .TRUE. # lsave
      



