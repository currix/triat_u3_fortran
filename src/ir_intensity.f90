PROGRAM IR_INTENSITY
  !
  ! Program to compute IR intensities between U(3) model eigenstates.
  ! Relevant formulae:
  ! PRA 77 032115 Eqs. (45,46,47,48)
  ! Always computing |<STATE_I^(l+1) | T_+ | STATE_J^(l)>|^2
  !
  ! T_+ = [t/sqrt(N)] D_+
  !
  ! by Currix TM.
  !
  IMPLICIT NONE
  !
  !Variable Definition
  INTEGER, PARAMETER :: Long=SELECTED_REAL_KIND(18,310) ! High precision real numbers
  !
  CHARACTER(LEN=64) FILE_1, FILE_2
  !
  INTEGER :: N1 = 0, N2 = 0 ! N values
  INTEGER :: L1 = 0, L2 = 0 ! Angular Momentum 
  INTEGER :: DIM_1 = 0, DIM_2 = 0 ! Eigenstate dimension
  !
  REAL (Long), ALLOCATABLE , DIMENSION(:,:) :: Tplus ! Local basis T+ matrix 
  REAL (Long), ALLOCATABLE , DIMENSION(:) :: STATE_1, STATE_2
  REAL (Long), ALLOCATABLE , DIMENSION(:) :: STATE_TEMP
  !
  !
  INTEGER :: DIM_TEMP, L_TEMP
  !
  REAL (Long) :: IRINTEN ! FINAL RESULT
  !
  INTERFACE
     SUBROUTINE READ_STATE(FILE,N,L,DIM,STATE)
       IMPLICIT NONE
       !
       INTEGER, PARAMETER :: Long=SELECTED_REAL_KIND(18,310) ! High precision real numbers
       !
       CHARACTER(LEN=64) FILE
       INTEGER, INTENT(OUT) :: N, L, DIM
       REAL (Long), ALLOCATABLE, DIMENSION (:), INTENT(OUT) :: STATE
     END SUBROUTINE READ_STATE
     SUBROUTINE TPMATRIX(TP,DIM,Nn,L)
       IMPLICIT NONE
       !
       INTEGER, PARAMETER :: Long=SELECTED_REAL_KIND(18,310) ! High precision real numbers
       !
       INTEGER, INTENT(IN) :: Nn, L, DIM
       REAL (Long), DIMENSION (:,:), INTENT(OUT) :: TP
     END SUBROUTINE TPMATRIX
  END INTERFACE
  !
  ! READ FILENAMES
  READ*, FILE_1, FILE_2
  !
  ! READ STATES
  CALL READ_STATE(FILE_1,N1,L1,DIM_1,STATE_1)
  !!PRINT *, FILE_1, N1, L1, DIM_1
  !!PRINT *, STATE_1
  CALL READ_STATE(FILE_2,N2,L2,DIM_2,STATE_2)
  !!PRINT *, FILE_2, N2, L2, DIM_2
  !!PRINT *, STATE_2
  !
  ! Tests
  IF (N1/=N2) STOP "N1 /= N2. PROGRAM ABORTED"
  !IF (MOD(N1,2)/=0) STOP "N1 = N2 SHOULD BE EVEN. PROGRAM ABORTED"
  IF (ABS(L1-L2)/=1) STOP "IR SELECTION RULE VIOLATED. PROGRAM ABORTED"
  !
  ! SWITCH VECTORS IF NECESSARY
  IF (L1 < L2) THEN
     ALLOCATE(STATE_TEMP(1:DIM_1))
     DIM_TEMP = DIM_1
     L_TEMP = L1
     STATE_TEMP=STATE_1
     DEALLOCATE(STATE_1)
     ALLOCATE(STATE_1(1:DIM_2))
     STATE_1 = STATE_2
     DIM_1 = DIM_2
     L1 = L2
     DEALLOCATE(STATE_2)
     DIM_2 = DIM_TEMP
     L2 = L_TEMP
     ALLOCATE(STATE_2(1:DIM_2))
     STATE_2 = STATE_TEMP
     DEALLOCATE(STATE_TEMP)
  ENDIF
  !
  ! BUILD T+ MATRIX
  ALLOCATE(Tplus(1:DIM_1,1:DIM_2))
  CALL TPMATRIX(Tplus,DIM_2,N2,L2)
  !
  !
  IRINTEN = DOT_PRODUCT(STATE_1,MATMUL(Tplus,STATE_2))
  !
  PRINT*, IRINTEN**2
  !
END PROGRAM IR_INTENSITY
!
!
!
SUBROUTINE READ_STATE(FILENAME,N,L,DIM,STATE)
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: Long=SELECTED_REAL_KIND(18,310) ! High precision real numbers
  !
  CHARACTER(LEN=64) FILENAME
  INTEGER, INTENT(OUT) :: N, L, DIM
  REAL (Long), ALLOCATABLE, DIMENSION (:), INTENT(OUT) :: STATE
  !
  INTEGER I
  !
  ! OPEN FILE
  OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME)
  !
  READ (10,*) N, L, DIM
  !
  ALLOCATE(STATE(1:DIM))
  DO I = 1, DIM
     READ(10,*) STATE(I) ! Cannot do in matrix form due to the form of the state output file
  ENDDO
! 
END SUBROUTINE READ_STATE
!
SUBROUTINE TPMATRIX(TP,DIM,Nn,L)
  !
  ! Subroutine to build the T_+ = 
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: Long=SELECTED_REAL_KIND(18,310) ! High precision real numbers
  !
  INTEGER, INTENT(IN) :: Nn, L, DIM
  REAL (Long), DIMENSION (:,:), INTENT(OUT) :: TP
  !
  INTEGER :: n ! U(2) label
  INTEGER :: COLIN = 0! COLUMN INDEX
  !
  ! BUILDING MATRIX
  !
  TP = 0.0_Long
  !
  IF (MOD(Nn-L,2) == 0) THEN
     ! N and L EVEN OR N and L ODD
     COLIN = 1
     !
     ! FIRST COLUMN
     n = Nn
     TP(COLIN,COLIN) = -SQRT(1.0_Long*(n-L)*(1.0_Long-(n-1.0_Long)/Nn))    ! <N-1^(l+1)| T_+ | N^l>
     ! INTERMEDIATE COLUMNS
     DO n = Nn-2, L+2, -2 ! notice reversed order
        COLIN = COLIN + 1
        TP(COLIN,COLIN) = -SQRT(1.0_Long*(n-L)*(1.0_Long-(n-1.0_Long)/Nn))    ! <n-1^(l+1)| T_+ | n^l>
        TP(COLIN-1,COLIN) = SQRT(1.0_Long*(n+L+2)*(1.0_Long-(1.0_Long*n)/Nn)) ! <n+1^(l+1)| T_+ | n^l>
     ENDDO
     !
     ! LAST COLUMN
     COLIN = DIM
     n = L
     TP(COLIN-1,COLIN) = SQRT(1.0_Long*(n+L+2)*(1.0_Long-(1.0_Long*n)/Nn)) ! <n+1^(l+1)| T_+ | n^l>
  ELSE
     ! N EVEN and L ODD OR N ODD and L EVEN
     COLIN = 1
     !
     ! INTERMEDIATE ROWS
     DO n = Nn-1, L+2, -2 ! notice reversed order
        TP(COLIN,COLIN) = SQRT(1.0_Long*(n+L+2)*(1.0_Long-(1.0_Long*n)/Nn)) ! <n+1^(l+1)| T_+ | n^l>
        TP(COLIN+1,COLIN) = -SQRT(1.0_Long*(n-L)*(1.0_Long-(n-1.0_Long)/Nn))    ! <n-1^(l+1)| T_+ | n^l>
        COLIN = COLIN + 1
     ENDDO
     !
     ! LAST COLUMN
     COLIN = DIM
     n = L
     TP(COLIN,COLIN) = SQRT(1.0_Long*(n+L+2)*(1.0_Long-(1.0_Long*n)/Nn)) ! <n+1^(l+1)| T_+ | n^l>
  ENDIF
  !
END SUBROUTINE TPMATRIX
