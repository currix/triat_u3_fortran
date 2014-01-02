      PROGRAM MINUIT_U3
C     
      IMPLICIT NONE
C     
      INTEGER IRD, IWR, ISAV
C     
      EXTERNAL FCN
C     
      DOUBLE PRECISION CHISQRE
      EXTERNAL CHISQRE
C
c
c     -----------------------------------------------
c     Initializing Minuit
c     Establishing unit numbers for data input/output
c     -----------------------------------------------
      ird=5                     ! Fortran Unit number for reading
      iwr=6                     ! Fortran Unit number for writing
      isav=7                    ! Fortran Unit number for  saving
c     call mninit(ird,iwr,isav)
      CALL MINTIO(IRD,IWR,ISAV)
c     
c     
c     *****************************************************************
c     
c     Calling to MINUIT
      CALL  MINUIT(FCN,CHISQRE)
c     
c     *****************************************************************
c     
c     
c     testing results....
c     call chisquare(NPAR,DUM,CHISQ,XV,IFLAG,FUTIL)
c     write(iwr,*)'chisquare: ',chisq
      STOP
      END
c          ! END MAIN PROGRAM 
c
c     *****************************************************************
