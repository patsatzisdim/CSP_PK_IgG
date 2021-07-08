cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS THE INITIAL CONDITIONS FOR THE PROBLEM TO RUN
c
c
c            WATCH OUT! IT GIVES YOU THE OPPORTUNITY TO READ THE INITIAL CONDITIONS FROM
C                       THE END OF THE CURRENT SOLUTION FILE. THIS HELPS FOR MULTPLE CALLS
C                       FROM THE SCRIPT FILE!!!!!!
C
C            (COMMENTED OUT CURRENTLY)
c

      subroutine InitCond(n,t,y)
      implicit none
      include 'paramet.i'
      integer :: i
      integer, intent(in) :: n
      double precision, intent(out) :: y(n)
      double precision, intent(inout) :: t
      integer :: IOState,iflag


c      iflag=2
      do i=1,n
       y(i)=1.0d0
      enddo

      y(1)=DoseA/Vp       ! A                   
      y(2)=0.0d0          ! B
      y(3)=C0             ! C
      y(4)=0.0d0          ! D

c      write(*,*) 'Initial Conditions Taken for steady-state'
c      else
ccpd   In case of initial conditions from other file
ccpd   read from the end of file
c      open(3,file='Asol.dat',action='read')
c 
c      do i=1,100000000
c       read(3,*,IOSTAT=IOstate) t,y(1:n)
c      if(IOstate.ne.0) exit
c      enddo
c      close(3)
c      write(*,*) 'Initial Conditions Taken from Asol.dat' 
c
c      y(2)=ICM*y(2) 
c    
c      endif
  
      return 
      end
cpd------------------------------------------------------------------
