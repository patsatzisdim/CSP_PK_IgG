cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS THE CODE TO COMPARE THE SOLUTION OBTAINED
c                 WITH MF=21 (USER SUPPLIED JACOBIAN) AND MF=22 (INTERNAL JACOBIAN CALCULATION)
c

      program Sol_comparison
      implicit none

      integer, parameter :: N=10
      double precision :: T, asolm1(N), asolm2(N), dum1, dum2, rel(N)
      integer :: i


      open(11,file='Asol_MF=21.dat',form='formatted',action='read')
      open(12,file='Asol_MF=22.dat',form='formatted',action='read')
      open(13,file='RelAsol.dat',form='formatted',status='unknown')

      T=0.0d0
      do while(T.lt.3.50d+3) 
       read(11,*) T, (asolm1(i),i=1,N), dum1, dum2
       read(12,*) T, (asolm2(i),i=1,N), dum1, dum2
c       write(*,*) T
       do i=1,N
       	rel(i)=dabs(asolm2(i))/asolm1(i)
       enddo
       write(13,31) T, (rel(i),i=1,N)
       if(t.gt.14.0d+2.and.t.lt.1401.0d0) then
        write(*,*) T,asolm2(3), asolm1(3)
       endif 
      enddo


   31 format(E18.11,10(E25.16))
      end