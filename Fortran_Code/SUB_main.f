C      Author: Dimitris Patsatzis
c  
c
c      last update: 12092019
c
c      run the script as sh script.sh
c     
#include "defs.h"
c
c     preprocessing for optional compiling
c
      program main
      implicit none

      EXTERNAL FEX, JEX 
      include 'paramet.i'

      character sp      
      integer,parameter :: n=4  ! 4 species
      integer,parameter :: m=5  ! 1 reversible and 3 unidirectional reactions      
      integer :: i,k,j
c      INTEGER :: I,J,INFO,ITEMP      

cpd---------------------------------------LSODE work parameters 
      integer :: MF,ITOL,ISTATE,ITASK,IOPT,LRW,LIW
      double precision :: RTOL,ATOL
      double precision, allocatable, dimension(:):: RWORK
      integer, allocatable, dimension(:):: IWORK
cpd---------------------------------------
      double precision, dimension(n) :: yy,ydot,rlmod
      double precision :: t,dtStart,dt,tout    
      double precision, dimension(n,n) :: pd,rjac,alpha,betta,po
c      double precision, dimension(:,:), allocatable :: Ar1,Br1,As1,Bs1,Ar2,Br2,As2,Bs2
      integer :: ML, MU, NRPD
      double precision, dimension(n,2) :: reim,tau
      double precision :: FiEig(n),tpi(n,m),api(n,m),fip(n),RR(m),impi(n,m),dpn(n), BS(n,m)
      double precision :: CSPrtol,CSPatol(n),temp,sum1
      integer :: noEM,index,count,mdum
      double precision :: PoS(n,n), tpiS(n,m), apiS(n,m), impiS(n,m)
      integer :: PoI(n,n), tpiI(n,m), apiI(n,m), impiI(n,m)
      double precision :: T_start,T_end                               ! file units for printing
      integer :: indP1,indP2,indCharPo, indCharTPI, indCharAPI
      common k
      k=m ! the number of Reactions Rates
cpd
cpd
cpd   CPU TIME for ANAL vs NUM Jacobi
cpd
      call cpu_time(T_start)
cpd
cpd--------------------------------------------------
      sp=char(9) ! necessary gap used for write	  		  
cpd--------------------------------------------------  
      t=0.0d0    
cpd--------------------------------------------------
cpd   connect subroutine to get the initial conditions
      call InitCond(n,t,yy)
      write(*,*) t 
cpd--------------------------------------------------
cpd   time and inital step for LSODE to have a good start
      dtStart=1.0d-8   
      dt=dtStart
      tout=t+dt
      write(*,*) "Running till time t=",tend
cpd--------------------------------------------------      
cpd   Lsode parameters, check lsode documentation. 
      MF=21                     ! internal(22) or user(21) supplied of Jacobian	  
      LRW=22+9*n+n**2
      LIW=20+n
      ALLOCATE(RWORK(LRW))
      ALLOCATE(IWORK(LIW))
      ITOL=1
      RTOL=5.0D-16
      RTOL=3.0D-16
      ATOL=1.0D-22
      ITASK=1
      ISTATE=1
      IOPT=1                      ! optional inputs (0.0d0 or 0 means the default)
      RWORK(5)=0.0D0              ! initial stepsize H0
      RWORK(6)=0.0D0              ! maximum stepsize allowed
      RWORK(7)=0.0D0              ! minimum stepsize allowed
      RWORK(8)=0.0D0              !
      RWORK(9)=0.0D0              !
      RWORK(10)=0.0D0             ! 
      IWORK(5)=0                  ! maximum order (integration method)
      IWORK(6)=500000000           ! maximum number of internal steps 
      IWORK(7)=10                 ! maximum number of problems printed
      IWORK(8)=0                  !
      IWORK(9)=0                  !
      IWORK(10)=0                 !   
cpd--------------------------------------------------
cpd   tollerances (relative and absolute) for CSP criterion of fast timescales
cpd   currently working with rtol
      CSPrtol=1.0d-3
      CSPatol(1)=1.0d-7
      CSPatol(2)=1.0d-7
      CSPatol(3)=1.0d-7
      CSPatol(4)=1.0d-7
c      CSPatol(5)=1.0d-3
c      CSPatol(6)=1.0d-5
c      CSPatol(7)=1.0d-7
cpd--------------------------------------------------
cpd   open .dat files to write results
      open(101,file='Asol.dat',form='formatted',access='append') ! solution in concentration
      open(102,file='Arhs.dat',form='formatted',access='append') ! g(y)
      open(103,file='Arates.dat',form='formatted',access='append') ! rates
      open(104,file='Aeigen.dat',form='formatted',access='append') ! eigen
cc     JACOBI EVALUATION ONLY
c      open(105,file='AJac.dat',form='formatted',access='append') ! tau      
      open(106,file='Atmscl.dat',form='formatted',access='append') ! tau
      open(107,file='AFi_PosNeg.dat',form='formatted',access='append') ! Fi calculated by eigenvectors     
      open(108,file='ANofExhMod.dat',form='formatted',access='append') ! No. exhausted modes
c       
      open(110,file='APointers.dat',form='formatted',access='append') ! 
      open(111,file='ATPI.dat',form='formatted',access='append') ! 
      open(112,file='AAPI.dat',form='formatted',access='append') ! 
      open(113,file='AII.dat',form='formatted',access='append')
      open(114,file='AFi_Pos.dat',form='formatted',access='append') ! 
c     
c
c     for smart print in different files for every mode
      open(121,file='ADiag/Mod1.dat',form='formatted',access='append')
      open(122,file='ADiag/Mod2.dat',form='formatted',access='append')
      open(123,file='ADiag/Mod3.dat',form='formatted',access='append')
      open(124,file='ADiag/Mod4.dat',form='formatted',access='append')
c      open(125,file='ADiag/Mod5.dat',form='formatted',access='append')
c      open(126,file='ADiag/Mod6.dat',form='formatted',access='append')
c      open(127,file='ADiag/Mod7.dat',form='formatted',access='append')
c
c
c     for printing of each tool in different files
c
c     TPIs
      open(131,file='ADiag/ATPI1.dat',form='formatted',access='append')
      open(132,file='ADiag/ATPI2.dat',form='formatted',access='append')
      open(133,file='ADiag/ATPI3.dat',form='formatted',access='append')
      open(134,file='ADiag/ATPI4.dat',form='formatted',access='append')
c      open(135,file='ADiag/ATPI5.dat',form='formatted',access='append')
c      open(136,file='ADiag/ATPI6.dat',form='formatted',access='append')
c      open(137,file='ADiag/ATPI7.dat',form='formatted',access='append')
c     
c     APIs      
      open(141,file='ADiag/AAPI1.dat',form='formatted',access='append')
      open(142,file='ADiag/AAPI2.dat',form='formatted',access='append')
      open(143,file='ADiag/AAPI3.dat',form='formatted',access='append')
      open(144,file='ADiag/AAPI4.dat',form='formatted',access='append')
c      open(145,file='ADiag/AAPI5.dat',form='formatted',access='append')
c      open(146,file='ADiag/AAPI6.dat',form='formatted',access='append')
c      open(147,file='ADiag/AAPI7.dat',form='formatted',access='append')
c
c
c     IIs
      open(151,file='ADiag/AII1.dat',form='formatted',access='append')
      open(152,file='ADiag/AII2.dat',form='formatted',access='append')
      open(153,file='ADiag/AII3.dat',form='formatted',access='append')
      open(154,file='ADiag/AII4.dat',form='formatted',access='append')
c      open(155,file='ADiag/AII5.dat',form='formatted',access='append')
c      open(156,file='ADiag/AII6.dat',form='formatted',access='append')
c      open(157,file='ADiag/AII7.dat',form='formatted',access='append')
c
c     Pointers
      open(161,file='ADiag/APo1.dat',form='formatted',access='append')
      open(162,file='ADiag/APo2.dat',form='formatted',access='append')
      open(163,file='ADiag/APo3.dat',form='formatted',access='append')
      open(164,file='ADiag/APo4.dat',form='formatted',access='append')
c      open(165,file='ADiag/APo5.dat',form='formatted',access='append')
c      open(166,file='ADiag/APo6.dat',form='formatted',access='append')
c      open(167,file='ADiag/APo7.dat',form='formatted',access='append')



cpd--------------------------------------------------
cpd   Compute eigenvalues-eigenvectors-timescales-Fi from eigs
      call jac(n,t,yy,rjac)   
      call eigen(rjac,alpha,betta,rlmod,reim)
      do i=1,n
        if(reim(i,1).eq.0.0d0.and.reim(i,2).eq.0.0d0) then
          tau(i,1)=1.0d+40
        else
          tau(i,1)=1.0d0/dsqrt(reim(i,1)*reim(i,1)+reim(i,2)*reim(i,2))
        endif
       if (reim(i,1).lt.0.0d0) then
        tau(i,2)=-1.0d0
       else
        tau(i,2)=1.0d0
       endif
      enddo
      call fex(n,t,yy,ydot)
      call smult(n,n,1,betta,ydot,FiEig,n,n,1)
      call rates(n,t,yy,k,RR)
cpd--------------------------------------------------
cpd   Print solution, rates, rhs, eigenvalues and timescales   
      write(101,30) t,yy(1:n)
      write(102,30) t,ydot(1:n)
      write(103,31) t,RR(1:k)
      write(104,32) t,(reim(i,1:2),i=1,n)
      write(106,30) t,tau(1:n,1)
      write(107,30) t,FiEig(1:n)  
cpd--------------------------------------------------                       Eigevectors for the criterion
cpd   number of exhausted modes with eigenvectors.
      call NEM_Eigen(n,k,t,yy,CSPrtol,CSPatol,rlmod,alpha,FiEig,noEM)
#ifdef CSPtools       
cpd--------------------------------------------------                         
cpd   Diagnostics from eigenvectors 
      call diag(n,k,NoEM,t,yy,alpha,betta,reim,po,tpi,api,impi,fip,BS)
      temp=1.0d-40
      do i=1,n
       if(temp.le.dabs(fip(i))) then
        temp=dabs(fip(i))
        index=i
       endif
      enddo
      write(108,35) t,noEM,index



C cpd   ONLY FOR SMOOTHENINGs APIs FOR FIGURES!!!!! 
        if (api(1,1).lt.0.0d0) api(1,1:k)=-api(1,1:k)
C        if (api(2,6).lt.0.0d0) api(2,1:k)=-api(2,1:k)
C        if (api(3,3).lt.0.0d0) api(3,1:k)=-api(3,1:k)
C        if (api(4,2).lt.0.0d0) api(4,1:k)=-api(4,1:k) 
C        if (api(5,2).lt.0.0d0) api(5,1:k)=-api(5,1:k)  
C        if (api(6,3).lt.0.0d0) api(6,1:k)=-api(6,1:k) 
C        if (api(7,10).lt.0.0d0) api(7,1:k)=-api(7,1:k)                              
cpd--------------------------------------------------
cpd   Print diagnostics Pointer, TPI, API, II and Fi 
      call Print_Diag(n,k,t,yy,po,tpi,api,impi,fip,110,111,112,113,114)
cpd   
cpd   and print them in seperate files. First TPI, then API, then II, then PO and last the diagonal of PO
      call Print_DiagS(n,k,t,po,tpi,api,impi,131,141,151,161)
ccpd--------------------------------------------------
cc     Diagnostics for only the characteristic mode
c      call charDiags_print(n,k,t,yy,NoEM,po,tpi,api,impi,indCharPo,indCharTPI,indCharAPI)
c      if (t.eq.0.0d0) then
c      call sortDiag(n,k,Po,tpi,api,impi,PoS,tpiS,apiS,impiS,PoI,tpiI,apiI,impiI)
c      call smart_print(n,k,t,yy,PoS,tpiS,apiS,impiS,PoI,tpiI,apiI,impiI,BS,121,127)
c      endif
ccpd--------------------------------------------------
#endif
cpd   START LSODE LOOP
   10 continue 
    
      call DLSODE (FEX, n, yy, t, tout, ITOL, RTOL, ATOL, ITASK,
     1               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
     
      if (ISTATE.lt.0)  go to 98
#ifdef CSPtools          
cpd--------------------------------------------------
cpd   Skip some time points that NoEM increases due to artifacts
c      if (t.gt.0.0709d0.and.t.lt.0.0719) go to 201
#endif      
cpd--------------------------------------------------
cpd   Compute eigenvalues-eigenvectors-timescales
      call jac(n,t,yy,rjac)      
      call eigen(rjac,alpha,betta,rlmod,reim)


      do i=1,n
        if(reim(i,1).eq.0.0d0.and.reim(i,2).eq.0.0d0) then
          tau(i,1)=1.0d+40
        else
          tau(i,1)=1.0d0/dsqrt(reim(i,1)*reim(i,1)+reim(i,2)*reim(i,2))
        endif
       if (reim(i,1).lt.0.0d0) then
        tau(i,2)=-1.0d0
       else
        tau(i,2)=1.0d0
       endif
      enddo
      call fex(n,t,yy,ydot)
      call smult(n,n,1,betta,ydot,FiEig,n,n,1)
      call rates(n,t,yy,k,RR)
cpd--------------------------------------------------
cpd   Print solution, rates, rhs, eigenvalues and timescales   
      write(101,30) t,yy(1:n)
      write(102,30) t,ydot(1:n)
      write(103,31) t,RR(1:k)
      write(104,32) t,(reim(i,1:2),i=1,n)
      write(106,30) t,tau(1:n,1)
      write(107,30) t, FiEig(1:n)
cpd   number of exhausted modes with Eigenvectors
      call NEM_Eigen(n,k,t,yy,CSPrtol,CSPatol,rlmod,alpha,FiEig,noEM)
#ifdef CSPtools      
cpd--------------------------------------------------                            
cpd   Diagnostics from eigenvectors 
      call diag(n,k,NoEM,t,yy,alpha,betta,reim,po,tpi,api,impi,fip,BS)
      temp=1.0d-40
      do i=1,n
       if(temp.lt.dabs(FiEig(i))) then
        temp=dabs(FiEig(i))
        index=i
       endif
      enddo
      write(108,35) t,noEM,index

C cpd   ONLY FOR SMOOTHENINGs APIs FOR FIGURES!!!!! 
        if (api(1,1).lt.0.0d0) api(1,1:k)=-api(1,1:k)
C        if (api(2,6).lt.0.0d0) api(2,1:k)=-api(2,1:k)
C        if (api(3,3).lt.0.0d0) api(3,1:k)=-api(3,1:k)
C        if (api(4,2).lt.0.0d0) api(4,1:k)=-api(4,1:k) 
C        if (api(5,2).lt.0.0d0) api(5,1:k)=-api(5,1:k)  
C        if (api(6,3).lt.0.0d0) api(6,1:k)=-api(6,1:k) 
C        if (api(7,10).lt.0.0d0) api(7,1:k)=-api(7,1:k)                   
cpd--------------------------------------------------
cpd   Print diagnostics Pointer, TPI, API, II and Fi 
      call Print_Diag(n,k,t,yy,po,tpi,api,impi,fip,110,111,112,113,114)
cpd   
cpd   and print them in seperate files. First TPI, then API, then II, then PO and last the diagonal of PO
      call Print_DiagS(n,k,t,po,tpi,api,impi,131,141,151,161)      
ccpd--------------------------------------------------
cc     Diagnostics for only the characteristic mode
c      call charDiags_print(n,k,t,yy,NoEM,po,tpi,api,impi,indCharPo,indCharTPI,indCharAPI)
cpd--------------------------------------------------
cpd   Diagnostics per point-Choose points first
      if (tout.gt.0.1d-1.and.tout.lt.0.1001d-1) go to 200
      if (tout.gt.0.1d+2.and.tout.lt.0.1001d+2) go to 200
      if (tout.gt.0.75d+2.and.tout.lt.0.7501d+2) go to 200            
      if (tout.gt.0.1d+3.and.tout.lt.0.1001d+3) go to 200
      if (tout.gt.0.1d+4.and.tout.lt.0.1001d+4) go to 200
      go to 201
cpd   Sort and print the diagnostics in selected time
  200 call sortDiag(n,k,Po,tpi,api,impi,PoS,tpiS,apiS,impiS,PoI,tpiI,apiI,impiI)
      call smart_print(n,k,t,yy,PoS,tpiS,apiS,impiS,PoI,tpiI,apiI,impiI,BS,121,124)
  201 continue      
cpd--------------------------------------------------
#endif
cpd   Define timestep
      dt=1.0d0
      if (tout.lt.1.0d+3) DT=1.0d-1
      if (tout.lt.1.0d+2) DT=1.0d-2
      if (tout.lt.1.0d+1) DT=1.0d-3 ! define time step
      if (tout.lt.1.0d+0) DT=1.0d-4 ! define time step
      if (tout.lt.1.0d-1) DT=1.0d-5 ! define time step
      if (tout.lt.1.0d-2) DT=1.0d-6 ! define time step
      if (tout.lt.1.0d-3) DT=1.0d-7 ! define time step
      if (tout.lt.1.0d-4) DT=1.0d-8 ! define time step
cpd--------------------------------------------------
cpd   Define next iteration   
      t=tout
      tout=tout+dt
      if(t.lt.tend) go to 10  ! Define tend above
cpd   END LSODE LOOP
c
cpd--------------------------------------------------  
      write(*,*) 'You may now close this window'   
c
cpd------------------------------------------------------------------
cpd   formats for all the .dat
   30 format(E18.11,4('    ',E25.16))   
   31 format(E18.11,5('    ',E25.16))
   32 format(E18.11,4('    ',E25.16,'   ',E10.3))
c   33 format(i18,29('    ',E25.16))     
   35 format(E18.11,2('    ',i3))
   36 format(E18.11,4('    ',E25.16))
   37 format(E18.8,i7,67(E25.10))
 
c   
   99 format(///'t=',D15.8,'    Error halt.. ISTATE =',I3)      
      go to 100
   98 write(*,99)  t, ISTATE
  100 continue   
cpd--------------------------------------------------   
cpd   close all the files
      close(101)
      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)      
      close(110)  
      close(111)
      close(112)
      close(113)
      close(114)
c      close(115)
c      close(116)
c      close(117)
c      close(118)
c      close(119)
c      
      close(121)
      close(122)
      close(123)
      close(124)
c      close(125)
c      close(126)
c      close(127)
      close(131)
      close(132)
      close(133)
      close(134)
c      close(135)
c      close(136)
c      close(137)
      close(141)
      close(142)
      close(143)
      close(144)
c      close(145)
c      close(146)
c      close(147)
      close(151)
      close(152)
      close(153)
      close(154)
c      close(155)
c      close(156)
c      close(157)
      close(161)
      close(162)
      close(163)
      close(164)
c      close(165)
c      close(166)
c      close(167)     
c      close(128)
c      close(129)
c      close(130)
c      close(131)      
c      close(132)
c      close(133)
c      close(134)
c      close(135)
c      close(136)
c      close(137)
c      close(138)
c      close(139)
c      close(140)
c      close(141)
c      close(142)
c      close(143)
c      close(144)
c      close(145)
c      close(146)
c      close(147)
c      close(148)
c      close(149)
c      close(150)
c      close(151)
c      close(152)
c      close(153)
c      close(154)
c
c
c
c
cpd     CPU TIME for ANAL vs NUM Jacobi
cpd
      call cpu_time(T_end)
c
c
      write(*,*) T_end-T_start
c
c    
      STOP
      END
ccpd------------------------------------------------------------------
ccpd   END of main routine
cc
ccpd------------------------------------------------------------------
ccpd   Subroutine for sorting diagnostics
ccpd   Inputs:  n: number of equations
ccpd            k: number of reactions
ccpd            Po, tpi, api: matrices to be sorted
cc
ccpd   Outputs: PoS, tpiS, apiS: sorted per rows diagnostic matrices 
ccpd            PoI, tpiI, apiI: keeping trace of the reactions. Index matrices  
      subroutine sortDiag(n,k,Po,tpi,api,impi,PoS,tpiS,apiS,impiS,PoI,tpiI,apiI,impiI)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: Po(n,n), tpi(n,k), api(n,k),impi(n,k)
      double precision, intent(out) :: PoS(n,n), tpiS(n,k), apiS(n,k),impiS(n,k)
      integer, intent(out) :: PoI(n,n), tpiI(n,k), apiI(n,k),impiI(n,k)
      integer :: i
      double precision :: PoRow(n), tpiRow(k), apiRow(k), impiRow(k)
      integer :: PoInd(n),tpiInd(k),apiInd(k), impiInd(k)
cpd--------------------------------------------------
cpd   Slice the matrices and sort over rows
      do i=1,n      
       PoRow(:)=Po(i,:)
       tpiRow(:)=tpi(i,:)
       apiRow(:)=api(i,:)
       impiRow(:)=impi(i,:)
cpd--------------------------------------------------
cpd    Sort over each slice
       call SortAbs(PoRow,n,PoInd)
       call SortAbs(tpiRow,k,tpiInd)
       call SortAbs(apiRow,k,apiInd)
       call SortAbs(impiRow,k,impiInd)
cpd--------------------------------------------------
cpd    Bind the slices into matrices and their indices        
       PoS(i,:)=PoRow(:)
       PoI(i,:)=PoInd(:)
       tpiS(i,:)=tpiRow(:)
       tpiI(i,:)=tpiInd(:)
       apiS(i,:)=apiRow(:)
       apiI(i,:)=apiInd(:)
       impiS(i,:)=impiRow(:)
       impiI(i,:)=impiInd(:)
      enddo 
cpd-------------------------------------------------- 
      return
      end    
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------           
cpd------------------------------------------------------------------
c
c
c
c
c                     RHS
c
c
c
c
c      
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------           
cpd------------------------------------------------------------------
cpd   SUBROUT used by lsode to solve the system
cpd   contains the RHS of the ODE system dy/dt=g(y) 
      subroutine FEX(n,t,yy,ydot)
      implicit none
      integer,intent(in) :: n
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: ydot(n)
cpd
      integer :: k            
      common k            ! number of reactions passed with common
cpd            
      double precision :: RR(k), st(n,k)

      include 'paramet.i'      
cpd
cpd   SUB_ProbDef contains these subs
cpd   generating fex through 
      call rates(n,t,yy,k,RR)
      call stoic(n,k,st)
      call smult(n,k,1,st,RR,ydot,n,k,1)

      return
      end
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------           
cpd------------------------------------------------------------------
c
c
c
c
c                       JACOBIAN                       
c
c
c
c      
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------           
cpd------------------------------------------------------------------
cpd   SUBROUT used by lsode to solve the system, 
cpd   contains the jacobian of the RHS of the ODE system dy/dt=g(y) 
      subroutine JEX(n,t,yy,ML,MU,pd,NRPD)
      implicit none
      integer, intent(in) :: n, ML, MU, NRPD
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: pd(n,n)
cpd
      integer :: k            
      common k 
      double precision :: st(n,k),gR(k,n)
cpd            
#ifdef NUM_JAC
      double precision :: ydot(n),epsfcn,wa1(n),wa2(n)
      integer ml1,mu1
      external fcn
#endif
      pd(:,:)=0.0d0
#ifdef ANAL_JAC
cpd
cpd   SUB_ProbDef contains these subs
      call stoic(n,k,st)
      call gradR(n,t,yy,k,gR)
      call smult(n,k,n,st,gR,pd,n,k,n)
#endif

#ifdef NUM_JAC
      epsfcn=0.0d0
      ml1=n                    ! use in banded Jacobian only
      mu1=n                    ! use in banded Jacobian only
      call fex(n,t,yy,ydot)
      call fdjac1(fcn,n,yy,ydot,pd,n,ml1,mu1,epsfcn,wa1,wa2)
#endif    
      return
      end
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------           
cpd------------------------------------------------------------------
c
c
c
c
c                       JACOBIAN                        
c
c
c
c      
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------      	 	  
cpd------------------------------------------------------------------
cpd   SUBROUT contains the jacobian of the RHS of the ODE system dy/dt=g(y) 
cpd   same with jex, user-use because jex is altered after calculations
      subroutine  jac(n,t,yy,pd)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: pd(n,n)
cpd
      integer :: k            
      common k 
      double precision :: st(n,k),gR(k,n)
cpd            
#ifdef NUM_JAC
      double precision :: ydot(n),epsfcn,wa1(n),wa2(n)
      integer ml1,mu1
      external fcn
#endif
      pd(:,:)=0.0d0
#ifdef ANAL_JAC
cpd
cpd   SUB_ProbDef contains these subs
      call stoic(n,k,st)
      call gradR(n,t,yy,k,gR)
      call smult(n,k,n,st,gR,pd,n,k,n)
#endif
c
#ifdef NUM_JAC
      epsfcn=0.0d0
      ml1=n                    ! use in banded Jacobian only
      mu1=n                    ! use in banded Jacobian only
      call fex(n,t,yy,ydot)
c      write(*,*) ydot
      call fdjac1(fcn,n,yy,ydot,pd,n,ml1,mu1,epsfcn,wa1,wa2)
#endif   
      return
      end
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------           
cpd------------------------------------------------------------------
c
c
c
c
c                       FEX in the appropriate form for the numerical Jacobi (fdjac2.f) to use it
c
c
c
c      
cpd------------------------------------------------------------------      
cpd------------------------------------------------------------------  
      subroutine fcn(n,x,fvec,iflag)
      integer n,iflag
      double precision x(n),fvec(n)
cpd            
      double precision :: t_dum   !time here is dummy, as an autonomous system restricts
c     
      t_dum=0.0d0
      call fex(n,t_dum,x,fvec)
c
      return
      end

