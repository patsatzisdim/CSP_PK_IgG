cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS THE SUBROUTINES THAT DEFINE THE PROBLEM
c
c                     REACTION RATES (rates)
c                     STOICHIOMETRY (stoic)     
c                     GRADIENT OF REACTION RATES (gradR)
c                     TIME DERIVATIVE OF JACOBIAN (Djac_dt)                  
c                 
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd   SUBROUT used to generate the Rates
      subroutine rates(n,t,yy,k,RR)
      implicit none 
      integer, intent(in) :: n,k
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: RR(k)
      double precision :: y1,y2,y3,y4,y5,y6,y7
      include 'paramet.i'

cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
  
cpd
cpd   Reaction Rates from mathematica
      RR(1)=(CLup*y1)/Ve      ! R1
      RR(2)=kdeg*y2           ! R2
      RR(3)=kon*y2*y3         ! R3f
      RR(4)=(CLup*y4)/Ve      ! R4
      RR(5)=koff*y4           ! R3b  
cpd          
      return
      end
cpd
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
c
c
c
c         Stoichiometry
c
c
c
c
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd   SUBROUT used for stoichiometry
      subroutine stoic(n,k,st)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(out) :: st(n,k)
      include 'paramet.i'
cpd
cpd   initialize first      
      st(:,:)=0.0d0
cpd
cpd   values from mathematica
      st(1,1)=-Ve/Vp
      st(1,4)=Ve/Vp
      st(2,1)=+1.0d0
      st(2,2)=-1.0d0
      st(2,3)=-1.0d0
      st(2,5)=+1.0d0
      st(3,3)=-1.0d0
      st(3,4)=+1.0d0
      st(3,5)=+1.0d0
      st(4,3)=+1.0d0
      st(4,4)=-1.0d0
      st(4,5)=-1.0d0      

cpd      
      return
      end      
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
c
c
c
c
c
c
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd   SUBROUT used to calulate grad(R) analytically from mathematica
      subroutine gradR(n,t,yy,k,gR) 
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: gR(k,n)
      double precision :: y1,y2,y3,y4,y5,y6,y7
      include 'paramet.i'
cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
cpd
cpd   initialize first
      gR(:,:)=0.0d0

cpd
cpd   values grad(R) from mathematica  
c     
      gR(1,1)=CLup/Ve
      gR(2,2)=kdeg
      gR(3,2)=kon*y3
      gR(3,3)=kon*y2
      gR(4,4)=CLup/Ve
      gR(5,4)=koff 
   
c
      return
      end


ccpd-----------------------------------------------
ccpd
ccpd------------------------------------------------------------------
ccpd   SUBROUT to calculate dJac/dt
      subroutine Djac_dt(n,t,yy,dJacdt)
      integer, intent(in) :: n
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: dJacdt(n,n)
      double precision :: y1,y2,y3,y4,y5,y6,y7,ydot(n)
      double precision :: dJ1(n,n),dJ2(n,n),dJ3(n,n),dJ4(n,n),dJ5(n,n),dJ6(n,n),dJ7(n,n)
      include 'paramet.i'
c
c     Right hand side
      call FEX(n,t,yy,ydot)
c      
cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
      y5=yy(5)
      y6=yy(6)
      y7=yy(7)
c
c
c
      dJ1(:,:)=0.0d0
      dJ1(:,:)=dJ1(:,:)*ydot(1)
cc
cc
      dJ2(:,:)=0.0d0
      dJ2(:,:)=dJ2(:,:)*ydot(2)
cc      
cc
      dJ3(:,:)=0.0d0
      dJ3(:,:)=dJ3(:,:)*ydot(3)
cc
cc     
      dJ4(:,:)=0.0d0 
      dJ4(:,:)=dJ4(:,:)*ydot(4)
cc
cc     
      dJ5(:,:)=0.0d0 
      dJ5(:,:)=dJ5(:,:)*ydot(5)
cc
cc     
      dJ6(:,:)=0.0d0 
      dJ6(:,:)=dJ6(:,:)*ydot(6)
cc
cc 
      dJ7(:,:)=0.0d0 
      dJ7(:,:)=dJ7(:,:)*ydot(7)
c            
c
      dJacdt(:,:)=dJ1(:,:)+dJ2(:,:)+dJ3(:,:)+dJ4(:,:)+dJ5(:,:)+dJ6(:,:)  
c
c           
      return
      end
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------