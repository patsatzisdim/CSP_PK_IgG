cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS THE CODE TO CALCULATE THE CSP BASIS VECTORS
c                                    AND TO CALCULATE THE NUMBER OF EXHAUSTED MODES
c                 
c     Inputs:     n: No. species
c                 k: No. reactions
c                 t: time point
c                 y: concentrations
c                 Ms: No. fast timescales
c                 PO: CSP Pointer based on eigenvectors
c                 NoPh: Flag for the number of PHASES
c
c     Outputs:    Ar1, As1, Br1, Bs1: CSP basis vectors after PHASE 1 : 1 Br and 1 Ar refinement
c                 Ar2, As2, Br2, Bs2: CSP basis vectors after PHASE 2 : 2 Br and 1 Ar refinement
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
c               NEEDS DEBUGGING!!!!!!!!!!!!!!!!!!!!!
C
C
C
C
      subroutine CSP_kern(n,k,t,y,Ms,PO,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2,NoPh)
      implicit none
      integer, intent(in) :: n,k,Ms,NoPh
      double precision, intent(in) :: t, y(n), PO(n,n)
      double precision :: Ar0(n,Ms),As0(n,n-Ms),Br0(Ms,n),Bs0(n-Ms,n),dBrdt(Ms,n)
      double precision :: Ar(n,Ms),As(n,n-Ms),Br(Ms,n),Bs(n-Ms,n)
      double precision, intent(out) :: Ar1(n,Ms),As1(n,n-Ms),Br1(Ms,n),Bs1(n-Ms,n)
      double precision, intent(out) :: Ar2(n,Ms),As2(n,n-Ms),Br2(Ms,n),Bs2(n-Ms,n)
      integer :: i
      logical :: Br1call

cpd-----------------------------------------------  
c      Prelimineries of CSP: initial Basis Vectors
cpd-----------------------------------------------  
      call Init_BV(n,Ms,PO,Ar0,As0,Br0,Bs0)
cpd-----------------------------------------------  
c     Phase 1: dBr0/dt=0     Br and Ar refs
cpd-----------------------------------------------  
c     Step 1: B^r refinement
      Br1call=.true.
      call Br_ref(n,k,Ms,t,y,Br1call,Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,dBrdt)    
c     Step 2: A^r refinement
      call Ar_ref(n,k,Ms,t,y,Ar,As,Br,Bs,Ar1,As1,Br1,Bs1)
      if(NoPh.eq.1) go to 10
cpd-----------------------------------------------  
c     Phase 2: dBr0/dt=/0    Br ref ONLY!!!                                       
cpd-----------------------------------------------  
c     Step 1: B^r refinement
      Br1call=.false.
      call Br_ref(n,k,Ms,t,y,Br1call,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2,dBrdt)
cpd-----------------------------------------------  
c     Cut the small values in Basis Vectors   !LEAVE THIS FOR NOW
cpd----------------------------------------------- 
   10 continue     
      return 
      end

cpd-------------------------------------------------------------------  
c
c     Subrout for Br refinement 
c
cpd-------------------------------------------------------------------
      subroutine Br_ref(n,k,m,t,y,Br1call,Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,dBrdt)
      implicit none	
      integer, intent(in) :: n,m,k
      double precision, intent(in) :: t,y(n),Ar0(n,m),As0(n,n-m),Br0(m,n),Bs0(n-m,n)
      double precision, intent(out) :: Ar(n,m),As(n,n-m),Br(m,n),Bs(n-m,n)
      double precision, intent(inout) :: dBrdt(m,n)   !Output at phase 1, Input at phase 2
      logical, intent(in) :: Br1call
      double precision :: st(n,k),gR(k,n),djac(n,n),dBr(m,n),Lrr(m,m),Trr(m,m),dJacdt(n,n)
      double precision :: temp1(m,n),Un(n,n),temp2(n,n),temp3(m,n)
      integer :: IPVT(n),INFO,i
cpd-----------------------------------------------
cpd   Calculate Lamba^r_r
      call stoic(n,k,st)
      call gradR(n,t,y,k,gR)
c      djac=matmul(st,gR)
      call smult(n,k,n,st,gR,djac,n,k,n)      ! Jacobian           
cpd-----------------------------------------------
      if (Br1call) then
       dBr(:,:)=0.0d0
      else
       dBr(:,:)=dBrdt(:,:)
      endif                               ! dBr0/dt
cpd-----------------------------------------------
c      temp1=matmul(Br0,djac)
      call smult(m,n,n,Br0,djac,temp1,m,n,n)  ! Br0.J
      temp1(:,:)=temp1(:,:)+dBr(:,:)          ! Br0.J+dBr0/dt
c      Lrr=matmul(temp1,Ar0)
      call smult(m,n,m,temp1,Ar0,Lrr,m,n,m)   ! (Br0.J+dBr0/dt).Ar0=Lrr   

      call sinve(m,Lrr,Trr,m,IPVT,INFO)       ! Trr=(Lrr)^-1
      if(info.ne.0) then
        write(*,*) info
        stop
      endif                   
cpd----------------------------------------------- 
c      Br=matmul(Trr,temp1)
      call smult(m,m,n,Trr,temp1,Br,m,m,n)    ! Br=Trr.Br0.J   
      Ar(:,:)=Ar0(:,:)					      ! Ar=Ar0
      Bs(:,:)=Bs0(:,:)					      ! Bs=Bs0   
      call unitary(n,Un,n)                    
      call smult(n,m,n,Ar,Br,temp2,n,m,n)	  ! Ar.Br
      temp2(:,:)=Un(:,:)-temp2(:,:)           ! I-Ar.Br
c      As=matmul(temp2,As0)
      call smult(n,n,n-m,temp2,As0,As,n,n,n-m)! (I-Ar.Br).As0 
cpd-----------------------------------------------
      if(Br1call) then                                                               
       temp1(:,:)=0.0d0
       call smult(m,m,n,Trr,Br0,temp1,m,m,n)  !Trr.Br0
       call Djac_dt(n,t,y,dJacdt)
       call smult(m,n,n,temp1,dJacdt,temp3,m,n,n) ! (Trr.Br0).dJ/dt
       call smult(m,n,n,temp3,temp2,dBrdt,m,n,n)       ! (Trr.Br0).dJ/dt.(I-Ar.Br)
      endif
cpd-----------------------------------------------
      return
      end
c
c
cpd-------------------------------------------------------------------  
c
c     Subrout for Ar refinement 
c
cpd-------------------------------------------------------------------
      subroutine Ar_ref(n,k,m,t,y,Ar0,As0,Br0,Bs0,Ar,As,Br,Bs)
      implicit none	
      integer, intent(in) :: n,m,k
      double precision, intent(in) :: t,y(n),Ar0(n,m),As0(n,n-m),
     - Br0(m,n),Bs0(n-m,n)
      double precision, intent(out) :: Ar(n,m),As(n,n-m),Br(m,n),
     - Bs(n-m,n)
      double precision :: st(n,k),gR(k,n),djac(n,n),dBr(m,n),Lrr(m,m),
     - Trr(m,m),dJacdt(n,n)
      double precision :: temp1(m,n),Un(n,n),temp2(n,m),temp3(n,n)
      integer :: IPVT(n),INFO
cpd-----------------------------------------------
cpd   Calculate Lamba^r_r
      call stoic(n,k,st)
      call gradR(n,t,y,k,gR)
c      djac=matmul(st,gR)
      call smult(n,k,n,st,gR,djac,n,k,n)      ! Jacobian              
cpd-----------------------------------------------
      call smult(m,n,n,Br0,djac,temp1,m,n,n)  ! Br0.J
      call smult(m,n,m,temp1,Ar0,Lrr,m,n,m)   ! (Br0.J).Ar0=Lrr
      call sinve(m,Lrr,Trr,m,IPVT,INFO)       ! Trr=(Lrr)^-1  
cpd----------------------------------------------- 
      call smult(n,n,m,djac,Ar0,temp2,n,n,m)   
      call smult(n,m,m,temp2,Trr,Ar,n,m,m)     ! Ar=(J.Ar0).Trr
      As(:,:)=As0(:,:)					       ! As=As0
      Br(:,:)=Br0(:,:)					       ! Br=Br0   
      call unitary(n,Un,n)                    
      call smult(n,m,n,Ar,Br,temp3,n,m,n)	   ! Ar.Br
      temp3(:,:)=Un(:,:)-temp3(:,:)            ! I-Ar.Br

      call smult(n-m,n,n,Bs0,temp3,Bs,n-m,n,n) ! Bs0.(I-Ar.Br) 


cpd-----------------------------------------------
      return
      end
c
c      
cpd-------------------------------------------------------------------  
c
c     Subrout for initial set of Basis vectors. 
c     According to the greatest specie of the PO in each mode.
c
cpd-------------------------------------------------------------------  
cpd-------------------------------------------------------------------  
cpd-------------------------------------------------------------------  
cpd------------------------------------------------------------------- 
      subroutine Init_BV(n,m,PO,Ar,As,Br,Bs)
      implicit none
      integer, intent(in) :: n,m
      double precision, intent(in) :: PO(n,n)
      double precision, intent(out) :: Ar(n,m),As(n,n-m),Br(m,n),
     - Bs(n-m,n)
      integer :: i,j,maxLocPo(n)
      double precision :: PO2(n)
cpd-----------------------------------------------  
      Ar(1:n,1:m)=0.0d0
      As(1:n,1:n-m)=0.0d0
      Br(1:m,1:n)=0.0d0
      Bs(1:n-m,1:n)=0.0d0
cpd-----------------------------------------------        
cpd   read CSP pointer and put 1's and 0's for the species of every mode
cpd   IN THE CASE WHERE IN A MODE THE FAST SPECIES IS FAST IN A PREVIOUS MODE, PUT 1 TO THE NEXT FASTER.
      do i=1,m
       maxLocPo(i)=maxloc(PO(i,:),n)          ! the greater  
        do j=1,i-1
         if(maxLocPo(j).eq.maxLocPo(i)) then            ! if the same species is pointed in another mode
          PO2(:)=PO(i,:)                                ! transfer the pointer in a new array
          PO2(maxLocPo(i))=0.0d0                        ! zero the fastest species
          maxLocPo(i)=maxloc(PO2,n)                     ! and find the second faster one
         endif
        enddo
       Ar(maxLocPo(i),i)=1.0d0                 
       Br(i,maxLocPo(i))=1.0d0
      enddo
c
cpd   for the slow ones check the same.
c      
      do i=1,n-m
       maxLocPo(i+m)=maxloc(PO(i+m,:),n)          ! the greater
       do j=1,i-1
         if(maxLocPo(j+m).eq.maxLocPo(i+m)) then            ! if the same species is pointed in another mode
          PO2(:)=PO(i+m,:)                                  ! transfer the pointer in a new array
          PO2(maxLocPo(i+m))=0.0d0                          ! zero the fastest species
          maxLocPo(i+m)=maxloc(PO2,n)                       ! and find the second faster one
         endif
        enddo                 
       As(maxLocPo(i+m),i)=1.0d0                 
       Bs(i,maxLocPo(i+m))=1.0d0
      enddo
cpd-----------------------------------------------  
      return
      end
cpd-----------------------------------------------     
c
c      
cpd-------------------------------------------------------------------  
c
c     Subrout for the calculation of the number of exhausted 
c     modes through the CSP basis vectors 
c
cpd-------------------------------------------------------------------  
cpd-------------------------------------------------------------------  
cpd-------------------------------------------------------------------  
cpd------------------------------------------------------------------- 
      subroutine NEM_CSPbv(n,k,t,yy,CSPrtol,CSPatol,rlmod,po,noEM)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: t, yy(n),CSPrtol,CSPatol(n),po(n,n),rlmod(n)
      integer, intent(out) :: noEM
      integer :: i,j,m1,m2,m,NoPh
      double precision :: yyer(n),temp(n),fi(n),ydot(n)
      double precision, dimension(:,:), allocatable :: Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,dummy,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2
      double precision, dimension(:), allocatable :: fr
cpd-----------------------------------------------
cpd-----------------------------------------------
      noEM=0
  200 continue
      m1=noEM+1
      m2=noEM+2                               !Next two modes
      do i=1,n
       temp(i)=0.0d0
      enddo
cpd-----------------------------------------------      
c     calculate CSP basis vectors (1-ref) with m1 exhausted modes
      allocate(Ar0(n,m1),As0(n,n-m1),Br0(m1,n),Bs0(n-m1,n),Ar(n,m1),As(n,n-m1),Br(m1,n),Bs(n-m1,n),dummy(m1,n),
     - Ar1(n,m1),As1(n,n-m1),Br1(m1,n),Bs1(n-m1,n),fr(m1),Ar2(n,m1),As2(n,n-m1),Br2(m1,n),Bs2(n-m1,n))
c      
      call CSP_kern(n,k,t,yy,m1,po,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2,2)  ! CSP bv
      call fex(n,t,yy,ydot)
      call smult(m1,N,1,Br2,ydot,fr,m1,N,1)                         ! fr
      call smult(n,m1,1,Ar2,fr,temp,n,m1,1)                         ! ar.fr

      deallocate(Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,dummy,Ar1,As1,Br1,Bs1,fr,Ar2,As2,Br2,Bs2) 


cpd-----------------------------------------------

      do j=1,n
       temp(j)=temp(j)/dabs(rlmod(m2))              ! multiply with tau+1
      enddo  
c    
      do i=1,n
       yyer(i)=dabs(rlmod(m2)/rlmod(m1))*yy(i)+CSPatol(i)        !Error
      enddo
      m=0
      do j=1,n
       if(dabs(temp(j)).lt.yyer(j)) m=m+1    !|Sum alpha*f|<rtol*y+atol
      enddo

      if(m.eq.n) then
       noEM=noEM+1
       go to 200
      endif


      return
      end


cpd-------------------------------------------------------------------  
cpd-------------------------------------------------------------------  
cpd-------------------------------------------------------------------  
cpd-------------------------------------------------------------------  
c
c     Subrout for the calculation of the number of exhausted 
c     modes through the eigenvectors 
c
cpd------------------------------------------------------------------- 
      subroutine NEM_Eigen(n,k,t,yy,CSPrtol,CSPatol,rlmod,alpha,FiEig,noEM)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: t, yy(n),CSPrtol,CSPatol(n),rlmod(n),alpha(n,n),FiEig(1,n)
      integer, intent(out) :: noEM
      integer :: i,j,m1,m2,m,NoPh
      double precision :: yyer(n),temp(n),fi(n),ydot(n)

cpd-----------------------------------------------
      noEM=0
  200 continue
      m1=noEM+1
      m2=noEM+2                               !Next two modes
      do i=1,n
       temp(i)=0.0d0
      enddo
cpd-----------------------------------------------      ! construct Ar. Fr for the given m1
      if(rlmod(m1).ne.rlmod(m2)) then         ! REAL EIGENVALUES
       do i=1,m1
        do j=1,n
         temp(j)=temp(j)+alpha(j,i)*FiEig(1,i)    ! XTEMP(J)+ALPHA(J,I)*FI(1,I)
        enddo
       enddo
      else                                    ! COMPLEX EIGENVALUES
       do i=1,m1+1
        do j=1,n
         temp(j)=temp(j)+alpha(j,i)*FiEig(1,i)    ! XTEMP(J)+ALPHA(J,I)*FI(1,I)
        enddo       
       enddo
      endif
cpd-----------------------------------------------

      do j=1,n
       temp(j)=temp(j)/rlmod(m2)              ! multiply with tau+1
      enddo  
c    
      do i=1,n
!       yyer(i)=dabs(rlmod(m2)/rlmod(m1))*yy(i)+CSPatol        !Error
       yyer(i)=CSPrtol*yy(i)+CSPatol(i)        !Error
      enddo
      m=0
      do j=1,n
       if(dabs(temp(j)).lt.yyer(j)) m=m+1    !|tau alpha*f|<rtol*y+atol
      enddo

      if(m.eq.n) then
       if(rlmod(m1).ne.rlmod(m2)) then         ! REAL EIGENVALUES
        noEM=noEM+1
       else                                    ! COMPLEX EIGENVALUES
        noEM=noEM+2
       endif
       go to 200
      endif

c      noEM=1

      return
      end
