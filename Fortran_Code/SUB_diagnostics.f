cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS THE SUBROUTINES TO CALCULATE THE CSP DIAGNOSTICS
c
c                     DIAGNOSTICS (Diag) 
c                     CSP Pointer separately (CSP_po)
c                     Printing diagnostics (Print_Diag) and also sorted (smart_print)  
c                     Printing the characteristic mode diagnostics only (charDiags_print)                    
c                 
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cpd------------------------------------------------------------------
cpd   Diagnostics subroutine. Printings are passed out at the main code.
cpd------------------------------------------------------------------
c     Inputs:     n: No. species
c                 k: No. reactions
c                 NoEM: No. Exhausted Modes
c                 t: time point
c                 y: concentrations
c                 alpha: right eigenvectors     (practically inputs but inout intent is for doing calculations with them)
c                 betta: left eigenvectors
c                 reim: real and imaginary part of eigenvalues
c
c     Outputs:    po: CSP Pointer (whole matrix)
c                 tpi: Timescale Participation Index
c                 bsr: Amplitude Participation Index
c                 impi: Importance Index
c                 fip: Amplitudes (only positive with possible change on the signs of eigenvectors)
c                 BS: Betta dot Stoichiometry matrix, needed as additional for API
c
c
      subroutine Diag(n,k,NoEM,t,yy,alpha,betta,reim,po,tpi,bsr,impi,fip,BS)
      implicit none
      integer, intent(in) :: n,k,NoEM
      integer :: i,j,ii,count,countSP,count1
      double precision, intent(in) :: t, yy(n), reim(n,2)
      double precision, intent(inout) ::  alpha(n,n), betta(n,n)
      double precision, intent(out) :: po(n,n),tpi(n,k),bsr(n,k),fip(n),impi(n,k), BS(n,k)
      double precision :: dpn(n), st(n,k), gR(k,n), gradSRk(n,n),tempS(n,1), tempGR(1,n), tempB(1,n), BgSRk(1,n), sum1, 
     - BgSRkA(n,k), tsum1, tsum2, sum2(n), RR(k), FI(4,n),tempSR(n,k), dumBSR(n,k), alphaS(n,n-NoEM), bettaS(n-NoEM,n),
     - asbs(n,n), asbsh(n,k), sumPo(n), sumdum
      
      count1=0

cpd
cpd   CSP Po through eigencvectors
      call CSP_po(n,reim,alpha,betta,po,dpn)     

cpd-----------------------------------------------------
cpd
cpd   B.S through eigenvectors
      call stoic(n,k,st)
      call smult(n,n,k,betta,st,BS,n,n,k)

cpd-----------------------------------------------------
cpd
cpd   CSP TPI through eigencvectors
      call stoic(n,k,st)
      call gradR(n,t,yy,k,gR)
cpd------------------------------
cpd   Loop over the reactions
      do j=1,k
cpd   Construct the grad(s_i R^i)
cpd   First s_i 
       do i=1,n
        tempS(i,1)=st(i,j)
       enddo       
cpd   Second grad(R^i)
       do i=1,n
        tempGR(1,i)=gR(j,i)
       enddo
c       gradSRk=matmul(tempS,tempGR)
       call smult(n,1,n,tempS,tempGR,gradSRk,n,1,n)  !  grad(S_i R^i) 
      countSP=0     ! mode counter for imaginary parts
cpd------------------------------
cpd   Loop over the species          
       do i=1,n    
        countSP=countSP+1
        if (reim(countSP,2).eq.0.0d0) then                  ! REAL EIGENVALUES
         do ii=1,n
          tempB(1,ii)=betta(countSP,ii)
         enddo
         call smult(1,n,n,tempB,gradSRk,BgSRk,1,n,n)  ! betta^j S_i gradR^i
         sum1=0.0d0  
         do ii=1,n
          sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP)
         enddo
         BgSRkA(countSP,j)=sum1                       ! betta^j S_i gradR^i alpha_j= c^n_k
cpd------------------------------
        else                                          ! COMPLEX EIGENVALUES
         tsum1=0.0d0
         tsum2=0.0d0                 
c                 
          do ii=1,n
           tempB(1,ii)=betta(countSP,ii)              ! b^1
          enddo
          call smult(1,n,n,tempB,gradSRk,BgSRk,1,n,n) ! b^1 J 
          sum1=0.0D0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP)    ! b^1 J a_1
          enddo
          tsum1=tsum1+sum1
          sum1=0.0d0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP+1)  ! b^1 J a_2
          enddo
          tsum2=tsum2+sum1
c         
          BgSRk(:,:)=0.0d0 
          do ii=1,n
           tempB(1,ii)=betta(countSP+1,ii)
          enddo
          call smult(1,n,n,tempB,gradSRk,BgSRk,1,n,n)
          sum1=0.0D0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP+1) ! b^2 J a_1
          enddo
          tsum2=tsum2-sum1                           ! b^1 J a_2 - b^2 J a_1
c          
          sum1=0.0D0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP)   ! b^2 J a_2
          enddo
          tsum1=tsum1+sum1                           ! b^1 J a_1 + b^2 J a_2

          BgSRkA(countSP,j)=0.5d0*tsum1              ! betta^j S_i gradR^i alpha_j= c^1_k
          BgSRkA(countSP+1,j)=0.5d0*tsum2            ! betta^j S_i gradR^i alpha_j= c^2_k
c
          countSP=countSP+1      ! if this is the first complex eigen from the pair, just pass, if not continue    
        endif 
        if(countSP.eq.n) exit
       enddo
cpd   End of loop over the species
      enddo
cpd------------------------------      
cpd   End of loop over the reactions   	
cpd   Now continue in construction of TPI     
      do i=1,n
       sum2(i)=0.0d0
       do j=1,k
        sum2(i)=sum2(i)+dabs(BgSRkA(i,j))
       enddo
      enddo
      do i=1,n
       do j=1,k
       	tpi(i,j)=BgSRkA(i,j)/sum2(i)
       enddo
      enddo

cpd--------------------------------------------------
cpd
cpd   CSP II through eigencvectors
      call stoic(n,k,st)
      call rates(n,t,yy,k,RR)
      do i=1,n-NoEM
       alphaS(:,i)=alpha(:,i+NoEM)
       bettaS(i,:)=betta(i+NoEM,:)
      enddo
      asbs(:,:)=0.0d0  
      asbsh(:,:)=0.0d0   
cpd------------------------------
      call smult(n,n-NoEM,n,alphaS,bettaS,asbs,n,n-NoEM,n)
      call smult(n,n,k,asbs,st,asbsh,n,n,k)
      do i=1,n
       do j=1,k
        impi(i,j)=asbsh(i,j)*RR(j)
       enddo
       tsum1=0.0d0
       do j=1,k
        tsum1=tsum1+dabs(impi(i,j))+1.0d-20 ! handle with care
       enddo
       do j=1,k
        impi(i,j)=impi(i,j)/tsum1
       enddo
      enddo
cpd--------------------------------------------------
cpd            
cpd   API construction through eigenvectors 
cpd   Amplitude construction
      call stoic(n,k,st)
      call rates(n,t,yy,k,RR)

      do i=1,n
       do j=1,k
       	tempSR(i,j)=st(i,j)*RR(j)     ! matrix of s^j_i R^j
       enddo
      enddo 
      call smult(n,n,k,betta,tempSR,BSR,n,n,k)    ! b^i s^j_i R^j 
cpd   Amplitudes to be positive
      do i=1,n                       ! Loop over species
       FI(1,i)=0.0d0
       FI(2,i)=0.0d0
       FI(3,i)=0.0d0
       FI(4,i)=0.0d0
       do j=1,k                            
       	FI(1,i)=FI(1,i)+BSR(i,j)             ! Sum over the reactions     
       	FI(2,i)=FI(2,i)+dabs(BSR(i,j))       ! abs
       	if(BSR(i,j).gt.0.0d0) FI(3,i)=FI(3,i)+BSR(i,j)      ! positive
       	if(BSR(i,j).lt.0.0d0) FI(4,i)=FI(4,i)+BSR(i,j)      ! negative
       enddo
       sum1=1.0d0/FI(2,i)                       ! 1/sum(abs)

       do j=1,k
       	dumBSR(i,j)=BSR(i,j)                       ! dummy to save the values
        BSR(i,j)=BSR(i,j)*sum1
       enddo	
      enddo   						 ! End of loop over species
c      fip(1:n)=FI(1,1:n)
      do i=1,n                       ! New loop over species after calculations
       do j=1,k
        BSR(i,j)=BSR(i,j)*dsign(1.d0,FI(1,i))
        BS(i,j)=BS(i,j)*dsign(1.0d0,FI(1,i))
        dumBSR(i,j)=dumBSR(i,j)*dsign(1.d0,FI(1,i))   ! take the sign of the sum of BSR
       enddo
       do j=1,n
       	alpha(j,i)=alpha(j,i)*dsign(1.d0,FI(1,i))     ! change the sign of eigenvectors
        betta(i,j)=betta(i,j)*dsign(1.d0,FI(1,i))
       enddo
       if(FI(1,i).lt.0.0d0) then
       	sum1=FI(3,i)
		    FI(3,i)=-FI(4,i)
		    FI(4,i)=-sum1
       endif
	     FI(1,i)=FI(1,i)*dsign(1.d0,FI(1,i))            ! F^i positive
       fip(i)=FI(1,i)
      enddo                          ! End of loop over species
      return 
      end
cpd--------------------------------------------------
c
c      Pointer is separate because it is called on SUB_main 
c
cpd--------------------------------------------------
c     Inputs:     n: No. species
c                 reim: real and imaginary part of eigenvalues (different calculation if complex)
c                 alpha: right eigenvectors
c                 betta: left eigenvectors
c
c     Outputs:    po: CSP pointer (the whole matrix)
c                 dpn: CSP pointer (only the diagonal entries of po)
c
      subroutine CSP_po(n,reim,alpha,betta,po,dpn)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: reim(n,2), alpha(n,n), betta(n,n)
      double precision, intent(out) :: po(n,n), dpn(n)
      integer :: i,count,j
      double precision :: sum1,sum2

cpd--------------------------------------------------
cpd   CSP Pointer through eigenvectors
cpd   Watch out for complex eigenvalues
      i=1
      do while (i.le.n)
       if (reim(i,2).eq.0.0d0) then     ! REAL and COMPLEX EIGENVALUES both of them are fine by this version.
        do j=1,n
         po(i,j)=betta(i,j)*alpha(j,i)  ! Pointer through eigenvectors
         if(i.eq.j) dpn(i)=po(i,j)      ! Diagonal element
        enddo
       else                             ! COMPLEX EIGENVALUES
        do j=1,n
         po(i,j)=0.5d0*(betta(i,j)*alpha(j,i)+betta(i+1,j)*alpha(j,i+1))
         po(i+1,j)=0.5d0*(betta(i+1,j)*alpha(j,i)-betta(i,j)*alpha(j,i+1))
        enddo
        i=i+1                           ! if this is the first complex eigen from the pair, continue without taking its conjucate pair        
       endif   
       i=i+1                   
      enddo
cpd--------------------------------------------------
      return 
      end

cpd--------------------------------------------------
c
c      Print the original diagnostics in different files
c
cpd--------------------------------------------------
cpd--------------------------------------------------
      subroutine Print_Diag(n,k,t,yy,po,tpi,api,impi,fip,unNo1,unNo2,
     - unNo3,unNo4,unNo5)
      implicit none
      integer,intent(in) :: n, k, unNo1, unNo2, unNo3, unNo4, unNo5
      double precision, intent(in) :: t, yy(n), po(n,n), api(n,k),
     - tpi(n,k), impi(n,k), fip(n)
      integer :: i

      do i=1,n
       write(unNo1,36) t, 'Spec  ',i, po(i,1:n) 
       write(unNo2,33) t, 'Mode  ',i, tpi(i,1:k)
       write(unNo3,33) t, 'Mode  ',i, api(i,1:k)
       write(unNo4,33) t, 'Spec  ',i, impi(i,1:k)
      enddo
  
      write(unNo5,30) t, fip(1:n)+1.0d-15

   30 format(E18.11,4('    ',E25.16))
   33 format(E18.11,a8,i3,5('    ',E25.16))
   36 format(E18.11,a8,i3,4('    ',E25.16))         

      return
      end 
cpd--------------------------------------------------
c
c      Print diagnostics seperately
c
cpd--------------------------------------------------
cpd--------------------------------------------------
      subroutine Print_DiagS(n,k,t,po,tpi,api,impi,unNo1,
     - unNo2,unNo3,unNo4)
      implicit none
      integer,intent(in) :: n, k, unNo1, unNo2, unNo3, unNo4
      double precision, intent(in) :: t, po(n,n), api(n,k),
     - tpi(n,k), impi(n,k)
      integer :: i

      do i=1,n
       write(unNo1+i-1,34) t, tpi(i,1:k)
       write(unNo2+i-1,34) t, api(i,1:k)
       write(unNo3+i-1,34) t, impi(i,1:k)
       write(unNo4+i-1,37) t, po(i,1:n) 
      enddo

   33 format(E18.11,a8,i3,5('    ',E25.16))
   36 format(E18.11,a8,i3,4('    ',E25.16))
   34 format(E18.11,5('    ',E25.16))            
   37 format(E18.11,4('    ',E25.16))
      return
      end         
cpd--------------------------------------------------
c
c      Print the sorted diagnostics in excel format
c
cpd--------------------------------------------------
cpd--------------------------------------------------
      subroutine smart_print(n,k,t,yy,PoS,tpiS,apiS,impiS,PoI,tpiI,
     - apiI,impiI,BS,unNoStart,unNoEnd)
      implicit none
      integer,intent(in) :: n, k, unNoStart, unNoEnd, PoI(n,n), 
     - tpiI(n,k), apiI(n,k), impiI(n,k)
      double precision, intent(in) :: t, yy(n), PoS(n,n), apiS(n,k),
     - tpiS(n,k), impiS(n,k), BS (n,k) 
      integer :: i,ii,j

      do ii=unNoStart,unNoEnd         !number of modes
       j=ii-unNoStart+1 
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP pointer'
       do i=1,n
        write(ii,40) t,i,PoI(j,i),PoS(j,i)
       enddo
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP API'
       do i=1,k 
        write(ii,41) t,i,apiI(j,i),apiS(j,i),BS(j,apiI(j,i))
       enddo
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP II'       
       do i=1,k
        write(ii,40) t,i,impiI(j,i),impiS(j,i)
       enddo
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP TPI'       
       do i=1,k
        write(ii,41) t,i,tpiI(j,i),tpiS(j,i),BS(j,tpiI(j,i))
       enddo
      enddo


   40 format(E18.11,i5,i5,F12.3)     
   41 format(E18.11,i5,i5,F12.3,F17.8)

      return
      end 
cpd--------------------------------------------------
c
c      for printing several of the diagnostics of the charactristic mode
c
cpd--------------------------------------------------
cpd--------------------------------------------------
      subroutine charDiags_print(n,k,t,yy,NoEM,po,tpi,api,impi,indW1,indW2,indW3)
      integer,intent(in) :: n,k,NoEM,indW1,indW2,indW3
      double precision, intent(in) :: t, yy(n),po(n,n),tpi(n,k),api(n,k),impi(n,k)


      write(indW1,51) T,po(NoEM+1,1:n)
      write(indW2,52) T,tpi(NoEM+1,1:k)
      write(indW3,52) T,api(NoEM+1,1:k)

   
   51 format(E18.11,4('    ',E15.6))  
   52 format(E18.11,5('    ',E15.6))

      return
      end
