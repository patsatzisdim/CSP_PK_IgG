cpd   Sub for sorting array. Absolute values considered here.
cpd   Inputs:         x: array to be sorted
cpd                   m: size of the array
cpd   Outputs:        x: sorted array
cpd            indexSor: index trace of original vector  
      subroutine  SortAbs(x, m, indexSor)
      implicit none
      integer, intent(in) :: m
      double precision, intent(inout) :: x (m)
      integer, intent (out) :: indexSor (m)
      double precision :: y (m)
      integer :: i, Loc
      

      do i = 1, m      
       Loc = minloc(dabs(x), m)   ! find minimum location
       y(m+1-i)= x(Loc)           ! store in the new output array. From the last
       indexSor(m+1-i)=Loc        ! keep the index
       x(Loc)=1.0d+40             ! large value for not taking it again
      enddo
      x(1:m)=y(1:m)               ! return it in a new array

      return
      end
