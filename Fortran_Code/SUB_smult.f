
      SUBROUTINE SMULT(N1,N2,N3,A,B,C,K1,K2,K3)
*
*     ==================================================================
*     - PURPOSE:
*      Perform a multiplication between two arrays A(N1,N2) and B(N2,N3). 
*      The resulting product is stored in a third  array C(N1,N3).
*     - INPUT:
*      N1   = the first leading DIMENSION of arrays A and C.
*              -data type: INTEGER scalar.
*      N2   = the second leading DIMENSION of array A and the first of B.
*              -data type: INTEGER scalar.
*      N3   = the second leading DIMENSION of arrays B and C.
*              -data type: INTEGER scalar.
*      A    = the left array of the product.
*              -data type: REAL array
*              -DIMENSION: A(N1,N2)
*      B    = the right array of the product.
*              -data type: REAL array
*              -DIMENSION: B(N2,N3)
*      K1   = the no of rows used of arrays A and C.
*              -data type: INTEGER scalar.
*      K2   = the no of A-columns and the no of B-rows used.
*              -data type: INTEGER scalar.
*      K3   = the no of columns used of arrays B and C.
*              -data type: INTEGER scalar.
*     - OUTPUT:
*      C    = the resulting array of the product.
*              -data type: REAL array
*              -DIMENSION: C(N1,N3)
*     
*     - INPUT/OUTPUT:
*     - CALLS:
*     - GENERAL DESCRIPTION:
*     - OTHER COMMENTS:
*      This code should be MIL-STD-1753 extended ANSI-77 compliant.
*     ==================================================================
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER(I-N) 
      DIMENSION A(N1,N2),B(N2,N3),C(N1,N3)
*
*     ==================================================================
*
      DO I=1,K1
         DO J=1,K3
            SUM=0.0D0
            DO K=1,K2
               SUM = SUM + A(I,K)*B(K,J)
            END DO
            C(I,J) = SUM
         END DO
      END DO
*
*     *----------------------------------------------------------------*
*     *----------------------------------------------------------------*
       RETURN
       END
