*     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*     ||                                                              ||
*     ||                    S-STEP Code Project.                      ||
*     ||                                                              ||
*     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*
      SUBROUTINE SINVE(N,A,B,M,IPVT,INFO)
*
*     ==================================================================
*     - PURPOSE:
*         Find the inverse of A (MxM) by solving A*B = I 
*         using LINPACK routines
*     
*         Note that the matrix A is replaced by its LU factorization.
*     
*     - HISTORY: 
*        ddmmyy  Programmer  Task
*        ------  ----------  ----
*        141099  Th71        V3-Public release
*                
*     - ARGUMENTS:
*     - INPUT:
*     - OUTPUT:
*     - INPUT/OUTPUT:
*     - CALLS:
*     - GENERAL DESCRIPTION:
*     - OTHER COMMENTS:
*      This code should be MIL-STD-1753 extended ANSI-77 compliant.
*     ==================================================================
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),  INTEGER (I-N)
      DIMENSION A(N,N),B(N,N),IPVT(N)
*
*     ==================================================================
*
C 
      INFO=0
C     ---------------------------------------------------- set B = I --- 
      CALL UNITARY(N,B,M) 
      CALL DGEFA(A,N,M,IPVT,INFO)
      IF (INFO .NE. 0) THEN
         WRITE(6, 37) INFO
cpd
c           write(*,*) a(INFO,INFO)
c           stop
cpd
      END IF
C 
      DO I=1,M
         CALL DGESL(A,N,M,IPVT,B(1,I),0)
      END DO
C 
*
*     *----------------------------------------------------------------*
 37   FORMAT('Singular matrix in SUBROUTINE DGEFA, row ', I3)
*
*     *----------------------------------------------------------------*
       RETURN
       END
