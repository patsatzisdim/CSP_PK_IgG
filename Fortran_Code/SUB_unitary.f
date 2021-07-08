*     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*     ||                                                              ||
*     ||                    S-STEP Code Project.                      ||
*     ||                                                              ||
*     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*
      SUBROUTINE UNITARY(N1,In,K1) 
*
*     ==================================================================
*     - PURPOSE:
*      This subroutine forms a unitary matrix < In > of dimensions n*n.
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
*       This code should be MIL-STD-1753 extended ANSI-77 compliant.
*     ==================================================================
      double precision In(N1,N1) 
*
*     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*
C     ---------------------------------------------- I(i,j)=0 --- 
      do 100 i=1,K1 
           do 70 j=1,K1 
               In(i,j)=0.0d0 
 70        continue 
C          ----------------------------------------- I(i,i)=1 --- 
           In(i,i)=1.0d0 
 100  continue
*
*     *----------------------------------------------------------------*
*     *----------------------------------------------------------------*
       RETURN
       END
