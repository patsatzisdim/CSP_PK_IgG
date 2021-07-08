*     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*     ||                                                              ||
*     ||                    S-STEP Code Project.                      ||
*     ||                                                              ||
*     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*
      SUBROUTINE EIGEN(RJAC,ALPHA,BETTA,RLMOD,REIM)
*
*     ==================================================================
*     - PURPOSE:
*       Given the jacobian  [RJAC(NSPEC,NSPEC)] we evaluate:
*       moduli of eigenvalues         [ RLMOD(NSPEC)       ]
*       real      part of eigenvalues [  REIM(NSPEC,1)     ]
*       imaginary part of eigenvalues [  REIM(NSPEC,2)     ]
*       left  eigenvectors            [ BETTA(NSPEC,NSPEC) ]
*       right eigenvectors            [ ALPHA(NSPEC,NSPEC) ]
*                                                       J=1,NSPEC
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N) 
      PARAMETER (NSPECMAX=4,NSPEC=4)
*
      DIMENSION RJAC(NSPECMAX,NSPECMAX)   ! JACOBIAN OF G(Z) WITH RESPECT TO Z
      DIMENSION ALPHA(NSPECMAX,NSPECMAX)  ! RIGHT EIGENVECTORS
      DIMENSION BETTA(NSPECMAX,NSPECMAX)  ! LEFT  EIGENVECTORS
      DIMENSION REIM(NSPECMAX,2)          ! REAL AND IMAGINARY
                                          !   PARTS OF THE EIGENVALUES
      DIMENSION RLMOD(NSPECMAX)           ! MODULI OF THE EIGENVALUES
*
      DIMENSION EIGR(NSPECMAX)            ! REAL     
                                          !   PARTS OF THE EIGENVALUES
      DIMENSION EIGI(NSPECMAX)            ! IMAGINARY
                                          !   PARTS OF THE EIGENVALUES
      DIMENSION AC(NSPECMAX,NSPECMAX)     ! REAL AND IMAGINARY
                                          !   PARTS OF THE EIGENVECTORS
                                          !   Warning: it gets modified
                                          !   in the process.
      DIMENSION RWKN(NSPECMAX)            ! WORKING ARRAY FOR rg
      DIMENSION IWKN(NSPECMAX)            ! WORKING ARRAY FOR rg
*
*     ==================================================================
      IERR=0
C
C     ------------------------------------------------------------------
C     Find eigenvalues and right eigenvectors.
C     ------------------------------------------------------------------
C
C
      CALL RG(NSPECMAX,NSPEC,RJAC,EIGR,EIGI,1,AC,IWKN,RWKN,IERR)

C
C      WRITE(*,*) EIGI(1),EIGI(2),EIGI(3)
      IF (IERR.NE.0) THEN
        WRITE(*,*) '(IERR.NE.0)'
        STOP
      ENDIF
C
C     ------------------------------------------------------------------
C     Transfer eigenvalues into REIM.
C     ------------------------------------------------------------------
C
      DO I=1,NSPEC
         REIM(I,1)=EIGR(I)
         REIM(I,2)=EIGI(I)
      END DO
C
C     ------------------------------------------------------------------
C     Transfer eigenvalue modulus into real RLMOD
C     ------------------------------------------------------------------
C
      DO I=1,NSPEC
         RLMOD(I)=DSQRT(REIM(I,1)**2+REIM(I,2)**2)
      END DO
C
C     ------------------------------------------------------------------
C     Order eigenvalues and eigenvectors.
C     1. Large magnitude first.
C     ------------------------------------------------------------------
C
      DO LA=1,NSPEC-1
C
C        Find LA'th largest element in list
C
         XX=-1.0D0
         IMOVE=LA
         DO I=LA,NSPEC
            XY=RLMOD(I)
            IF (XY.GT.XX) THEN   
               XX=XY
               IMOVE=I
            END IF
         END DO
C
C        Swap with the element indexed LA, if necessary
C
         IF (LA.NE.IMOVE) THEN
            DO K=1,2
               COM=REIM(LA,K)
               REIM(LA,K)=REIM(IMOVE,K)
               REIM(IMOVE,K)=COM
            END DO
            COM = RLMOD(LA)
            RLMOD(LA) = RLMOD(IMOVE)
            RLMOD(IMOVE) = COM
            DO J=1,NSPEC
               COM=AC(J,LA)
               AC(J,LA)=AC(J,IMOVE)
               AC(J,IMOVE)=COM
            END DO
         END IF
      END DO
C CPatsatzis 03112019
C c
C c     For complex eigenvalues sort first the positive imaginary part and then the negative
C c
C       DO 14 I=1,NSPEC-1
C          XY1=REIM(I,2)
C          XY2=REIM(I+1,2)
C          IF ((XY1.EQ.0.0D0).or.(XY2.EQ.0.0D0)) GO TO 14
C          XZ1=DABS(XY1)
C          XZ2=DABS(XY2)
C          XY12=DABS((XZ1-XZ2)/XZ1)
C          IF (XY12.GT.1.0D-7) GO TO 14
C          REIM(I,2)=DABS(REIM(I,2))
C          REIM(I+1,2)=-DABS(REIM(I+1,2))
C    14 CONTINUE

c          if(t.gt.0.0235d0.and.t.lt.0.0237d0) then
c           write(*,*) t, reim(1:7,2)
c          endif      


C
C
C     ------------------------------------------------------------------
C     2. If complex pair then two ALPHA's are formed from the
C        real & imaginary parts of eigenvectors
C     ------------------------------------------------------------------
C
      DO 4 I=1,NSPEC-1
         XY1=REIM(I,2)
         XY2=REIM(I+1,2)
         IF ((XY1.EQ.0.0D0).or.(XY2.EQ.0.0D0)) GO TO 4
         XZ1=DABS(XY1)
         XZ2=DABS(XY2)
         XY12=DABS((XZ1-XZ2)/XZ1)
         IF (XY12.GT.1.0D-7) GO TO 4
         DO J=1,NSPEC
            X1=AC(J,I)
            X2=AC(J,I+1)
            AC(J,I  )=X1
            AC(J,I+1)=X2			
         END DO
    4 CONTINUE
C
C     ------------------------------------------------------------------
C     3. Transfer right eigenvectors into real ALPHA.
C     ------------------------------------------------------------------
C
      DO I=1,NSPEC
         DO J=1,NSPEC
            ALPHA(I,J)=AC(I,J)
         END DO
      END DO

    
C
C     ------------------------------------------------------------------
C     4. Normalize the ALPHA's with the largest element.
C     ------------------------------------------------------------------
C

      DO I=1,NSPEC
        XY1=REIM(I,2)
        IF (I.EQ.NSPEC) THEN
          IF (DABS(XY1).LT.1.0D-7) GO TO 7    !  the last eigenvalue is real
          GO TO 6                             !  the last eigenvalue is comlex
        ENDIF
        XY2=REIM(I+1,2)
        IF ((XY1.EQ.0.0D0).or.(XY2.EQ.0.0D0)) GO TO 5
        XZ1=DABS(XY1)
        XZ2=DABS(XY2)
        XY12=DABS((XZ1-XZ2)/XZ1)
        IF (XY12.GT.1.0D-7) GO TO 6
        LPHAMAX=0
        ALPHAMAX=0.0D0
          DO J=1,NSPEC
            ALPHATMP=DABS(dsqrt(ALPHA(J,I)**2+ALPHA(J,I+1)**2))
            IF (ALPHATMP.GT.ALPHAMAX) LPHAMAX=J
            IF (ALPHATMP.GT.ALPHAMAX) ALPHAMAX=ALPHATMP
          END DO
          coe1=1.0d0/(ALPHA(LPHAMAX,I)**2+ALPHA(LPHAMAX,I+1)**2)
          coe2=ALPHA(LPHAMAX,I)
          coe3=ALPHA(LPHAMAX,I+1)
        DO J=1,NSPEC
           coe4=coe1*(coe2*ALPHA(J,I)+coe3*ALPHA(J,I+1))
           coe5=coe1*(coe2*ALPHA(J,I+1)-coe3*ALPHA(J,I))
           ALPHA(J,I)=coe4
           ALPHA(J,I+1)=coe5
        END DO
        GO TO 6
    5   CONTINUE
           IF ((XY1.NE.0.0D0).and.(XY2.EQ.0.0D0)) GO TO 6         ! if the next one is real 
    7   CONTINUE       
           ALPHAMAX=0.0D0
         DO J=1,NSPEC
c           ALPHATMP=DABS(ALPHA(J,I))
           ALPHATMP=ALPHA(J,I)
c           IF (ALPHATMP.GT.ALPHAMAX) ALPHAMAX=ALPHATMP
           IF (dabs(ALPHATMP).GT.dabs(ALPHAMAX)) ALPHAMAX=ALPHATMP
         END DO
         DO J=1,NSPEC
           ALPHA(J,I)=ALPHA(J,I)/ALPHAMAX
         END DO
    6   CONTINUE
        END DO

C
C     ------------------------------------------------------------------
C     Normalize the ALPHA's.
C     ------------------------------------------------------------------
C
c      DO I=1,NSPEC
c         ALPHAMAX=0.0D0
c         DO J=1,NSPEC
c            ALPHATMP=DABS(ALPHA(J,I))
c            IF (ALPHATMP.GT.ALPHAMAX) ALPHAMAX=ALPHATMP
c         END DO
c         DO J=1,NSPEC
c            ALPHA(J,I)=ALPHA(J,I)/ALPHAMAX
c         END DO
c      END DO
CC
c      DO I=1,NSPEC
c         XY1=REIM(I,2)
c         XY2=REIM(I+1,2)
c         IF ((XY1.EQ.0.0D0).or.(XY2.EQ.0.0D0)) GO TO 5
c         XZ1=DABS(XY1)
c         XZ2=DABS(XY2)
c         XY12=DABS((XZ1-XZ2)/XZ1)
c         IF (XY12.GT.1.0D-7) GO TO 6
c         ALPHAMAX=0.0D0
c         DO J=1,NSPEC
c              ALPHATMP=DABS(dsqrt(ALPHA(J,I)**2+ALPHA(J,I+1)**2))
c              IF (ALPHATMP.GT.ALPHAMAX) ALPHAMAX=ALPHATMP
c         END DO
c         DO J=1,NSPEC
c         ALPHA(J,I)=ALPHA(J,I)/ALPHAMAX
c         ALPHA(J,I+1)=ALPHA(J,I+1)/ALPHAMAX
c         END DO
c         GO TO 6
c    5    ALPHAMAX=0.0D0
c         DO J=1,NSPEC
c            ALPHATMP=DABS(ALPHA(J,I))
c            IF (ALPHATMP.GT.ALPHAMAX) ALPHAMAX=ALPHATMP
c         END DO
c         DO J=1,NSPEC
c            ALPHA(J,I)=ALPHA(J,I)/ALPHAMAX
c         END DO
c    6 Continue
c      END DO
C     ------------------------------------------------------------------
C     Compute the real BETTA's.  (BETTA=1/ALPHA)
C     ------------------------------------------------------------------
C
      DO I=1,NSPEC
         DO J=1,NSPEC
           AC(I,J) = ALPHA(I,J)
         END DO
      END DO
C
      IINFO=0
      CALL SINVE(NSPECMAX,AC,BETTA,NSPEC,IWKN,IINFO)
C
      IF (IINFO.NE.0) THEN
        WRITE (6,*) 'Execution stops in SUBROUTINE EIGEN'
        STOP
      END IF


C
*
*     *----------------------------------------------------------------*
*     *----------------------------------------------------------------*
       RETURN
       END

