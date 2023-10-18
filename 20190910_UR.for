C**********************************************************************
C     THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C     DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C**********************************************************************
      SUBROUTINE SPECTRAL(A,D,V)
C
C     THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C     A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C     THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C     EIGENVALUES IN ASCENDING ORDER, AND
C     A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C     EIGENVECTORS.
C**********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NP=3)
      DIMENSION D(NP),V(NP,NP)
      DIMENSION A(3,3),E(NP,NP)

      DO 2 I = 1,3
        DO 1 J= 1,3
          E(I,J) = A(I,J)
 1      CONTINUE
 2    CONTINUE

      CALL JACOBI(E,3,NP,D,V,NROT)
      CALL EIGSRT(D,V,3,NP)

      RETURN
      END

C**********************************************************************
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C     COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C     NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C     ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C     AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C     VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C     [V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C     [A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C     EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C     WHICH WERE REQUIRED.

C     THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C**********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX =100)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C     INITIALIZE [V] TO THE IDENTITY MATRIX

      DO 12 IP = 1,N
        DO 11 IQ = 1,N
          V(IP,IQ) = 0.D0
 11     CONTINUE
          V(IP,IP) = 1.D0
 12   CONTINUE

C     INITIALIZE [B] AND [D] TO THE DIAGONAL OF [A], AND Z TO ZERO.
C     THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C     IN EQUATION (11.1.14)

      DO 13 IP = 1,N
        B(IP) = A(IP,IP)
        D(IP) = B(IP)
        Z(IP) = 0.D0
 13   CONTINUE
C
      NROT = 0
      DO 24 I = 1,50

C     SUM OFF-DIAGONAL ELEMENTS

        SM = 0.D0
        DO 15 IP = 1, N-1
          DO 14 IQ = IP + 1, N
            SM = SM + DABS ( A(IP,IQ ))
 14       CONTINUE
 15     CONTINUE

C     IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C     WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C     UNDERFLOW.

        IF ( SM .EQ. 0.D0) RETURN
C
C     IF ( SM .LT. 1.0D-15) RETURN

C     IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C     |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, 
C     SEE EQUATION (11.1.25). THEREAFTER TRESH = 0.

        IF ( I .LT. 4) THEN
          TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          ENDIF
C
        DO 22 IP = 1, N-1
          DO 21 IQ = IP+1,N
            G = 100.D0*DABS(A(IP,IQ))

C     AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C     OFF-DIAGONAL ELEMENT IS SMALL.

            IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C     T = 1./(2.*THETA), EQUATION(11.1.10)

                T = A(IP,IQ)/H
           ELSE
                THETA = 0.5D0*H/A(IP,IQ)
                T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
                IF (THETA .LT. 0.D0) T = -T
           ENDIF
           C = 1.D0/DSQRT(1.D0 + T**2)
           S = T*C
           TAU = S/(1.D0 + C)
           H = T*A(IP,IQ)
           Z(IP) = Z(IP) - H
           Z(IQ) = Z(IQ) + H
           D(IP) = D(IP) - H
           D(IQ) = D(IQ) + H
           A(IP,IQ) = 0.D0

C     CASE OF ROTATIONS 1 <= J < P

            DO 16 J = 1, IP-1
              G = A(J,IP)
              H = A(J,IQ)
              A(J,IP) = G - S*(H + G*TAU)
              A(J,IQ) = H + S*(G - H*TAU)
 16         CONTINUE

C     CASE OF ROTATIONS P < J < Q

            DO 17 J = IP+1, IQ-1
              G = A(IP,J)
              H = A(J,IQ)
              A(IP,J) = G - S*(H + G*TAU)
              A(J,IQ) = H + S*(G - H*TAU)
 17         CONTINUE

C     CASE OF ROTATIONS Q < J <= N

            DO 18 J = IQ+1, N
              G = A(IP,J)
              H = A(IQ,J)
              A(IP,J) = G - S*(H + G*TAU)
              A(IQ,J) = H + S*(G - H*TAU)
 18         CONTINUE
            DO 19 J = 1,N
              G = V(J,IP)
              H = V(J,IQ)
              V(J,IP) = G - S*(H + G*TAU)
              V(J,IQ) = H + S*(G - H*TAU)
 19         CONTINUE
            NROT = NROT + 1
              ENDIF
 21       CONTINUE
 22     CONTINUE

C     UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

        DO 23 IP = 1, N
          B(IP) = B(IP) + Z(IP)
          D(IP) = B(IP)
          Z(IP) = 0.D0
 23     CONTINUE
 24   CONTINUE

C     IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C     THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C     AND STOP.

      WRITE (80,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

      RETURN
      END

C**********************************************************************
      SUBROUTINE EIGSRT(D,V,N,NP)

C     GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C     JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
C     AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.

C     THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C**********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(NP),V(NP,NP)

      DO 13 I = 1,N-1
        K = I
        P = D(I)
        DO 11 J = I+1,N
          IF (D(J) .GE. P) THEN
          K = J
          P = D(J)
          END IF
 11     CONTINUE
        IF (K .NE. I) THEN
          D(K) = D(I)
          D(I) = P
          DO 12 J = 1,N
            P = V(J,I)
            V(J,I) = V(J,K)
            V(J,K) = P
 12       CONTINUE
        ENDIF
 13   CONTINUE

      RETURN
      END

C**********************************************************************
C	    THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C**********************************************************************
      SUBROUTINE ZEROV(V,SIZE)

C	    THIS SUBROUTINE STORES THE ZERO VECTOR IN A VECTOR V
C	    OF SIZE SIZE.
C**********************************************************************

      INTEGER SIZE
      REAL*8 V(0:SIZE-1)

      DO 1 I=0,SIZE
        V(I) = 0.0D0
 1    CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE ZEROM(A)
C
C	    THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

      REAL*8 A(3,3)

      DO 2 I=1,3
        DO 1 J=1,3
          A(I,J) = 0.0D0
 1      CONTINUE
 2    CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE ONEM(A)

C	    THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	    3 BY 3 MATRIX [A]
C**********************************************************************

      REAL*8 A(3,3)

      DO 2 I=1,3
      	DO 1 J=1,3
          IF (I .EQ. J) THEN
            A(I,J) = 1.0D0
          ELSE
            A(I,J) = 0.0D0
          ENDIF
 1      CONTINUE   
 2    CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MTRANS(A,ATRANS)
 
C	    THIS SUBROUTINE CALCULATES THE TRANSPOSE OF A 3 BY 3 
C	    MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

      REAL*8 A(3,3),ATRANS(3,3)

      DO 2 I=1,3
        DO 1 J=1,3
          ATRANS(J,I) = A(I,J)
 1      CONTINUE
 2    CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MPROD(A,B,C)
 
C 	  THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	  AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

      REAL*8 A(3,3),B(3,3),C(3,3)

      DO 3 I = 1, 3
        DO 2 J = 1, 3
          C(I,J) = 0.D0
          DO 1 K = 1, 3
            C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
 1        CONTINUE
 2      CONTINUE
 3    CONTINUE
      
      RETURN
      END

C**********************************************************************
      SUBROUTINE MPROD4(A,B,C)
 
C	    THIS SUBROUTINE MULTIPLIES TWO 4 BY 4 MATRICES [A] AND [B],
C 	  AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

      REAL*8 A(4,4),B(4,4),C(4,4)

      DO 3 I = 1, 4
        DO 2 J = 1, 4
          C(I,J) = 0.D0
          DO 1 K = 1, 4
            C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
 1        CONTINUE
 2      CONTINUE
 3    CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE DOTPM(A,B,C)

C	    THIS SUBROUTINE CALCULATES THE SCALAR PRODUCT OF TWO
C	    3 BY 3 MATRICES [A] AND [B] AND STORES THE RESULT IN THE
C	    SCALAR C.
C**********************************************************************

      REAL*8 A(3,3),B(3,3),C

      C = 0.D0
      DO 2 I = 1,3
        DO 1 J = 1,3
          C = C + A(I,J)*B(I,J)
 1      CONTINUE
 2    CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MDET(A,DET)
 
C 	  THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	  OF A 3 BY 3 MATRIX [A].
C**********************************************************************

      REAL*8  A(3,3), DET

      DET = A(1,1)*A(2,2)*A(3,3) 
     +      + A(1,2)*A(2,3)*A(3,1)
     +      + A(1,3)*A(2,1)*A(3,2)
     +      - A(3,1)*A(2,2)*A(1,3)
     +      - A(3,2)*A(2,3)*A(1,1)
     +      - A(3,3)*A(2,1)*A(1,2)

      RETURN
      END

C**********************************************************************
      SUBROUTINE M3INV(A,AINV)

C 	  THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	      [A] AND PLACES THE RESULT IN [AINV]. 
C 	  IF DET(A) IS ZERO, THE CALCULATION
C 	  IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

      REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	    A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	    DET		-- THE COMPUTED DETERMINANT OF [A].
C	    ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C	    		   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C	    		   IS CALLED THE COFACTOR OF A(I,J).
C	    AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C	    		   OBTAINED BY REPLACING EACH ELEMENT OF
C	    		   [A] BY ITS COFACTOR, AND THEN TAKING
C	    		   TRANSPOSE OF THE RESULTING MATRIX.
C	    AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C	    		   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

      CALL MDET(A,DET)
      IF ( DET .EQ. 0.D0 ) THEN
        WRITE(80,100)
        WRITE(*,100)
        STOP
      ENDIF
      CALL MCOFAC(A,ACOFAC)
      CALL MTRANS(ACOFAC,AADJ)
      DO 2 I = 1,3
        DO 1 J = 1,3
          AINV(I,J) = AADJ(I,J)/DET
 1      CONTINUE
 2    CONTINUE

 100  FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +       10X,'PROGRAM TERMINATED')

      RETURN
      END

C**********************************************************************
      SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	  THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	  AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

      REAL*8  A(3,3), ACOFAC(3,3)

      ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
      ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
      ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
      ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
      ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
      ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

      RETURN
      END

C**********************************************************************
      SUBROUTINE INVAR(A,IA,IIA,IIIA)

C	    THIS SUBROUTINE CALCULATES THE PRINCIPAL INVARIANTS 
C	    IA, IIA, IIIA OF A TENSOR [A].
C**********************************************************************

      REAL*8 A(3,3), AD(3,3),AD2(3,3), DETA, IA,IIA,IIIA

      DO 2 I=1,3
        DO 1 J=1,3
          AD(I,J) = A(I,J)
 1      CONTINUE
 2    CONTINUE
      
      IA = AD(1,1) + AD(2,2) + AD(3,3)

C	    CALCULATE THE SQUARE OF [AD]

      CALL MPROD(AD,AD,AD2)
      IIA =0.5D0 * ( IA*IA - ( AD2(1,1) + AD2(2,2) + AD2(3,3) ) )

      CALL  MDET(AD,DETA)
      IIIA = DETA

      RETURN
      END

C**********************************************************************
      SUBROUTINE TRACEM(A,TRA)

C	    THIS SUBROUTINE CALCULATES THE TRACE OF A 3 BY 3 MATRIX [A]
C	    AND STORES THE RESULT IN THE SCALAR TRA
C**********************************************************************

      REAL*8 A(3,3),TRA

      TRA = A(1,1) + A(2,2) + A(3,3)

      RETURN 
      END

C**********************************************************************
      SUBROUTINE DEVM(A,ADEV)

C	    THIS SUBROUTINE CALCULATES THE DEVIATORIC PART OF A
C	    3 BY 3 MATRIX [A]
C**********************************************************************

      REAL*8 A(3,3),TRA,ADEV(3,3),IDEN(3,3)

      CALL TRACEM(A,TRA)
      CALL ONEM(IDEN)
      CALL ZEROM(ADEV)

      DO 2 I = 1,3
        DO 1 J = 1,3
          ADEV(I,J) = A(I,J) - (1.0D0/3.0D0)*TRA*IDEN(I,J)
 1      CONTINUE
 2    CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE EQUIVS(S,SB)

C	    THIS SUBROUTINE CALCULATES THE EQUIVALENT SHEAR STRESS SB
C	    CORRESPONDING TO A 3 BY 3 STRESS MATRIX [S]
C**********************************************************************

      REAL*8 S(3,3),SDEV(3,3),SDOTS,SB

      SB = 0.D0
      SDOTS = 0.D0

      CALL DEVM(S,SDEV)
      CALL DOTPM(SDEV,SDEV,SDOTS)
C      SB = DSQRT(1.5D0*SDOTS)
      SB = DSQRT(0.5D0*SDOTS)
      RETURN
      END
C***********************************************************************
      SUBROUTINE PRESSURE(A,PRA)
C
C     THIS SUBROUTINE CALCULATES THE MEAN NORMAL PRESSURE
C     OF A 3 BY 3 MATRIX [A]
C     AND STORES THE RESULT IN THE SCALAR PRA
C ----------------------------------------------------------------------
C       VARIABLES
C
      REAL*8 A(3,3),PRA

      PRA = -(1.D0 / 3.D0)*( A(1,1) + A(2,2) + A(3,3) )

      RETURN 
      END

C***********************************************************************
      SUBROUTINE MAGM(A,MAGA)
C
C     THIS SUBROUTINE CALCULATES THE MAGNITUDE OF A 3 BY 3 MATRIX [A]
C     THE MAGNITUDE IS DEFINED AS MAGA = SQRT(A(i,j)A(j,i))
C     SUM ON REPEATED INDICIES
C
C***********************************************************************
      REAL*8 A(3,3),MAGA,ADDA
      ADDA = 0.0D0
      DO 2 I = 1,3
        DO 1 J = 1,3
        ADDA = ADDA + A(I,J)*A(J,I) 
 1      CONTINUE
 2    CONTINUE
      MAGA = DSQRT(ADDA)
      RETURN
      END

C**********************************************************************
      SUBROUTINE PRTMAT(A,M,N)
C
C     THIS SUBROUTINE PRINTS A MATRIX [A] WITH DIMENSION M BY N
C       TO IO 80 and TERMINAL
C**********************************************************************

      INTEGER M,N
      REAL*8 A(M,N)  

      DO 10 K=1,M
C        WRITE(80,'(2X,6E12.4,2X)') (A(K,L), L=1,N)
        WRITE(6,'(2X,6E12.4,2X)') (A(K,L), L=1,N)
 10   CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE PRTVEC(A,M)
C
C     THIS SUBROUTINE PRINTS A VECTOR [A] WITH DIMENSION M
C       TO IO 80 and TERMINAL
C**********************************************************************

      INTEGER M
      REAL*8 A(M)  

      WRITE(80,'(2X,6E12.4,2X)') (A(K), K=1,M)
      WRITE(*,'(2X,6E12.4,2X)') (A(K), K=1,M)
      RETURN
      END

C**********************************************************************
      SUBROUTINE M3INV_SPFAST(A,AINV)

      REAL*4 A(3,3), AINV(3,3)
      REAL*4 DET

      DET = A(1,1)*(A(3,3)*A(2,2) - A(3,2)*A(2,3)) -
     1      A(2,1)*(A(3,3)*A(1,2) - A(3,2)*A(1,3)) + 
     2      A(3,1)*(A(2,3)*A(1,2) - A(2,2)*A(1,3))

      AINV(1,1) =  (1.0/DET)*(A(3,3)*A(2,2) - A(3,2)*A(2,3))
      AINV(1,2) = -(1.0/DET)*(A(3,3)*A(1,2) - A(3,2)*A(1,3))
      AINV(1,3) =  (1.0/DET)*(A(2,3)*A(1,2) - A(2,2)*A(1,3))

      AINV(2,1) = -(1.0/DET)*(A(3,3)*A(2,1) - A(3,1)*A(2,3))
      AINV(2,2) =  (1.0/DET)*(A(3,3)*A(1,1) - A(3,1)*A(1,3))
      AINV(2,3) = -(1.0/DET)*(A(2,3)*A(1,1) - A(2,1)*A(1,3))

      AINV(3,1) =  (1.0/DET)*(A(3,2)*A(2,1) - A(3,1)*A(2,2))
      AINV(3,2) = -(1.0/DET)*(A(3,2)*A(1,1) - A(3,1)*A(1,2))
      AINV(3,3) =  (1.0/DET)*(A(2,2)*A(1,1) - A(2,1)*A(1,2))

      RETURN
      END

C**********************************************************************










