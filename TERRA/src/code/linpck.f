      SUBROUTINE SPOCO(A,LDA,N,RCOND,Z,INFO)
C***BEGIN PROLOGUE  SPOCO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=SINGLE PRECISION(SPOCO-S DPOCO-D CPOCO-C),
C             CONDITION NUMBER,LINEAR ALGEBRA,MATRIX,
C             MATRIX FACTORIZATION,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a real SYMMETRIC POSITIVE DEFINITE MATRIX
C            and estimates the condition number of the matrix.
C***DESCRIPTION
C
C     SPOCO factors a real symmetric positive definite matrix
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, SPOFA is slightly faster.
C     To solve  A*X = B , follow SPOCO by SPOSL.
C     To compute  INVERSE(A)*C , follow SPOCO by SPOSL.
C     To compute  DETERMINANT(A) , follow SPOCO by SPODI.
C     To compute  INVERSE(A) , follow SPOCO by SPODI.
C
C     On Entry
C
C        A       REAL(LDA, N)
C                the symmetric matrix to be factored.  Only the
C                diagonal and upper triangle are used.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix  R  so that  A = TRANS(R)*R
C                where  TRANS(R)  is the transpose.
C                The strict lower triangle is unaltered.
C                If  INFO .NE. 0 , the factorization is not complete.
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.  If INFO .NE. 0 , RCOND is unchanged.
C
C        Z       REAL(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                If  INFO .NE. 0 , Z  is unchanged.
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     LINPACK SPOFA
C     BLAS SAXPY,SDOT,SSCAL,SASUM
C     Fortran ABS,MAX,REAL,SIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SASUM,SAXPY,SDOT,SPOFA,SSCAL
C***END PROLOGUE  SPOCO
      INTEGER LDA,N,INFO
      REAL A(LDA,1),Z(1)
      REAL RCOND
C
      REAL SDOT,EK,T,WK,WKM
      REAL ANORM,S,SASUM,SM,YNORM
      INTEGER I,J,JM1,K,KB,KP1
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
C***FIRST EXECUTABLE STATEMENT  SPOCO
      DO 30 J = 1, N
         Z(J) = SASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + ABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = MAX(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL SPOFA(A,LDA,N,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE TRANS(R)*W = E
C
         EK = 1.0E0
         DO 50 J = 1, N
            Z(J) = 0.0E0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
            IF (ABS(EK-Z(K)) .LE. A(K,K)) GO TO 60
               S = A(K,K)/ABS(EK-Z(K))
               CALL SSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = ABS(WK)
            SM = ABS(WKM)
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
            KP1 = K + 1
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + ABS(Z(J)+WKM*A(K,J))
                  Z(J) = Z(J) + WK*A(K,J)
                  S = S + ABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*A(K,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
C
C        SOLVE R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 120
               S = A(K,K)/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL SAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
C
         YNORM = 1.0E0
C
C        SOLVE TRANS(R)*V = Y
C
         DO 150 K = 1, N
            Z(K) = Z(K) - SDOT(K-1,A(1,K),1,Z(1),1)
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 140
               S = A(K,K)/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/A(K,K)
  150    CONTINUE
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE R*Z = V
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 160
               S = A(K,K)/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL SAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
  180 CONTINUE
      RETURN
      END
      SUBROUTINE SPODI(A,LDA,N,DET,JOB)
C***BEGIN PROLOGUE  SPODI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2B1B,D3B1B
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=SINGLE PRECISION(SPODI-S DPODI-D CPODI-C),
C             DETERMINANT,INVERSE,LINEAR ALGEBRA,MATRIX,
C             MATRIX FACTORIZATION,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Computes the determinant and inverse of a certain
C            REAL SYMMETRIC POSITIVE DEFINITE matrix (see abstract)
C            using the factors computed by SPOCO, SPOFA or SQRDC.
C***DESCRIPTION
C
C     SPODI computes the determinant and inverse of a certain
C     real symmetric positive definite matrix (see below)
C     using the factors computed by SPOCO, SPOFA or SQRDC.
C
C     On Entry
C
C        A       REAL(LDA, N)
C                the output  A  from SPOCO or SPOFA
C                or the output  X  from SQRDC.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       If SPOCO or SPOFA was used to factor  A , then
C                SPODI produces the upper half of INVERSE(A) .
C                If SQRDC was used to decompose  X , then
C                SPODI produces the upper half of INVERSE(TRANS(X)*X),
C                where TRANS(X) is the transpose.
C                Elements of  A  below the diagonal are unchanged.
C                If the units digit of JOB is zero,  A  is unchanged.
C
C        DET     REAL(2)
C                determinant of  A  or of  TRANS(X)*X  if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DET(1) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if SPOCO or SPOFA has set INFO .EQ. 0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS SAXPY,SSCAL
C     Fortran MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SAXPY,SSCAL
C***END PROLOGUE  SPODI
      INTEGER LDA,N,JOB
      REAL A(LDA,1)
      REAL DET(2)
C
      REAL T
      REAL S
      INTEGER I,J,JM1,K,KP1
C
C     COMPUTE DETERMINANT
C
C***FIRST EXECUTABLE STATEMENT  SPODI
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0
         DET(2) = 0.0E0
         S = 10.0E0
         DO 50 I = 1, N
            DET(1) = A(I,I)**2*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0E0) GO TO 60
   10       IF (DET(1) .GE. 1.0E0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0E0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0E0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(R)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         DO 100 K = 1, N
            A(K,K) = 1.0E0/A(K,K)
            T = -A(K,K)
            CALL SSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0E0
               CALL SAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM  INVERSE(R) * TRANS(INVERSE(R))
C
         DO 130 J = 1, N
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = A(K,J)
               CALL SAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
            T = A(J,J)
            CALL SSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
      SUBROUTINE SPOFA(A,LDA,N,INFO)
C***BEGIN PROLOGUE  SPOFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=SINGLE PRECISION(SPOFA-S DPOFA-D CPOFA-C),
C             LINEAR ALGEBRA,MATRIX,MATRIX FACTORIZATION,
C             POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a real SYMMETRIC POSITIVE DEFINITE matrix.
C***DESCRIPTION
C
C     SPOFA factors a real symmetric positive definite matrix.
C
C     SPOFA is usually called by SPOCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for SPOCO) = (1 + 18/N)*(Time for SPOFA) .
C
C     On Entry
C
C        A       REAL(LDA, N)
C                the symmetric matrix to be factored.  Only the
C                diagonal and upper triangle are used.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix  R  so that  A = TRANS(R)*R
C                where  TRANS(R)  is the transpose.
C                The strict lower triangle is unaltered.
C                If  INFO .NE. 0 , the factorization is not complete.
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS SDOT
C     Fortran SQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SDOT
C***END PROLOGUE  SPOFA
      INTEGER LDA,N,INFO
      REAL A(LDA,1)
C
      REAL SDOT,T
      REAL S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
C***FIRST EXECUTABLE STATEMENT  SPOFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - SDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
C     ......EXIT
            IF (S .LE. 0.0E0) GO TO 40
            A(J,J) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SPOSL(A,LDA,N,B)
C***BEGIN PROLOGUE  SPOSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=SINGLE PRECISION(SPOSL-S DPOSL-D CPOSL-C),
C             LINEAR ALGEBRA,MATRIX,POSITIVE DEFINITE,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the real SYMMETRIC POSITIVE DEFINITE system A*X=B
C            using the factors computed by SPOCO or SPOFA.
C***DESCRIPTION
C
C     SPOSL solves the real symmetric positive definite system
C     A * X = B
C     using the factors computed by SPOCO or SPOFA.
C
C     On Entry
C
C        A       REAL(LDA, N)
C                the output from SPOCO or SPOFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        B       REAL(N)
C                the right hand side vector.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal.  Technically, this indicates
C        singularity, but it is usually caused by improper subroutine
C        arguments.  It will not occur if the subroutines are called
C        correctly and  INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL SPOCO(A,LDA,N,RCOND,Z,INFO)
C           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL SPOSL(A,LDA,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS SAXPY,SDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SAXPY,SDOT
C***END PROLOGUE  SPOSL
      INTEGER LDA,N
      REAL A(LDA,1),B(1)
C
      REAL SDOT,T
      INTEGER K,KB
C
C     SOLVE TRANS(R)*Y = B
C
C***FIRST EXECUTABLE STATEMENT  SPOSL
      DO 10 K = 1, N
         T = SDOT(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL SAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
      END
      REAL FUNCTION SASUM(N,SX,INCX)
C***BEGIN PROLOGUE  SASUM
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A3A
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=SINGLE PRECISION(SASUM-S DASUM-D SCASUM-C),ADD,
C             LINEAR ALGEBRA,MAGNITUDE,SUM,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Sum of magnitudes of s.p vector components
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(S)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C
C     --Output--
C    SASUM  single precision result (zero if N .LE. 0)
C
C     Returns sum of magnitudes of single precision SX.
C     SASUM = sum from 0 to N-1 of  ABS(SX(1+I*INCX))
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SASUM
C
      REAL SX(1)
C***FIRST EXECUTABLE STATEMENT  SASUM
      SASUM = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        SASUM = SASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     1  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SAXPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A7
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=SINGLE PRECISION(SAXPY-S DAXPY-D CAXPY-C),
C             LINEAR ALGEBRA,TRIAD,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  S.P. computation y = a*x + y
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SA  single precision scalar multiplier
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SY  single precision result (unchanged if N .LE. 0)
C
C     Overwrite single precision SY with single precision SA*SX +SY.
C     For I = 0 to N-1, replace  SY(LY+I*INCY) with SA*SX(LX+I*INCX) +
C       SY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
C       and LY is defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SAXPY
C
      REAL SX(1),SY(1),SA
C***FIRST EXECUTABLE STATEMENT  SAXPY
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A4
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=SINGLE PRECISION(SDOT-S DDOT-D CDOTU-C),
C             INNER PRODUCT,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  S.P. inner product of s.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C     SDOT  single precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of single precision SX and SY.
C     SDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SDOT
C
      REAL SX(1),SY(1)
C***FIRST EXECUTABLE STATEMENT  SDOT
      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     1   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
C
C        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
      SUBROUTINE SSCAL(N,SA,SX,INCX)
C***BEGIN PROLOGUE  SSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A6
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=SINGLE PRECISION(SSCAL-S DSCAL-D CSCAL-C),
C             LINEAR ALGEBRA,SCALE,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  S.P. vector scale x = a*x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SA  single precision scale factor
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C
C     --Output--
C       SX  single precision result (unchanged if N .LE. 0)
C
C     Replace single precision SX by single precision SA*SX.
C     For I = 0 to N-1, replace SX(1+I*INCX) with  SA * SX(1+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSCAL
C
      REAL SA,SX(1)
C***FIRST EXECUTABLE STATEMENT  SSCAL
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END
      subroutine spbfa(abd,lda,n,m,info)
      integer lda,n,m,info
      real abd(lda,1)
c
c     spbfa factors a real symmetric positive definite matrix
c     stored in band form.
c
c     spbfa is usually called by spbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                the matrix to be factored.  the columns of the upper
c                triangle are stored in the columns of abd and the
c                diagonals of the upper triangle are stored in the
c                rows of abd .  see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. m + 1 .
c
c        n       integer
c                the order of the matrix  a .
c
c        m       integer
c                the number of diagonals above the main diagonal.
c                0 .le. m .lt. n .
c
c     on return
c
c        abd     an upper triangular matrix  r , stored in band
c                form, so that  a = trans(r)*r .
c
c        info    integer
c                = 0  for normal return.
c                = k  if the leading minor of order  k  is not
c                     positive definite.
c
c     band storage
c
c           if  a  is a symmetric positive definite band matrix,
c           the following program segment will set up the input.
c
c                   m = (band width above diagonal)
c                   do 20 j = 1, n
c                      i1 = max0(1, j-m)
c                      do 10 i = i1, j
c                         k = i-j+m+1
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas sdot
c     fortran max0,sqrt
c
c     internal variables
c
      real sdot,t
      real s
      integer ik,j,jk,k,mu
c     begin block with ...exits to 40
c
c
         do 30 j = 1, n
            info = j
            s = 0.0e0
            ik = m + 1
            jk = max0(j-m,1)
            mu = max0(m+2-j,1)
            if (m .lt. mu) go to 20
            do 10 k = mu, m
               t = abd(k,j) - sdot(k-mu,abd(ik,jk),1,abd(mu,j),1)
               t = t/abd(m+1,jk)
               abd(k,j) = t
               s = s + t*t
               ik = ik - 1
               jk = jk + 1
   10       continue
   20       continue
            s = abd(m+1,j) - s
c     ......exit
            if (s .le. 0.0e0) go to 40
            abd(m+1,j) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
      subroutine scopy(nn,vin,iin,vout,iout)
 
c     This routine copies array vin into array vout.
 
      real vout(*), vin(*)
 
      do 10 ii=1,nn-2,3
      vout(ii)   = vin(ii)
      vout(ii+1) = vin(ii+1)
      vout(ii+2) = vin(ii+2)
 10   continue
 
      do 20 ii=3*(nn/3)+1,nn
      vout(ii)   = vin(ii)
 20   continue
 
      end
