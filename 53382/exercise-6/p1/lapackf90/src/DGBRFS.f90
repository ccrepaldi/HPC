      SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, &
                         IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, &
                         INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), &
                         BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBRFS improves the computed solution to a system of linear
!  equations when the coefficient matrix is banded, and provides
!  error bounds and backward error estimates for the solution.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X.  NRHS >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The original band matrix A, stored in rows 1 to KL+KU+1.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KL+KU+1.
!
!  AFB     (input) DOUBLE PRECISION array, dimension (LDAFB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by DGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAFB   (input) INTEGER
!          The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGBTRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          The right hand side matrix B.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS)
!          On entry, the solution matrix X, as computed by DGBTRS.
!          On exit, the improved solution matrix X.
!
!  LDX     (input) INTEGER
!          The leading dimension of the array X.  LDX >= max(1,N).
!
!  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The estimated forward error bound for each solution vector
!          X(j) (the j-th column of the solution matrix X).
!          If XTRUE is the true solution corresponding to X(j), FERR(j)
!          is an estimated upper bound for the magnitude of the largest
!          element in (X(j) - XTRUE) divided by the magnitude of the
!          largest element in X(j).  The estimate is as reliable as
!          the estimate for RCOND, and is almost always a slight
!          overestimate of the true error.
!
!  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The componentwise relative backward error of each solution
!          vector X(j) (i.e., the smallest relative change in
!          any element of A or B that makes X(j) an exact solution).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Internal Parameters
!  ===================
!
!  ITMAX is the maximum number of steps of iterative refinement.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
      DOUBLE PRECISION   THREE
      PARAMETER          ( THREE = 3.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
      CHARACTER          TRANST
      INTEGER            COUNT, I, J, K, KASE, KK, NZ
      DOUBLE PRECISION   EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGBMV, DGBTRS, DLACON, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
          LSAME( TRANS, 'C' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( KL < 0 ) then
         INFO = -3
      else if ( KU < 0 ) then
         INFO = -4
      else if ( NRHS < 0 ) then
         INFO = -5
      else if ( LDAB < KL+KU+1 ) then
         INFO = -7
      else if ( LDAFB < 2*KL+KU+1 ) then
         INFO = -9
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -12
      else if ( LDX < MAX( 1, N ) ) then
         INFO = -14
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGBRFS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 .OR. NRHS.EQ.0 ) then
         FERR(1:nrhs) = ZERO
         BERR(1:nrhs) = ZERO
         RETURN
      end if
!
      if ( NOTRAN ) then
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      end if
!
!     NZ = maximum number of nonzero elements in each row of A, plus 1
!
      NZ = MIN( KL+KU+2, N+1 )
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
!
!     Do for each right hand side
!
      DO 140 J = 1, NRHS
!
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
!
!        Loop until stopping criterion is satisfied.
!
!        Compute residual R = B - op(A) * X,
!        where op(A) = A, A**T, or A**H, depending on TRANS.
!
         CALL DCOPY( N, B( 1, J ), 1, WORK( N+1 ), 1 )
         CALL DGBMV( TRANS, N, N, KL, KU, -ONE, AB, LDAB, X( 1, J ), 1, &
                     ONE, WORK( N+1 ), 1 )
!
!        Compute componentwise relative backward error from formula
!
!        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
!
!        where abs(Z) is the componentwise absolute value of the matrix
!        or vector Z.  If the i-th component of the denominator is less
!        than SAFE2, then SAFE1 is added to the i-th components of the
!        numerator and denominator before dividing.
!
         DO 30 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   30    CONTINUE
!
!        Compute abs(op(A))*abs(X) + abs(B).
!
         if ( NOTRAN ) then
            DO 50 K = 1, N
               KK = KU + 1 - K
               XK = ABS( X( K, J ) )
               DO 40 I = MAX( 1, K-KU ), MIN( N, K+KL )
                  WORK( I ) = WORK( I ) + ABS( AB( KK+I, K ) )*XK
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 70 K = 1, N
               S = ZERO
               KK = KU + 1 - K
               DO 60 I = MAX( 1, K-KU ), MIN( N, K+KL )
                  S = S + ABS( AB( KK+I, K ) )*ABS( X( I, J ) )
   60          CONTINUE
               WORK( K ) = WORK( K ) + S
   70       CONTINUE
         end if
         S = ZERO
         DO 80 I = 1, N
            if ( WORK( I ).GT.SAFE2 ) then
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / &
                   ( WORK( I )+SAFE1 ) )
            end if
   80    CONTINUE
         BERR( J ) = S
!
!        Test stopping criterion. Continue iterating if
!           1) The residual BERR(J) is larger than machine epsilon, and
!           2) BERR(J) decreased by at least a factor of 2 during the
!              last iteration, and
!           3) At most ITMAX iterations tried.
!
         if ( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND. &
             COUNT.LE.ITMAX ) then
!
!           Update solution and try again.
!
            CALL DGBTRS( TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, &
                         WORK( N+1 ), N, INFO )
            CALL DAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         end if
!
!        Bound error from formula
!
!        norm(X - XTRUE) / norm(X) .le. FERR =
!        norm( abs(inv(op(A)))*
!           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
!
!        where
!          norm(Z) is the magnitude of the largest component of Z
!          inv(op(A)) is the inverse of op(A)
!          abs(Z) is the componentwise absolute value of the matrix or
!             vector Z
!          NZ is the maximum number of nonzeros in any row of A, plus 1
!          EPS is machine epsilon
!
!        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
!        is incremented by SAFE1 if the i-th component of
!        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
!
!        Use DLACON to estimate the infinity-norm of the matrix
!           inv(op(A)) * diag(W),
!        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
!
         DO 90 I = 1, N
            if ( WORK( I ).GT.SAFE2 ) then
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            end if
   90    CONTINUE
!
         KASE = 0
  100    CONTINUE
         CALL DLACON( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), &
                      KASE )
         if ( KASE.NE.0 ) then
            if ( KASE == 1 ) then
!
!              Multiply by diag(W)*inv(op(A)**T).
!
               CALL DGBTRS( TRANST, N, KL, KU, 1, AFB, LDAFB, IPIV, &
                            WORK( N+1 ), N, INFO )
               DO 110 I = 1, N
                  WORK( N+I ) = WORK( N+I )*WORK( I )
  110          CONTINUE
            ELSE
!
!              Multiply by inv(op(A))*diag(W).
!
               DO 120 I = 1, N
                  WORK( N+I ) = WORK( N+I )*WORK( I )
  120          CONTINUE
               CALL DGBTRS( TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, &
                            WORK( N+1 ), N, INFO )
            end if
            GO TO 100
         end if
!
!        Normalize error.
!
         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  130    CONTINUE
         if ( LSTRES.NE.ZERO ) &
            FERR( J ) = FERR( J ) / LSTRES
!
  140 CONTINUE
!
      RETURN
!
!     End of DGBRFS
!
      END