      SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, &
                         INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
!     ..
!
!  Purpose
!  =======
!
!  DGTTRS solves one of the systems of equations
!     A*X = B  or  A'*X = B,
!  with a tridiagonal matrix A using the LU factorization computed
!  by DGTTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  DL      (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) multipliers that define the matrix L from the
!          LU factorization of A.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the upper triangular matrix U from
!          the LU factorization of A.
!
!  DU      (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) elements of the first super-diagonal of U.
!
!  DU2     (input) DOUBLE PRECISION array, dimension (N-2)
!          The (n-2) elements of the second super-diagonal of U.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= n, row i of the matrix was
!          interchanged with row IPIV(i).  IPIV(i) will always be either
!          i or i+1; IPIV(i) = i indicates a row interchange was not
!          required.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the matrix of right hand side vectors B.
!          On exit, B is overwritten by the solution vectors X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            ITRANS, J, JB, NB
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGTTS2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      NOTRAN = ( TRANS == 'N' .OR. TRANS.EQ.'n' )
      if ( .NOT.NOTRAN .AND. .NOT.( TRANS == 'T' .OR. TRANS.EQ. &
          't' ) .AND. .NOT.( TRANS == 'C' .OR. TRANS.EQ.'c' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDB < MAX( N, 1 ) ) then
         INFO = -10
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGTTRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 .OR. NRHS.EQ.0 ) &
         RETURN
!
!     Decode TRANS
!
      if ( NOTRAN ) then
         ITRANS = 0
      ELSE
         ITRANS = 1
      end if
!
!     Determine the number of right-hand sides to solve at a time.
!
      if ( NRHS == 1 ) then
         NB = 1
      ELSE
         NB = MAX( 1, ILAENV( 1, 'DGTTRS', TRANS, N, NRHS, -1, -1 ) )
      end if
!
      if ( NB.GE.NRHS ) then
         CALL DGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL DGTTS2( ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ), &
                         LDB )
   10    CONTINUE
      end if
!
!     End of DGTTRS
!
      END
