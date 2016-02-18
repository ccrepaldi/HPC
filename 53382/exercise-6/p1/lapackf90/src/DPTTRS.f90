      SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTTRS solves a tridiagonal system of the form
!     A * X = B
!  using the L*D*L' factorization of A computed by DPTTRF.  D is a
!  diagonal matrix specified in the vector D, L is a unit bidiagonal
!  matrix whose subdiagonal is specified in the vector E, and X and B
!  are N by NRHS matrices.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the diagonal matrix D from the
!          L*D*L' factorization of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) subdiagonal elements of the unit bidiagonal factor
!          L from the L*D*L' factorization of A.  E can also be regarded
!          as the superdiagonal of the unit bidiagonal factor U from the
!          factorization A = U'*D*U.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            J, JB, NB
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DPTTS2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      if ( N < 0 ) then
         INFO = -1
      else if ( NRHS < 0 ) then
         INFO = -2
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -6
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPTTRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 .OR. NRHS.EQ.0 ) &
         RETURN
!
!     Determine the number of right-hand sides to solve at a time.
!
      if ( NRHS == 1 ) then
         NB = 1
      ELSE
         NB = MAX( 1, ILAENV( 1, 'DPTTRS', ' ', N, NRHS, -1, -1 ) )
      end if
!
      if ( NB.GE.NRHS ) then
         CALL DPTTS2( N, NRHS, D, E, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL DPTTS2( N, JB, D, E, B( 1, J ), LDB )
   10    CONTINUE
      end if
!
      RETURN
!
!     End of DPTTRS
!
      END
