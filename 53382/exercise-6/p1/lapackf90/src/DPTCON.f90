      SUBROUTINE DPTCON( N, D, E, ANORM, RCOND, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
      DOUBLE PRECISION   ANORM, RCOND
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTCON computes the reciprocal of the condition number (in the
!  1-norm) of a real symmetric positive definite tridiagonal matrix
!  using the factorization A = L*D*L**T or A = U**T*D*U computed by
!  DPTTRF.
!
!  Norm(inv(A)) is computed by a direct method, and the reciprocal of
!  the condition number is computed as
!               RCOND = 1 / (ANORM * norm(inv(A))).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the diagonal matrix D from the
!          factorization of A, as computed by DPTTRF.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) off-diagonal elements of the unit bidiagonal factor
!          U or L from the factorization of A,  as computed by DPTTRF.
!
!  ANORM   (input) DOUBLE PRECISION
!          The 1-norm of the original matrix A.
!
!  RCOND   (output) DOUBLE PRECISION
!          The reciprocal of the condition number of the matrix A,
!          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the
!          1-norm of inv(A) computed in this routine.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The method used is described in Nicholas J. Higham, "Efficient
!  Algorithms for Computing the Condition Number of a Tridiagonal
!  Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IX
      DOUBLE PRECISION   AINVNM
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      if ( N < 0 ) then
         INFO = -1
      else if ( ANORM < ZERO ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPTCON', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      RCOND = ZERO
      if ( N == 0 ) then
         RCOND = ONE
         RETURN
      else if ( ANORM == ZERO ) then
         RETURN
      end if
!
!     Check that D(1:N) is positive.
!
      DO 10 I = 1, N
         if ( D( I ).LE.ZERO ) &
            RETURN
   10 CONTINUE
!
!     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
!
!        m(i,j) =  abs(A(i,j)), i = j,
!        m(i,j) = -abs(A(i,j)), i .ne. j,
!
!     and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.
!
!     Solve M(L) * x = e.
!
      WORK( 1 ) = ONE
      DO 20 I = 2, N
         WORK( I ) = ONE + WORK( I-1 )*ABS( E( I-1 ) )
   20 CONTINUE
!
!     Solve D * M(L)' * x = b.
!
      WORK( N ) = WORK( N ) / D( N )
      DO 30 I = N - 1, 1, -1
         WORK( I ) = WORK( I ) / D( I ) + WORK( I+1 )*ABS( E( I ) )
   30 CONTINUE
!
!     Compute AINVNM = max(x(i)), 1<=i<=n.
!
      IX = IDAMAX( N, WORK, 1 )
      AINVNM = ABS( WORK( IX ) )
!
!     Compute the reciprocal condition number.
!
      if ( AINVNM.NE.ZERO ) &
         RCOND = ( ONE / AINVNM ) / ANORM
!
      RETURN
!
!     End of DPTCON
!
      END
