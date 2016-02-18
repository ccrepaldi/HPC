      SUBROUTINE DDISNA( JOB, M, N, D, SEP, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), SEP( * )
!     ..
!
!  Purpose
!  =======
!
!  DDISNA computes the reciprocal condition numbers for the eigenvectors
!  of a real symmetric or complex Hermitian matrix or for the left or
!  right singular vectors of a general m-by-n matrix. The reciprocal
!  condition number is the 'gap' between the corresponding eigenvalue or
!  singular value and the nearest other one.
!
!  The bound on the error, measured by angle in radians, in the I-th
!  computed vector is given by
!
!         DLAMCH( 'E' ) * ( ANORM / SEP( I ) )
!
!  where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed
!  to be smaller than DLAMCH( 'E' )*ANORM in order to limit the size of
!  the error bound.
!
!  DDISNA may also be used to compute error bounds for eigenvectors of
!  the generalized symmetric definite eigenproblem.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies for which problem the reciprocal condition numbers
!          should be computed:
!          = 'E':  the eigenvectors of a symmetric/Hermitian matrix;
!          = 'L':  the left singular vectors of a general matrix;
!          = 'R':  the right singular vectors of a general matrix.
!
!  M       (input) INTEGER
!          The number of rows of the matrix. M >= 0.
!
!  N       (input) INTEGER
!          If JOB = 'L' or 'R', the number of columns of the matrix,
!          in which case N >= 0. Ignored if JOB = 'E'.
!
!  D       (input) DOUBLE PRECISION array, dimension (M) if JOB = 'E'
!                              dimension (min(M,N)) if JOB = 'L' or 'R'
!          The eigenvalues (if JOB = 'E') or singular values (if JOB =
!          'L' or 'R') of the matrix, in either increasing or decreasing
!          order. If singular values, they must be non-negative.
!
!  SEP     (output) DOUBLE PRECISION array, dimension (M) if JOB = 'E'
!                               dimension (min(M,N)) if JOB = 'L' or 'R'
!          The reciprocal condition numbers of the vectors.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DECR, EIGEN, INCR, LEFT, RIGHT, SING
      INTEGER            I, K
      DOUBLE PRECISION   ANORM, EPS, NEWGAP, OLDGAP, SAFMIN, THRESH
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      EIGEN = LSAME( JOB, 'E' )
      LEFT = LSAME( JOB, 'L' )
      RIGHT = LSAME( JOB, 'R' )
      SING = LEFT .OR. RIGHT
      if ( EIGEN ) then
         K = M
      else if ( SING ) then
         K = MIN( M, N )
      end if
      if ( .NOT.EIGEN .AND. .NOT.SING ) then
         INFO = -1
      else if ( M < 0 ) then
         INFO = -2
      else if ( K < 0 ) then
         INFO = -3
      ELSE
         INCR = .TRUE.
         DECR = .TRUE.
         DO I = 1, K - 1
            if ( INCR ) &
               INCR = INCR .AND. D( I ).LE.D( I+1 )
            if ( DECR ) &
               DECR = DECR .AND. D( I ).GE.D( I+1 )
         end do
         if ( SING .AND. K.GT.0 ) then
            if ( INCR ) &
               INCR = INCR .AND. ZERO.LE.D( 1 )
            if ( DECR ) &
               DECR = DECR .AND. D( K ).GE.ZERO
         end if
         if ( .NOT.( INCR .OR. DECR ) ) &
            INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DDISNA', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( K == 0 ) &
         RETURN
!
!     Compute reciprocal condition numbers
!
      if ( K == 1 ) then
         SEP( 1 ) = DLAMCH( 'O' )
      ELSE
         OLDGAP = ABS( D( 2 )-D( 1 ) )
         SEP( 1 ) = OLDGAP
         DO 20 I = 2, K - 1
            NEWGAP = ABS( D( I+1 )-D( I ) )
            SEP( I ) = MIN( OLDGAP, NEWGAP )
            OLDGAP = NEWGAP
   20    CONTINUE
         SEP( K ) = OLDGAP
      end if
      if ( SING ) then
         if ( ( LEFT .AND. M.GT.N ) .OR. ( RIGHT .AND. M < N ) ) then
            if ( INCR ) &
               SEP( 1 ) = MIN( SEP( 1 ), D( 1 ) )
            if ( DECR ) &
               SEP( K ) = MIN( SEP( K ), D( K ) )
         end if
      end if
!
!     Ensure that reciprocal condition numbers are not less than
!     threshold, in order to limit the size of the error bound
!
      EPS = DLAMCH( 'E' )
      SAFMIN = DLAMCH( 'S' )
      ANORM = MAX( ABS( D( 1 ) ), ABS( D( K ) ) )
      if ( ANORM == ZERO ) then
         THRESH = EPS
      ELSE
         THRESH = MAX( EPS*ANORM, SAFMIN )
      end if
      DO 30 I = 1, K
         SEP( I ) = MAX( SEP( I ), THRESH )
   30 CONTINUE
!
      RETURN
!
!     End of DDISNA
!
      END
