      SUBROUTINE DLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, JOB, N
      DOUBLE PRECISION   TOL
!     ..
!     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAGTS may be used to solve one of the systems of equations
!
!     (T - lambda*I)*x = y   or   (T - lambda*I)'*x = y,
!
!  where T is an n by n tridiagonal matrix, for x, following the
!  factorization of (T - lambda*I) as
!
!     (T - lambda*I) = P*L*U ,
!
!  by routine DLAGTF. The choice of equation to be solved is
!  controlled by the argument JOB, and in each case there is an option
!  to perturb zero or very small diagonal elements of U, this option
!  being intended for use in applications such as inverse iteration.
!
!  Arguments
!  =========
!
!  JOB     (input) INTEGER
!          Specifies the job to be performed by DLAGTS as follows:
!          =  1: The equations  (T - lambda*I)x = y  are to be solved,
!                but diagonal elements of U are not to be perturbed.
!          = -1: The equations  (T - lambda*I)x = y  are to be solved
!                and, if overflow would otherwise occur, the diagonal
!                elements of U are to be perturbed. See argument TOL
!                below.
!          =  2: The equations  (T - lambda*I)'x = y  are to be solved,
!                but diagonal elements of U are not to be perturbed.
!          = -2: The equations  (T - lambda*I)'x = y  are to be solved
!                and, if overflow would otherwise occur, the diagonal
!                elements of U are to be perturbed. See argument TOL
!                below.
!
!  N       (input) INTEGER
!          The order of the matrix T.
!
!  A       (input) DOUBLE PRECISION array, dimension (N)
!          On entry, A must contain the diagonal elements of U as
!          returned from DLAGTF.
!
!  B       (input) DOUBLE PRECISION array, dimension (N-1)
!          On entry, B must contain the first super-diagonal elements of
!          U as returned from DLAGTF.
!
!  C       (input) DOUBLE PRECISION array, dimension (N-1)
!          On entry, C must contain the sub-diagonal elements of L as
!          returned from DLAGTF.
!
!  D       (input) DOUBLE PRECISION array, dimension (N-2)
!          On entry, D must contain the second super-diagonal elements
!          of U as returned from DLAGTF.
!
!  IN      (input) INTEGER array, dimension (N)
!          On entry, IN must contain details of the matrix P as returned
!          from DLAGTF.
!
!  Y       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the right hand side vector y.
!          On exit, Y is overwritten by the solution vector x.
!
!  TOL     (input/output) DOUBLE PRECISION
!          On entry, with  JOB .lt. 0, TOL should be the minimum
!          perturbation to be made to very small diagonal elements of U.
!          TOL should normally be chosen as about eps*norm(U), where eps
!          is the relative machine precision, but if TOL is supplied as
!          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).
!          If  JOB .gt. 0  then TOL is not referenced.
!
!          On exit, TOL is changed as described above, only if TOL is
!          non-positive on entry. Otherwise TOL is unchanged.
!
!  INFO    (output) INTEGER
!          = 0   : successful exit
!          .lt. 0: if INFO = -i, the i-th argument had an illegal value
!          .gt. 0: overflow would occur when computing the INFO(th)
!                  element of the solution vector x. This can only occur
!                  when JOB is supplied as positive and either means
!                  that a diagonal element of U is very small, or that
!                  the elements of the right-hand side vector y are very
!                  large.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            K
      DOUBLE PRECISION   ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      if ( ( ABS( JOB ).GT.2 ) .OR. ( JOB == 0 ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DLAGTS', -INFO )
         RETURN
      end if
!
      if ( N == 0 ) &
         RETURN
!
      EPS = DLAMCH( 'Epsilon' )
      SFMIN = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SFMIN
!
      if ( JOB < 0 ) then
         if ( TOL.LE.ZERO ) then
            TOL = ABS( A( 1 ) )
            if ( N.GT.1 ) &
               TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) )
            DO 10 K = 3, N
               TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ), &
                     ABS( D( K-2 ) ) )
   10       CONTINUE
            TOL = TOL*EPS
            if ( TOL == ZERO ) &
               TOL = EPS
         end if
      end if
!
      if ( ABS( JOB ) == 1 ) then
         DO 20 K = 2, N
            if ( IN( K-1 ) == 0 ) then
               Y( K ) = Y( K ) - C( K-1 )*Y( K-1 )
            ELSE
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            end if
   20    CONTINUE
         if ( JOB == 1 ) then
            DO 30 K = N, 1, -1
               if ( K.LE.N-2 ) then
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               else if ( K == N-1 ) then
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               ELSE
                  TEMP = Y( K )
               end if
               AK = A( K )
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) then
                  if ( ABSAK < SFMIN ) then
                     if ( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     end if
                  else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) then
                     INFO = K
                     RETURN
                  end if
               end if
               Y( K ) = TEMP / AK
   30       CONTINUE
         ELSE
            DO 50 K = N, 1, -1
               if ( K.LE.N-2 ) then
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               else if ( K == N-1 ) then
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               ELSE
                  TEMP = Y( K )
               end if
               AK = A( K )
               PERT = SIGN( TOL, AK )
   40          CONTINUE
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) then
                  if ( ABSAK < SFMIN ) then
                     if ( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 40
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     end if
                  else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) then
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 40
                  end if
               end if
               Y( K ) = TEMP / AK
   50       CONTINUE
         end if
      ELSE
!
!        Come to here if  JOB = 2 or -2
!
         if ( JOB == 2 ) then
            DO 60 K = 1, N
               if ( K.GE.3 ) then
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               else if ( K == 2 ) then
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               ELSE
                  TEMP = Y( K )
               end if
               AK = A( K )
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) then
                  if ( ABSAK < SFMIN ) then
                     if ( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     end if
                  else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) then
                     INFO = K
                     RETURN
                  end if
               end if
               Y( K ) = TEMP / AK
   60       CONTINUE
         ELSE
            DO 80 K = 1, N
               if ( K.GE.3 ) then
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               else if ( K == 2 ) then
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               ELSE
                  TEMP = Y( K )
               end if
               AK = A( K )
               PERT = SIGN( TOL, AK )
   70          CONTINUE
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) then
                  if ( ABSAK < SFMIN ) then
                     if ( ABSAK == ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 70
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     end if
                  else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) then
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 70
                  end if
               end if
               Y( K ) = TEMP / AK
   80       CONTINUE
         end if
!
         DO 90 K = N, 2, -1
            if ( IN( K-1 ) == 0 ) then
               Y( K-1 ) = Y( K-1 ) - C( K-1 )*Y( K )
            ELSE
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            end if
   90    CONTINUE
      end if
!
!     End of DLAGTS
!
      END
