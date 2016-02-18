      SUBROUTINE DSTERF( N, D, E, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm failed to find all of the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, &
                         NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC, &
                         OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN, &
                         SIGMA, SSFMAX, SSFMIN
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
!     Quick return if possible
!
      if ( N < 0 ) then
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      end if
      if ( N.LE.1 ) &
         RETURN
!
!     Determine the unit roundoff for this environment.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
!
   10 CONTINUE
      if ( L1.GT.N ) &
         GO TO 170
      if ( L1.GT.1 ) &
         E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         if ( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
             1 ) ) ) )*EPS ) then
            E( M ) = ZERO
            GO TO 30
         end if
   20 CONTINUE
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      if ( LEND == L ) &
         GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      if ( ANORM.GT.SSFMAX ) then
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, &
                      INFO )
      else if ( ANORM < SSFMIN ) then
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, &
                      INFO )
      end if
!
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
!
!     Choose between QL and QR iteration
!
      if ( ABS( D( LEND ) ) < ABS( D( L ) ) ) then
         LEND = LSV
         L = LENDSV
      end if
!
      if ( LEND.GE.L ) then
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   50    CONTINUE
         if ( L.NE.LEND ) then
            DO 60 M = L, LEND - 1
               if ( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) ) &
                  GO TO 70
   60       CONTINUE
         end if
         M = LEND
!
   70    CONTINUE
         if ( M < LEND ) &
            E( M ) = ZERO
         P = D( L )
         if ( M == L ) &
            GO TO 90
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         if ( M == L+1 ) then
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            if ( L.LE.LEND ) &
               GO TO 50
            GO TO 150
         end if
!
         if ( JTOT == NMAXIT ) &
            GO TO 150
         JTOT = JTOT + 1
!
!        Form shift.
!
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
!
!        Inner loop
!
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            if ( I.NE.M-1 ) &
               E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            if ( C.NE.ZERO ) then
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            end if
   80    CONTINUE
!
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
!
!        Eigenvalue found.
!
   90    CONTINUE
         D( L ) = P
!
         L = L + 1
         if ( L.LE.LEND ) &
            GO TO 50
         GO TO 150
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            if ( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) ) &
               GO TO 120
  110    CONTINUE
         M = LEND
!
  120    CONTINUE
         if ( M.GT.LEND ) &
            E( M-1 ) = ZERO
         P = D( L )
         if ( M == L ) &
            GO TO 140
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         if ( M == L-1 ) then
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            if ( L.GE.LEND ) &
               GO TO 100
            GO TO 150
         end if
!
         if ( JTOT == NMAXIT ) &
            GO TO 150
         JTOT = JTOT + 1
!
!        Form shift.
!
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
!
!        Inner loop
!
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            if ( I.NE.M ) &
               E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            if ( C.NE.ZERO ) then
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            end if
  130    CONTINUE
!
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
!
!        Eigenvalue found.
!
  140    CONTINUE
         D( L ) = P
!
         L = L - 1
         if ( L.GE.LEND ) &
            GO TO 100
         GO TO 150
!
      end if
!
!     Undo scaling if necessary
!
  150 CONTINUE
      if ( ISCALE == 1 ) &
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
      if ( ISCALE == 2 ) &
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      if ( JTOT < NMAXIT ) &
         GO TO 10
      DO 160 I = 1, N - 1
         if ( E( I ).NE.ZERO ) &
            INFO = INFO + 1
  160 CONTINUE
      GO TO 180
!
!     Sort eigenvalues in increasing order.
!
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
!
  180 CONTINUE
      RETURN
!
!     End of DSTERF
!
      END
