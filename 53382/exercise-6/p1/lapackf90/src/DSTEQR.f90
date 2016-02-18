      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band symmetric matrix can also be found
!  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
!  tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  symmetric matrix.  On entry, Z must contain the
!                  orthogonal matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero; on exit, D
!                and E contain the elements of a symmetric tridiagonal
!                matrix which is orthogonally similar to the original
!                matrix.
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
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                         LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, &
                         NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                         S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR, &
                         DLASRT, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      if ( LSAME( COMPZ, 'N' ) ) then
         ICOMPZ = 0
      else if ( LSAME( COMPZ, 'V' ) ) then
         ICOMPZ = 1
      else if ( LSAME( COMPZ, 'I' ) ) then
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      end if
      if ( ICOMPZ < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( ( LDZ < 1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, &
               N ) ) ) then
         INFO = -6
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
      if ( N == 1 ) then
         if ( ICOMPZ == 2 ) &
            Z( 1, 1 ) = ONE
         RETURN
      end if
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      if ( ICOMPZ == 2 ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      if ( L1.GT.N ) &
         GO TO 160
      if ( L1.GT.1 ) &
         E( L1-1 ) = ZERO
      if ( L1.LE.NM1 ) then
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            if ( TST == ZERO ) &
               GO TO 30
            if ( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
                1 ) ) ) )*EPS ) then
               E( M ) = ZERO
               GO TO 30
            end if
   20    CONTINUE
      end if
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
      if ( ANORM == ZERO ) &
         GO TO 10
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
!     Choose between QL and QR iteration
!
      if ( ABS( D( LEND ) ) < ABS( D( L ) ) ) then
         LEND = LSV
         L = LENDSV
      end if
!
      if ( LEND.GT.L ) then
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE
         if ( L.NE.LEND ) then
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               if ( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ &
                   SAFMIN )GO TO 60
   50       CONTINUE
         end if
!
         M = LEND
!
   60    CONTINUE
         if ( M < LEND ) &
            E( M ) = ZERO
         P = D( L )
         if ( M == L ) &
            GO TO 80
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         if ( M == L+1 ) then
            if ( ICOMPZ.GT.0 ) then
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ), &
                           WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            end if
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            if ( L.LE.LEND ) &
               GO TO 40
            GO TO 140
         end if
!
         if ( JTOT == NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            if ( I.NE.M-1 ) &
               E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            if ( ICOMPZ.GT.0 ) then
               WORK( I ) = C
               WORK( N-1+I ) = -S
            end if
!
   70    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         if ( ICOMPZ.GT.0 ) then
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), &
                        Z( 1, L ), LDZ )
         end if
!
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
!
!        Eigenvalue found.
!
   80    CONTINUE
         D( L ) = P
!
         L = L + 1
         if ( L.LE.LEND ) &
            GO TO 40
         GO TO 140
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE
         if ( L.NE.LEND ) then
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               if ( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ &
                   SAFMIN )GO TO 110
  100       CONTINUE
         end if
!
         M = LEND
!
  110    CONTINUE
         if ( M.GT.LEND ) &
            E( M-1 ) = ZERO
         P = D( L )
         if ( M == L ) &
            GO TO 130
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         if ( M == L-1 ) then
            if ( ICOMPZ.GT.0 ) then
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ), &
                           WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            end if
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            if ( L.GE.LEND ) &
               GO TO 90
            GO TO 140
         end if
!
         if ( JTOT == NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            if ( I.NE.M ) &
               E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            if ( ICOMPZ.GT.0 ) then
               WORK( I ) = C
               WORK( N-1+I ) = S
            end if
!
  120    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         if ( ICOMPZ.GT.0 ) then
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), &
                        Z( 1, M ), LDZ )
         end if
!
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
!
!        Eigenvalue found.
!
  130    CONTINUE
         D( L ) = P
!
         L = L - 1
         if ( L.GE.LEND ) &
            GO TO 90
         GO TO 140
!
      end if
!
!     Undo scaling if necessary
!
  140 CONTINUE
      if ( ISCALE == 1 ) then
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      else if ( ISCALE == 2 ) then
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      end if
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      if ( JTOT < NMAXIT ) &
         GO TO 10
      DO 150 I = 1, N - 1
         if ( E( I ).NE.ZERO ) &
            INFO = INFO + 1
  150 CONTINUE
      GO TO 190
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      if ( ICOMPZ == 0 ) then
!
!        Use Quick Sort
!
         CALL DLASRT( 'I', N, D, INFO )
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               if ( D( J ) < P ) then
                  K = J
                  P = D( J )
               end if
  170       CONTINUE
            if ( K.NE.I ) then
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            end if
  180    CONTINUE
      end if
!
  190 CONTINUE
      RETURN
!
!     End of DSTEQR
!
      END
