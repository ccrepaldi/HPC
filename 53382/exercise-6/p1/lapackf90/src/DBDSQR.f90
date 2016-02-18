      SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
                         LDU, C, LDC, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ), &
                         VT( LDVT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DBDSQR computes the singular value decomposition (SVD) of a real
!  N-by-N (upper or lower) bidiagonal matrix B:  B = Q * S * P' (P'
!  denotes the transpose of P), where S is a diagonal matrix with
!  non-negative diagonal elements (the singular values of B), and Q
!  and P are orthogonal matrices.
!
!  The routine computes S, and optionally computes U * Q, P' * VT,
!  or Q' * C, for given real input matrices U, VT, and C.
!
!  See "Computing  Small Singular Values of Bidiagonal Matrices With
!  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
!  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
!  no. 5, pp. 873-912, Sept 1990) and
!  "Accurate singular values and differential qd algorithms," by
!  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
!  Department, University of California at Berkeley, July 1992
!  for a detailed description of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  B is upper bidiagonal;
!          = 'L':  B is lower bidiagonal.
!
!  N       (input) INTEGER
!          The order of the matrix B.  N >= 0.
!
!  NCVT    (input) INTEGER
!          The number of columns of the matrix VT. NCVT >= 0.
!
!  NRU     (input) INTEGER
!          The number of rows of the matrix U. NRU >= 0.
!
!  NCC     (input) INTEGER
!          The number of columns of the matrix C. NCC >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the bidiagonal matrix B.
!          On exit, if INFO=0, the singular values of B in decreasing
!          order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the elements of E contain the
!          offdiagonal elements of the bidiagonal matrix whose SVD
!          is desired. On normal exit (INFO = 0), E is destroyed.
!          If the algorithm does not converge (INFO > 0), D and E
!          will contain the diagonal and superdiagonal elements of a
!          bidiagonal matrix orthogonally equivalent to the one given
!          as input. E(N) is used for workspace.
!
!  VT      (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)
!          On entry, an N-by-NCVT matrix VT.
!          On exit, VT is overwritten by P' * VT.
!          VT is not referenced if NCVT = 0.
!
!  LDVT    (input) INTEGER
!          The leading dimension of the array VT.
!          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
!
!  U       (input/output) DOUBLE PRECISION array, dimension (LDU, N)
!          On entry, an NRU-by-N matrix U.
!          On exit, U is overwritten by U * Q.
!          U is not referenced if NRU = 0.
!
!  LDU     (input) INTEGER
!          The leading dimension of the array U.  LDU >= max(1,NRU).
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)
!          On entry, an N-by-NCC matrix C.
!          On exit, C is overwritten by Q' * C.
!          C is not referenced if NCC = 0.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C.
!          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  If INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm did not converge; D and E contain the
!                elements of a bidiagonal matrix which is orthogonally
!                similar to the input matrix B;  if INFO = i, i
!                elements of E have not converged to zero.
!
!  Internal Parameters
!  ===================
!
!  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
!          TOLMUL controls the convergence criterion of the QR loop.
!          If it is positive, TOLMUL*EPS is the desired relative
!             precision in the computed singular values.
!          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
!             desired absolute accuracy in the computed singular
!             values (corresponds to relative accuracy
!             abs(TOLMUL*EPS) in the largest singular value.
!          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
!             between 10 (for fast convergence) and .1/EPS
!             (for there to be some accuracy in the results).
!          Default is to lose at either one eighth or 2 of the
!             available decimal digits in each computed singular value
!             (whichever is smaller).
!
!  MAXITR  INTEGER, default = 6
!          MAXITR controls the maximum number of passes of the
!          algorithm through its inner loop. The algorithms stops
!          (and so fails to converge) if the number of passes
!          through the inner loop exceeds MAXITR*N**2.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   NEGONE
      PARAMETER          ( NEGONE = -1.0D0 )
      DOUBLE PRECISION   HNDRTH
      PARAMETER          ( HNDRTH = 0.01D0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 10.0D0 )
      DOUBLE PRECISION   HNDRD
      PARAMETER          ( HNDRD = 100.0D0 )
      DOUBLE PRECISION   MEIGTH
      PARAMETER          ( MEIGTH = -0.125D0 )
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 6 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, ROTATE
      INTEGER            I, IDIR, ISUB, ITER, J, LL, LLL, M, MAXIT, NM1, &
                         NM12, NM13, OLDLL, OLDM
      DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU, &
                         OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL, &
                         SINR, SLL, SMAX, SMIN, SMINL
!     DOUBLE PRECISION   SMINLO
      DOUBLE PRECISION   SMINOA, &
                         SN, THRESH, TOL, TOLMUL, UNFL
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARTG, DLAS2, DLASQ1, DLASR, DLASV2, DROT, &
                         DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NCVT < 0 ) then
         INFO = -3
      else if ( NRU < 0 ) then
         INFO = -4
      else if ( NCC < 0 ) then
         INFO = -5
      else if ( ( NCVT == 0 .AND. LDVT < 1 ) .OR. &
               ( NCVT.GT.0 .AND. LDVT < MAX( 1, N ) ) ) then
         INFO = -9
      else if ( LDU < MAX( 1, NRU ) ) then
         INFO = -11
      else if ( ( NCC == 0 .AND. LDC < 1 ) .OR. &
               ( NCC.GT.0 .AND. LDC < MAX( 1, N ) ) ) then
         INFO = -13
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DBDSQR', -INFO )
         RETURN
      end if
      if ( N == 0 ) &
         RETURN
      if ( N == 1 ) &
         GO TO 160
!
!     ROTATE is true if any singular vectors desired, false otherwise
!
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
!
!     If no singular vectors desired, use qd algorithm
!
      if ( .NOT.ROTATE ) then
         CALL DLASQ1( N, D, E, WORK, INFO )
         RETURN
      end if
!
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
      IDIR = 0
!
!     Get machine constants
!
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
      if ( LOWER ) then
         DO I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            WORK( I ) = CS
            WORK( NM1+I ) = SN
         end do
!
!  Update singular vectors if desired
!
         if ( NRU.GT.0 ) &
            CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U, &
                        LDU )
         if ( NCC.GT.0 ) &
            CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C, &
                        LDC )
      end if
!
!     Compute singular values to relative accuracy TOL
!     (By setting TOL to be negative, algorithm will compute
!     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
!
      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS
!
!     Compute approximate maximum, minimum singular values
!
      SMAX = ZERO
      DO 20 I = 1, N
         SMAX = MAX( SMAX, ABS( D( I ) ) )
   20 CONTINUE
      DO 30 I = 1, N - 1
         SMAX = MAX( SMAX, ABS( E( I ) ) )
   30 CONTINUE
      SMINL = ZERO
      if ( TOL.GE.ZERO ) then
!
!        Relative accuracy desired
!
         SMINOA = ABS( D( 1 ) )
         if ( SMINOA == ZERO ) &
            GO TO 50
         MU = SMINOA
         DO 40 I = 2, N
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
            SMINOA = MIN( SMINOA, MU )
            if ( SMINOA == ZERO ) &
               GO TO 50
   40    CONTINUE
   50    CONTINUE
         SMINOA = SMINOA / SQRT( DBLE( N ) )
         THRESH = MAX( TOL*SMINOA, MAXITR*N*N*UNFL )
      ELSE
!
!        Absolute accuracy desired
!
         THRESH = MAX( ABS( TOL )*SMAX, MAXITR*N*N*UNFL )
      end if
!
!     Prepare for main iteration loop for the singular values
!     (MAXIT is the maximum number of passes through the inner
!     loop permitted before nonconvergence signalled.)
!
      MAXIT = MAXITR*N*N
      ITER = 0
      OLDLL = -1
      OLDM = -1
!
!     M points to last element of unconverged part of matrix
!
      M = N
!
!     Begin main iteration loop
!
   60 CONTINUE
!
!     Check for convergence or exceeding iteration count
!
      if ( M.LE.1 ) &
         GO TO 160
      if ( ITER.GT.MAXIT ) &
         GO TO 200
!
!     Find diagonal block of matrix to work on
!
      if ( TOL < ZERO .AND. ABS( D( M ) ).LE.THRESH ) &
         D( M ) = ZERO
      SMAX = ABS( D( M ) )
      SMIN = SMAX
      DO 70 LLL = 1, M - 1
         LL = M - LLL
         ABSS = ABS( D( LL ) )
         ABSE = ABS( E( LL ) )
         if ( TOL < ZERO .AND. ABSS.LE.THRESH ) &
            D( LL ) = ZERO
         if ( ABSE.LE.THRESH ) &
            GO TO 80
         SMIN = MIN( SMIN, ABSS )
         SMAX = MAX( SMAX, ABSS, ABSE )
   70 CONTINUE
      LL = 0
      GO TO 90
   80 CONTINUE
      E( LL ) = ZERO
!
!     Matrix splits since E(LL) = 0
!
      if ( LL == M-1 ) then
!
!        Convergence of bottom singular value, return to top of loop
!
         M = M - 1
         GO TO 60
      end if
   90 CONTINUE
      LL = LL + 1
!
!     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
!
      if ( LL == M-1 ) then
!
!        2 by 2 block, handle separately
!
         CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR, &
                      COSR, SINL, COSL )
         D( M-1 ) = SIGMX
         E( M-1 ) = ZERO
         D( M ) = SIGMN
!
!        Compute singular vectors, if desired
!
         if ( NCVT.GT.0 ) &
            CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR, &
                       SINR )
         if ( NRU.GT.0 ) &
            CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
         if ( NCC.GT.0 ) &
            CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL, &
                       SINL )
         M = M - 2
         GO TO 60
      end if
!
!     If working on new submatrix, choose shift direction
!     (from larger end diagonal element towards smaller)
!
      if ( LL.GT.OLDM .OR. M < OLDLL ) then
         if ( ABS( D( LL ) ).GE.ABS( D( M ) ) ) then
!
!           Chase bulge from top (big end) to bottom (small end)
!
            IDIR = 1
         ELSE
!
!           Chase bulge from bottom (big end) to top (small end)
!
            IDIR = 2
         end if
      end if
!
!     Apply convergence tests
!
      if ( IDIR == 1 ) then
!
!        Run convergence test in forward direction
!        First apply standard test to bottom of matrix
!
         if ( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR. &
             ( TOL < ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) then
            E( M-1 ) = ZERO
            GO TO 60
         end if

         if ( TOL.GE.ZERO ) then
!
!  If relative accuracy desired, apply convergence criterion forward
!
            MU = ABS( D( LL ) )
            SMINL = MU
            DO LLL = LL, M - 1
               if ( ABS( E( LLL ) ).LE.TOL*MU ) then
                  E( LLL ) = ZERO
                  GO TO 60
               end if
!              SMINLO = SMINL
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
            end do
         end if

      ELSE
!
!        Run convergence test in backward direction
!        First apply standard test to top of matrix
!
         if ( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR. &
             ( TOL < ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) then
            E( LL ) = ZERO
            GO TO 60
         end if
!
         if ( TOL.GE.ZERO ) then
!
!           If relative accuracy desired,
!           apply convergence criterion backward
!
            MU = ABS( D( M ) )
            SMINL = MU
            DO LLL = M - 1, LL, -1
               if ( ABS( E( LLL ) ).LE.TOL*MU ) then
                  E( LLL ) = ZERO
                  GO TO 60
               end if
!              SMINLO = SMINL
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
            end do
         end if
      end if
      OLDLL = LL
      OLDM = M
!
!     Compute shift.  First, test if shifting would ruin relative
!     accuracy, and if so set the shift to zero.
!
      if ( TOL.GE.ZERO .AND. N*TOL*( SMINL / SMAX ).LE. &
          MAX( EPS, HNDRTH*TOL ) ) then
!
!        Use a zero shift to avoid loss of relative accuracy
!
         SHIFT = ZERO
      ELSE
!
!        Compute the shift from 2-by-2 block at end of matrix
!
         if ( IDIR == 1 ) then
            SLL = ABS( D( LL ) )
            CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
         ELSE
            SLL = ABS( D( M ) )
            CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
         end if
!
!        Test if shift negligible, and if so set to zero
!
         if ( SLL.GT.ZERO ) then
            if ( ( SHIFT / SLL )**2 < EPS ) &
               SHIFT = ZERO
         end if
      end if
!
!     Increment iteration count
!
      ITER = ITER + M - LL
!
!     If SHIFT = 0, do simplified QR iteration
!
      if ( SHIFT == ZERO ) then
         if ( IDIR == 1 ) then
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
            CS = ONE
            OLDCS = ONE
            DO I = LL, M - 1
               CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
               if ( I.GT.LL ) &
                  E( I-1 ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
               WORK( I-LL+1 ) = CS
               WORK( I-LL+1+NM1 ) = SN
               WORK( I-LL+1+NM12 ) = OLDCS
               WORK( I-LL+1+NM13 ) = OLDSN
            end do
            H = D( M )*CS
            D( M ) = H*OLDCS
            E( M-1 ) = H*OLDSN
!
!           Update singular vectors
!
            if ( NCVT.GT.0 ) &
               CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), &
                           WORK( N ), VT( LL, 1 ), LDVT )
            if ( NRU.GT.0 ) &
               CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), &
                           WORK( NM13+1 ), U( 1, LL ), LDU )
            if ( NCC.GT.0 ) &
               CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), &
                           WORK( NM13+1 ), C( LL, 1 ), LDC )
!
!           Test convergence
!
            if ( ABS( E( M-1 ) ).LE.THRESH ) &
               E( M-1 ) = ZERO
!
         ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
            CS = ONE
            OLDCS = ONE
            DO I = M, LL + 1, -1
               CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
               if ( I < M ) &
                  E( I ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
               WORK( I-LL ) = CS
               WORK( I-LL+NM1 ) = -SN
               WORK( I-LL+NM12 ) = OLDCS
               WORK( I-LL+NM13 ) = -OLDSN
            end do
            H = D( LL )*CS
            D( LL ) = H*OLDCS
            E( LL ) = H*OLDSN
!
!           Update singular vectors
!
            if ( NCVT.GT.0 ) &
               CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                           WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            if ( NRU.GT.0 ) &
               CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), &
                           WORK( N ), U( 1, LL ), LDU )
            if ( NCC.GT.0 ) &
               CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                           WORK( N ), C( LL, 1 ), LDC )
!
!           Test convergence
!
            if ( ABS( E( LL ) ).LE.THRESH ) &
               E( LL ) = ZERO
         end if
      ELSE
!
!        Use nonzero shift
!
         if ( IDIR == 1 ) then
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
            F = ( ABS( D( LL ) )-SHIFT )* &
                ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
            G = E( LL )
            DO I = LL, M - 1
               CALL DLARTG( F, G, COSR, SINR, R )
               if ( I.GT.LL ) &
                  E( I-1 ) = R
               F = COSR*D( I ) + SINR*E( I )
               E( I ) = COSR*E( I ) - SINR*D( I )
               G = SINR*D( I+1 )
               D( I+1 ) = COSR*D( I+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I ) + SINL*D( I+1 )
               D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
               if ( I < M-1 ) then
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
               end if
               WORK( I-LL+1 ) = COSR
               WORK( I-LL+1+NM1 ) = SINR
               WORK( I-LL+1+NM12 ) = COSL
               WORK( I-LL+1+NM13 ) = SINL
            end do
            E( M-1 ) = F
!
!           Update singular vectors
!
            if ( NCVT.GT.0 ) &
               CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), &
                           WORK( N ), VT( LL, 1 ), LDVT )
            if ( NRU.GT.0 ) &
               CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), &
                           WORK( NM13+1 ), U( 1, LL ), LDU )
            if ( NCC.GT.0 ) &
               CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), &
                           WORK( NM13+1 ), C( LL, 1 ), LDC )
!
!           Test convergence
!
            if ( ABS( E( M-1 ) ).LE.THRESH ) &
               E( M-1 ) = ZERO
!
         ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT / &
                D( M ) )
            G = E( M-1 )
            DO I = M, LL + 1, -1
               CALL DLARTG( F, G, COSR, SINR, R )
               if ( I < M ) &
                  E( I ) = R
               F = COSR*D( I ) + SINR*E( I-1 )
               E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
               G = SINR*D( I-1 )
               D( I-1 ) = COSR*D( I-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I-1 ) + SINL*D( I-1 )
               D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
               if ( I.GT.LL+1 ) then
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
               end if
               WORK( I-LL ) = COSR
               WORK( I-LL+NM1 ) = -SINR
               WORK( I-LL+NM12 ) = COSL
               WORK( I-LL+NM13 ) = -SINL
            end do
            E( LL ) = F
!
!           Test convergence
!
            if ( ABS( E( LL ) ).LE.THRESH ) &
               E( LL ) = ZERO
!
!           Update singular vectors if desired
!
            if ( NCVT.GT.0 ) &
               CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                           WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            if ( NRU.GT.0 ) &
               CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), &
                           WORK( N ), U( 1, LL ), LDU )
            if ( NCC.GT.0 ) &
               CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                           WORK( N ), C( LL, 1 ), LDC )
         end if
      end if
!
!     QR iteration finished, go back and check convergence
!
      GO TO 60
!
!     All singular values converged, so make them positive
!
  160 CONTINUE
      DO I = 1, N
         if ( D( I ) < ZERO ) then
            D( I ) = -D( I )
!
!           Change sign of singular vectors, if desired
!
            if ( NCVT.GT.0 ) &
               CALL DSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
         end if
      end do
!
!     Sort the singular values into decreasing order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
      DO I = 1, N - 1
!
!        Scan for smallest D(I)
!
         ISUB = 1
         SMIN = D( 1 )
         DO J = 2, N + 1 - I
            if ( D( J ).LE.SMIN ) then
               ISUB = J
               SMIN = D( J )
            end if
         end do
         if ( ISUB.NE.N+1-I ) then
!
!           Swap singular values and vectors
!
            D( ISUB ) = D( N+1-I )
            D( N+1-I ) = SMIN
            if ( NCVT.GT.0 ) &
               CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ), &
                           LDVT )
            if ( NRU.GT.0 ) &
               CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
            if ( NCC.GT.0 ) &
               CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
         end if
      end do
      GO TO 220
!
!     Maximum number of iterations exceeded, failure to converge
!
  200 CONTINUE
      INFO = 0
      DO I = 1, N - 1
         if ( E( I ).NE.ZERO ) &
            INFO = INFO + 1
      end do
  220 CONTINUE
      RETURN
!
!     End of DBDSQR
!
      END
