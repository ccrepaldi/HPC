      SUBROUTINE DTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, &
                         LDZ, J1, N1, N2, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, LWORK, N, N1, N2
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), &
                         WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DTGEX2 swaps adjacent diagonal blocks (A11, B11) and (A22, B22)
!  of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
!  (A, B) by an orthogonal equivalence transformation.
!
!  (A, B) must be in generalized real Schur canonical form (as returned
!  by DGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
!  diagonal blocks. B is upper triangular.
!
!  Optionally, the matrices Q and Z of generalized Schur vectors are
!  updated.
!
!         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'
!         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'
!
!
!  Arguments
!  =========
!
!  WANTQ   (input) LOGICAL
!          .TRUE. : update the left transformation matrix Q;
!          .FALSE.: do not update Q.
!
!  WANTZ   (input) LOGICAL
!          .TRUE. : update the right transformation matrix Z;
!          .FALSE.: do not update Z.
!
!  N       (input) INTEGER
!          The order of the matrices A and B. N >= 0.
!
!  A      (input/output) DOUBLE PRECISION arrays, dimensions (LDA,N)
!          On entry, the matrix A in the pair (A, B).
!          On exit, the updated matrix A.
!
!  LDA     (input)  INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  B      (input/output) DOUBLE PRECISION arrays, dimensions (LDB,N)
!          On entry, the matrix B in the pair (A, B).
!          On exit, the updated matrix B.
!
!  LDB     (input)  INTEGER
!          The leading dimension of the array B. LDB >= max(1,N).
!
!  Q       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
!          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.
!          On exit, the updated matrix Q.
!          Not referenced if WANTQ = .FALSE..
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q. LDQ >= 1.
!          If WANTQ = .TRUE., LDQ >= N.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
!          On entry, if WANTZ =.TRUE., the orthogonal matrix Z.
!          On exit, the updated matrix Z.
!          Not referenced if WANTZ = .FALSE..
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z. LDZ >= 1.
!          If WANTZ = .TRUE., LDZ >= N.
!
!  J1      (input) INTEGER
!          The index to the first block (A11, B11). 1 <= J1 <= N.
!
!  N1      (input) INTEGER
!          The order of the first block (A11, B11). N1 = 0, 1 or 2.
!
!  N2      (input) INTEGER
!          The order of the second block (A22, B22). N2 = 0, 1 or 2.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK).
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          LWORK >=  MAX( N*(N2+N1), (N2+N1)*(N2+N1)*2 )
!
!  INFO    (output) INTEGER
!            =0: Successful exit
!            >0: If INFO = 1, the transformed matrix (A, B) would be
!                too far from generalized Schur form; the blocks are
!                not swapped and (A, B) and (Q, Z) are unchanged.
!                The problem of swapping is too ill-conditioned.
!            <0: If INFO = -16: LWORK is too small. Appropriate value
!                for LWORK is returned in WORK(1).
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!     Umea University, S-901 87 Umea, Sweden.
!
!  In the current code both weak and strong stability tests are
!  performed. The user can omit the strong stability test by changing
!  the internal logical parameter WANDS to .FALSE.. See ref. [2] for
!  details.
!
!  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
!      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
!      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
!      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
!
!  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
!      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
!      Estimation: Theory, Algorithms and Software,
!      Report UMINF - 94.04, Department of Computing Science, Umea
!      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working
!      Note 87. To appear in Numerical Algorithms, 1996.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 1.0D+01 )
      INTEGER            LDST
      PARAMETER          ( LDST = 4 )
      LOGICAL            WANDS
      PARAMETER          ( WANDS = .TRUE. )
!     ..
!     .. Local Scalars ..
      LOGICAL            DTRONG, WEAK
      INTEGER            I, IDUM, LINFO, M
      DOUBLE PRECISION   BQRA21, BRQA21, DDUM, DNORM, DSCALE, DSUM, EPS, &
                         F, G, SA, SB, SCALE, SMLNUM, SS, THRESH, WS
!     ..
!     .. Local Arrays ..
      INTEGER            IWORK( LDST )
      DOUBLE PRECISION   AI( 2 ), AR( 2 ), BE( 2 ), IR( LDST, LDST ), &
                         IRCOP( LDST, LDST ), LI( LDST, LDST ), &
                         LICOP( LDST, LDST ), S( LDST, LDST ), &
                         SCPY( LDST, LDST ), T( LDST, LDST ), &
                         TAUL( LDST ), TAUR( LDST ), TCPY( LDST, LDST )
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGEQR2, DGERQ2, DLACPY, DLAGV2, &
                         DLARTG, DLASSQ, DORG2R, DORGR2, DORM2R, DORMR2, &
                         DROT, DSCAL, DTGSY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      if ( N.LE.1 .OR. N1.LE.0 .OR. N2.LE.0 ) &
         RETURN
      if ( N1.GT.N .OR. ( J1+N1 ).GT.N ) &
         RETURN
      M = N1 + N2
      if ( LWORK < MAX( N*M, M*M*2 ) ) then
         INFO = -16
         WORK( 1 ) = MAX( N*M, M*M*2 )
         RETURN
      end if
!
      WEAK = .FALSE.
      DTRONG = .FALSE.
!
!     Make a local copy of selected block
!
      CALL DCOPY( LDST*LDST, ZERO, 0, LI, 1 )
      CALL DCOPY( LDST*LDST, ZERO, 0, IR, 1 )
      CALL DLACPY( 'Full', M, M, A( J1, J1 ), LDA, S, LDST )
      CALL DLACPY( 'Full', M, M, B( J1, J1 ), LDB, T, LDST )
!
!     Compute threshold for testing acceptance of swapping.
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      DSCALE = ZERO
      DSUM = ONE
      CALL DLACPY( 'Full', M, M, S, LDST, WORK, M )
      CALL DLASSQ( M*M, WORK, 1, DSCALE, DSUM )
      CALL DLACPY( 'Full', M, M, T, LDST, WORK, M )
      CALL DLASSQ( M*M, WORK, 1, DSCALE, DSUM )
      DNORM = DSCALE*SQRT( DSUM )
      THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
!
      if ( M == 2 ) then
!
!        CASE 1: Swap 1-by-1 and 1-by-1 blocks.
!
!        Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks
!        using Givens rotations and perform the swap tentatively.
!
         F = S( 2, 2 )*T( 1, 1 ) - T( 2, 2 )*S( 1, 1 )
         G = S( 2, 2 )*T( 1, 2 ) - T( 2, 2 )*S( 1, 2 )
         SB = ABS( T( 2, 2 ) )
         SA = ABS( S( 2, 2 ) )
         CALL DLARTG( F, G, IR( 1, 2 ), IR( 1, 1 ), DDUM )
         IR( 2, 1 ) = -IR( 1, 2 )
         IR( 2, 2 ) = IR( 1, 1 )
         CALL DROT( 2, S( 1, 1 ), 1, S( 1, 2 ), 1, IR( 1, 1 ), &
                    IR( 2, 1 ) )
         CALL DROT( 2, T( 1, 1 ), 1, T( 1, 2 ), 1, IR( 1, 1 ), &
                    IR( 2, 1 ) )
         if ( SA.GE.SB ) then
            CALL DLARTG( S( 1, 1 ), S( 2, 1 ), LI( 1, 1 ), LI( 2, 1 ), &
                         DDUM )
         ELSE
            CALL DLARTG( T( 1, 1 ), T( 2, 1 ), LI( 1, 1 ), LI( 2, 1 ), &
                         DDUM )
         end if
         CALL DROT( 2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, LI( 1, 1 ), &
                    LI( 2, 1 ) )
         CALL DROT( 2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, LI( 1, 1 ), &
                    LI( 2, 1 ) )
         LI( 2, 2 ) = LI( 1, 1 )
         LI( 1, 2 ) = -LI( 2, 1 )
!
!        Weak stability test:
!           |S21| + |T21| <= O(EPS * F-norm((S, T)))
!
         WS = ABS( S( 2, 1 ) ) + ABS( T( 2, 1 ) )
         WEAK = WS.LE.THRESH
         if ( .NOT.WEAK ) &
            GO TO 70
!
         if ( WANDS ) then
!
!           Strong stability test:
!             F-norm((A-QL'*S*QR, B-QL'*T*QR)) <= O(EPS*F-norm((A,B)))
!
            CALL DLACPY( 'Full', M, M, A( J1, J1 ), LDA, WORK( M*M+1 ), &
                         M )
            CALL DGEMM( 'N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, &
                        WORK, M )
            CALL DGEMM( 'N', 'T', M, M, M, -ONE, WORK, M, IR, LDST, ONE, &
                        WORK( M*M+1 ), M )
            DSCALE = ZERO
            DSUM = ONE
            CALL DLASSQ( M*M, WORK( M*M+1 ), 1, DSCALE, DSUM )
!
            CALL DLACPY( 'Full', M, M, B( J1, J1 ), LDB, WORK( M*M+1 ), &
                         M )
            CALL DGEMM( 'N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, &
                        WORK, M )
            CALL DGEMM( 'N', 'T', M, M, M, -ONE, WORK, M, IR, LDST, ONE, &
                        WORK( M*M+1 ), M )
            CALL DLASSQ( M*M, WORK( M*M+1 ), 1, DSCALE, DSUM )
            SS = DSCALE*SQRT( DSUM )
            DTRONG = SS.LE.THRESH
            if ( .NOT.DTRONG ) &
               GO TO 70
         end if
!
!        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
!               (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
!
         CALL DROT( J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, IR( 1, 1 ), &
                    IR( 2, 1 ) )
         CALL DROT( J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, IR( 1, 1 ), &
                    IR( 2, 1 ) )
         CALL DROT( N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, &
                    LI( 1, 1 ), LI( 2, 1 ) )
         CALL DROT( N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, &
                    LI( 1, 1 ), LI( 2, 1 ) )
!
!        Set  N1-by-N2 (2,1) - blocks to ZERO.
!
         A( J1+1, J1 ) = ZERO
         B( J1+1, J1 ) = ZERO
!
!        Accumulate transformations into Q and Z if requested.
!
         if ( WANTZ ) &
            CALL DROT( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, IR( 1, 1 ), &
                       IR( 2, 1 ) )
         if ( WANTQ ) &
            CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, LI( 1, 1 ), &
                       LI( 2, 1 ) )
!
!        Exit with INFO = 0 if swap was successfully performed.
!
         RETURN
!
      ELSE
!
!        CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
!                and 2-by-2 blocks.
!
!        Solve the generalized Sylvester equation
!                 S11 * R - L * S22 = SCALE * S12
!                 T11 * R - L * T22 = SCALE * T12
!        for R and L. Solutions in LI and IR.
!
         CALL DLACPY( 'Full', N1, N2, T( 1, N1+1 ), LDST, LI, LDST )
         CALL DLACPY( 'Full', N1, N2, S( 1, N1+1 ), LDST, &
                      IR( N2+1, N1+1 ), LDST )
         CALL DTGSY2( 'N', 0, N1, N2, S, LDST, S( N1+1, N1+1 ), LDST, &
                      IR( N2+1, N1+1 ), LDST, T, LDST, T( N1+1, N1+1 ), &
                      LDST, LI, LDST, SCALE, DSUM, DSCALE, IWORK, IDUM, &
                      LINFO )
!
!        Compute orthogonal matrix QL:
!
!                    QL' * LI = [ TL ]
!                               [ 0  ]
!        where
!                    LI =  [      -L              ]
!                          [ SCALE * identity(N2) ]
!
         DO 10 I = 1, N2
            CALL DSCAL( N1, -ONE, LI( 1, I ), 1 )
            LI( N1+I, I ) = SCALE
   10    CONTINUE
         CALL DGEQR2( M, N2, LI, LDST, TAUL, WORK, LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
         CALL DORG2R( M, M, N2, LI, LDST, TAUL, WORK, LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
!
!        Compute orthogonal matrix RQ:
!
!                    IR * RQ' =   [ 0  TR],
!
!         where IR = [ SCALE * identity(N1), R ]
!
         DO 20 I = 1, N1
            IR( N2+I, I ) = SCALE
   20    CONTINUE
         CALL DGERQ2( N1, M, IR( N2+1, 1 ), LDST, TAUR, WORK, LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
         CALL DORGR2( M, M, N1, IR, LDST, TAUR, WORK, LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
!
!        Perform the swapping tentatively:
!
         CALL DGEMM( 'T', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, &
                     WORK, M )
         CALL DGEMM( 'N', 'T', M, M, M, ONE, WORK, M, IR, LDST, ZERO, S, &
                     LDST )
         CALL DGEMM( 'T', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, &
                     WORK, M )
         CALL DGEMM( 'N', 'T', M, M, M, ONE, WORK, M, IR, LDST, ZERO, T, &
                     LDST )
         CALL DLACPY( 'F', M, M, S, LDST, SCPY, LDST )
         CALL DLACPY( 'F', M, M, T, LDST, TCPY, LDST )
         CALL DLACPY( 'F', M, M, IR, LDST, IRCOP, LDST )
         CALL DLACPY( 'F', M, M, LI, LDST, LICOP, LDST )
!
!        Triangularize the B-part by an RQ factorization.
!        Apply transformation (from left) to A-part, giving S.
!
         CALL DGERQ2( M, M, T, LDST, TAUR, WORK, LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
         CALL DORMR2( 'R', 'T', M, M, M, T, LDST, TAUR, S, LDST, WORK, &
                      LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
         CALL DORMR2( 'L', 'N', M, M, M, T, LDST, TAUR, IR, LDST, WORK, &
                      LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
!
!        Compute F-norm(S21) in BRQA21. (T21 is 0.)
!
         DSCALE = ZERO
         DSUM = ONE
         DO 30 I = 1, N2
            CALL DLASSQ( N1, S( N2+1, I ), 1, DSCALE, DSUM )
   30    CONTINUE
         BRQA21 = DSCALE*SQRT( DSUM )
!
!        Triangularize the B-part by a QR factorization.
!        Apply transformation (from right) to A-part, giving S.
!
         CALL DGEQR2( M, M, TCPY, LDST, TAUL, WORK, LINFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
         CALL DORM2R( 'L', 'T', M, M, M, TCPY, LDST, TAUL, SCPY, LDST, &
                      WORK, INFO )
         CALL DORM2R( 'R', 'N', M, M, M, TCPY, LDST, TAUL, LICOP, LDST, &
                      WORK, INFO )
         if ( LINFO.NE.0 ) &
            GO TO 70
!
!        Compute F-norm(S21) in BQRA21. (T21 is 0.)
!
         DSCALE = ZERO
         DSUM = ONE
         DO 40 I = 1, N2
            CALL DLASSQ( N1, SCPY( N2+1, I ), 1, DSCALE, DSUM )
   40    CONTINUE
         BQRA21 = DSCALE*SQRT( DSUM )
!
!        Decide which method to use.
!          Weak stability test:
!             F-norm(S21) <= O(EPS * F-norm((S, T)))
!
         if ( BQRA21.LE.BRQA21 .AND. BQRA21.LE.THRESH ) then
            CALL DLACPY( 'F', M, M, SCPY, LDST, S, LDST )
            CALL DLACPY( 'F', M, M, TCPY, LDST, T, LDST )
            CALL DLACPY( 'F', M, M, IRCOP, LDST, IR, LDST )
            CALL DLACPY( 'F', M, M, LICOP, LDST, LI, LDST )
         else if ( BRQA21.GE.THRESH ) then
            GO TO 70
         end if
!
!        Set lower triangle of B-part to zero
!
         DO 50 I = 2, M
            CALL DCOPY( M-I+1, ZERO, 0, T( I, I-1 ), 1 )
   50    CONTINUE
!
         if ( WANDS ) then
!
!           Strong stability test:
!              F-norm((A-QL*S*QR', B-QL*T*QR')) <= O(EPS*F-norm((A,B)))
!
            CALL DLACPY( 'Full', M, M, A( J1, J1 ), LDA, WORK( M*M+1 ), &
                         M )
            CALL DGEMM( 'N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, &
                        WORK, M )
            CALL DGEMM( 'N', 'N', M, M, M, -ONE, WORK, M, IR, LDST, ONE, &
                        WORK( M*M+1 ), M )
            DSCALE = ZERO
            DSUM = ONE
            CALL DLASSQ( M*M, WORK( M*M+1 ), 1, DSCALE, DSUM )
!
            CALL DLACPY( 'Full', M, M, B( J1, J1 ), LDB, WORK( M*M+1 ), &
                         M )
            CALL DGEMM( 'N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, &
                        WORK, M )
            CALL DGEMM( 'N', 'N', M, M, M, -ONE, WORK, M, IR, LDST, ONE, &
                        WORK( M*M+1 ), M )
            CALL DLASSQ( M*M, WORK( M*M+1 ), 1, DSCALE, DSUM )
            SS = DSCALE*SQRT( DSUM )
            DTRONG = ( SS.LE.THRESH )
            if ( .NOT.DTRONG ) &
               GO TO 70
!
         end if
!
!        If the swap is accepted ("weakly" and "strongly"), apply the
!        transformations and set N1-by-N2 (2,1)-block to zero.
!
         DO 60 I = 1, N2
            CALL DCOPY( N1, ZERO, 0, S( N2+1, I ), 1 )
   60    CONTINUE
!
!        copy back M-by-M diagonal block starting at index J1 of (A, B)
!
         CALL DLACPY( 'F', M, M, S, LDST, A( J1, J1 ), LDA )
         CALL DLACPY( 'F', M, M, T, LDST, B( J1, J1 ), LDB )
         CALL DCOPY( LDST*LDST, ZERO, 0, T, 1 )
!
!        Standardize existing 2-by-2 blocks.
!
         CALL DCOPY( M*M, ZERO, 0, WORK, 1 )
         WORK( 1 ) = ONE
         T( 1, 1 ) = ONE
         IDUM = LWORK - M*M - 2
         if ( N2.GT.1 ) then
            CALL DLAGV2( A( J1, J1 ), LDA, B( J1, J1 ), LDB, AR, AI, BE, &
                         WORK( 1 ), WORK( 2 ), T( 1, 1 ), T( 2, 1 ) )
            WORK( M+1 ) = -WORK( 2 )
            WORK( M+2 ) = WORK( 1 )
            T( N2, N2 ) = T( 1, 1 )
            T( 1, 2 ) = -T( 2, 1 )
         end if
         WORK( M*M ) = ONE
         T( M, M ) = ONE
!
         if ( N1.GT.1 ) then
            CALL DLAGV2( A( J1+N2, J1+N2 ), LDA, B( J1+N2, J1+N2 ), LDB, &
                         TAUR, TAUL, WORK( M*M+1 ), WORK( N2*M+N2+1 ), &
                         WORK( N2*M+N2+2 ), T( N2+1, N2+1 ), &
                         T( M, M-1 ) )
            WORK( M*M ) = WORK( N2*M+N2+1 )
            WORK( M*M-1 ) = -WORK( N2*M+N2+2 )
            T( M, M ) = T( N2+1, N2+1 )
            T( M-1, M ) = -T( M, M-1 )
         end if
         CALL DGEMM( 'T', 'N', N2, N1, N2, ONE, WORK, M, A( J1, J1+N2 ), &
                     LDA, ZERO, WORK( M*M+1 ), N2 )
         CALL DLACPY( 'Full', N2, N1, WORK( M*M+1 ), N2, A( J1, J1+N2 ), &
                      LDA )
         CALL DGEMM( 'T', 'N', N2, N1, N2, ONE, WORK, M, B( J1, J1+N2 ), &
                     LDB, ZERO, WORK( M*M+1 ), N2 )
         CALL DLACPY( 'Full', N2, N1, WORK( M*M+1 ), N2, B( J1, J1+N2 ), &
                      LDB )
         CALL DGEMM( 'N', 'N', M, M, M, ONE, LI, LDST, WORK, M, ZERO, &
                     WORK( M*M+1 ), M )
         CALL DLACPY( 'Full', M, M, WORK( M*M+1 ), M, LI, LDST )
         CALL DGEMM( 'N', 'N', N2, N1, N1, ONE, A( J1, J1+N2 ), LDA, &
                     T( N2+1, N2+1 ), LDST, ZERO, WORK, N2 )
         CALL DLACPY( 'Full', N2, N1, WORK, N2, A( J1, J1+N2 ), LDA )
         CALL DGEMM( 'N', 'N', N2, N1, N1, ONE, B( J1, J1+N2 ), LDA, &
                     T( N2+1, N2+1 ), LDST, ZERO, WORK, N2 )
         CALL DLACPY( 'Full', N2, N1, WORK, N2, B( J1, J1+N2 ), LDB )
         CALL DGEMM( 'T', 'N', M, M, M, ONE, IR, LDST, T, LDST, ZERO, &
                     WORK, M )
         CALL DLACPY( 'Full', M, M, WORK, M, IR, LDST )
!
!        Accumulate transformations into Q and Z if requested.
!
         if ( WANTQ ) then
            CALL DGEMM( 'N', 'N', N, M, M, ONE, Q( 1, J1 ), LDQ, LI, &
                        LDST, ZERO, WORK, N )
            CALL DLACPY( 'Full', N, M, WORK, N, Q( 1, J1 ), LDQ )
!
         end if
!
         if ( WANTZ ) then
            CALL DGEMM( 'N', 'N', N, M, M, ONE, Z( 1, J1 ), LDZ, IR, &
                        LDST, ZERO, WORK, N )
            CALL DLACPY( 'Full', N, M, WORK, N, Z( 1, J1 ), LDZ )
!
         end if
!
!        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
!                (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
!
         I = J1 + M
         if ( I.LE.N ) then
            CALL DGEMM( 'T', 'N', M, N-I+1, M, ONE, LI, LDST, &
                        A( J1, I ), LDA, ZERO, WORK, M )
            CALL DLACPY( 'Full', M, N-I+1, WORK, M, A( J1, I ), LDA )
            CALL DGEMM( 'T', 'N', M, N-I+1, M, ONE, LI, LDST, &
                        B( J1, I ), LDA, ZERO, WORK, M )
            CALL DLACPY( 'Full', M, N-I+1, WORK, M, B( J1, I ), LDA )
         end if
         I = J1 - 1
         if ( I.GT.0 ) then
            CALL DGEMM( 'N', 'N', I, M, M, ONE, A( 1, J1 ), LDA, IR, &
                        LDST, ZERO, WORK, I )
            CALL DLACPY( 'Full', I, M, WORK, I, A( 1, J1 ), LDA )
            CALL DGEMM( 'N', 'N', I, M, M, ONE, B( 1, J1 ), LDB, IR, &
                        LDST, ZERO, WORK, I )
            CALL DLACPY( 'Full', I, M, WORK, I, B( 1, J1 ), LDB )
         end if
!
!        Exit with INFO = 0 if swap was successfully performed.
!
         RETURN
!
      end if
!
!     Exit with INFO = 1 if swap was rejected.
!
   70 CONTINUE
!
      INFO = 1
      RETURN
!
!     End of DTGEX2
!
      END
