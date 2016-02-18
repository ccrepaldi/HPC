      SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, &
                         RSCALE, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ), &
                         RSCALE( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGGBAL balances a pair of general real matrices (A,B).  This
!  involves, first, permuting A and B by similarity transformations to
!  isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N
!  elements on the diagonal; and second, applying a diagonal similarity
!  transformation to rows and columns ILO to IHI to make the rows
!  and columns as close in norm as possible. Both steps are optional.
!
!  Balancing may reduce the 1-norm of the matrices, and improve the
!  accuracy of the computed eigenvalues and/or eigenvectors in the
!  generalized eigenvalue problem A*x = lambda*B*x.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the operations to be performed on A and B:
!          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
!                  and RSCALE(I) = 1.0 for i = 1,...,N.
!          = 'P':  permute only;
!          = 'S':  scale only;
!          = 'B':  both permute and scale.
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the input matrix A.
!          On exit,  A is overwritten by the balanced matrix.
!          If JOB = 'N', A is not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
!          On entry, the input matrix B.
!          On exit,  B is overwritten by the balanced matrix.
!          If JOB = 'N', B is not referenced.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= max(1,N).
!
!  ILO     (output) INTEGER
!  IHI     (output) INTEGER
!          ILO and IHI are set to integers such that on exit
!          A(i,j) = 0 and B(i,j) = 0 if i > j and
!          j = 1,...,ILO-1 or i = IHI+1,...,N.
!          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!
!  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied
!          to the left side of A and B.  If P(j) is the index of the
!          row interchanged with row j, and D(j)
!          is the scaling factor applied to row j, then
!            LSCALE(j) = P(j)    for J = 1,...,ILO-1
!                      = D(j)    for J = ILO,...,IHI
!                      = P(j)    for J = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied
!          to the right side of A and B.  If P(j) is the index of the
!          column interchanged with column j, and D(j)
!          is the scaling factor applied to column j, then
!            LSCALE(j) = P(j)    for J = 1,...,ILO-1
!                      = D(j)    for J = ILO,...,IHI
!                      = P(j)    for J = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  See R.C. WARD, Balancing the generalized eigenvalue problem,
!                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   THREE, SCLFAC
      PARAMETER          ( THREE = 3.0D+0, SCLFAC = 1.0D+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1, &
                         K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN, &
                         M, NR, NRP2
      DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2, &
                         COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX, &
                         SFMIN, SUM, T, TA, TB, TC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DDOT, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG10, MAX, MIN, SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      if ( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -4
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -5
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGGBAL', -INFO )
         RETURN
      end if
!
      K = 1
      L = N
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
      if ( LSAME( JOB, 'N' ) ) then
         ILO = 1
         IHI = N
         DO 10 I = 1, N
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
   10    CONTINUE
         RETURN
      end if
!
      if ( K == L ) then
         ILO = 1
         IHI = 1
         LSCALE( 1 ) = ONE
         RSCALE( 1 ) = ONE
         RETURN
      end if
!
      if ( LSAME( JOB, 'S' ) ) &
         GO TO 190
!
      GO TO 30
!
!     Permute the matrices A and B to isolate the eigenvalues.
!
!     Find row with one nonzero in columns 1 through L
!
   20 CONTINUE
      L = LM1
      if ( L.NE.1 ) &
         GO TO 30
!
      RSCALE( 1 ) = 1
      LSCALE( 1 ) = 1
      GO TO 190
!
   30 CONTINUE
      LM1 = L - 1
      DO 80 I = L, 1, -1
         DO 40 J = 1, LM1
            JP1 = J + 1
            if ( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO ) &
               GO TO 50
   40    CONTINUE
         J = L
         GO TO 70
!
   50    CONTINUE
         DO 60 J = JP1, L
            if ( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO ) &
               GO TO 80
   60    CONTINUE
         J = JP1 - 1
!
   70    CONTINUE
         M = L
         IFLOW = 1
         GO TO 160
   80 CONTINUE
      GO TO 100
!
!     Find column with one nonzero in rows K through N
!
   90 CONTINUE
      K = K + 1
!
  100 CONTINUE
      DO 150 J = K, L
         DO 110 I = K, LM1
            IP1 = I + 1
            if ( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO ) &
               GO TO 120
  110    CONTINUE
         I = L
         GO TO 140
  120    CONTINUE
         DO 130 I = IP1, L
            if ( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO ) &
               GO TO 150
  130    CONTINUE
         I = IP1 - 1
  140    CONTINUE
         M = K
         IFLOW = 2
         GO TO 160
  150 CONTINUE
      GO TO 190
!
!     Permute rows M and I
!
  160 CONTINUE
      LSCALE( M ) = I
      if ( I == M ) &
         GO TO 170
      CALL DSWAP( N-K+1, A( I, K ), LDA, A( M, K ), LDA )
      CALL DSWAP( N-K+1, B( I, K ), LDB, B( M, K ), LDB )
!
!     Permute columns M and J
!
  170 CONTINUE
      RSCALE( M ) = J
      if ( J == M ) &
         GO TO 180
      CALL DSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL DSWAP( L, B( 1, J ), 1, B( 1, M ), 1 )
!
  180 CONTINUE
      GO TO ( 20, 90 )IFLOW
!
  190 CONTINUE
      ILO = K
      IHI = L
!
      if ( ILO == IHI ) &
         RETURN
!
      if ( LSAME( JOB, 'P' ) ) &
         RETURN
!
!     Balance the submatrix in rows ILO to IHI.
!
      NR = IHI - ILO + 1
      DO 200 I = ILO, IHI
         RSCALE( I ) = ZERO
         LSCALE( I ) = ZERO
!
         WORK( I ) = ZERO
         WORK( I+N ) = ZERO
         WORK( I+2*N ) = ZERO
         WORK( I+3*N ) = ZERO
         WORK( I+4*N ) = ZERO
         WORK( I+5*N ) = ZERO
  200 CONTINUE
!
!     Compute right side vector in resulting linear equations
!
      BASL = LOG10( SCLFAC )
      DO 240 I = ILO, IHI
         DO 230 J = ILO, IHI
            TB = B( I, J )
            TA = A( I, J )
            if ( TA == ZERO ) &
               GO TO 210
            TA = LOG10( ABS( TA ) ) / BASL
  210       CONTINUE
            if ( TB == ZERO ) &
               GO TO 220
            TB = LOG10( ABS( TB ) ) / BASL
  220       CONTINUE
            WORK( I+4*N ) = WORK( I+4*N ) - TA - TB
            WORK( J+5*N ) = WORK( J+5*N ) - TA - TB
  230    CONTINUE
  240 CONTINUE
!
      COEF = ONE / DBLE( 2*NR )
      COEF2 = COEF*COEF
      COEF5 = HALF*COEF2
      NRP2 = NR + 2
      BETA = ZERO
      IT = 1
!
!     Start generalized conjugate gradient iteration
!
  250 CONTINUE
!
      GAMMA = DDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) + &
              DDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 )
!
      EW = ZERO
      EWC = ZERO
      DO 260 I = ILO, IHI
         EW = EW + WORK( I+4*N )
         EWC = EWC + WORK( I+5*N )
  260 CONTINUE
!
      GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2
      if ( GAMMA == ZERO ) &
         GO TO 350
      if ( IT.NE.1 ) &
         BETA = GAMMA / PGAMMA
      T = COEF5*( EWC-THREE*EW )
      TC = COEF5*( EW-THREE*EWC )
!
      CALL DSCAL( NR, BETA, WORK( ILO ), 1 )
      CALL DSCAL( NR, BETA, WORK( ILO+N ), 1 )
!
      CALL DAXPY( NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 )
      CALL DAXPY( NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 )
!
      DO 270 I = ILO, IHI
         WORK( I ) = WORK( I ) + TC
         WORK( I+N ) = WORK( I+N ) + T
  270 CONTINUE
!
!     Apply matrix to vector
!
      DO 300 I = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 290 J = ILO, IHI
            if ( A( I, J ) == ZERO ) &
               GO TO 280
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  280       CONTINUE
            if ( B( I, J ) == ZERO ) &
               GO TO 290
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  290    CONTINUE
         WORK( I+2*N ) = DBLE( KOUNT )*WORK( I+N ) + SUM
  300 CONTINUE
!
      DO 330 J = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 320 I = ILO, IHI
            if ( A( I, J ) == ZERO ) &
               GO TO 310
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  310       CONTINUE
            if ( B( I, J ) == ZERO ) &
               GO TO 320
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  320    CONTINUE
         WORK( J+3*N ) = DBLE( KOUNT )*WORK( J ) + SUM
  330 CONTINUE
!
      SUM = DDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) + &
            DDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 )
      ALPHA = GAMMA / SUM
!
!     Determine correction to current iteration
!
      CMAX = ZERO
      DO 340 I = ILO, IHI
         COR = ALPHA*WORK( I+N )
         if ( ABS( COR ).GT.CMAX ) &
            CMAX = ABS( COR )
         LSCALE( I ) = LSCALE( I ) + COR
         COR = ALPHA*WORK( I )
         if ( ABS( COR ).GT.CMAX ) &
            CMAX = ABS( COR )
         RSCALE( I ) = RSCALE( I ) + COR
  340 CONTINUE
      if ( CMAX < HALF ) &
         GO TO 350
!
      CALL DAXPY( NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 )
      CALL DAXPY( NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 )
!
      PGAMMA = GAMMA
      IT = IT + 1
      if ( IT.LE.NRP2 ) &
         GO TO 250
!
!     End generalized conjugate gradient iteration
!
  350 CONTINUE
      SFMIN = DLAMCH( 'S' )
      SFMAX = ONE / SFMIN
      LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE )
      LSFMAX = INT( LOG10( SFMAX ) / BASL )
      DO 360 I = ILO, IHI
         IRAB = IDAMAX( N-ILO+1, A( I, ILO ), LDA )
         RAB = ABS( A( I, IRAB+ILO-1 ) )
         IRAB = IDAMAX( N-ILO+1, B( I, ILO ), LDA )
         RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) )
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE )
         IR = LSCALE( I ) + SIGN( HALF, LSCALE( I ) )
         IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
         LSCALE( I ) = SCLFAC**IR
         ICAB = IDAMAX( IHI, A( 1, I ), 1 )
         CAB = ABS( A( ICAB, I ) )
         ICAB = IDAMAX( IHI, B( 1, I ), 1 )
         CAB = MAX( CAB, ABS( B( ICAB, I ) ) )
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE )
         JC = RSCALE( I ) + SIGN( HALF, RSCALE( I ) )
         JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
         RSCALE( I ) = SCLFAC**JC
  360 CONTINUE
!
!     Row scaling of matrices A and B
!
      DO 370 I = ILO, IHI
         CALL DSCAL( N-ILO+1, LSCALE( I ), A( I, ILO ), LDA )
         CALL DSCAL( N-ILO+1, LSCALE( I ), B( I, ILO ), LDB )
  370 CONTINUE
!
!     Column scaling of matrices A and B
!
      DO 380 J = ILO, IHI
         CALL DSCAL( IHI, RSCALE( J ), A( 1, J ), 1 )
         CALL DSCAL( IHI, RSCALE( J ), B( 1, J ), 1 )
  380 CONTINUE
!
      RETURN
!
!     End of DGGBAL
!
      END
