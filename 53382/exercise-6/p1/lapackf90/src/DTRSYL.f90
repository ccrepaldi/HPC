      SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, &
                         LDC, SCALE, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSYL solves the real Sylvester matrix equation:
!
!     op(A)*X + X*op(B) = scale*C or
!     op(A)*X - X*op(B) = scale*C,
!
!  where op(A) = A or A**T, and  A and B are both upper quasi-
!  triangular. A is M-by-M and B is N-by-N; the right hand side C and
!  the solution X are M-by-N; and scale is an output scale factor, set
!  <= 1 to avoid overflow in X.
!
!  A and B must be in Schur canonical form (as returned by DHSEQR), that
!  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!  each 2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  TRANA   (input) CHARACTER*1
!          Specifies the option op(A):
!          = 'N': op(A) = A    (No transpose)
!          = 'T': op(A) = A**T (Transpose)
!          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
!
!  TRANB   (input) CHARACTER*1
!          Specifies the option op(B):
!          = 'N': op(B) = B    (No transpose)
!          = 'T': op(B) = B**T (Transpose)
!          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
!
!  ISGN    (input) INTEGER
!          Specifies the sign in the equation:
!          = +1: solve op(A)*X + X*op(B) = scale*C
!          = -1: solve op(A)*X - X*op(B) = scale*C
!
!  M       (input) INTEGER
!          The order of the matrix A, and the number of rows in the
!          matrices X and C. M >= 0.
!
!  N       (input) INTEGER
!          The order of the matrix B, and the number of columns in the
!          matrices X and C. N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,M)
!          The upper quasi-triangular matrix A, in Schur canonical form.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
!          The upper quasi-triangular matrix B, in Schur canonical form.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= max(1,N).
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N right hand side matrix C.
!          On exit, C is overwritten by the solution matrix X.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M)
!
!  SCALE   (output) DOUBLE PRECISION
!          The scale factor, scale, set <= 1 to avoid overflow in X.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          = 1: A and B have common or very close eigenvalues; perturbed
!               values were used to solve the equation (but the matrices
!               A and B are unchanged).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT
      DOUBLE PRECISION   A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, &
                         SMLNUM, SUML, SUMR, XNORM
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 ), VEC( 2, 2 ), X( 2, 2 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANGE
      EXTERNAL           LSAME, DDOT, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLALN2, DLASY2, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and Test input parameters
!
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
!
      INFO = 0
      if ( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT. &
          LSAME( TRANA, 'C' ) ) then
         INFO = -1
      else if ( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT. &
               LSAME( TRANB, 'C' ) ) then
         INFO = -2
      else if ( ISGN.NE.1 .AND. ISGN.NE.-1 ) then
         INFO = -3
      else if ( M < 0 ) then
         INFO = -4
      else if ( N < 0 ) then
         INFO = -5
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -7
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -9
      else if ( LDC < MAX( 1, M ) ) then
         INFO = -11
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DTRSYL', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 ) &
         RETURN
!
!     Set constants to control overflow
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
!
      SMIN = MAX( SMLNUM, EPS*DLANGE( 'M', M, M, A, LDA, DUM ), &
             EPS*DLANGE( 'M', N, N, B, LDB, DUM ) )
!
      SCALE = ONE
      SGN = ISGN
!
      if ( NOTRNA .AND. NOTRNB ) then
!
!        Solve    A*X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                  M                         L-1
!        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
!                I=K+1                       J=1
!
!        Start column loop (index = L)
!        L1 (L2) : column index of the first (first) row of X(K,L).
!
         LNEXT = 1
         DO 60 L = 1, N
            if ( L < LNEXT ) &
               GO TO 60
            if ( L == N ) then
               L1 = L
               L2 = L
            ELSE
               if ( B( L+1, L ).NE.ZERO ) then
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               end if
            end if
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L).
!
            KNEXT = M
            DO 50 K = M, 1, -1
               if ( K.GT.KNEXT ) &
                  GO TO 50
               if ( K == 1 ) then
                  K1 = K
                  K2 = K
               ELSE
                  if ( A( K, K-1 ).NE.ZERO ) then
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  end if
               end if
!
               if ( L1 == L2 .AND. K1.EQ.K2 ) then
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) then
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  end if
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11 < ONE .AND. DB.GT.ONE ) then
                     if ( DB.GT.BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  end if
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  if ( SCALOC.NE.ONE ) then
                     DO 10 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
!
               else if ( L1 == L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 20 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1 == K2 ) then
!
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                         C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 30 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   30                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL DLASY2( .FALSE., .FALSE., ISGN, 2, 2, &
                               A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, &
                               2, SCALOC, X, 2, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 40 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               end if
!
   50       CONTINUE
!
   60    CONTINUE
!
      else if ( .NOT.NOTRNA .AND. NOTRNB ) then
!
!        Solve    A' *X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        upper-left corner column by column by
!
!          A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                   K-1                        L-1
!          R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
!                   I=1                        J=1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         LNEXT = 1
         DO 120 L = 1, N
            if ( L < LNEXT ) &
               GO TO 120
            if ( L == N ) then
               L1 = L
               L2 = L
            ELSE
               if ( B( L+1, L ).NE.ZERO ) then
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               end if
            end if
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
            KNEXT = 1
            DO 110 K = 1, M
               if ( K < KNEXT ) &
                  GO TO 110
               if ( K == M ) then
                  K1 = K
                  K2 = K
               ELSE
                  if ( A( K+1, K ).NE.ZERO ) then
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  end if
               end if
!
               if ( L1 == L2 .AND. K1.EQ.K2 ) then
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) then
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  end if
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11 < ONE .AND. DB.GT.ONE ) then
                     if ( DB.GT.BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  end if
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  if ( SCALOC.NE.ONE ) then
                     DO 70 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   70                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
!
               else if ( L1 == L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 80 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1 == K2 ) then
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 90 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL DLASY2( .TRUE., .FALSE., ISGN, 2, 2, A( K1, K1 ), &
                               LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, &
                               2, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 100 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               end if
!
  110       CONTINUE
  120    CONTINUE
!
      else if ( .NOT.NOTRNA .AND. .NOT.NOTRNB ) then
!
!        Solve    A'*X + ISGN*X*B' = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        top-right corner column by column by
!
!           A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
!
!        Where
!                     K-1                          N
!            R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
!                     I=1                        J=L+1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         LNEXT = N
         DO 180 L = N, 1, -1
            if ( L.GT.LNEXT ) &
               GO TO 180
            if ( L == 1 ) then
               L1 = L
               L2 = L
            ELSE
               if ( B( L, L-1 ).NE.ZERO ) then
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               end if
            end if
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
            KNEXT = 1
            DO 170 K = 1, M
               if ( K < KNEXT ) &
                  GO TO 170
               if ( K == M ) then
                  K1 = K
                  K2 = K
               ELSE
                  if ( A( K+1, K ).NE.ZERO ) then
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  end if
               end if
!
               if ( L1 == L2 .AND. K1.EQ.K2 ) then
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC, &
                         B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) then
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  end if
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11 < ONE .AND. DB.GT.ONE ) then
                     if ( DB.GT.BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  end if
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  if ( SCALOC.NE.ONE ) then
                     DO 130 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  130                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
!
               else if ( L1 == L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 140 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  140                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1 == K2 ) then
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 150 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  150                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL DLASY2( .TRUE., .TRUE., ISGN, 2, 2, A( K1, K1 ), &
                               LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, &
                               2, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 160 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  160                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               end if
!
  170       CONTINUE
  180    CONTINUE
!
      else if ( NOTRNA .AND. .NOT.NOTRNB ) then
!
!        Solve    A*X + ISGN*X*B' = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-right corner column by column by
!
!            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
!
!        Where
!                      M                          N
!            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
!                    I=K+1                      J=L+1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         LNEXT = N
         DO 240 L = N, 1, -1
            if ( L.GT.LNEXT ) &
               GO TO 240
            if ( L == 1 ) then
               L1 = L
               L2 = L
            ELSE
               if ( B( L, L-1 ).NE.ZERO ) then
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               end if
            end if
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
            KNEXT = M
            DO 230 K = M, 1, -1
               if ( K.GT.KNEXT ) &
                  GO TO 230
               if ( K == 1 ) then
                  K1 = K
                  K2 = K
               ELSE
                  if ( A( K, K-1 ).NE.ZERO ) then
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  end if
               end if
!
               if ( L1 == L2 .AND. K1.EQ.K2 ) then
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC, &
                         B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) then
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  end if
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11 < ONE .AND. DB.GT.ONE ) then
                     if ( DB.GT.BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  end if
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  if ( SCALOC.NE.ONE ) then
                     DO 190 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  190                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
!
               else if ( L1 == L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 200 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  200                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1 == K2 ) then
!
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                         C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 210 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  210                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               else if ( L1.NE.L2 .AND. K1.NE.K2 ) then
!
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL DLASY2( .FALSE., .TRUE., ISGN, 2, 2, A( K1, K1 ), &
                               LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, &
                               2, XNORM, IERR )
                  if ( IERR.NE.0 ) &
                     INFO = 1
!
                  if ( SCALOC.NE.ONE ) then
                     DO 220 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  220                CONTINUE
                     SCALE = SCALE*SCALOC
                  end if
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               end if
!
  230       CONTINUE
  240    CONTINUE
!
      end if
!
      RETURN
!
!     End of DTRSYL
!
      END
