      SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          VECT
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGBR generates one of the real orthogonal matrices Q or P**T
!  determined by DGEBRD when reducing a real matrix A to bidiagonal
!  form: A = Q * B * P**T.  Q and P**T are defined as products of
!  elementary reflectors H(i) or G(i) respectively.
!
!  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
!  is of order M:
!  if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
!  columns of Q, where m >= n >= k;
!  if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
!  M-by-M matrix.
!
!  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
!  is of order N:
!  if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m
!  rows of P**T, where n >= m >= k;
!  if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as
!  an N-by-N matrix.
!
!  Arguments
!  =========
!
!  VECT    (input) CHARACTER*1
!          Specifies whether the matrix Q or the matrix P**T is
!          required, as defined in the transformation applied by DGEBRD:
!          = 'Q':  generate Q;
!          = 'P':  generate P**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q or P**T to be returned.
!          M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q or P**T to be returned.
!          N >= 0.
!          If VECT = 'Q', M >= N >= min(M,K);
!          if VECT = 'P', N >= M >= min(N,K).
!
!  K       (input) INTEGER
!          If VECT = 'Q', the number of columns in the original M-by-K
!          matrix reduced by DGEBRD.
!          If VECT = 'P', the number of rows in the original K-by-N
!          matrix reduced by DGEBRD.
!          K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by DGEBRD.
!          On exit, the M-by-N matrix Q or P**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension
!                                (min(M,K)) if VECT = 'Q'
!                                (min(N,K)) if VECT = 'P'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i) or G(i), which determines Q or P**T, as
!          returned by DGEBRD in its array argument TAUQ or TAUP.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
!          For optimum performance LWORK >= min(M,N)*NB, where NB
!          is the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, WANTQ
      INTEGER            I, IINFO, J, LWKOPT, MN, NB
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORGLQ, DORGQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) then
         INFO = -1
      else if ( M < 0 ) then
         INFO = -2
      else if ( N < 0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M, &
               K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M <  &
               MIN( N, K ) ) ) ) then
         INFO = -3
      else if ( K < 0 ) then
         INFO = -4
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -6
      else if ( LWORK < MAX( 1, MN ) .AND. .NOT.LQUERY ) then
         INFO = -9
      end if
!
      if ( INFO == 0 ) then
         if ( WANTQ ) then
            NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
         ELSE
            NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
         end if
         LWKOPT = MAX( 1, MN )*NB
         WORK( 1 ) = LWKOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DORGBR', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 ) then
         WORK( 1 ) = 1
         RETURN
      end if
!
      if ( WANTQ ) then
!
!        Form Q, determined by a call to DGEBRD to reduce an m-by-k
!        matrix
!
         if ( M.GE.K ) then
!
!           If m >= k, assume m >= n >= k
!
            CALL DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
!
         ELSE
!
!           If m < k, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           column to the right, and set the first row and column of Q
!           to those of the unit matrix
!
            DO 20 J = M, 2, -1
               A( 1, J ) = ZERO
               DO 10 I = J + 1, M
                  A( I, J ) = A( I, J-1 )
   10          CONTINUE
   20       CONTINUE
            A( 1, 1 ) = ONE
            DO 30 I = 2, M
               A( I, 1 ) = ZERO
   30       CONTINUE
            if ( M.GT.1 ) then
!
!              Form Q(2:m,2:m)
!
               CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, &
                            LWORK, IINFO )
            end if
         end if
      ELSE
!
!        Form P', determined by a call to DGEBRD to reduce a k-by-n
!        matrix
!
         if ( K < N ) then
!
!           If k < n, assume k <= m <= n
!
            CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
!
         ELSE
!
!           If k >= n, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           row downward, and set the first row and column of P' to
!           those of the unit matrix
!
            A( 1, 1 ) = ONE
            DO 40 I = 2, N
               A( I, 1 ) = ZERO
   40       CONTINUE
            DO 60 J = 2, N
               DO 50 I = J - 1, 2, -1
                  A( I, J ) = A( I-1, J )
   50          CONTINUE
               A( 1, J ) = ZERO
   60       CONTINUE
            if ( N.GT.1 ) then
!
!              Form P'(2:n,2:n)
!
               CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, &
                            LWORK, IINFO )
            end if
         end if
      end if
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORGBR
!
      END
