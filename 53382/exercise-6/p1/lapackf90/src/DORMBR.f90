      SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
                         LDC, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, VECT
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C
!  with
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
!  with
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      P * C          C * P
!  TRANS = 'T':      P**T * C       C * P**T
!
!  Here Q and P**T are the orthogonal matrices determined by DGEBRD when
!  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
!  P**T are defined as products of elementary reflectors H(i) and G(i)
!  respectively.
!
!  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
!  order of the orthogonal matrix Q or P**T that is applied.
!
!  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
!  if nq >= k, Q = H(1) H(2) . . . H(k);
!  if nq < k, Q = H(1) H(2) . . . H(nq-1).
!
!  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
!  if k < nq, P = G(1) G(2) . . . G(k);
!  if k >= nq, P = G(1) G(2) . . . G(nq-1).
!
!  Arguments
!  =========
!
!  VECT    (input) CHARACTER*1
!          = 'Q': apply Q or Q**T;
!          = 'P': apply P or P**T.
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q, Q**T, P or P**T from the Left;
!          = 'R': apply Q, Q**T, P or P**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q  or P;
!          = 'T':  Transpose, apply Q**T or P**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          If VECT = 'Q', the number of columns in the original
!          matrix reduced by DGEBRD.
!          If VECT = 'P', the number of rows in the original
!          matrix reduced by DGEBRD.
!          K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                                (LDA,min(nq,K)) if VECT = 'Q'
!                                (LDA,nq)        if VECT = 'P'
!          The vectors which define the elementary reflectors H(i) and
!          G(i), whose products determine the matrices Q and P, as
!          returned by DGEBRD.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If VECT = 'Q', LDA >= max(1,nq);
!          if VECT = 'P', LDA >= max(1,min(nq,K)).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (min(nq,K))
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i) or G(i) which determines Q or P, as returned
!          by DGEBRD in the array argument TAUQ or TAUP.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
!          or P*C or P**T*C or C*P or C*P**T.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
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
!     .. Local Scalars ..
      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORMLQ, DORMQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      APPLYQ = LSAME( VECT, 'Q' )
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK == -1 )
!
!     NQ is the order of Q or P and NW is the minimum dimension of WORK
!
      if ( LEFT ) then
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      end if
      if ( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) then
         INFO = -1
      else if ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) then
         INFO = -2
      else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) then
         INFO = -3
      else if ( M < 0 ) then
         INFO = -4
      else if ( N < 0 ) then
         INFO = -5
      else if ( K < 0 ) then
         INFO = -6
      else if ( ( APPLYQ .AND. LDA < MAX( 1, NQ ) ) .OR. &
               ( .NOT.APPLYQ .AND. LDA < MAX( 1, MIN( NQ, K ) ) ) ) &
                THEN
         INFO = -8
      else if ( LDC < MAX( 1, M ) ) then
         INFO = -11
      else if ( LWORK < MAX( 1, NW ) .AND. .NOT.LQUERY ) then
         INFO = -13
      end if
!
      if ( INFO == 0 ) then
         if ( APPLYQ ) then
            if ( LEFT ) then
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M-1, N, M-1, &
                    -1 )
            ELSE
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N-1, N-1, &
                    -1 )
            end if
         ELSE
            if ( LEFT ) then
               NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M-1, N, M-1, &
                    -1 )
            ELSE
               NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N-1, N-1, &
                    -1 )
            end if
         end if
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DORMBR', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      WORK( 1 ) = 1
      if ( M == 0 .OR. N.EQ.0 ) &
         RETURN
!
      if ( APPLYQ ) then
!
!        Apply Q
!
         if ( NQ.GE.K ) then
!
!           Q was determined by a call to DGEBRD with nq >= k
!
            CALL DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, IINFO )
         else if ( NQ.GT.1 ) then
!
!           Q was determined by a call to DGEBRD with nq < k
!
            if ( LEFT ) then
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            end if
            CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, &
                         C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         end if
      ELSE
!
!        Apply P
!
         if ( NOTRAN ) then
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         end if
         if ( NQ.GT.K ) then
!
!           P was determined by a call to DGEBRD with nq > k
!
            CALL DORMLQ( SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, IINFO )
         else if ( NQ.GT.1 ) then
!
!           P was determined by a call to DGEBRD with nq <= k
!
            if ( LEFT ) then
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            end if
            CALL DORMLQ( SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA, &
                         TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         end if
      end if
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMBR
!
      END
