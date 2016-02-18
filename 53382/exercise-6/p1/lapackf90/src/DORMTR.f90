      SUBROUTINE DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMTR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by DSYTRD:
!
!  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from DSYTRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from DSYTRD.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                               (LDA,M) if SIDE = 'L'
!                               (LDA,N) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by DSYTRD.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!
!  TAU     (input) DOUBLE PRECISION array, dimension
!                               (M-1) if SIDE = 'L'
!                               (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSYTRD.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
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
      LOGICAL            LEFT, LQUERY, UPPER
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORMQL, DORMQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      if ( LEFT ) then
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      end if
      if ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) then
         INFO = -1
      else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -2
      else if ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) ) &
                THEN
         INFO = -3
      else if ( M < 0 ) then
         INFO = -4
      else if ( N < 0 ) then
         INFO = -5
      else if ( LDA < MAX( 1, NQ ) ) then
         INFO = -7
      else if ( LDC < MAX( 1, M ) ) then
         INFO = -10
      else if ( LWORK < MAX( 1, NW ) .AND. .NOT.LQUERY ) then
         INFO = -12
      end if
!
      if ( INFO == 0 ) then
         if ( UPPER ) then
            if ( LEFT ) then
               NB = ILAENV( 1, 'DORMQL', SIDE // TRANS, M-1, N, M-1, &
                    -1 )
            ELSE
               NB = ILAENV( 1, 'DORMQL', SIDE // TRANS, M, N-1, N-1, &
                    -1 )
            end if
         ELSE
            if ( LEFT ) then
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M-1, N, M-1, &
                    -1 )
            ELSE
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N-1, N-1, &
                    -1 )
            end if
         end if
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DORMTR', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) then
         WORK( 1 ) = 1
         RETURN
      end if
!
      if ( LEFT ) then
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      end if
!
      if ( UPPER ) then
!
!        Q was determined by a call to DSYTRD with UPLO = 'U'
!
         CALL DORMQL( SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C, &
                      LDC, WORK, LWORK, IINFO )
      ELSE
!
!        Q was determined by a call to DSYTRD with UPLO = 'L'
!
         if ( LEFT ) then
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         end if
         CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, &
                      C( I1, I2 ), LDC, WORK, LWORK, IINFO )
      end if
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMTR
!
      END
