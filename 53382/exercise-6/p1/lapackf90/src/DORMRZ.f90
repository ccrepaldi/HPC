      SUBROUTINE DORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMRZ overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
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
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  L       (input) INTEGER
!          The number of columns of the matrix A containing
!          the meaningful part of the Householder reflectors.
!          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DTZRZF in the last k rows of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DTZRZF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
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
!  Further Details
!  ===============
!
!  Based on contributions by
!    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JA, JC, &
                         LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARZB, DLARZT, DORMR3, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
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
      else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) then
         INFO = -2
      else if ( M < 0 ) then
         INFO = -3
      else if ( N < 0 ) then
         INFO = -4
      else if ( K < 0 .OR. K.GT.NQ ) then
         INFO = -5
      else if ( L < 0 .OR. ( LEFT .AND. ( L.GT.M ) ) .OR. &
               ( .NOT.LEFT .AND. ( L.GT.N ) ) ) then
         INFO = -6
      else if ( LDA < MAX( 1, K ) ) then
         INFO = -8
      else if ( LDC < MAX( 1, M ) ) then
         INFO = -11
      else if ( LWORK < MAX( 1, NW ) .AND. .NOT.LQUERY ) then
         INFO = -13
      end if
!
      if ( INFO == 0 ) then
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         NB = MIN( NBMAX, ILAENV( 1, 'DORMRQ', SIDE // TRANS, M, N, K, &
              -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DORMRZ', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 .OR. K.EQ.0 ) then
         WORK( 1 ) = 1
         RETURN
      end if
!
      NBMIN = 2
      LDWORK = NW
      if ( NB.GT.1 .AND. NB < K ) then
         IWS = NW*NB
         if ( LWORK < IWS ) then
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMRQ', SIDE // TRANS, M, N, K, &
                    -1 ) )
         end if
      ELSE
         IWS = NW
      end if
!
      if ( NB < NBMIN .OR. NB.GE.K ) then
!
!        Use unblocked code
!
         CALL DORMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, &
                      WORK, IINFO )
      ELSE
!
!        Use blocked code
!
         if ( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. NOTRAN ) ) then
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         end if
!
         if ( LEFT ) then
            NI = N
            JC = 1
            JA = M - L + 1
         ELSE
            MI = M
            IC = 1
            JA = N - L + 1
         end if
!
         if ( NOTRAN ) then
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         end if
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL DLARZT( 'Backward', 'Rowwise', L, IB, A( I, JA ), LDA, &
                         TAU( I ), T, LDT )
!
            if ( LEFT ) then
!
!              H or H' is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H' is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            end if
!
!           Apply H or H'
!
            CALL DLARZB( SIDE, TRANST, 'Backward', 'Rowwise', MI, NI, &
                         IB, L, A( I, JA ), LDA, T, LDT, C( IC, JC ), &
                         LDC, WORK, LDWORK )
   10    CONTINUE
!
      end if
!
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of DORMRZ
!
      END
