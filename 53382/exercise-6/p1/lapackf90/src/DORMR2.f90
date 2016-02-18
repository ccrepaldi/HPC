      SUBROUTINE DORMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMR2 overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGERQF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
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
!  A       (input) DOUBLE PRECISION array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGERQF in the last k rows of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGERQF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, MI, NI, NQ
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
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
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      if ( LEFT ) then
         NQ = M
      ELSE
         NQ = N
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
      else if ( LDA < MAX( 1, K ) ) then
         INFO = -7
      else if ( LDC < MAX( 1, M ) ) then
         INFO = -10
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DORMR2', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      if ( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) &
           THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      end if
!
      if ( LEFT ) then
         NI = N
      ELSE
         MI = M
      end if
!
      DO 10 I = I1, I2, I3
         if ( LEFT ) then
!
!           H(i) is applied to C(1:m-k+i,1:n)
!
            MI = M - K + I
         ELSE
!
!           H(i) is applied to C(1:m,1:n-k+i)
!
            NI = N - K + I
         end if
!
!        Apply H(i)
!
         AII = A( I, NQ-K+I )
         A( I, NQ-K+I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, 1 ), LDA, TAU( I ), C, LDC, &
                     WORK )
         A( I, NQ-K+I ) = AII
   10 CONTINUE
      RETURN
!
!     End of DORMR2
!
      END
