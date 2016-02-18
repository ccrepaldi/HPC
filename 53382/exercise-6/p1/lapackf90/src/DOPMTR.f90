      SUBROUTINE DOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, &
                         INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DOPMTR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by DSPTRD using packed
!  storage:
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
!          = 'U': Upper triangular packed storage used in previous
!                 call to DSPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to DSPTRD.
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
!  AP      (input) DOUBLE PRECISION array, dimension
!                               (M*(M+1)/2) if SIDE = 'L'
!                               (N*(N+1)/2) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by DSPTRD.  AP is modified by the routine but
!          restored on exit.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L'
!                                     or (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSPTRD.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L'
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FORWRD, LEFT, NOTRAN, UPPER
      INTEGER            I, I1, I2, I3, IC, II, JC, MI, NI, NQ
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
      UPPER = LSAME( UPLO, 'U' )
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
      else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -2
      else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) then
         INFO = -3
      else if ( M < 0 ) then
         INFO = -4
      else if ( N < 0 ) then
         INFO = -5
      else if ( LDC < MAX( 1, M ) ) then
         INFO = -9
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DOPMTR', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 ) &
         RETURN
!
      if ( UPPER ) then
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
         FORWRD = ( LEFT .AND. NOTRAN ) .OR. &
                  ( .NOT.LEFT .AND. .NOT.NOTRAN )
!
         if ( FORWRD ) then
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
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
!              H(i) is applied to C(1:i,1:n)
!
               MI = I
            ELSE
!
!              H(i) is applied to C(1:m,1:i)
!
               NI = I
            end if
!
!           Apply H(i)
!
            AII = AP( II )
            AP( II ) = ONE
            CALL DLARF( SIDE, MI, NI, AP( II-I+1 ), 1, TAU( I ), C, LDC, &
                        WORK )
            AP( II ) = AII
!
            if ( FORWRD ) then
               II = II + I + 2
            ELSE
               II = II - I - 1
            end if
   10    CONTINUE
      ELSE
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
         FORWRD = ( LEFT .AND. .NOT.NOTRAN ) .OR. &
                  ( .NOT.LEFT .AND. NOTRAN )
!
         if ( FORWRD ) then
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
         end if
!
         if ( LEFT ) then
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         end if
!
         DO 20 I = I1, I2, I3
            AII = AP( II )
            AP( II ) = ONE
            if ( LEFT ) then
!
!              H(i) is applied to C(i+1:m,1:n)
!
               MI = M - I
               IC = I + 1
            ELSE
!
!              H(i) is applied to C(1:m,i+1:n)
!
               NI = N - I
               JC = I + 1
            end if
!
!           Apply H(i)
!
            CALL DLARF( SIDE, MI, NI, AP( II ), 1, TAU( I ), &
                        C( IC, JC ), LDC, WORK )
            AP( II ) = AII
!
            if ( FORWRD ) then
               II = II + NQ - I + 1
            ELSE
               II = II - NQ + I - 2
            end if
   20    CONTINUE
      end if
      RETURN
!
!     End of DOPMTR
!
      END
