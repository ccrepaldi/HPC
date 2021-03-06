      SUBROUTINE DORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGRQ generates an M-by-N real matrix Q with orthonormal rows,
!  which is defined as the last M rows of a product of K elementary
!  reflectors of order N
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGERQF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. N >= M.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. M >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (m-k+i)-th row must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGERQF in the last k rows of its array argument
!          A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGERQF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,M).
!          For optimum performance LWORK >= M*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, II, IINFO, IWS, J, KK, L, LDWORK, &
                         LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORGR2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'DORGRQ', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, M )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) then
         INFO = -1
      else if ( N < M ) then
         INFO = -2
      else if ( K < 0 .OR. K.GT.M ) then
         INFO = -3
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -5
      else if ( LWORK < MAX( 1, M ) .AND. .NOT.LQUERY ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DORGRQ', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M.LE.0 ) then
         WORK( 1 ) = 1
         RETURN
      end if
!
      NBMIN = 2
      NX = 0
      IWS = M
      if ( NB.GT.1 .AND. NB < K ) then
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'DORGRQ', ' ', M, N, K, -1 ) )
         if ( NX < K ) then
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = M
            IWS = LDWORK*NB
            if ( LWORK < IWS ) then
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGRQ', ' ', M, N, K, -1 ) )
            end if
         end if
      end if
!
      if ( NB.GE.NBMIN .AND. NB < K .AND. NX.LT.K ) then
!
!        Use blocked code after the first block.
!        The last kk rows are handled by the block method.
!
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
!
!        Set A(1:m-kk,n-kk+1:n) to zero.
!
         DO 20 J = N - KK + 1, N
            DO 10 I = 1, M - KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      end if
!
!     Use unblocked code for the first or only block.
!
      CALL DORGR2( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
!
      if ( KK.GT.0 ) then
!
!        Use blocked code
!
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            II = M - K + I
            if ( II.GT.1 ) then
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL DLARFT( 'Backward', 'Rowwise', N-K+I+IB-1, IB, &
                            A( II, 1 ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H' to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
!
               CALL DLARFB( 'Right', 'Transpose', 'Backward', 'Rowwise', &
                            II-1, N-K+I+IB-1, IB, A( II, 1 ), LDA, WORK, &
                            LDWORK, A, LDA, WORK( IB+1 ), LDWORK )
            end if
!
!           Apply H' to columns 1:n-k+i+ib-1 of current block
!
            CALL DORGR2( IB, N-K+I+IB-1, IB, A( II, 1 ), LDA, TAU( I ), &
                         WORK, IINFO )
!
!           Set columns n-k+i+ib:n of current block to zero
!
            DO 40 L = N - K + I + IB, N
               DO 30 J = II, II + IB - 1
                  A( J, L ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      end if
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORGRQ
!
      END
