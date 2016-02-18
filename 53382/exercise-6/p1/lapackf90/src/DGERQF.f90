      SUBROUTINE DGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGERQF computes an RQ factorization of a real M-by-N matrix A:
!  A = R * Q.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!          if m <= n, the upper triangle of the subarray
!          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
!          if m >= n, the elements on and above the (m-n)-th subdiagonal
!          contain the M-by-N upper trapezoidal matrix R;
!          the remaining elements, with the array TAU, represent the
!          orthogonal matrix Q as a product of min(m,n) elementary
!          reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,M).
!          For optimum performance LWORK >= M*NB, where NB is
!          the optimal blocksize.
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
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
!  A(m-k+i,1:n-k+i-1), and tau in TAU(i).
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, KI, KK, LDWORK, LWKOPT, &
                         MU, NB, NBMIN, NU, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGERQ2, DLARFB, DLARFT, XERBLA
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
      NB = ILAENV( 1, 'DGERQF', ' ', M, N, -1, -1 )
      LWKOPT = M*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -4
      else if ( LWORK < MAX( 1, M ) .AND. .NOT.LQUERY ) then
         INFO = -7
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGERQF', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      K = MIN( M, N )
      if ( K == 0 ) then
         WORK( 1 ) = 1
         RETURN
      end if
!
      NBMIN = 2
      NX = 1
      IWS = M
      if ( NB.GT.1 .AND. NB < K ) then
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'DGERQF', ' ', M, N, -1, -1 ) )
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
               NBMIN = MAX( 2, ILAENV( 2, 'DGERQF', ' ', M, N, -1, &
                       -1 ) )
            end if
         end if
      end if
!
      if ( NB.GE.NBMIN .AND. NB < K .AND. NX.LT.K ) then
!
!        Use blocked code initially.
!        The last kk rows are handled by the block method.
!
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!
         DO 10 I = K - KK + KI + 1, K - KK + 1, -NB
            IB = MIN( K-I+1, NB )
!
!           Compute the RQ factorization of the current block
!           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)
!
            CALL DGERQ2( IB, N-K+I+IB-1, A( M-K+I, 1 ), LDA, TAU( I ), &
                         WORK, IINFO )
            if ( M-K+I.GT.1 ) then
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL DLARFT( 'Backward', 'Rowwise', N-K+I+IB-1, IB, &
                            A( M-K+I, 1 ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
!
               CALL DLARFB( 'Right', 'No transpose', 'Backward', &
                            'Rowwise', M-K+I-1, N-K+I+IB-1, IB, &
                            A( M-K+I, 1 ), LDA, WORK, LDWORK, A, LDA, &
                            WORK( IB+1 ), LDWORK )
            end if
   10    CONTINUE
         MU = M - K + I + NB - 1
         NU = N - K + I + NB - 1
      ELSE
         MU = M
         NU = N
      end if
!
!     Use unblocked code to factor the last or only block
!
      if ( MU.GT.0 .AND. NU.GT.0 ) &
         CALL DGERQ2( MU, NU, A, LDA, TAU, WORK, IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DGERQF
!
      END