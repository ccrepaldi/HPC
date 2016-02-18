      SUBROUTINE DTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
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
!  DTZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
!  to upper triangular form by means of orthogonal transformations.
!
!  The upper trapezoidal matrix A is factored as
!
!     A = ( R  0 ) * Z,
!
!  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
!  triangular matrix.
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
!          On entry, the leading M-by-N upper trapezoidal part of the
!          array A must contain the matrix to be factorized.
!          On exit, the leading M-by-M upper triangular part of A
!          contains the upper triangular matrix R, and elements M+1 to
!          N of the first M rows of A, with the array TAU, represent the
!          orthogonal matrix Z as a product of M elementary reflectors.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (M)
!          The scalar factors of the elementary reflectors.
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
!  Based on contributions by
!    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!  The factorization is obtained by Householder's method.  The kth
!  transformation matrix, Z( k ), which is used to introduce zeros into
!  the ( m - k + 1 )th row of A, is given in the form
!
!     Z( k ) = ( I     0   ),
!              ( 0  T( k ) )
!
!  where
!
!     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
!                                                 (   0    )
!                                                 ( z( k ) )
!
!  tau is a scalar and z( k ) is an ( n - m ) element vector.
!  tau and z( k ) are chosen to annihilate the elements of the kth row
!  of X.
!
!  The scalar tau is returned in the kth element of TAU and the vector
!  u( k ) in the kth row of A, such that the elements of z( k ) are
!  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in
!  the upper triangular part of A.
!
!  Z is given by
!
!     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IWS, KI, KK, LDWORK, LWKOPT, M1, MU, NB, &
                         NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARZB, DLARZT, DLATRZ, XERBLA
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
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) then
         INFO = -1
      else if ( N < M ) then
         INFO = -2
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -4
      else if ( LWORK < MAX( 1, M ) .AND. .NOT.LQUERY ) then
         INFO = -7
      end if
!
      if ( INFO == 0 ) then
!
!        Determine the block size.
!
         NB = ILAENV( 1, 'DGERQF', ' ', M, N, -1, -1 )
         LWKOPT = M*NB
         WORK( 1 ) = LWKOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DTZRZF', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 ) then
         WORK( 1 ) = 1
         RETURN
      else if ( M == N ) then
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
         WORK( 1 ) = 1
         RETURN
      end if
!
      NBMIN = 2
      NX = 1
      IWS = M
      if ( NB.GT.1 .AND. NB < M ) then
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'DGERQF', ' ', M, N, -1, -1 ) )
         if ( NX < M ) then
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
      if ( NB.GE.NBMIN .AND. NB < M .AND. NX.LT.M ) then
!
!        Use blocked code initially.
!        The last kk rows are handled by the block method.
!
         M1 = MIN( M+1, N )
         KI = ( ( M-NX-1 ) / NB )*NB
         KK = MIN( M, KI+NB )
!
         DO 20 I = M - KK + KI + 1, M - KK + 1, -NB
            IB = MIN( M-I+1, NB )
!
!           Compute the TZ factorization of the current block
!           A(i:i+ib-1,i:n)
!
            CALL DLATRZ( IB, N-I+1, N-M, A( I, I ), LDA, TAU( I ), &
                         WORK )
            if ( I.GT.1 ) then
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL DLARZT( 'Backward', 'Rowwise', N-M, IB, A( I, M1 ), &
                            LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:i-1,i:n) from the right
!
               CALL DLARZB( 'Right', 'No transpose', 'Backward', &
                            'Rowwise', I-1, N-I+1, IB, N-M, A( I, M1 ), &
                            LDA, WORK, LDWORK, A( 1, I ), LDA, &
                            WORK( IB+1 ), LDWORK )
            end if
   20    CONTINUE
         MU = I + NB - 1
      ELSE
         MU = M
      end if
!
!     Use unblocked code to factor the last or only block
!
      if ( MU.GT.0 ) &
         CALL DLATRZ( MU, N, N-M, A, LDA, TAU, WORK )
!
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of DTZRZF
!
      END
