      SUBROUTINE DGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )
!
!  -- LAPACK test routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  This routine is deprecated and has been replaced by routine DGEQP3.
!
!  DGEQPF computes a QR factorization with column pivoting of a
!  real M-by-N matrix A: A*P = Q*R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A. N >= 0
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the upper triangle of the array contains the
!          min(M,N)-by-N upper triangular matrix R; the elements
!          below the diagonal, together with the array TAU,
!          represent the orthogonal matrix Q as a product of
!          min(m,n) elementary reflectors.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  JPVT    (input/output) INTEGER array, dimension (N)
!          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
!          to the front of A*P (a leading column); if JPVT(i) = 0,
!          the i-th column of A is a free column.
!          On exit, if JPVT(i) = k, then the i-th column of A*P
!          was the k-th column of A.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
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
!     Q = H(1) H(2) . . . H(n)
!
!  Each H(i) has the form
!
!     H = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
!
!  The matrix P is represented in jpvt as follows: If
!     jpvt(j) = i
!  then the jth column of P is the ith canonical unit vector.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITEMP, J, MA, MN, PVT
      DOUBLE PRECISION   AII, TEMP, TEMP2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQR2, DLARF, DLARFG, DORM2R, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DNRM2
      EXTERNAL           IDAMAX, DNRM2
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEQPF', -INFO )
         RETURN
      end if
!
      MN = MIN( M, N )
!
!     Move initial columns up front
!
      ITEMP = 1
      DO 10 I = 1, N
         if ( JPVT( I ).NE.0 ) then
            if ( I.NE.ITEMP ) then
               CALL DSWAP( M, A( 1, I ), 1, A( 1, ITEMP ), 1 )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            ELSE
               JPVT( I ) = I
            end if
            ITEMP = ITEMP + 1
         ELSE
            JPVT( I ) = I
         end if
   10 CONTINUE
      ITEMP = ITEMP - 1
!
!     Compute the QR factorization and update remaining columns
!
      if ( ITEMP.GT.0 ) then
         MA = MIN( ITEMP, M )
         CALL DGEQR2( M, MA, A, LDA, TAU, WORK, INFO )
         if ( MA < N ) then
            CALL DORM2R( 'Left', 'Transpose', M, N-MA, MA, A, LDA, TAU, &
                         A( 1, MA+1 ), LDA, WORK, INFO )
         end if
      end if
!
      if ( ITEMP < MN ) then
!
!        Initialize partial column norms. The first n elements of
!        work store the exact column norms.
!
         DO 20 I = ITEMP + 1, N
            WORK( I ) = DNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            WORK( N+I ) = WORK( I )
   20    CONTINUE
!
!        Compute factorization
!
         DO 40 I = ITEMP + 1, MN
!
!           Determine ith pivot column and swap if necessary
!
            PVT = ( I-1 ) + IDAMAX( N-I+1, WORK( I ), 1 )
!
            if ( PVT.NE.I ) then
               CALL DSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               WORK( PVT ) = WORK( I )
               WORK( N+PVT ) = WORK( N+I )
            end if
!
!           Generate elementary reflector H(i)
!
            if ( I < M ) then
               CALL DLARFG( M-I+1, A( I, I ), A( I+1, I ), 1, TAU( I ) )
            ELSE
               CALL DLARFG( 1, A( M, M ), A( M, M ), 1, TAU( M ) )
            end if
!
            if ( I < N ) then
!
!              Apply H(i) to A(i:m,i+1:n) from the left
!
               AII = A( I, I )
               A( I, I ) = ONE
               CALL DLARF( 'LEFT', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                           A( I, I+1 ), LDA, WORK( 2*N+1 ) )
               A( I, I ) = AII
            end if
!
!           Update partial column norms
!
            DO 30 J = I + 1, N
               if ( WORK( J ).NE.ZERO ) then
                  TEMP = ONE - ( ABS( A( I, J ) ) / WORK( J ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = ONE + 0.05D0*TEMP* &
                          ( WORK( J ) / WORK( N+J ) )**2
                  if ( TEMP2 == ONE ) then
                     if ( M-I.GT.0 ) then
                        WORK( J ) = DNRM2( M-I, A( I+1, J ), 1 )
                        WORK( N+J ) = WORK( J )
                     ELSE
                        WORK( J ) = ZERO
                        WORK( N+J ) = ZERO
                     end if
                  ELSE
                     WORK( J ) = WORK( J )*SQRT( TEMP )
                  end if
               end if
   30       CONTINUE
!
   40    CONTINUE
      end if
      RETURN
!
!     End of DGEQPF
!
      END
