      SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGL2 generates an m by n real matrix Q with orthonormal rows,
!  which is defined as the first m rows of a product of k elementary
!  reflectors of order n
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGELQF.
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
!          On entry, the i-th row must contain the vector which defines
!          the elementary reflector H(i), for i = 1,2,...,k, as returned
!          by DGELQF in the first k rows of its array argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      if ( M < 0 ) then
         INFO = -1
      else if ( N < M ) then
         INFO = -2
      else if ( K < 0 .OR. K.GT.M ) then
         INFO = -3
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -5
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DORGL2', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M.LE.0 ) &
         RETURN
!
      if ( K < M ) then
!
!        Initialise rows k+1:m to rows of the unit matrix
!
         DO 20 J = 1, N
            DO 10 L = K + 1, M
               A( L, J ) = ZERO
   10       CONTINUE
            if ( J.GT.K .AND. J.LE.M ) &
               A( J, J ) = ONE
   20    CONTINUE
      end if
!
      DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the right
!
         if ( I < N ) then
            if ( I < M ) then
               A( I, I ) = ONE
               CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, &
                           TAU( I ), A( I+1, I ), LDA, WORK )
            end if
            CALL DSCAL( N-I, -TAU( I ), A( I, I+1 ), LDA )
         end if
         A( I, I ) = ONE - TAU( I )
!
!        Set A(i,1:i-1) to zero
!
         DO 30 L = 1, I - 1
            A( I, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of DORGL2
!
      END
