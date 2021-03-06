      SUBROUTINE DPPTRI( UPLO, N, AP, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * )
!     ..
!
!  Purpose
!  =======
!
!  DPPTRI computes the inverse of a real symmetric positive definite
!  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
!  computed by DPPTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangular factor is stored in AP;
!          = 'L':  Lower triangular factor is stored in AP.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the triangular factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T, packed columnwise as
!          a linear array.  The j-th column of U or L is stored in the
!          array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
!
!          On exit, the upper or lower triangle of the (symmetric)
!          inverse of A, overwriting the input factor U or L.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the (i,i) element of the factor U or L is
!                zero, and the inverse could not be computed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ, JJN
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSPR, DTPMV, DTPTRI, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPPTRI', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Invert the triangular Cholesky factor U or L.
!
      CALL DTPTRI( UPLO, 'Non-unit', N, AP, INFO )
      if ( INFO.GT.0 ) &
         RETURN
!
      if ( UPPER ) then
!
!        Compute the product inv(U) * inv(U)'.
!
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
            if ( J.GT.1 ) &
               CALL DSPR( 'Upper', J-1, ONE, AP( JC ), 1, AP )
            AJJ = AP( JJ )
            CALL DSCAL( J, AJJ, AP( JC ), 1 )
   10    CONTINUE
!
      ELSE
!
!        Compute the product inv(L)' * inv(L).
!
         JJ = 1
         DO 20 J = 1, N
            JJN = JJ + N - J + 1
            AP( JJ ) = DDOT( N-J+1, AP( JJ ), 1, AP( JJ ), 1 )
            if ( J < N ) &
               CALL DTPMV( 'Lower', 'Transpose', 'Non-unit', N-J, &
                           AP( JJN ), AP( JJ+1 ), 1 )
            JJ = JJN
   20    CONTINUE
      end if
!
      RETURN
!
!     End of DPPTRI
!
      END
