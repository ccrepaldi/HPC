      SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AP( * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DSPTRS solves a system of linear equations A*X = B with a real
!  symmetric matrix A stored in packed format using the factorization
!  A = U*D*U**T or A = L*D*L**T computed by DSPTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular, form is A = U*D*U**T;
!          = 'L':  Lower triangular, form is A = L*D*L**T.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The block diagonal matrix D and the multipliers used to
!          obtain the factor U or L as computed by DSPTRF, stored as a
!          packed triangular matrix.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSPTRF.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KC, KP
      DOUBLE PRECISION   AK, AKM1, AKM1K, BK, BKM1, DENOM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -7
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSPTRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 .OR. NRHS.EQ.0 ) &
         RETURN
!
      if ( UPPER ) then
!
!        Solve A*X = B, where A = U*D*U'.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
         KC = N*( N+1 ) / 2 + 1
   10    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 ) &
            GO TO 30
!
         KC = KC - K
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
            CALL DGER( K-1, NRHS, -ONE, AP( KC ), 1, B( K, 1 ), LDB, &
                       B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            CALL DSCAL( NRHS, ONE / AP( KC+K-1 ), B( K, 1 ), LDB )
            K = K - 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K-1 ) &
               CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
            CALL DGER( K-2, NRHS, -ONE, AP( KC ), 1, B( K, 1 ), LDB, &
                       B( 1, 1 ), LDB )
            CALL DGER( K-2, NRHS, -ONE, AP( KC-( K-1 ) ), 1, &
                       B( K-1, 1 ), LDB, B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = AP( KC+K-2 )
            AKM1 = AP( KC-1 ) / AKM1K
            AK = AP( KC+K-1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / AKM1K
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            KC = KC - K + 1
            K = K - 2
         end if
!
         GO TO 10
   30    CONTINUE
!
!        Next solve U'*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
         KC = 1
   40    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K.GT.N ) &
            GO TO 50
!
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U'(K)), where U(K) is the transformation
!           stored in column K of A.
!
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, AP( KC ), &
                        1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC + K
            K = K + 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U'(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, AP( KC ), &
                        1, ONE, B( K, 1 ), LDB )
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, &
                        AP( KC+K ), 1, ONE, B( K+1, 1 ), LDB )
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K ) &
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC + 2*K + 1
            K = K + 2
         end if
!
         GO TO 40
   50    CONTINUE
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L'.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
         KC = 1
   60    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K.GT.N ) &
            GO TO 80
!
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
            if ( K < N ) &
               CALL DGER( N-K, NRHS, -ONE, AP( KC+1 ), 1, B( K, 1 ), &
                          LDB, B( K+1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            CALL DSCAL( NRHS, ONE / AP( KC ), B( K, 1 ), LDB )
            KC = KC + N - K + 1
            K = K + 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K+1 ) &
               CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
            if ( K < N-1 ) then
               CALL DGER( N-K-1, NRHS, -ONE, AP( KC+2 ), 1, B( K, 1 ), &
                          LDB, B( K+2, 1 ), LDB )
               CALL DGER( N-K-1, NRHS, -ONE, AP( KC+N-K+2 ), 1, &
                          B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            end if
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = AP( KC+1 )
            AKM1 = AP( KC ) / AKM1K
            AK = AP( KC+N-K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / AKM1K
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            KC = KC + 2*( N-K ) + 1
            K = K + 2
         end if
!
         GO TO 60
   80    CONTINUE
!
!        Next solve L'*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
         KC = N*( N+1 ) / 2 + 1
   90    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 ) &
            GO TO 100
!
         KC = KC - ( N-K+1 )
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Multiply by inv(L'(K)), where L(K) is the transformation
!           stored in column K of A.
!
            if ( K < N ) &
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), &
                           LDB, AP( KC+1 ), 1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L'(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
            if ( K < N ) then
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), &
                           LDB, AP( KC+1 ), 1, ONE, B( K, 1 ), LDB )
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), &
                           LDB, AP( KC-( N-K ) ), 1, ONE, B( K-1, 1 ), &
                           LDB )
            end if
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K ) &
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC - ( N-K+2 )
            K = K - 2
         end if
!
         GO TO 90
  100    CONTINUE
      end if
!
      RETURN
!
!     End of DSPTRS
!
      END
