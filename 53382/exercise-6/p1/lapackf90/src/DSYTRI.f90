      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRI computes the inverse of a real symmetric indefinite matrix
!  A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!  DSYTRF.
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
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the block diagonal matrix D and the multipliers
!          used to obtain the factor U or L as computed by DSYTRF.
!
!          On exit, if INFO = 0, the (symmetric) inverse of the original
!          matrix.  If UPLO = 'U', the upper triangular part of the
!          inverse is formed and the part of A below the diagonal is not
!          referenced; if UPLO = 'L' the lower triangular part of the
!          inverse is formed and the part of A above the diagonal is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSYTRF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!               inverse could not be computed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KP, KSTEP
      DOUBLE PRECISION   AK, AKKP1, AKP1, D, T, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DSWAP, DSYMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
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
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSYTRI', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Check that the diagonal matrix D is nonsingular.
!
      if ( UPPER ) then
!
!        Upper triangular storage: examine D from bottom to top
!
         DO 10 INFO = N, 1, -1
            if ( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ) == ZERO ) &
               RETURN
   10    CONTINUE
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO 20 INFO = 1, N
            if ( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ) == ZERO ) &
               RETURN
   20    CONTINUE
      end if
      INFO = 0
!
      if ( UPPER ) then
!
!        Compute inv(A) from the factorization A = U*D*U'.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   30    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K.GT.N ) &
            GO TO 40
!
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
            A( K, K ) = ONE / A( K, K )
!
!           Compute column K of the inverse.
!
            if ( K.GT.1 ) then
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, &
                           A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ), &
                           1 )
            end if
            KSTEP = 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
            T = ABS( A( K, K+1 ) )
            AK = A( K, K ) / T
            AKP1 = A( K+1, K+1 ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D
!
!           Compute columns K and K+1 of the inverse.
!
            if ( K.GT.1 ) then
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, &
                           A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ), &
                           1 )
               A( K, K+1 ) = A( K, K+1 ) - &
                             DDOT( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL DCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, &
                           A( 1, K+1 ), 1 )
               A( K+1, K+1 ) = A( K+1, K+1 ) - &
                               DDOT( K-1, WORK, 1, A( 1, K+1 ), 1 )
            end if
            KSTEP = 2
         end if
!
         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) then
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
            CALL DSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL DSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            if ( KSTEP == 2 ) then
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            end if
         end if
!
         K = K + KSTEP
         GO TO 30
   40    CONTINUE
!
      ELSE
!
!        Compute inv(A) from the factorization A = L*D*L'.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   50    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 ) &
            GO TO 60
!
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
            A( K, K ) = ONE / A( K, K )
!
!           Compute column K of the inverse.
!
            if ( K < N ) then
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, &
                           ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ), &
                           1 )
            end if
            KSTEP = 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
            T = ABS( A( K, K-1 ) )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
!
!           Compute columns K-1 and K of the inverse.
!
            if ( K < N ) then
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, &
                           ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ), &
                           1 )
               A( K, K-1 ) = A( K, K-1 ) - &
                             DDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ), &
                             1 )
               CALL DCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, &
                           ZERO, A( K+1, K-1 ), 1 )
               A( K-1, K-1 ) = A( K-1, K-1 ) - &
                               DDOT( N-K, WORK, 1, A( K+1, K-1 ), 1 )
            end if
            KSTEP = 2
         end if
!
         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) then
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
            if ( KP < N ) &
               CALL DSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL DSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            if ( KSTEP == 2 ) then
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            end if
         end if
!
         K = K - KSTEP
         GO TO 50
   60    CONTINUE
      end if
!
      RETURN
!
!     End of DSYTRI
!
      END
