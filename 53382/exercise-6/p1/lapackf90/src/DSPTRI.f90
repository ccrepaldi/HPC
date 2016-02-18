      SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
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
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AP( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSPTRI computes the inverse of a real symmetric indefinite matrix
!  A in packed storage using the factorization A = U*D*U**T or
!  A = L*D*L**T computed by DSPTRF.
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
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the block diagonal matrix D and the multipliers
!          used to obtain the factor U or L as computed by DSPTRF,
!          stored as a packed triangular matrix.
!
!          On exit, if INFO = 0, the (symmetric) inverse of the original
!          matrix, stored as a packed triangular matrix. The j-th column
!          of inv(A) is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
!          if UPLO = 'L',
!             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSPTRF.
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
      INTEGER            J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP
      DOUBLE PRECISION   AK, AKKP1, AKP1, D, T, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DSPMV, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
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
         CALL XERBLA( 'DSPTRI', -INFO )
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
         KP = N*( N+1 ) / 2
         DO 10 INFO = N, 1, -1
            if ( IPIV( INFO ).GT.0 .AND. AP( KP ) == ZERO ) &
               RETURN
            KP = KP - INFO
   10    CONTINUE
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         KP = 1
         DO 20 INFO = 1, N
            if ( IPIV( INFO ).GT.0 .AND. AP( KP ) == ZERO ) &
               RETURN
            KP = KP + N - INFO + 1
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
         KC = 1
   30    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K.GT.N ) &
            GO TO 50
!
         KCNEXT = KC + K
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
            AP( KC+K-1 ) = ONE / AP( KC+K-1 )
!
!           Compute column K of the inverse.
!
            if ( K.GT.1 ) then
               CALL DCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL DSPMV( UPLO, K-1, -ONE, AP, WORK, 1, ZERO, AP( KC ), &
                           1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) - &
                              DDOT( K-1, WORK, 1, AP( KC ), 1 )
            end if
            KSTEP = 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
            T = ABS( AP( KCNEXT+K-1 ) )
            AK = AP( KC+K-1 ) / T
            AKP1 = AP( KCNEXT+K ) / T
            AKKP1 = AP( KCNEXT+K-1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KC+K-1 ) = AKP1 / D
            AP( KCNEXT+K ) = AK / D
            AP( KCNEXT+K-1 ) = -AKKP1 / D
!
!           Compute columns K and K+1 of the inverse.
!
            if ( K.GT.1 ) then
               CALL DCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL DSPMV( UPLO, K-1, -ONE, AP, WORK, 1, ZERO, AP( KC ), &
                           1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) - &
                              DDOT( K-1, WORK, 1, AP( KC ), 1 )
               AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) - &
                                  DDOT( K-1, AP( KC ), 1, AP( KCNEXT ), &
                                  1 )
               CALL DCOPY( K-1, AP( KCNEXT ), 1, WORK, 1 )
               CALL DSPMV( UPLO, K-1, -ONE, AP, WORK, 1, ZERO, &
                           AP( KCNEXT ), 1 )
               AP( KCNEXT+K ) = AP( KCNEXT+K ) - &
                                DDOT( K-1, WORK, 1, AP( KCNEXT ), 1 )
            end if
            KSTEP = 2
            KCNEXT = KCNEXT + K + 1
         end if
!
         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) then
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
            KPC = ( KP-1 )*KP / 2 + 1
            CALL DSWAP( KP-1, AP( KC ), 1, AP( KPC ), 1 )
            KX = KPC + KP - 1
            DO 40 J = KP + 1, K - 1
               KX = KX + J - 1
               TEMP = AP( KC+J-1 )
               AP( KC+J-1 ) = AP( KX )
               AP( KX ) = TEMP
   40       CONTINUE
            TEMP = AP( KC+K-1 )
            AP( KC+K-1 ) = AP( KPC+KP-1 )
            AP( KPC+KP-1 ) = TEMP
            if ( KSTEP == 2 ) then
               TEMP = AP( KC+K+K-1 )
               AP( KC+K+K-1 ) = AP( KC+K+KP-1 )
               AP( KC+K+KP-1 ) = TEMP
            end if
         end if
!
         K = K + KSTEP
         KC = KCNEXT
         GO TO 30
   50    CONTINUE
!
      ELSE
!
!        Compute inv(A) from the factorization A = L*D*L'.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         NPP = N*( N+1 ) / 2
         K = N
         KC = NPP
   60    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 ) &
            GO TO 80
!
         KCNEXT = KC - ( N-K+2 )
         if ( IPIV( K ).GT.0 ) then
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
            AP( KC ) = ONE / AP( KC )
!
!           Compute column K of the inverse.
!
            if ( K < N ) then
               CALL DCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL DSPMV( UPLO, N-K, -ONE, AP( KC+N-K+1 ), WORK, 1, &
                           ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1, AP( KC+1 ), 1 )
            end if
            KSTEP = 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
            T = ABS( AP( KCNEXT+1 ) )
            AK = AP( KCNEXT ) / T
            AKP1 = AP( KC ) / T
            AKKP1 = AP( KCNEXT+1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KCNEXT ) = AKP1 / D
            AP( KC ) = AK / D
            AP( KCNEXT+1 ) = -AKKP1 / D
!
!           Compute columns K-1 and K of the inverse.
!
            if ( K < N ) then
               CALL DCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL DSPMV( UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1, &
                           ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1, AP( KC+1 ), 1 )
               AP( KCNEXT+1 ) = AP( KCNEXT+1 ) - &
                                DDOT( N-K, AP( KC+1 ), 1, &
                                AP( KCNEXT+2 ), 1 )
               CALL DCOPY( N-K, AP( KCNEXT+2 ), 1, WORK, 1 )
               CALL DSPMV( UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1, &
                           ZERO, AP( KCNEXT+2 ), 1 )
               AP( KCNEXT ) = AP( KCNEXT ) - &
                              DDOT( N-K, WORK, 1, AP( KCNEXT+2 ), 1 )
            end if
            KSTEP = 2
            KCNEXT = KCNEXT - ( N-K+3 )
         end if
!
         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) then
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
            KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1
            if ( KP < N ) &
               CALL DSWAP( N-KP, AP( KC+KP-K+1 ), 1, AP( KPC+1 ), 1 )
            KX = KC + KP - K
            DO 70 J = K + 1, KP - 1
               KX = KX + N - J + 1
               TEMP = AP( KC+J-K )
               AP( KC+J-K ) = AP( KX )
               AP( KX ) = TEMP
   70       CONTINUE
            TEMP = AP( KC )
            AP( KC ) = AP( KPC )
            AP( KPC ) = TEMP
            if ( KSTEP == 2 ) then
               TEMP = AP( KC-N+K-1 )
               AP( KC-N+K-1 ) = AP( KC-N+KP-1 )
               AP( KC-N+KP-1 ) = TEMP
            end if
         end if
!
         K = K - KSTEP
         KC = KCNEXT
         GO TO 60
   80    CONTINUE
      end if
!
      RETURN
!
!     End of DSPTRI
!
      END
