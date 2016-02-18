      SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTF2 computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method:
!
!     A = U*D*U'  or  A = L*D*L'
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, U' is the transpose of U, and D is symmetric and
!  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, and division by zero will occur if it
!               is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
!         Company
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1, &
                         ROWMAX, T, WK, WKM1, WKP1
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, DSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
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
         CALL XERBLA( 'DSYTF2', -INFO )
         RETURN
      end if
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
      if ( UPPER ) then
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop
!
         if ( K < 1 ) &
            GO TO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         if ( K.GT.1 ) then
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         end if
!
         if ( MAX( ABSAKK, COLMAX ) == ZERO ) then
!
!           Column K is zero: set INFO and continue
!
            if ( INFO == 0 ) &
               INFO = K
            KP = K
         ELSE
            if ( ABSAKK.GE.ALPHA*COLMAX ) then
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               if ( IMAX.GT.1 ) then
                  JMAX = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               end if
!
               if ( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               else if ( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!
            KK = K - KSTEP + 1
            if ( KP.NE.KK ) then
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
               CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), &
                           LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) then
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               end if
            end if
!
!           Update the leading submatrix
!
            if ( KSTEP == 1 ) then
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!
!              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
!
               R1 = ONE / A( K, K )
               CALL DSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
!
!              Store U(k) in column k
!
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
!
               if ( K.GT.2 ) then
!
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )
                  D12 = T / D12
!
                  DO 30 J = K - 2, 1, -1
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - A( I, K )*WK - &
                                    A( I, K-1 )*WKM1
   20                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
   30             CONTINUE
!
               end if
!
            end if
         end if
!
!        Store details of the interchanges in IPIV
!
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         end if
!
!        Decrease K and return to the start of the main loop
!
         K = K - KSTEP
         GO TO 10
!
      ELSE
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
         K = 1
   40    CONTINUE
!
!        If K > N, exit from loop
!
         if ( K.GT.N ) &
            GO TO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         if ( K < N ) then
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         end if
!
         if ( MAX( ABSAKK, COLMAX ) == ZERO ) then
!
!           Column K is zero: set INFO and continue
!
            if ( INFO == 0 ) &
               INFO = K
            KP = K
         ELSE
            if ( ABSAKK.GE.ALPHA*COLMAX ) then
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               if ( IMAX < N ) then
                  JMAX = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               end if
!
               if ( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               else if ( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!
            KK = K + KSTEP - 1
            if ( KP.NE.KK ) then
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
               if ( KP < N ) &
                  CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), &
                           LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) then
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               end if
            end if
!
!           Update the trailing submatrix
!
            if ( KSTEP == 1 ) then
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
               if ( K < N ) then
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!
!                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
!
                  D11 = ONE / A( K, K )
                  CALL DSYR( UPLO, N-K, -D11, A( K+1, K ), 1, &
                             A( K+1, K+1 ), LDA )
!
!                 Store L(k) in column K
!
                  CALL DSCAL( N-K, D11, A( K+1, K ), 1 )
               end if
            ELSE
!
!              2-by-2 pivot block D(k)
!
               if ( K < N-1 ) then
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
!
!                 where L(k) and L(k+1) are the k-th and (k+1)-th
!                 columns of L
!
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
!
                  DO 60 J = K + 2, N
!
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
!
                     DO 50 I = J, N
                        A( I, J ) = A( I, J ) - A( I, K )*WK - &
                                    A( I, K+1 )*WKP1
   50                CONTINUE
!
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
!
   60             CONTINUE
               end if
            end if
         end if
!
!        Store details of the interchanges in IPIV
!
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         end if
!
!        Increase K and return to the start of the main loop
!
         K = K + KSTEP
         GO TO 40
!
      end if
!
   70 CONTINUE
!
      RETURN
!
!     End of DSYTF2
!
      END
