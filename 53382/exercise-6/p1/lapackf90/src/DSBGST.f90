      SUBROUTINE DSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, &
                         LDX, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO, VECT
      INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), WORK( * ), &
                         X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  DSBGST reduces a real symmetric-definite banded generalized
!  eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,
!  such that C has the same bandwidth as A.
!
!  B must have been previously factorized as S**T*S by DPBSTF, using a
!  split Cholesky factorization. A is overwritten by C = X**T*A*X, where
!  X = S**(-1)*Q and Q is an orthogonal matrix chosen to preserve the
!  bandwidth of A.
!
!  Arguments
!  =========
!
!  VECT    (input) CHARACTER*1
!          = 'N':  do not form the transformation matrix X;
!          = 'V':  form X.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  KA      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.
!
!  KB      (input) INTEGER
!          The number of superdiagonals of the matrix B if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first ka+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
!
!          On exit, the transformed matrix X**T*A*X, stored in the same
!          format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KA+1.
!
!  BB      (input) DOUBLE PRECISION array, dimension (LDBB,N)
!          The banded factor S from the split Cholesky factorization of
!          B, as returned by DPBSTF, stored in the first KB+1 rows of
!          the array.
!
!  LDBB    (input) INTEGER
!          The leading dimension of the array BB.  LDBB >= KB+1.
!
!  X       (output) DOUBLE PRECISION array, dimension (LDX,N)
!          If VECT = 'V', the n-by-n matrix X.
!          If VECT = 'N', the array X is not referenced.
!
!  LDX     (input) INTEGER
!          The leading dimension of the array X.
!          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPDATE, UPPER, WANTX
      INTEGER            I, I0, I1, I2, INCA, J, J1, J1T, J2, J2T, K, &
                         KA1, KB1, KBT, L, M, NR, NRT, NX
      DOUBLE PRECISION   BII, RA, RA1, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGER, DLAR2V, DLARGV, DLARTG, DLARTV, DLASET, &
                         DROT, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      WANTX = LSAME( VECT, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      KA1 = KA + 1
      KB1 = KB + 1
      INFO = 0
      if ( .NOT.WANTX .AND. .NOT.LSAME( VECT, 'N' ) ) then
         INFO = -1
      else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( KA < 0 ) then
         INFO = -4
      else if ( KB < 0 ) then
         INFO = -5
      else if ( LDAB < KA+1 ) then
         INFO = -7
      else if ( LDBB < KB+1 ) then
         INFO = -9
      else if ( LDX < 1 .OR. WANTX .AND. LDX.LT.MAX( 1, N ) ) then
         INFO = -11
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSBGST', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
      INCA = LDAB*KA1
!
!     Initialize X to the unit matrix, if needed
!
      if ( WANTX ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, X, LDX )
!
!     Set M to the splitting point m. It must be the same value as is
!     used in DPBSTF. The chosen value allows the arrays WORK and RWORK
!     to be of dimension (N).
!
      M = ( N+KB ) / 2
!
!     The routine works in two phases, corresponding to the two halves
!     of the split Cholesky factorization of B as S**T*S where
!
!     S = ( U    )
!         ( M  L )
!
!     with U upper triangular of order m, and L lower triangular of
!     order n-m. S has the same bandwidth as B.
!
!     S is treated as a product of elementary matrices:
!
!     S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n)
!
!     where S(i) is determined by the i-th row of S.
!
!     In phase 1, the index i takes the values n, n-1, ... , m+1;
!     in phase 2, it takes the values 1, 2, ... , m.
!
!     For each value of i, the current matrix A is updated by forming
!     inv(S(i))**T*A*inv(S(i)). This creates a triangular bulge outside
!     the band of A. The bulge is then pushed down toward the bottom of
!     A in phase 1, and up toward the top of A in phase 2, by applying
!     plane rotations.
!
!     There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1
!     of them are linearly independent, so annihilating a bulge requires
!     only 2*kb-1 plane rotations. The rotations are divided into a 1st
!     set of kb-1 rotations, and a 2nd set of kb rotations.
!
!     Wherever possible, rotations are generated and applied in vector
!     operations of length NR between the indices J1 and J2 (sometimes
!     replaced by modified values NRT, J1T or J2T).
!
!     The cosines and sines of the rotations are stored in the array
!     WORK. The cosines of the 1st set of rotations are stored in
!     elements n+2:n+m-kb-1 and the sines of the 1st set in elements
!     2:m-kb-1; the cosines of the 2nd set are stored in elements
!     n+m-kb+1:2*n and the sines of the second set in elements m-kb+1:n.
!
!     The bulges are not formed explicitly; nonzero elements outside the
!     band are created only when they are required for generating new
!     rotations; they are stored in the array WORK, in positions where
!     they are later overwritten by the sines of the rotations which
!     annihilate them.
!
!     **************************** Phase 1 *****************************
!
!     The logical structure of this phase is:
!
!     UPDATE = .TRUE.
!     DO I = N, M + 1, -1
!        use S(i) to update A and create a new bulge
!        apply rotations to push all bulges KA positions downward
!     END DO
!     UPDATE = .FALSE.
!     DO I = M + KA + 1, N - 1
!        apply rotations to push all bulges KA positions downward
!     END DO
!
!     To avoid duplicating code, the two loops are merged.
!
      UPDATE = .TRUE.
      I = N + 1
   10 CONTINUE
      if ( UPDATE ) then
         I = I - 1
         KBT = MIN( KB, I-1 )
         I0 = I - 1
         I1 = MIN( N, I+KA )
         I2 = I - KBT + KA1
         if ( I < M+1 ) then
            UPDATE = .FALSE.
            I = I + 1
            I0 = M
            if ( KA == 0 ) &
               GO TO 480
            GO TO 10
         end if
      ELSE
         I = I + KA
         if ( I.GT.N-1 ) &
            GO TO 480
      end if
!
      if ( UPPER ) then
!
!        Transform A, working with the upper triangle
!
         if ( UPDATE ) then
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
            BII = BB( KB1, I )
            DO 20 J = I, I1
               AB( I-J+KA1, J ) = AB( I-J+KA1, J ) / BII
   20       CONTINUE
            DO 30 J = MAX( 1, I-KA ), I
               AB( J-I+KA1, I ) = AB( J-I+KA1, I ) / BII
   30       CONTINUE
            DO 60 K = I - KBT, I - 1
               DO 40 J = I - KBT, K
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) - &
                                     BB( J-I+KB1, I )*AB( K-I+KA1, I ) - &
                                     BB( K-I+KB1, I )*AB( J-I+KA1, I ) + &
                                     AB( KA1, I )*BB( J-I+KB1, I )* &
                                     BB( K-I+KB1, I )
   40          CONTINUE
               DO 50 J = MAX( 1, I-KA ), I - KBT - 1
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) - &
                                     BB( K-I+KB1, I )*AB( J-I+KA1, I )
   50          CONTINUE
   60       CONTINUE
            DO 80 J = I, I1
               DO 70 K = MAX( J-KA, I-KBT ), I - 1
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) - &
                                     BB( K-I+KB1, I )*AB( I-J+KA1, J )
   70          CONTINUE
   80       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by inv(S(i))
!
               CALL DSCAL( N-M, ONE / BII, X( M+1, I ), 1 )
               if ( KBT.GT.0 ) &
                  CALL DGER( N-M, KBT, -ONE, X( M+1, I ), 1, &
                             BB( KB1-KBT, I ), 1, X( M+1, I-KBT ), LDX )
            end if
!
!           store a(i,i1) in RA1 for use in next loop over K
!
            RA1 = AB( I-I1+KA1, I1 )
         end if
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions down toward the bottom of the
!        band
!
         DO 130 K = 1, KB - 1
            if ( UPDATE ) then
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
               if ( I-K+KA < N .AND. I-K.GT.1 ) then
!
!                 generate rotation to annihilate a(i,i-k+ka+1)
!
                  CALL DLARTG( AB( K+1, I-K+KA ), RA1, &
                               WORK( N+I-K+KA-M ), WORK( I-K+KA-M ), &
                               RA )
!
!                 create nonzero element a(i-k,i-k+ka+1) outside the
!                 band and store it in WORK(i-k)
!
                  T = -BB( KB1-K, I )*RA1
                  WORK( I-K ) = WORK( N+I-K+KA-M )*T - &
                                WORK( I-K+KA-M )*AB( 1, I-K+KA )
                  AB( 1, I-K+KA ) = WORK( I-K+KA-M )*T + &
                                    WORK( N+I-K+KA-M )*AB( 1, I-K+KA )
                  RA1 = RA
               end if
            end if
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            if ( UPDATE ) then
               J2T = MAX( J2, I+2*KA-K+1 )
            ELSE
               J2T = J2
            end if
            NRT = ( N-J2T+KA ) / KA1
            DO 90 J = J2T, J1, KA1
!
!              create nonzero element a(j-ka,j+1) outside the band
!              and store it in WORK(j-m)
!
               WORK( J-M ) = WORK( J-M )*AB( 1, J+1 )
               AB( 1, J+1 ) = WORK( N+J-M )*AB( 1, J+1 )
   90       CONTINUE
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
            if ( NRT.GT.0 ) &
               CALL DLARGV( NRT, AB( 1, J2T ), INCA, WORK( J2T-M ), KA1, &
                            WORK( N+J2T-M ), KA1 )
            if ( NR.GT.0 ) then
!
!              apply rotations in 1st set from the right
!
               DO 100 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J2 ), INCA, &
                               AB( KA-L, J2+1 ), INCA, WORK( N+J2-M ), &
                               WORK( J2-M ), KA1 )
  100          CONTINUE
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( KA1, J2 ), AB( KA1, J2+1 ), &
                            AB( KA, J2+1 ), INCA, WORK( N+J2-M ), &
                            WORK( J2-M ), KA1 )
!
            end if
!
!           start applying rotations in 1st set from the left
!
            DO 110 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J2+KA1-L ), INCA, &
                               AB( L+1, J2+KA1-L ), INCA, &
                               WORK( N+J2-M ), WORK( J2-M ), KA1 )
  110       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 1st set
!
               DO 120 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, &
                             WORK( N+J-M ), WORK( J-M ) )
  120          CONTINUE
            end if
  130    CONTINUE
!
         if ( UPDATE ) then
            if ( I2.LE.N .AND. KBT.GT.0 ) then
!
!              create nonzero element a(i-kbt,i-kbt+ka+1) outside the
!              band and store it in WORK(i-kbt)
!
               WORK( I-KBT ) = -BB( KB1-KBT, I )*RA1
            end if
         end if
!
         DO 170 K = KB, 1, -1
            if ( UPDATE ) then
               J2 = I - K - 1 + MAX( 2, K-I0+1 )*KA1
            ELSE
               J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            end if
!
!           finish applying rotations in 2nd set from the left
!
            DO 140 L = KB - K, 1, -1
               NRT = ( N-J2+KA+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J2-L+1 ), INCA, &
                               AB( L+1, J2-L+1 ), INCA, WORK( N+J2-KA ), &
                               WORK( J2-KA ), KA1 )
  140       CONTINUE
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            DO 150 J = J1, J2, -KA1
               WORK( J ) = WORK( J-KA )
               WORK( N+J ) = WORK( N+J-KA )
  150       CONTINUE
            DO 160 J = J2, J1, KA1
!
!              create nonzero element a(j-ka,j+1) outside the band
!              and store it in WORK(j)
!
               WORK( J ) = WORK( J )*AB( 1, J+1 )
               AB( 1, J+1 ) = WORK( N+J )*AB( 1, J+1 )
  160       CONTINUE
            if ( UPDATE ) then
               if ( I-K < N-KA .AND. K.LE.KBT ) &
                  WORK( I-K+KA ) = WORK( I-K )
            end if
  170    CONTINUE
!
         DO 210 K = KB, 1, -1
            J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            if ( NR.GT.0 ) then
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
               CALL DLARGV( NR, AB( 1, J2 ), INCA, WORK( J2 ), KA1, &
                            WORK( N+J2 ), KA1 )
!
!              apply rotations in 2nd set from the right
!
               DO 180 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J2 ), INCA, &
                               AB( KA-L, J2+1 ), INCA, WORK( N+J2 ), &
                               WORK( J2 ), KA1 )
  180          CONTINUE
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( KA1, J2 ), AB( KA1, J2+1 ), &
                            AB( KA, J2+1 ), INCA, WORK( N+J2 ), &
                            WORK( J2 ), KA1 )
!
            end if
!
!           start applying rotations in 2nd set from the left
!
            DO 190 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J2+KA1-L ), INCA, &
                               AB( L+1, J2+KA1-L ), INCA, WORK( N+J2 ), &
                               WORK( J2 ), KA1 )
  190       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 2nd set
!
               DO 200 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, &
                             WORK( N+J ), WORK( J ) )
  200          CONTINUE
            end if
  210    CONTINUE
!
         DO 230 K = 1, KB - 1
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
!
!           finish applying rotations in 1st set from the left
!
            DO 220 L = KB - K, 1, -1
               NRT = ( N-J2+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J2+KA1-L ), INCA, &
                               AB( L+1, J2+KA1-L ), INCA, &
                               WORK( N+J2-M ), WORK( J2-M ), KA1 )
  220       CONTINUE
  230    CONTINUE
!
         if ( KB.GT.1 ) then
            DO 240 J = N - 1, I - KB + 2*KA + 1, -1
               WORK( N+J-M ) = WORK( N+J-KA-M )
               WORK( J-M ) = WORK( J-KA-M )
  240       CONTINUE
         end if
!
      ELSE
!
!        Transform A, working with the lower triangle
!
         if ( UPDATE ) then
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
            BII = BB( 1, I )
            DO 250 J = I, I1
               AB( J-I+1, I ) = AB( J-I+1, I ) / BII
  250       CONTINUE
            DO 260 J = MAX( 1, I-KA ), I
               AB( I-J+1, J ) = AB( I-J+1, J ) / BII
  260       CONTINUE
            DO 290 K = I - KBT, I - 1
               DO 270 J = I - KBT, K
                  AB( K-J+1, J ) = AB( K-J+1, J ) - &
                                   BB( I-J+1, J )*AB( I-K+1, K ) - &
                                   BB( I-K+1, K )*AB( I-J+1, J ) + &
                                   AB( 1, I )*BB( I-J+1, J )* &
                                   BB( I-K+1, K )
  270          CONTINUE
               DO 280 J = MAX( 1, I-KA ), I - KBT - 1
                  AB( K-J+1, J ) = AB( K-J+1, J ) - &
                                   BB( I-K+1, K )*AB( I-J+1, J )
  280          CONTINUE
  290       CONTINUE
            DO 310 J = I, I1
               DO 300 K = MAX( J-KA, I-KBT ), I - 1
                  AB( J-K+1, K ) = AB( J-K+1, K ) - &
                                   BB( I-K+1, K )*AB( J-I+1, I )
  300          CONTINUE
  310       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by inv(S(i))
!
               CALL DSCAL( N-M, ONE / BII, X( M+1, I ), 1 )
               if ( KBT.GT.0 ) &
                  CALL DGER( N-M, KBT, -ONE, X( M+1, I ), 1, &
                             BB( KBT+1, I-KBT ), LDBB-1, &
                             X( M+1, I-KBT ), LDX )
            end if
!
!           store a(i1,i) in RA1 for use in next loop over K
!
            RA1 = AB( I1-I+1, I )
         end if
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions down toward the bottom of the
!        band
!
         DO 360 K = 1, KB - 1
            if ( UPDATE ) then
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
               if ( I-K+KA < N .AND. I-K.GT.1 ) then
!
!                 generate rotation to annihilate a(i-k+ka+1,i)
!
                  CALL DLARTG( AB( KA1-K, I ), RA1, WORK( N+I-K+KA-M ), &
                               WORK( I-K+KA-M ), RA )
!
!                 create nonzero element a(i-k+ka+1,i-k) outside the
!                 band and store it in WORK(i-k)
!
                  T = -BB( K+1, I-K )*RA1
                  WORK( I-K ) = WORK( N+I-K+KA-M )*T - &
                                WORK( I-K+KA-M )*AB( KA1, I-K )
                  AB( KA1, I-K ) = WORK( I-K+KA-M )*T + &
                                   WORK( N+I-K+KA-M )*AB( KA1, I-K )
                  RA1 = RA
               end if
            end if
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            if ( UPDATE ) then
               J2T = MAX( J2, I+2*KA-K+1 )
            ELSE
               J2T = J2
            end if
            NRT = ( N-J2T+KA ) / KA1
            DO 320 J = J2T, J1, KA1
!
!              create nonzero element a(j+1,j-ka) outside the band
!              and store it in WORK(j-m)
!
               WORK( J-M ) = WORK( J-M )*AB( KA1, J-KA+1 )
               AB( KA1, J-KA+1 ) = WORK( N+J-M )*AB( KA1, J-KA+1 )
  320       CONTINUE
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
            if ( NRT.GT.0 ) &
               CALL DLARGV( NRT, AB( KA1, J2T-KA ), INCA, WORK( J2T-M ), &
                            KA1, WORK( N+J2T-M ), KA1 )
            if ( NR.GT.0 ) then
!
!              apply rotations in 1st set from the left
!
               DO 330 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J2-L ), INCA, &
                               AB( L+2, J2-L ), INCA, WORK( N+J2-M ), &
                               WORK( J2-M ), KA1 )
  330          CONTINUE
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( 1, J2 ), AB( 1, J2+1 ), AB( 2, J2 ), &
                            INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 )
!
            end if
!
!           start applying rotations in 1st set from the right
!
            DO 340 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J2 ), INCA, &
                               AB( KA1-L, J2+1 ), INCA, WORK( N+J2-M ), &
                               WORK( J2-M ), KA1 )
  340       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 1st set
!
               DO 350 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, &
                             WORK( N+J-M ), WORK( J-M ) )
  350          CONTINUE
            end if
  360    CONTINUE
!
         if ( UPDATE ) then
            if ( I2.LE.N .AND. KBT.GT.0 ) then
!
!              create nonzero element a(i-kbt+ka+1,i-kbt) outside the
!              band and store it in WORK(i-kbt)
!
               WORK( I-KBT ) = -BB( KBT+1, I-KBT )*RA1
            end if
         end if
!
         DO 400 K = KB, 1, -1
            if ( UPDATE ) then
               J2 = I - K - 1 + MAX( 2, K-I0+1 )*KA1
            ELSE
               J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            end if
!
!           finish applying rotations in 2nd set from the right
!
            DO 370 L = KB - K, 1, -1
               NRT = ( N-J2+KA+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J2-KA ), INCA, &
                               AB( KA1-L, J2-KA+1 ), INCA, &
                               WORK( N+J2-KA ), WORK( J2-KA ), KA1 )
  370       CONTINUE
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            DO 380 J = J1, J2, -KA1
               WORK( J ) = WORK( J-KA )
               WORK( N+J ) = WORK( N+J-KA )
  380       CONTINUE
            DO 390 J = J2, J1, KA1
!
!              create nonzero element a(j+1,j-ka) outside the band
!              and store it in WORK(j)
!
               WORK( J ) = WORK( J )*AB( KA1, J-KA+1 )
               AB( KA1, J-KA+1 ) = WORK( N+J )*AB( KA1, J-KA+1 )
  390       CONTINUE
            if ( UPDATE ) then
               if ( I-K < N-KA .AND. K.LE.KBT ) &
                  WORK( I-K+KA ) = WORK( I-K )
            end if
  400    CONTINUE
!
         DO 440 K = KB, 1, -1
            J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            if ( NR.GT.0 ) then
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
               CALL DLARGV( NR, AB( KA1, J2-KA ), INCA, WORK( J2 ), KA1, &
                            WORK( N+J2 ), KA1 )
!
!              apply rotations in 2nd set from the left
!
               DO 410 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J2-L ), INCA, &
                               AB( L+2, J2-L ), INCA, WORK( N+J2 ), &
                               WORK( J2 ), KA1 )
  410          CONTINUE
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( 1, J2 ), AB( 1, J2+1 ), AB( 2, J2 ), &
                            INCA, WORK( N+J2 ), WORK( J2 ), KA1 )
!
            end if
!
!           start applying rotations in 2nd set from the right
!
            DO 420 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J2 ), INCA, &
                               AB( KA1-L, J2+1 ), INCA, WORK( N+J2 ), &
                               WORK( J2 ), KA1 )
  420       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 2nd set
!
               DO 430 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, &
                             WORK( N+J ), WORK( J ) )
  430          CONTINUE
            end if
  440    CONTINUE
!
         DO 460 K = 1, KB - 1
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
!
!           finish applying rotations in 1st set from the right
!
            DO 450 L = KB - K, 1, -1
               NRT = ( N-J2+L ) / KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J2 ), INCA, &
                               AB( KA1-L, J2+1 ), INCA, WORK( N+J2-M ), &
                               WORK( J2-M ), KA1 )
  450       CONTINUE
  460    CONTINUE
!
         if ( KB.GT.1 ) then
            DO 470 J = N - 1, I - KB + 2*KA + 1, -1
               WORK( N+J-M ) = WORK( N+J-KA-M )
               WORK( J-M ) = WORK( J-KA-M )
  470       CONTINUE
         end if
!
      end if
!
      GO TO 10
!
  480 CONTINUE
!
!     **************************** Phase 2 *****************************
!
!     The logical structure of this phase is:
!
!     UPDATE = .TRUE.
!     DO I = 1, M
!        use S(i) to update A and create a new bulge
!        apply rotations to push all bulges KA positions upward
!     END DO
!     UPDATE = .FALSE.
!     DO I = M - KA - 1, 2, -1
!        apply rotations to push all bulges KA positions upward
!     END DO
!
!     To avoid duplicating code, the two loops are merged.
!
      UPDATE = .TRUE.
      I = 0
  490 CONTINUE
      if ( UPDATE ) then
         I = I + 1
         KBT = MIN( KB, M-I )
         I0 = I + 1
         I1 = MAX( 1, I-KA )
         I2 = I + KBT - KA1
         if ( I.GT.M ) then
            UPDATE = .FALSE.
            I = I - 1
            I0 = M + 1
            if ( KA == 0 ) &
               RETURN
            GO TO 490
         end if
      ELSE
         I = I - KA
         if ( I < 2 ) &
            RETURN
      end if
!
      if ( I < M-KBT ) then
         NX = M
      ELSE
         NX = N
      end if
!
      if ( UPPER ) then
!
!        Transform A, working with the upper triangle
!
         if ( UPDATE ) then
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
            BII = BB( KB1, I )
            DO 500 J = I1, I
               AB( J-I+KA1, I ) = AB( J-I+KA1, I ) / BII
  500       CONTINUE
            DO 510 J = I, MIN( N, I+KA )
               AB( I-J+KA1, J ) = AB( I-J+KA1, J ) / BII
  510       CONTINUE
            DO 540 K = I + 1, I + KBT
               DO 520 J = K, I + KBT
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) - &
                                     BB( I-J+KB1, J )*AB( I-K+KA1, K ) - &
                                     BB( I-K+KB1, K )*AB( I-J+KA1, J ) + &
                                     AB( KA1, I )*BB( I-J+KB1, J )* &
                                     BB( I-K+KB1, K )
  520          CONTINUE
               DO 530 J = I + KBT + 1, MIN( N, I+KA )
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) - &
                                     BB( I-K+KB1, K )*AB( I-J+KA1, J )
  530          CONTINUE
  540       CONTINUE
            DO 560 J = I1, I
               DO 550 K = I + 1, MIN( J+KA, I+KBT )
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) - &
                                     BB( I-K+KB1, K )*AB( J-I+KA1, I )
  550          CONTINUE
  560       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by inv(S(i))
!
               CALL DSCAL( NX, ONE / BII, X( 1, I ), 1 )
               if ( KBT.GT.0 ) &
                  CALL DGER( NX, KBT, -ONE, X( 1, I ), 1, BB( KB, I+1 ), &
                             LDBB-1, X( 1, I+1 ), LDX )
            end if
!
!           store a(i1,i) in RA1 for use in next loop over K
!
            RA1 = AB( I1-I+KA1, I )
         end if
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions up toward the top of the band
!
         DO 610 K = 1, KB - 1
            if ( UPDATE ) then
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
               if ( I+K-KA1.GT.0 .AND. I+K < M ) then
!
!                 generate rotation to annihilate a(i+k-ka-1,i)
!
                  CALL DLARTG( AB( K+1, I ), RA1, WORK( N+I+K-KA ), &
                               WORK( I+K-KA ), RA )
!
!                 create nonzero element a(i+k-ka-1,i+k) outside the
!                 band and store it in WORK(m-kb+i+k)
!
                  T = -BB( KB1-K, I+K )*RA1
                  WORK( M-KB+I+K ) = WORK( N+I+K-KA )*T - &
                                     WORK( I+K-KA )*AB( 1, I+K )
                  AB( 1, I+K ) = WORK( I+K-KA )*T + &
                                 WORK( N+I+K-KA )*AB( 1, I+K )
                  RA1 = RA
               end if
            end if
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            if ( UPDATE ) then
               J2T = MIN( J2, I-2*KA+K-1 )
            ELSE
               J2T = J2
            end if
            NRT = ( J2T+KA-1 ) / KA1
            DO 570 J = J1, J2T, KA1
!
!              create nonzero element a(j-1,j+ka) outside the band
!              and store it in WORK(j)
!
               WORK( J ) = WORK( J )*AB( 1, J+KA-1 )
               AB( 1, J+KA-1 ) = WORK( N+J )*AB( 1, J+KA-1 )
  570       CONTINUE
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
            if ( NRT.GT.0 ) &
               CALL DLARGV( NRT, AB( 1, J1+KA ), INCA, WORK( J1 ), KA1, &
                            WORK( N+J1 ), KA1 )
            if ( NR.GT.0 ) then
!
!              apply rotations in 1st set from the left
!
               DO 580 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J1+L ), INCA, &
                               AB( KA-L, J1+L ), INCA, WORK( N+J1 ), &
                               WORK( J1 ), KA1 )
  580          CONTINUE
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( KA1, J1 ), AB( KA1, J1-1 ), &
                            AB( KA, J1 ), INCA, WORK( N+J1 ), &
                            WORK( J1 ), KA1 )
!
            end if
!
!           start applying rotations in 1st set from the right
!
            DO 590 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J1T ), INCA, &
                               AB( L+1, J1T-1 ), INCA, WORK( N+J1T ), &
                               WORK( J1T ), KA1 )
  590       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 1st set
!
               DO 600 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1, &
                             WORK( N+J ), WORK( J ) )
  600          CONTINUE
            end if
  610    CONTINUE
!
         if ( UPDATE ) then
            if ( I2.GT.0 .AND. KBT.GT.0 ) then
!
!              create nonzero element a(i+kbt-ka-1,i+kbt) outside the
!              band and store it in WORK(m-kb+i+kbt)
!
               WORK( M-KB+I+KBT ) = -BB( KB1-KBT, I+KBT )*RA1
            end if
         end if
!
         DO 650 K = KB, 1, -1
            if ( UPDATE ) then
               J2 = I + K + 1 - MAX( 2, K+I0-M )*KA1
            ELSE
               J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            end if
!
!           finish applying rotations in 2nd set from the right
!
            DO 620 L = KB - K, 1, -1
               NRT = ( J2+KA+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J1T+KA ), INCA, &
                               AB( L+1, J1T+KA-1 ), INCA, &
                               WORK( N+M-KB+J1T+KA ), &
                               WORK( M-KB+J1T+KA ), KA1 )
  620       CONTINUE
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            DO 630 J = J1, J2, KA1
               WORK( M-KB+J ) = WORK( M-KB+J+KA )
               WORK( N+M-KB+J ) = WORK( N+M-KB+J+KA )
  630       CONTINUE
            DO 640 J = J1, J2, KA1
!
!              create nonzero element a(j-1,j+ka) outside the band
!              and store it in WORK(m-kb+j)
!
               WORK( M-KB+J ) = WORK( M-KB+J )*AB( 1, J+KA-1 )
               AB( 1, J+KA-1 ) = WORK( N+M-KB+J )*AB( 1, J+KA-1 )
  640       CONTINUE
            if ( UPDATE ) then
               if ( I+K.GT.KA1 .AND. K.LE.KBT ) &
                  WORK( M-KB+I+K-KA ) = WORK( M-KB+I+K )
            end if
  650    CONTINUE
!
         DO 690 K = KB, 1, -1
            J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            if ( NR.GT.0 ) then
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
               CALL DLARGV( NR, AB( 1, J1+KA ), INCA, WORK( M-KB+J1 ), &
                            KA1, WORK( N+M-KB+J1 ), KA1 )
!
!              apply rotations in 2nd set from the left
!
               DO 660 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J1+L ), INCA, &
                               AB( KA-L, J1+L ), INCA, &
                               WORK( N+M-KB+J1 ), WORK( M-KB+J1 ), KA1 )
  660          CONTINUE
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( KA1, J1 ), AB( KA1, J1-1 ), &
                            AB( KA, J1 ), INCA, WORK( N+M-KB+J1 ), &
                            WORK( M-KB+J1 ), KA1 )
!
            end if
!
!           start applying rotations in 2nd set from the right
!
            DO 670 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J1T ), INCA, &
                               AB( L+1, J1T-1 ), INCA, &
                               WORK( N+M-KB+J1T ), WORK( M-KB+J1T ), &
                               KA1 )
  670       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 2nd set
!
               DO 680 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1, &
                             WORK( N+M-KB+J ), WORK( M-KB+J ) )
  680          CONTINUE
            end if
  690    CONTINUE
!
         DO 710 K = 1, KB - 1
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
!
!           finish applying rotations in 1st set from the right
!
            DO 700 L = KB - K, 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( L, J1T ), INCA, &
                               AB( L+1, J1T-1 ), INCA, WORK( N+J1T ), &
                               WORK( J1T ), KA1 )
  700       CONTINUE
  710    CONTINUE
!
         if ( KB.GT.1 ) then
            DO 720 J = 2, MIN( I+KB, M ) - 2*KA - 1
               WORK( N+J ) = WORK( N+J+KA )
               WORK( J ) = WORK( J+KA )
  720       CONTINUE
         end if
!
      ELSE
!
!        Transform A, working with the lower triangle
!
         if ( UPDATE ) then
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
            BII = BB( 1, I )
            DO 730 J = I1, I
               AB( I-J+1, J ) = AB( I-J+1, J ) / BII
  730       CONTINUE
            DO 740 J = I, MIN( N, I+KA )
               AB( J-I+1, I ) = AB( J-I+1, I ) / BII
  740       CONTINUE
            DO 770 K = I + 1, I + KBT
               DO 750 J = K, I + KBT
                  AB( J-K+1, K ) = AB( J-K+1, K ) - &
                                   BB( J-I+1, I )*AB( K-I+1, I ) - &
                                   BB( K-I+1, I )*AB( J-I+1, I ) + &
                                   AB( 1, I )*BB( J-I+1, I )* &
                                   BB( K-I+1, I )
  750          CONTINUE
               DO 760 J = I + KBT + 1, MIN( N, I+KA )
                  AB( J-K+1, K ) = AB( J-K+1, K ) - &
                                   BB( K-I+1, I )*AB( J-I+1, I )
  760          CONTINUE
  770       CONTINUE
            DO 790 J = I1, I
               DO 780 K = I + 1, MIN( J+KA, I+KBT )
                  AB( K-J+1, J ) = AB( K-J+1, J ) - &
                                   BB( K-I+1, I )*AB( I-J+1, J )
  780          CONTINUE
  790       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by inv(S(i))
!
               CALL DSCAL( NX, ONE / BII, X( 1, I ), 1 )
               if ( KBT.GT.0 ) &
                  CALL DGER( NX, KBT, -ONE, X( 1, I ), 1, BB( 2, I ), 1, &
                             X( 1, I+1 ), LDX )
            end if
!
!           store a(i,i1) in RA1 for use in next loop over K
!
            RA1 = AB( I-I1+1, I1 )
         end if
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions up toward the top of the band
!
         DO 840 K = 1, KB - 1
            if ( UPDATE ) then
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
               if ( I+K-KA1.GT.0 .AND. I+K < M ) then
!
!                 generate rotation to annihilate a(i,i+k-ka-1)
!
                  CALL DLARTG( AB( KA1-K, I+K-KA ), RA1, &
                               WORK( N+I+K-KA ), WORK( I+K-KA ), RA )
!
!                 create nonzero element a(i+k,i+k-ka-1) outside the
!                 band and store it in WORK(m-kb+i+k)
!
                  T = -BB( K+1, I )*RA1
                  WORK( M-KB+I+K ) = WORK( N+I+K-KA )*T - &
                                     WORK( I+K-KA )*AB( KA1, I+K-KA )
                  AB( KA1, I+K-KA ) = WORK( I+K-KA )*T + &
                                      WORK( N+I+K-KA )*AB( KA1, I+K-KA )
                  RA1 = RA
               end if
            end if
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            if ( UPDATE ) then
               J2T = MIN( J2, I-2*KA+K-1 )
            ELSE
               J2T = J2
            end if
            NRT = ( J2T+KA-1 ) / KA1
            DO 800 J = J1, J2T, KA1
!
!              create nonzero element a(j+ka,j-1) outside the band
!              and store it in WORK(j)
!
               WORK( J ) = WORK( J )*AB( KA1, J-1 )
               AB( KA1, J-1 ) = WORK( N+J )*AB( KA1, J-1 )
  800       CONTINUE
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
            if ( NRT.GT.0 ) &
               CALL DLARGV( NRT, AB( KA1, J1 ), INCA, WORK( J1 ), KA1, &
                            WORK( N+J1 ), KA1 )
            if ( NR.GT.0 ) then
!
!              apply rotations in 1st set from the right
!
               DO 810 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J1 ), INCA, AB( L+2, J1-1 ), &
                               INCA, WORK( N+J1 ), WORK( J1 ), KA1 )
  810          CONTINUE
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( 1, J1 ), AB( 1, J1-1 ), &
                            AB( 2, J1-1 ), INCA, WORK( N+J1 ), &
                            WORK( J1 ), KA1 )
!
            end if
!
!           start applying rotations in 1st set from the left
!
            DO 820 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA, &
                               AB( KA1-L, J1T-KA1+L ), INCA, &
                               WORK( N+J1T ), WORK( J1T ), KA1 )
  820       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 1st set
!
               DO 830 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1, &
                             WORK( N+J ), WORK( J ) )
  830          CONTINUE
            end if
  840    CONTINUE
!
         if ( UPDATE ) then
            if ( I2.GT.0 .AND. KBT.GT.0 ) then
!
!              create nonzero element a(i+kbt,i+kbt-ka-1) outside the
!              band and store it in WORK(m-kb+i+kbt)
!
               WORK( M-KB+I+KBT ) = -BB( KBT+1, I )*RA1
            end if
         end if
!
         DO 880 K = KB, 1, -1
            if ( UPDATE ) then
               J2 = I + K + 1 - MAX( 2, K+I0-M )*KA1
            ELSE
               J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            end if
!
!           finish applying rotations in 2nd set from the left
!
            DO 850 L = KB - K, 1, -1
               NRT = ( J2+KA+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J1T+L-1 ), INCA, &
                               AB( KA1-L, J1T+L-1 ), INCA, &
                               WORK( N+M-KB+J1T+KA ), &
                               WORK( M-KB+J1T+KA ), KA1 )
  850       CONTINUE
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            DO 860 J = J1, J2, KA1
               WORK( M-KB+J ) = WORK( M-KB+J+KA )
               WORK( N+M-KB+J ) = WORK( N+M-KB+J+KA )
  860       CONTINUE
            DO 870 J = J1, J2, KA1
!
!              create nonzero element a(j+ka,j-1) outside the band
!              and store it in WORK(m-kb+j)
!
               WORK( M-KB+J ) = WORK( M-KB+J )*AB( KA1, J-1 )
               AB( KA1, J-1 ) = WORK( N+M-KB+J )*AB( KA1, J-1 )
  870       CONTINUE
            if ( UPDATE ) then
               if ( I+K.GT.KA1 .AND. K.LE.KBT ) &
                  WORK( M-KB+I+K-KA ) = WORK( M-KB+I+K )
            end if
  880    CONTINUE
!
         DO 920 K = KB, 1, -1
            J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            if ( NR.GT.0 ) then
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
               CALL DLARGV( NR, AB( KA1, J1 ), INCA, WORK( M-KB+J1 ), &
                            KA1, WORK( N+M-KB+J1 ), KA1 )
!
!              apply rotations in 2nd set from the right
!
               DO 890 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J1 ), INCA, AB( L+2, J1-1 ), &
                               INCA, WORK( N+M-KB+J1 ), WORK( M-KB+J1 ), &
                               KA1 )
  890          CONTINUE
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
               CALL DLAR2V( NR, AB( 1, J1 ), AB( 1, J1-1 ), &
                            AB( 2, J1-1 ), INCA, WORK( N+M-KB+J1 ), &
                            WORK( M-KB+J1 ), KA1 )
!
            end if
!
!           start applying rotations in 2nd set from the left
!
            DO 900 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA, &
                               AB( KA1-L, J1T-KA1+L ), INCA, &
                               WORK( N+M-KB+J1T ), WORK( M-KB+J1T ), &
                               KA1 )
  900       CONTINUE
!
            if ( WANTX ) then
!
!              post-multiply X by product of rotations in 2nd set
!
               DO 910 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1, &
                             WORK( N+M-KB+J ), WORK( M-KB+J ) )
  910          CONTINUE
            end if
  920    CONTINUE
!
         DO 940 K = 1, KB - 1
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
!
!           finish applying rotations in 1st set from the left
!
            DO 930 L = KB - K, 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               if ( NRT.GT.0 ) &
                  CALL DLARTV( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA, &
                               AB( KA1-L, J1T-KA1+L ), INCA, &
                               WORK( N+J1T ), WORK( J1T ), KA1 )
  930       CONTINUE
  940    CONTINUE
!
         if ( KB.GT.1 ) then
            DO 950 J = 2, MIN( I+KB, M ) - 2*KA - 1
               WORK( N+J ) = WORK( N+J+KA )
               WORK( J ) = WORK( J+KA )
  950       CONTINUE
         end if
!
      end if
!
      GO TO 490
!
!     End of DSBGST
!
      END
