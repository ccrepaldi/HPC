      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTRF computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  Contributed by
!  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, IB, II, J, JJ, NB
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   WORK( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DPBTF2, DPOTF2, DSYRK, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      if ( ( .NOT.LSAME( UPLO, 'U' ) ) .AND. &
          ( .NOT.LSAME( UPLO, 'L' ) ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( KD < 0 ) then
         INFO = -3
      else if ( LDAB < KD+1 ) then
         INFO = -5
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPBTRF', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
      NB = MIN( NB, NBMAX )
!
      if ( NB.LE.1 .OR. NB.GT.KD ) then
!
!        Use unblocked code
!
         CALL DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
!
!        Use blocked code
!
         if ( LSAME( UPLO, 'U' ) ) then
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               if ( II.NE.0 ) then
                  INFO = I + II - 1
                  GO TO 150
               end if
               if ( I+IB.LE.N ) then
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  if ( I2.GT.0 ) then
!
!                    Update A12
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose', &
                                 'Non-unit', IB, I2, ONE, AB( KD+1, I ), &
                                 LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Upper', 'Transpose', I2, IB, -ONE, &
                                 AB( KD+1-IB, I+IB ), LDAB-1, ONE, &
                                 AB( KD+1, I+IB ), LDAB-1 )
                  end if
!
                  if ( I3.GT.0 ) then
!
!                    Copy the lower triangle of A13 into the work array.
!
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE
!
!                    Update A13 (in the work array).
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose', &
                                 'Non-unit', IB, I3, ONE, AB( KD+1, I ), &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A23
!
                     if ( I2.GT.0 ) &
                        CALL DGEMM( 'Transpose', 'No Transpose', I2, I3, &
                                    IB, -ONE, AB( KD+1-IB, I+IB ), &
                                    LDAB-1, WORK, LDWORK, ONE, &
                                    AB( 1+IB, I+KD ), LDAB-1 )
!
!                    Update A33
!
                     CALL DSYRK( 'Upper', 'Transpose', I3, IB, -ONE, &
                                 WORK, LDWORK, ONE, AB( KD+1, I+KD ), &
                                 LDAB-1 )
!
!                    Copy the lower triangle of A13 back into place.
!
                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  end if
               end if
   70       CONTINUE
         ELSE
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               if ( II.NE.0 ) then
                  INFO = I + II - 1
                  GO TO 150
               end if
               if ( I+IB.LE.N ) then
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  if ( I2.GT.0 ) then
!
!                    Update A21
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose', &
                                 'Non-unit', I2, IB, ONE, AB( 1, I ), &
                                 LDAB-1, AB( 1+IB, I ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Lower', 'No Transpose', I2, IB, -ONE, &
                                 AB( 1+IB, I ), LDAB-1, ONE, &
                                 AB( 1, I+IB ), LDAB-1 )
                  end if
!
                  if ( I3.GT.0 ) then
!
!                    Copy the upper triangle of A31 into the work array.
!
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE
!
!                    Update A31 (in the work array).
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose', &
                                 'Non-unit', I3, IB, ONE, AB( 1, I ), &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A32
!
                     if ( I2.GT.0 ) &
                        CALL DGEMM( 'No transpose', 'Transpose', I3, I2, &
                                    IB, -ONE, WORK, LDWORK, &
                                    AB( 1+IB, I ), LDAB-1, ONE, &
                                    AB( 1+KD-IB, I+IB ), LDAB-1 )
!
!                    Update A33
!
                     CALL DSYRK( 'Lower', 'No Transpose', I3, IB, -ONE, &
                                 WORK, LDWORK, ONE, AB( 1, I+KD ), &
                                 LDAB-1 )
!
!                    Copy the upper triangle of A31 back into place.
!
                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  end if
               end if
  140       CONTINUE
         end if
      end if
      RETURN
!
  150 CONTINUE
      RETURN
!
!     End of DPBTRF
!
      END