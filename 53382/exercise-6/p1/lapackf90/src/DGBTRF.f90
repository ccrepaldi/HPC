      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
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
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
                         JU, K2, KM, KV, NB, NW
      DOUBLE PRECISION   TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   WORK13( LDWORK, NBMAX ), &
                         WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX, ILAENV
      EXTERNAL           IDAMAX, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL, &
                         DSWAP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( KL < 0 ) then
         INFO = -3
      else if ( KU < 0 ) then
         INFO = -4
      else if ( LDAB < KL+KV+1 ) then
         INFO = -6
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGBTRF', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      if ( NB.LE.1 .OR. NB.GT.KL ) then
!
!        Use unblocked code
!
         CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            WORK13(1:j-1, J ) = ZERO
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               if ( JJ+KV.LE.N ) then
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               end if
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               if ( AB( KV+JP, JJ ).NE.ZERO ) then
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  if ( JP.NE.1 ) then
!
!                    Apply interchange to columns J to J+JB-1
!
                     if ( JP+JJ-1 < J+KL ) then
!
                        CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                                    AB( KV+JP, JJ ), LDAB-1 )
                     end if
                  end if
!
!                 Compute multipliers
!
                  CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
                              1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  if ( JM.GT.JJ ) &
                     CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
                                AB( KV, JJ+1 ), LDAB-1, &
                                AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  if ( INFO == 0 ) &
                     INFO = JJ
               end if
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               if ( NW.GT.0 ) &
                  CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
                              WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            if ( J+JB.LE.N ) then
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use DLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
                            IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     if ( IP.NE.II ) then
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     end if
                  end do
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               if ( J2.GT.0 ) then
!
!                 Update A12
!
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J2, ONE, AB( KV+1, J ), LDAB-1, &
                              AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  if ( I2.GT.0 ) then
!
!                    Update A22
!
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J2, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+1, J+JB ), LDAB-1 )
                  end if
!
                  if ( I3.GT.0 ) then
!
!                    Update A32
!
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J2, &
                                 JB, -ONE, WORK31, LDWORK, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  end if
               end if
!
               if ( J3.GT.0 ) then
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J3, ONE, AB( KV+1, J ), LDAB-1, &
                              WORK13, LDWORK )
!
                  if ( I2.GT.0 ) then
!
!                    Update A23
!
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J3, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
                                 LDAB-1 )
                  end if
!
                  if ( I3.GT.0 ) then
!
!                    Update A33
!
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J3, &
                                 JB, -ONE, WORK31, LDWORK, WORK13, &
                                 LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  end if
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               end if
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            end if
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               if ( JP.NE.1 ) then
!
!                 Apply interchange to columns J to JJ-1
!
                  if ( JP+JJ-1 < J+KL ) then
!
!                    The interchange does not affect A31
!
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  end if
               end if
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               if ( NW.GT.0 ) &
                  CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
                              AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      end if
!
      RETURN
!
!     End of DGBTRF
!
      END
