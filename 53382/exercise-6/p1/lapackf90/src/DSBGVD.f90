      SUBROUTINE DSBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, &
                         Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), W( * ), &
                         WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSBGVD computes all the eigenvalues, and optionally, the eigenvectors
!  of a real generalized symmetric-definite banded eigenproblem, of the
!  form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and
!  banded, and B is also positive definite.  If eigenvectors are
!  desired, it uses a divide and conquer algorithm.
!
!  The divide and conquer algorithm makes very mild assumptions about
!  floating point arithmetic. It will work on machines with a guard
!  digit in add/subtract, or on those binary machines without guard
!  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!  Cray-2. It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangles of A and B are stored;
!          = 'L':  Lower triangles of A and B are stored.
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
!          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first ka+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
!
!          On exit, the contents of AB are destroyed.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KA+1.
!
!  BB      (input/output) DOUBLE PRECISION array, dimension (LDBB, N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix B, stored in the first kb+1 rows of the array.  The
!          j-th column of B is stored in the j-th column of the array BB
!          as follows:
!          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
!          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).
!
!          On exit, the factor S from the split Cholesky factorization
!          B = S**T*S, as returned by DPBSTF.
!
!  LDBB    (input) INTEGER
!          The leading dimension of the array BB.  LDBB >= KB+1.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
!          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
!          eigenvectors, with the i-th column of Z holding the
!          eigenvector associated with W(i).  The eigenvectors are
!          normalized so Z**T*B*Z = I.
!          If JOBZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If N <= 1,               LWORK >= 1.
!          If JOBZ = 'N' and N > 1, LWORK >= 3*N.
!          If JOBZ = 'V' and N > 1, LWORK >= 1 + 5*N + 2*N**2.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
!          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.
!          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, and i is:
!             <= N:  the algorithm failed to converge:
!                    i off-diagonal elements of an intermediate
!                    tridiagonal form did not converge to zero;
!             > N:   if INFO = N + i, for 1 <= i <= N, then DPBSTF
!                    returned INFO = i: B is not positive definite.
!                    The factorization of B could not be completed and
!                    no eigenvalues or eigenvectors were computed.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, WANTZ
      CHARACTER          VECT
      INTEGER            IINFO, INDE, INDWK2, INDWRK, LIWMIN, LLWRK2, &
                         LWMIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DPBSTF, DSBGST, DSBTRD, DSTEDC, &
                         DSTERF, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 .OR. LIWORK.EQ.-1 )
!
      INFO = 0
      if ( N.LE.1 ) then
         LIWMIN = 1
         LWMIN = 1
      ELSE
         if ( WANTZ ) then
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 5*N + 2*N**2
         ELSE
            LIWMIN = 1
            LWMIN = 2*N
         end if
      end if
!
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( KA < 0 ) then
         INFO = -4
      else if ( KB < 0 .OR. KB.GT.KA ) then
         INFO = -5
      else if ( LDAB < KA+1 ) then
         INFO = -7
      else if ( LDBB < KB+1 ) then
         INFO = -9
      else if ( LDZ < 1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) then
         INFO = -12
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -14
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -16
      end if
!
      if ( INFO == 0 ) then
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSBGVD', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Form a split Cholesky factorization of B.
!
      CALL DPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      if ( INFO.NE.0 ) then
         INFO = N + INFO
         RETURN
      end if
!
!     Transform problem to standard eigenvalue problem.
!
      INDE = 1
      INDWRK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWRK2 = LWORK - INDWK2 + 1
      CALL DSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, &
                   WORK( INDWRK ), IINFO )
!
!     Reduce to tridiagonal form.
!
      if ( WANTZ ) then
         VECT = 'U'
      ELSE
         VECT = 'N'
      end if
      CALL DSBTRD( VECT, UPLO, N, KA, AB, LDAB, W, WORK( INDE ), Z, LDZ, &
                   WORK( INDWRK ), IINFO )
!
!     For eigenvalues only, call DSTERF. For eigenvectors, call SSTEDC.
!
      if ( .NOT.WANTZ ) then
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DSTEDC( 'I', N, W, WORK( INDE ), WORK( INDWRK ), N, &
                      WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO )
         CALL DGEMM( 'N', 'N', N, N, N, ONE, Z, LDZ, WORK( INDWRK ), N, &
                     ZERO, WORK( INDWK2 ), N )
         CALL DLACPY( 'A', N, N, WORK( INDWK2 ), N, Z, LDZ )
      end if
!
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of DSBGVD
!
      END