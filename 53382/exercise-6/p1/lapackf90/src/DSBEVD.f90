      SUBROUTINE DSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, &
                         LWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSBEVD computes all the eigenvalues and, optionally, eigenvectors of
!  a real symmetric band matrix A. If eigenvectors are desired, it uses
!  a divide and conquer algorithm.
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
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, AB is overwritten by values generated during the
!          reduction to tridiagonal form.  If UPLO = 'U', the first
!          superdiagonal and the diagonal of the tridiagonal matrix T
!          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
!          the diagonal and first subdiagonal of T are returned in the
!          first two rows of AB.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD + 1.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
!          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!          eigenvectors of the matrix A, with the i-th column of Z
!          holding the eigenvector associated with W(i).
!          If JOBZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array,
!                                         dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          IF N <= 1,                LWORK must be at least 1.
!          If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N.
!          If JOBZ  = 'V' and N > 2, LWORK must be at least
!                         ( 1 + 5*N + 2*N**2 ).
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
!          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array LIWORK.
!          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.
!          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N.
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, &
                         LLWRK2, LWMIN
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                         SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSB
      EXTERNAL           LSAME, DLAMCH, DLANSB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASCL, DSBTRD, DSCAL, DSTEDC, &
                         DSTERF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
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
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( KD < 0 ) then
         INFO = -4
      else if ( LDAB < KD+1 ) then
         INFO = -6
      else if ( LDZ < 1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) then
         INFO = -9
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -11
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -13
      end if
!
      if ( INFO == 0 ) then
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSBEVD', -INFO )
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
      if ( N == 1 ) then
         W( 1 ) = AB( 1, 1 )
         if ( WANTZ ) &
            Z( 1, 1 ) = ONE
         RETURN
      end if
!
!     Get machine constants.
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
      ANRM = DLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM < RMIN ) then
         ISCALE = 1
         SIGMA = RMIN / ANRM
      else if ( ANRM.GT.RMAX ) then
         ISCALE = 1
         SIGMA = RMAX / ANRM
      end if
      if ( ISCALE == 1 ) then
         if ( LOWER ) then
            CALL DLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         ELSE
            CALL DLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         end if
      end if
!
!     Call DSBTRD to reduce symmetric band matrix to tridiagonal form.
!
      INDE = 1
      INDWRK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWRK2 = LWORK - INDWK2 + 1
      CALL DSBTRD( JOBZ, UPLO, N, KD, AB, LDAB, W, WORK( INDE ), Z, LDZ, &
                   WORK( INDWRK ), IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC.
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
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      if ( ISCALE == 1 ) &
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
!
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
!
!     End of DSBEVD
!
      END
