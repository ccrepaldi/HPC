      SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
                         LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
!  real symmetric matrix A. If eigenvectors are desired, it uses a
!  divide and conquer algorithm.
!
!  The divide and conquer algorithm makes very mild assumptions about
!  floating point arithmetic. It will work on machines with a guard
!  digit in add/subtract, or on those binary machines without guard
!  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!  Cray-2. It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Because of large use of BLAS of level 3, DSYEVD needs N**2 more
!  workspace than DSYEVX.
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
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) DOUBLE PRECISION array,
!                                         dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If N <= 1,               LWORK must be at least 1.
!          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
!          If JOBZ = 'V' and N > 1, LWORK must be at least
!                                                1 + 6*N + 2*N**2.
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
!          The dimension of the array IWORK.
!          If N <= 1,                LIWORK must be at least 1.
!          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
!          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
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
!  Further Details
!  ===============
!
!  Based on contributions by
!     Jeff Rutter, Computer Science Division, University of California
!     at Berkeley, USA
!  Modified by Francoise Tisseur, University of Tennessee.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
!
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, INDE, INDTAU, INDWK2, INDWRK, ISCALE, &
                         LIOPT, LIWMIN, LLWORK, LLWRK2, LOPT, LWMIN
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                         SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, DLAMCH, DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLACPY, DLASCL, DORMTR, DSCAL, DSTEDC, DSTERF, &
                         DSYTRD, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
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
         LOPT = LWMIN
         LIOPT = LIWMIN
      ELSE
         if ( WANTZ ) then
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 6*N + 2*N**2
         ELSE
            LIWMIN = 1
            LWMIN = 2*N + 1
         end if
         LOPT = LWMIN
         LIOPT = LIWMIN
      end if
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -5
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -8
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -10
      end if
!
      if ( INFO == 0 ) then
         WORK( 1 ) = LOPT
         IWORK( 1 ) = LIOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSYEVD', -INFO )
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
         W( 1 ) = A( 1, 1 )
         if ( WANTZ ) &
            A( 1, 1 ) = ONE
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
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM < RMIN ) then
         ISCALE = 1
         SIGMA = RMIN / ANRM
      else if ( ANRM.GT.RMAX ) then
         ISCALE = 1
         SIGMA = RMAX / ANRM
      end if
      if ( ISCALE == 1 ) &
         CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
!
!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
!
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      INDWK2 = INDWRK + N*N
      LLWRK2 = LWORK - INDWK2 + 1
!
      CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), &
                   WORK( INDWRK ), LLWORK, IINFO )
      LOPT = 2*N + WORK( INDWRK )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
!     tridiagonal matrix, then call DORMTR to multiply it by the
!     Householder transformations stored in A.
!
      if ( .NOT.WANTZ ) then
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DSTEDC( 'I', N, W, WORK( INDE ), WORK( INDWRK ), N, &
                      WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO )
         CALL DORMTR( 'L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), &
                      WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO )
         CALL DLACPY( 'A', N, N, WORK( INDWRK ), N, A, LDA )
         LOPT = MAX( LOPT, 1+6*N+2*N**2 )
      end if
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      if ( ISCALE == 1 ) &
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
!
      WORK( 1 ) = LOPT
      IWORK( 1 ) = LIOPT
!
      RETURN
!
!     End of DSYEVD
!
      END
