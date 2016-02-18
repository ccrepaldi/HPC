      SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, &
                         IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSPEVD computes all the eigenvalues and, optionally, eigenvectors
!  of a real symmetric matrix A in packed storage. If eigenvectors are
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
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!
!          On exit, AP is overwritten by values generated during the
!          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!          and first superdiagonal of the tridiagonal matrix T overwrite
!          the corresponding elements of A, and if UPLO = 'L', the
!          diagonal and first subdiagonal of T overwrite the
!          corresponding elements of A.
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
!          If N <= 1,               LWORK must be at least 1.
!          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N.
!          If JOBZ = 'V' and N > 1, LWORK must be at least
!                                                 1 + 6*N + N**2.
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
!          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.
!          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
      LOGICAL            LQUERY, WANTZ
      INTEGER            IINFO, INDE, INDTAU, INDWRK, ISCALE, LIWMIN, &
                         LLWORK, LWMIN
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                         SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSP
      EXTERNAL           LSAME, DLAMCH, DLANSP
!     ..
!     .. External Subroutines ..
      EXTERNAL           DOPMTR, DSCAL, DSPTRD, DSTEDC, DSTERF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      LQUERY = ( LWORK == -1 .OR. LIWORK.EQ.-1 )
!
      INFO = 0
      if ( N.LE.1 ) then
         LIWMIN = 1
         LWMIN = 1
      ELSE
         if ( WANTZ ) then
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 6*N + N**2
         ELSE
            LIWMIN = 1
            LWMIN = 2*N
         end if
      end if
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( .NOT.( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) &
                THEN
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( LDZ < 1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) then
         INFO = -7
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -9
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -11
      end if
!
      if ( INFO == 0 ) then
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSPEVD', -INFO )
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
         W( 1 ) = AP( 1 )
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
      ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM < RMIN ) then
         ISCALE = 1
         SIGMA = RMIN / ANRM
      else if ( ANRM.GT.RMAX ) then
         ISCALE = 1
         SIGMA = RMAX / ANRM
      end if
      if ( ISCALE == 1 ) then
         CALL DSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
      end if
!
!     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
!
      INDE = 1
      INDTAU = INDE + N
      CALL DSPTRD( UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
!     tridiagonal matrix, then call DOPMTR to multiply it by the
!     Householder transformations represented in AP.
!
      if ( .NOT.WANTZ ) then
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         INDWRK = INDTAU + N
         LLWORK = LWORK - INDWRK + 1
         CALL DSTEDC( 'I', N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), &
                      LLWORK, IWORK, LIWORK, INFO )
         CALL DOPMTR( 'L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ, &
                      WORK( INDWRK ), IINFO )
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
!     End of DSPEVD
!
      END
