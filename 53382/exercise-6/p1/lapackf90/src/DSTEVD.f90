      SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, &
                         LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ
      INTEGER            INFO, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEVD computes all eigenvalues and, optionally, eigenvectors of a
!  real symmetric tridiagonal matrix. If eigenvectors are desired, it
!  uses a divide and conquer algorithm.
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
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix
!          A.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix A, stored in elements 1 to N-1 of E; E(N) need not
!          be set, but is used by the routine.
!          On exit, the contents of E are destroyed.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
!          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!          eigenvectors of the matrix A, with the i-th column of Z
!          holding the eigenvector associated with D(i).
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
!          If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1.
!          If JOBZ  = 'V' and N > 1 then LWORK must be at least
!                         ( 1 + 4*N + N**2 ).
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
!          If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1.
!          If JOBZ  = 'V' and N > 1 then LIWORK must be at least 3+5*N.
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
!                off-diagonal elements of E did not converge to zero.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, WANTZ
      INTEGER            ISCALE, LIWMIN, LWMIN
      DOUBLE PRECISION   BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, &
                         TNRM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTEDC, DSTERF, XERBLA
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
      LIWMIN = 1
      LWMIN = 1
      if ( N.GT.1 .AND. WANTZ ) then
         LWMIN = 1 + 4*N + N**2
         LIWMIN = 3 + 5*N
      end if
!
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDZ < 1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) then
         INFO = -6
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -8
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -10
      end if
!
      if ( INFO == 0 ) then
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSTEVD', -INFO )
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
      ISCALE = 0
      TNRM = DLANST( 'M', N, D, E )
      if ( TNRM.GT.ZERO .AND. TNRM < RMIN ) then
         ISCALE = 1
         SIGMA = RMIN / TNRM
      else if ( TNRM.GT.RMAX ) then
         ISCALE = 1
         SIGMA = RMAX / TNRM
      end if
      if ( ISCALE == 1 ) then
         CALL DSCAL( N, SIGMA, D, 1 )
         CALL DSCAL( N-1, SIGMA, E( 1 ), 1 )
      end if
!
!     For eigenvalues only, call DSTERF.  For eigenvalues and
!     eigenvectors, call DSTEDC.
!
      if ( .NOT.WANTZ ) then
         CALL DSTERF( N, D, E, INFO )
      ELSE
         CALL DSTEDC( 'I', N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, &
                      INFO )
      end if
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      if ( ISCALE == 1 ) &
         CALL DSCAL( N, ONE / SIGMA, D, 1 )
!
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of DSTEVD
!
      END
