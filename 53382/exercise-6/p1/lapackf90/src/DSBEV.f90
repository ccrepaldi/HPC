      SUBROUTINE DSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, &
                        INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KD, LDAB, LDZ, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSBEV computes all the eigenvalues and, optionally, eigenvectors of
!  a real symmetric band matrix A.
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
!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,3*N-2))
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
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDWRK, ISCALE
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                         SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSB
      EXTERNAL           LSAME, DLAMCH, DLANSB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASCL, DSBTRD, DSCAL, DSTEQR, DSTERF, XERBLA
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
!
      INFO = 0
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
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSBEV ', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
      if ( N == 1 ) then
         if ( LOWER ) then
            W( 1 ) = AB( 1, 1 )
         ELSE
            W( 1 ) = AB( KD+1, 1 )
         end if
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
      CALL DSBTRD( JOBZ, UPLO, N, KD, AB, LDAB, W, WORK( INDE ), Z, LDZ, &
                   WORK( INDWRK ), IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.
!
      if ( .NOT.WANTZ ) then
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), &
                      INFO )
      end if
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      if ( ISCALE == 1 ) then
         if ( INFO == 0 ) then
            IMAX = N
         ELSE
            IMAX = INFO - 1
         end if
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      end if
!
      RETURN
!
!     End of DSBEV
!
      END
