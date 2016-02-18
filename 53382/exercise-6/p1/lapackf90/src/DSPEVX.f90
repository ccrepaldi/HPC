      SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, &
                         ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, &
                         INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDZ, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSPEVX computes selected eigenvalues and, optionally, eigenvectors
!  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
!  can be selected by specifying either a range of values or a range of
!  indices for the desired eigenvalues.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  RANGE   (input) CHARACTER*1
!          = 'A': all eigenvalues will be found;
!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!                 will be found;
!          = 'I': the IL-th through IU-th eigenvalues will be found.
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
!  VL      (input) DOUBLE PRECISION
!  VU      (input) DOUBLE PRECISION
!          If RANGE='V', the lower and upper bounds of the interval to
!          be searched for eigenvalues. VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!
!  IL      (input) INTEGER
!  IU      (input) INTEGER
!          If RANGE='I', the indices (in ascending order) of the
!          smallest and largest eigenvalues to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The absolute error tolerance for the eigenvalues.
!          An approximate eigenvalue is accepted as converged
!          when it is determined to lie in an interval [a,b]
!          of width less than or equal to
!
!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!
!          where EPS is the machine precision.  If ABSTOL is less than
!          or equal to zero, then  EPS*|T|  will be used in its place,
!          where |T| is the 1-norm of the tridiagonal matrix obtained
!          by reducing AP to tridiagonal form.
!
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!          If this routine returns with INFO>0, indicating that some
!          eigenvectors did not converge, try setting ABSTOL to
!          2*DLAMCH('S').
!
!          See "Computing Small Singular Values of Bidiagonal Matrices
!          with Guaranteed High Relative Accuracy," by Demmel and
!          Kahan, LAPACK Working Note #3.
!
!  M       (output) INTEGER
!          The total number of eigenvalues found.  0 <= M <= N.
!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the selected eigenvalues in ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!          contain the orthonormal eigenvectors of the matrix A
!          corresponding to the selected eigenvalues, with the i-th
!          column of Z holding the eigenvector associated with W(i).
!          If an eigenvector fails to converge, then that column of Z
!          contains the latest approximation to the eigenvector, and the
!          index of the eigenvector is returned in IFAIL.
!          If JOBZ = 'N', then Z is not referenced.
!          Note: the user must ensure that at least max(1,M) columns are
!          supplied in the array Z; if RANGE = 'V', the exact value of M
!          is not known in advance and an upper bound must be used.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
!
!  IWORK   (workspace) INTEGER array, dimension (5*N)
!
!  IFAIL   (output) INTEGER array, dimension (N)
!          If JOBZ = 'V', then if INFO = 0, the first M elements of
!          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!          indices of the eigenvectors that failed to converge.
!          If JOBZ = 'N', then IFAIL is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, then i eigenvectors failed to converge.
!                Their indices are stored in array IFAIL.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, &
                         INDISP, INDIWO, INDTAU, INDWRK, ISCALE, ITMP1, &
                         J, JJ, NSPLIT
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, &
                         SIGMA, SMLNUM, TMP1, VLL, VUU
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSP
      EXTERNAL           LSAME, DLAMCH, DLANSP
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DOPGTR, DOPMTR, DSCAL, DSPTRD, DSTEBZ, &
                         DSTEIN, DSTEQR, DSTERF, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
!
      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) then
         INFO = -2
      else if ( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) ) &
                THEN
         INFO = -3
      else if ( N < 0 ) then
         INFO = -4
      ELSE
         if ( VALEIG ) then
            if ( N.GT.0 .AND. VU.LE.VL ) &
               INFO = -7
         else if ( INDEIG ) then
            if ( IL < 1 .OR. IL.GT.MAX( 1, N ) ) then
               INFO = -8
            else if ( IU < MIN( N, IL ) .OR. IU.GT.N ) then
               INFO = -9
            end if
         end if
      end if
      if ( INFO == 0 ) then
         if ( LDZ < 1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) &
            INFO = -14
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSPEVX', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      M = 0
      if ( N == 0 ) &
         RETURN
!
      if ( N == 1 ) then
         if ( ALLEIG .OR. INDEIG ) then
            M = 1
            W( 1 ) = AP( 1 )
         ELSE
            if ( VL < AP( 1 ) .AND. VU.GE.AP( 1 ) ) then
               M = 1
               W( 1 ) = AP( 1 )
            end if
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
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
!
!     Scale matrix to allowable range, if necessary.
!
      ISCALE = 0
      ABSTLL = ABSTOL
      if ( VALEIG ) then
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      end if
      ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
      if ( ANRM.GT.ZERO .AND. ANRM < RMIN ) then
         ISCALE = 1
         SIGMA = RMIN / ANRM
      else if ( ANRM.GT.RMAX ) then
         ISCALE = 1
         SIGMA = RMAX / ANRM
      end if
      if ( ISCALE == 1 ) then
         CALL DSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
         if ( ABSTOL.GT.0 ) &
            ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) then
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         end if
      end if
!
!     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
!
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDWRK = INDD + N
      CALL DSPTRD( UPLO, N, AP, WORK( INDD ), WORK( INDE ), &
                   WORK( INDTAU ), IINFO )
!
!     If all eigenvalues are desired and ABSTOL is less than or equal
!     to zero, then call DSTERF or DOPGTR and SSTEQR.  If this fails
!     for some eigenvalue, then try DSTEBZ.
!
      if ( ( ALLEIG .OR. ( INDEIG .AND. IL == 1 .AND. IU.EQ.N ) ) .AND. &
          ( ABSTOL.LE.ZERO ) ) then
         CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         if ( .NOT.WANTZ ) then
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL DOPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ, &
                         WORK( INDWRK ), IINFO )
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ, &
                         WORK( INDWRK ), INFO )
            if ( INFO == 0 ) then
               DO 10 I = 1, N
                  IFAIL( I ) = 0
   10          CONTINUE
            end if
         end if
         if ( INFO == 0 ) then
            M = N
            GO TO 20
         end if
         INFO = 0
      end if
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.
!
      if ( WANTZ ) then
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      end if
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, &
                   WORK( INDD ), WORK( INDE ), M, NSPLIT, W, &
                   IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ), &
                   IWORK( INDIWO ), INFO )
!
      if ( WANTZ ) then
         CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), M, W, &
                      IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, &
                      WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
         CALL DOPMTR( 'L', UPLO, 'N', N, M, AP, WORK( INDTAU ), Z, LDZ, &
                      WORK( INDWRK ), INFO )
      end if
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
   20 CONTINUE
      if ( ISCALE == 1 ) then
         if ( INFO == 0 ) then
            IMAX = M
         ELSE
            IMAX = INFO - 1
         end if
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      end if
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
      if ( WANTZ ) then
         DO 40 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 30 JJ = J + 1, M
               if ( W( JJ ) < TMP1 ) then
                  I = JJ
                  TMP1 = W( JJ )
               end if
   30       CONTINUE
!
            if ( I.NE.0 ) then
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               if ( INFO.NE.0 ) then
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               end if
            end if
   40    CONTINUE
      end if
!
      RETURN
!
!     End of DSPEVX
!
      END
