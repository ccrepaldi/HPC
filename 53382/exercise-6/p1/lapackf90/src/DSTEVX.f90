      SUBROUTINE DSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, &
                         M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            IL, INFO, IU, LDZ, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEVX computes selected eigenvalues and, optionally, eigenvectors
!  of a real symmetric tridiagonal matrix A.  Eigenvalues and
!  eigenvectors can be selected by specifying either a range of values
!  or a range of indices for the desired eigenvalues.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  RANGE   (input) CHARACTER*1
!          = 'A': all eigenvalues will be found.
!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!                 will be found.
!          = 'I': the IL-th through IU-th eigenvalues will be found.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix
!          A.
!          On exit, D may be multiplied by a constant factor chosen
!          to avoid over/underflow in computing the eigenvalues.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix A in elements 1 to N-1 of E; E(N) need not be set.
!          On exit, E may be multiplied by a constant factor chosen
!          to avoid over/underflow in computing the eigenvalues.
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
!          where EPS is the machine precision.  If ABSTOL is less
!          than or equal to zero, then  EPS*|T|  will be used in
!          its place, where |T| is the 1-norm of the tridiagonal
!          matrix.
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
!          The first M elements contain the selected eigenvalues in
!          ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!          contain the orthonormal eigenvectors of the matrix A
!          corresponding to the selected eigenvalues, with the i-th
!          column of Z holding the eigenvector associated with W(i).
!          If an eigenvector fails to converge (INFO > 0), then that
!          column of Z contains the latest approximation to the
!          eigenvector, and the index of the eigenvector is returned
!          in IFAIL.  If JOBZ = 'N', then Z is not referenced.
!          Note: the user must ensure that at least max(1,M) columns are
!          supplied in the array Z; if RANGE = 'V', the exact value of M
!          is not known in advance and an upper bound must be used.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
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
      INTEGER            I, IMAX, INDIBL, INDISP, INDIWO, INDWRK, &
                         ISCALE, ITMP1, J, JJ, NSPLIT
      DOUBLE PRECISION   BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, &
                         TMP1, TNRM, VLL, VUU
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DSCAL, DSTEBZ, DSTEIN, DSTEQR, DSTERF, &
                         DSWAP, XERBLA
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
      else if ( N < 0 ) then
         INFO = -3
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
         CALL XERBLA( 'DSTEVX', -INFO )
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
            W( 1 ) = D( 1 )
         ELSE
            if ( VL < D( 1 ) .AND. VU.GE.D( 1 ) ) then
               M = 1
               W( 1 ) = D( 1 )
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
      if ( VALEIG ) then
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      end if
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
         if ( VALEIG ) then
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         end if
      end if
!
!     If all eigenvalues are desired and ABSTOL is less than zero, then
!     call DSTERF or SSTEQR.  If this fails for some eigenvalue, then
!     try DSTEBZ.
!
      if ( ( ALLEIG .OR. ( INDEIG .AND. IL == 1 .AND. IU.EQ.N ) ) .AND. &
          ( ABSTOL.LE.ZERO ) ) then
         CALL DCOPY( N, D, 1, W, 1 )
         CALL DCOPY( N-1, E( 1 ), 1, WORK( 1 ), 1 )
         INDWRK = N + 1
         if ( .NOT.WANTZ ) then
            CALL DSTERF( N, W, WORK, INFO )
         ELSE
            CALL DSTEQR( 'I', N, W, WORK, Z, LDZ, WORK( INDWRK ), INFO )
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
      INDWRK = 1
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E, M, &
                   NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), &
                   WORK( INDWRK ), IWORK( INDIWO ), INFO )
!
      if ( WANTZ ) then
         CALL DSTEIN( N, D, E, M, W, IWORK( INDIBL ), IWORK( INDISP ), &
                      Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, &
                      INFO )
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
!     End of DSTEVX
!
      END
