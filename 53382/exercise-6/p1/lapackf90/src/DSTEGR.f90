      SUBROUTINE DSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, &
                         M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, &
                         LIWORK, INFO )
!
!  -- LAPACK computational routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
! DSTEGR computes selected eigenvalues and, optionally, eigenvectors
! of a real symmetric tridiagonal matrix T.  Eigenvalues and
! eigenvectors can be selected by specifying either a range of values
! or a range of indices for the desired eigenvalues. The eigenvalues
! are computed by the dqds algorithm, while orthogonal eigenvectors are
! computed from various ``good'' L D L^T representations (also known as
! Relatively Robust Representations). Gram-Schmidt orthogonalization is
! avoided as far as possible. More specifically, the various steps of
! the algorithm are as follows. For the i-th unreduced block of T,
!     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
!         is a relatively robust representation,
!     (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
!         relative accuracy by the dqds algorithm,
!     (c) If there is a cluster of close eigenvalues, "choose" sigma_i
!         close to the cluster, and go to step (a),
!     (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
!         compute the corresponding eigenvector by forming a
!         rank-revealing twisted factorization.
!  The desired accuracy of the output can be specified by the input
!  parameter ABSTOL.
!
!  For more details, see "A new O(n^2) algorithm for the symmetric
!  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
!  Computer Science Division Technical Report No. UCB/CSD-97-971,
!  UC Berkeley, May 1997.
!
!  Note 1 : Currently DSTEGR is only set up to find ALL the n
!  eigenvalues and eigenvectors of T in O(n^2) time
!  Note 2 : Currently the routine DSTEIN is called when an appropriate
!  sigma_i cannot be chosen in step (c) above. DSTEIN invokes modified
!  Gram-Schmidt when eigenvalues are close.
!  Note 3 : DSTEGR works only on machines which follow ieee-754
!  floating-point standard in their handling of infinities and NaNs.
!  Normal execution of DSTEGR may create NaNs and infinities and hence
!  may abort due to a floating point exception in environments which
!  do not conform to the ieee standard.
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
!********* Only RANGE = 'A' is currently supported *********************
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix
!          T. On exit, D is overwritten.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix T in elements 1 to N-1 of E; E(N) need not be set.
!          On exit, E is overwritten.
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
!          The absolute error tolerance for the
!          eigenvalues/eigenvectors. IF JOBZ = 'V', the eigenvalues and
!          eigenvectors output have residual norms bounded by ABSTOL,
!          and the dot products between different eigenvectors are
!          bounded by ABSTOL. If ABSTOL is less than N*EPS*|T|, then
!          N*EPS*|T| will be used in its place, where EPS is the
!          machine precision and |T| is the 1-norm of the tridiagonal
!          matrix. The eigenvalues are computed to an accuracy of
!          EPS*|T| irrespective of ABSTOL. If high relative accuracy
!          is important, set ABSTOL to DLAMCH( 'Safe minimum' ).
!          See Barlow and Demmel "Computing Accurate Eigensystems of
!          Scaled Diagonally Dominant Matrices", LAPACK Working Note #7
!          for a discussion of which matrices define their eigenvalues
!          to high relative accuracy.
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
!          contain the orthonormal eigenvectors of the matrix T
!          corresponding to the selected eigenvalues, with the i-th
!          column of Z holding the eigenvector associated with W(i).
!          If JOBZ = 'N', then Z is not referenced.
!          Note: the user must ensure that at least max(1,M) columns are
!          supplied in the array Z; if RANGE = 'V', the exact value of M
!          is not known in advance and an upper bound must be used.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )
!          The support of the eigenvectors in Z, i.e., the indices
!          indicating the nonzero elements in Z. The i-th eigenvector
!          is nonzero only in elements ISUPPZ( 2*i-1 ) through
!          ISUPPZ( 2*i ).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal
!          (and minimal) LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,18*N)
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
!          The dimension of the array IWORK.  LIWORK >= max(1,10*N)
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = 1, internal error in DLARRE,
!                if INFO = 2, internal error in DLARRV.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Inderjit Dhillon, IBM Almaden, USA
!     Osni Marques, LBNL/NERSC, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ
      INTEGER            I, IBEGIN, IEND, IINDBL, IINDWK, IINFO, IINSPL, &
                         INDGRS, INDWOF, INDWRK, ITMP, J, JJ, LIWMIN, &
                         LWMIN, NSPLIT
      DOUBLE PRECISION   BIGNUM, EPS, RMAX, RMIN, SAFMIN, SCALE, SMLNUM, &
                         THRESH, TMP, TNRM, TOL
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARRE, DLARRV, DLASET, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
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
      LQUERY = ( ( LWORK == -1 ) .OR. ( LIWORK.EQ.-1 ) )
      LWMIN = 18*N
      LIWMIN = 10*N
!
      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) then
         INFO = -2
!
!     The following two lines need to be removed once the
!     RANGE = 'V' and RANGE = 'I' options are provided.
!
      else if ( VALEIG .OR. INDEIG ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) then
         INFO = -7
      else if ( INDEIG .AND. IL < 1 ) then
         INFO = -8
!     The following change should be made in DSTEVX also, otherwise
!     IL can be specified as N+1 and IU as N.
!     else if ( INDEIG .AND. ( IU < MIN( N, IL ) .OR. IU.GT.N ) ) then
      else if ( INDEIG .AND. ( IU < IL .OR. IU.GT.N ) ) then
         INFO = -9
      else if ( LDZ < 1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) then
         INFO = -14
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -17
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -19
      end if
      if ( INFO == 0 ) then
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSTEGR', -INFO )
         RETURN
      else if ( LQUERY ) then
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
      SCALE = ONE
      TNRM = DLANST( 'M', N, D, E )
      if ( TNRM.GT.ZERO .AND. TNRM < RMIN ) then
         SCALE = RMIN / TNRM
      else if ( TNRM.GT.RMAX ) then
         SCALE = RMAX / TNRM
      end if
      if ( SCALE.NE.ONE ) then
         CALL DSCAL( N, SCALE, D, 1 )
         CALL DSCAL( N-1, SCALE, E, 1 )
         TNRM = TNRM*SCALE
      end if
      INDGRS = 1
      INDWOF = 2*N + 1
      INDWRK = 3*N + 1
!
      IINSPL = 1
      IINDBL = N + 1
      IINDWK = 2*N + 1
!
      CALL DLASET( 'Full', N, N, ZERO, ZERO, Z, LDZ )
!
!     Compute the desired eigenvalues of the tridiagonal after splitting
!     into smaller subblocks if the corresponding of-diagonal elements
!     are small
!
      THRESH = EPS*TNRM
      CALL DLARRE( N, D, E, THRESH, NSPLIT, IWORK( IINSPL ), M, W, &
                   WORK( INDWOF ), WORK( INDGRS ), WORK( INDWRK ), &
                   IINFO )
      if ( IINFO.NE.0 ) then
         INFO = 1
         RETURN
      end if
!
      if ( WANTZ ) then
!
!        Compute the desired eigenvectors corresponding to the computed
!        eigenvalues
!
         TOL = MAX( ABSTOL, DBLE( N )*THRESH )
         IBEGIN = 1
         DO 20 I = 1, NSPLIT
            IEND = IWORK( IINSPL+I-1 )
            DO 10 J = IBEGIN, IEND
               IWORK( IINDBL+J-1 ) = I
   10       CONTINUE
            IBEGIN = IEND + 1
   20    CONTINUE
!
         CALL DLARRV( N, D, E, IWORK( IINSPL ), M, W, IWORK( IINDBL ), &
                      WORK( INDGRS ), TOL, Z, LDZ, ISUPPZ, &
                      WORK( INDWRK ), IWORK( IINDWK ), IINFO )
         if ( IINFO.NE.0 ) then
            INFO = 2
            RETURN
         end if
!
      end if
!
      IBEGIN = 1
      DO 40 I = 1, NSPLIT
         IEND = IWORK( IINSPL+I-1 )
         DO 30 J = IBEGIN, IEND
            W( J ) = W( J ) + WORK( INDWOF+I-1 )
   30    CONTINUE
         IBEGIN = IEND + 1
   40 CONTINUE
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      if ( SCALE.NE.ONE ) then
         CALL DSCAL( M, ONE / SCALE, W, 1 )
      end if
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
      if ( NSPLIT.GT.1 ) then
         DO 60 J = 1, M - 1
            I = 0
            TMP = W( J )
            DO 50 JJ = J + 1, M
               if ( W( JJ ) < TMP ) then
                  I = JJ
                  TMP = W( JJ )
               end if
   50       CONTINUE
            if ( I.NE.0 ) then
               W( I ) = W( J )
               W( J ) = TMP
               if ( WANTZ ) then
                  CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
                  ITMP = ISUPPZ( 2*I-1 )
                  ISUPPZ( 2*I-1 ) = ISUPPZ( 2*J-1 )
                  ISUPPZ( 2*J-1 ) = ITMP
                  ITMP = ISUPPZ( 2*I )
                  ISUPPZ( 2*I ) = ISUPPZ( 2*J )
                  ISUPPZ( 2*J ) = ITMP
               end if
            end if
   60    CONTINUE
      end if
!
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
!
!     End of DSTEGR
!
      END
