      SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
                         ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, &
                         IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 20, 2000
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYEVR computes selected eigenvalues and, optionally, eigenvectors
!  of a real symmetric matrix T.  Eigenvalues and eigenvectors can be
!  selected by specifying either a range of values or a range of
!  indices for the desired eigenvalues.
!
!  Whenever possible, DSYEVR calls DSTEGR to compute the
!  eigenspectrum using Relatively Robust Representations.  DSTEGR
!  computes eigenvalues by the dqds algorithm, while orthogonal
!  eigenvectors are computed from various "good" L D L^T representations
!  (also known as Relatively Robust Representations). Gram-Schmidt
!  orthogonalization is avoided as far as possible. More specifically,
!  the various steps of the algorithm are as follows. For the i-th
!  unreduced block of T,
!     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
!          is a relatively robust representation,
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
!  Computer Science Division Technical Report No. UCB//CSD-97-971,
!  UC Berkeley, May 1997.
!
!
!  Note 1 : DSYEVR calls DSTEGR when the full spectrum is requested
!  on machines which conform to the ieee-754 floating point standard.
!  DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and
!  when partial spectrum requests are made.
!
!  Normal execution of DSTEGR may create NaNs and infinities and
!  hence may abort due to a floating point exception in environments
!  which do not handle NaNs and infinities in the ieee standard default
!  manner.
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
!********* For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
!********* DSTEIN are called
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
!          On exit, the lower triangle (if UPLO='L') or the upper
!          triangle (if UPLO='U') of A, including the diagonal, is
!          destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
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
!          by reducing A to tridiagonal form.
!
!          See "Computing Small Singular Values of Bidiagonal Matrices
!          with Guaranteed High Relative Accuracy," by Demmel and
!          Kahan, LAPACK Working Note #3.
!
!          If high relative accuracy is important, set ABSTOL to
!          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
!          eigenvalues are computed to high relative accuracy when
!          possible in future releases.  The current code does not
!          make any guarantees about high relative accuracy, but
!          furutre releases will. See J. Barlow and J. Demmel,
!          "Computing Accurate Eigensystems of Scaled Diagonally
!          Dominant Matrices", LAPACK Working Note #7, for a discussion
!          of which matrices define their eigenvalues to high relative
!          accuracy.
!
!  M       (output) INTEGER
!          The total number of eigenvalues found.  0 <= M <= N.
!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          The first M elements contain the selected eigenvalues in
!          ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!          contain the orthonormal eigenvectors of the matrix A
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
!  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
!          The support of the eigenvectors in Z, i.e., the indices
!          indicating the nonzero elements in Z. The i-th eigenvector
!          is nonzero only in elements ISUPPZ( 2*i-1 ) through
!          ISUPPZ( 2*i ).
!********* Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,26*N).
!          For optimal efficiency, LWORK >= (NB+6)*N,
!          where NB is the max of the blocksize for DSYTRD and DORMTR
!          returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
!          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  Internal error
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Inderjit Dhillon, IBM Almaden, USA
!     Osni Marques, LBNL/NERSC, USA
!     Ken Stanley, Computer Science Division, University of
!       California at Berkeley, USA
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, LQUERY, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IEEEOK, IINFO, IMAX, INDD, INDDD, INDE, &
                         INDEE, INDIBL, INDIFL, INDISP, INDIWO, INDTAU, &
                         INDWK, INDWKN, ISCALE, ITMP1, J, JJ, LIWMIN, &
                         LLWORK, LLWRKN, LWKOPT, LWMIN, NB, NSPLIT
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, &
                         SIGMA, SMLNUM, TMP1, VLL, VUU
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DORMTR, DSCAL, DSTEBZ, DSTEGR, DSTEIN, &
                         DSTERF, DSWAP, DSYTRD, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      IEEEOK = ILAENV( 10, 'DSYEVR', 'N', 1, 2, 3, 4 )
!
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
!
      LQUERY = ( ( LWORK == -1 ) .OR. ( LIWORK.EQ.-1 ) )
!
      LWMIN = MAX( 1, 26*N )
      LIWMIN = MAX( 1, 10*N )
!
      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -1
      else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) then
         INFO = -2
      else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) then
         INFO = -3
      else if ( N < 0 ) then
         INFO = -4
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -6
      ELSE
         if ( VALEIG ) then
            if ( N.GT.0 .AND. VU.LE.VL ) &
               INFO = -8
         else if ( INDEIG ) then
            if ( IL < 1 .OR. IL.GT.MAX( 1, N ) ) then
               INFO = -9
            else if ( IU < MIN( N, IL ) .OR. IU.GT.N ) then
               INFO = -10
            end if
         end if
      end if
      if ( INFO == 0 ) then
         if ( LDZ < 1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) then
            INFO = -15
         else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
            INFO = -18
         else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
            INFO = -20
         end if
      end if
!
      if ( INFO == 0 ) then
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         NB = MAX( NB, ILAENV( 1, 'ZUNMTR', UPLO, N, -1, -1, -1 ) )
         LWKOPT = MAX( ( NB+1 )*N, LWMIN )
         WORK( 1 ) = LWKOPT
         IWORK( 1 ) = LIWMIN
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSYEVR', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      M = 0
      if ( N == 0 ) then
         WORK( 1 ) = 1
         RETURN
      end if
!
      if ( N == 1 ) then
         WORK( 1 ) = 7
         if ( ALLEIG .OR. INDEIG ) then
            M = 1
            W( 1 ) = A( 1, 1 )
         ELSE
            if ( VL < A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) then
               M = 1
               W( 1 ) = A( 1, 1 )
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
      VLL = VL
      VUU = VU
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      if ( ANRM.GT.ZERO .AND. ANRM < RMIN ) then
         ISCALE = 1
         SIGMA = RMIN / ANRM
      else if ( ANRM.GT.RMAX ) then
         ISCALE = 1
         SIGMA = RMAX / ANRM
      end if
      if ( ISCALE == 1 ) then
         if ( LOWER ) then
            DO 10 J = 1, N
               CALL DSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL DSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         end if
         if ( ABSTOL.GT.0 ) &
            ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) then
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         end if
      end if
!
!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
!
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDEE = INDD + N
      INDDD = INDEE + N
      INDIFL = INDDD + N
      INDWK = INDIFL + N
      LLWORK = LWORK - INDWK + 1
      CALL DSYTRD( UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), &
                   WORK( INDTAU ), WORK( INDWK ), LLWORK, IINFO )
!
!     If all eigenvalues are desired
!     then call DSTERF or SSTEGR and DORMTR.
!
      if ( ( ALLEIG .OR. ( INDEIG .AND. IL == 1 .AND. IU.EQ.N ) ) .AND. &
          IEEEOK == 1 ) then
         if ( .NOT.WANTZ ) then
            CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DCOPY( N, WORK( INDD ), 1, WORK( INDDD ), 1 )
!
            CALL DSTEGR( JOBZ, 'A', N, WORK( INDDD ), WORK( INDEE ), &
                         VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, &
                         WORK( INDWK ), LWORK, IWORK, LIWORK, INFO )
!
!
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
            if ( WANTZ .AND. INFO == 0 ) then
               INDWKN = INDE
               LLWRKN = LWORK - INDWKN + 1
               CALL DORMTR( 'L', UPLO, 'N', N, M, A, LDA, &
                            WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), &
                            LLWRKN, IINFO )
            end if
         end if
!
!
         if ( INFO == 0 ) then
            M = N
            GO TO 30
         end if
         INFO = 0
      end if
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.
!     Also call DSTEBZ and SSTEIN if SSTEGR fails.
!
      if ( WANTZ ) then
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      end if
      INDIFL = 1
      INDIBL = INDIFL + N
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, &
                   WORK( INDD ), WORK( INDE ), M, NSPLIT, W, &
                   IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWK ), &
                   IWORK( INDIWO ), INFO )
!
      if ( WANTZ ) then
         CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), M, W, &
                      IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, &
                      WORK( INDWK ), IWORK( INDIWO ), IWORK( INDIFL ), &
                      INFO )
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         CALL DORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, &
                      LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      end if
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
   30 CONTINUE
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
         DO 50 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, M
               if ( W( JJ ) < TMP1 ) then
                  I = JJ
                  TMP1 = W( JJ )
               end if
   40       CONTINUE
!
            if ( I.NE.0 ) then
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            end if
   50    CONTINUE
      end if
!
!     Set WORK(1) to optimal workspace size.
!
      WORK( 1 ) = LWKOPT
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of DSYEVR
!
      END
