      SUBROUTINE DGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, &
                         WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, &
                         IWORK, LIWORK, BWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVS, SENSE, SORT
      INTEGER            INFO, LDA, LDVS, LIWORK, LWORK, N, SDIM
      DOUBLE PRECISION   RCONDE, RCONDV
!     ..
!     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ), &
                         WR( * )
!     ..
!     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
!     ..
!
!  Purpose
!  =======
!
!  DGEESX computes for an N-by-N real nonsymmetric matrix A, the
!  eigenvalues, the real Schur form T, and, optionally, the matrix of
!  Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
!
!  Optionally, it also orders the eigenvalues on the diagonal of the
!  real Schur form so that selected eigenvalues are at the top left;
!  computes a reciprocal condition number for the average of the
!  selected eigenvalues (RCONDE); and computes a reciprocal condition
!  number for the right invariant subspace corresponding to the
!  selected eigenvalues (RCONDV).  The leading columns of Z form an
!  orthonormal basis for this invariant subspace.
!
!  For further explanation of the reciprocal condition numbers RCONDE
!  and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where
!  these quantities are called s and sep respectively).
!
!  A real matrix is in real Schur form if it is upper quasi-triangular
!  with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in
!  the form
!            [  a  b  ]
!            [  c  a  ]
!
!  where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
!
!  Arguments
!  =========
!
!  JOBVS   (input) CHARACTER*1
!          = 'N': Schur vectors are not computed;
!          = 'V': Schur vectors are computed.
!
!  SORT    (input) CHARACTER*1
!          Specifies whether or not to order the eigenvalues on the
!          diagonal of the Schur form.
!          = 'N': Eigenvalues are not ordered;
!          = 'S': Eigenvalues are ordered (see SELECT).
!
!  SELECT  (input) LOGICAL FUNCTION of two DOUBLE PRECISION arguments
!          SELECT must be declared EXTERNAL in the calling subroutine.
!          If SORT = 'S', SELECT is used to select eigenvalues to sort
!          to the top left of the Schur form.
!          If SORT = 'N', SELECT is not referenced.
!          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
!          SELECT(WR(j),WI(j)) is true; i.e., if either one of a
!          complex conjugate pair of eigenvalues is selected, then both
!          are.  Note that a selected complex eigenvalue may no longer
!          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
!          ordering may change the value of complex eigenvalues
!          (especially if the eigenvalue is ill-conditioned); in this
!          case INFO may be set to N+3 (see INFO below).
!
!  SENSE   (input) CHARACTER*1
!          Determines which reciprocal condition numbers are computed.
!          = 'N': None are computed;
!          = 'E': Computed for average of selected eigenvalues only;
!          = 'V': Computed for selected right invariant subspace only;
!          = 'B': Computed for both.
!          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'.
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the N-by-N matrix A.
!          On exit, A is overwritten by its real Schur form T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  SDIM    (output) INTEGER
!          If SORT = 'N', SDIM = 0.
!          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
!                         for which SELECT is true. (Complex conjugate
!                         pairs for which SELECT is true for either
!                         eigenvalue count as 2.)
!
!  WR      (output) DOUBLE PRECISION array, dimension (N)
!  WI      (output) DOUBLE PRECISION array, dimension (N)
!          WR and WI contain the real and imaginary parts, respectively,
!          of the computed eigenvalues, in the same order that they
!          appear on the diagonal of the output Schur form T.  Complex
!          conjugate pairs of eigenvalues appear consecutively with the
!          eigenvalue having the positive imaginary part first.
!
!  VS      (output) DOUBLE PRECISION array, dimension (LDVS,N)
!          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
!          vectors.
!          If JOBVS = 'N', VS is not referenced.
!
!  LDVS    (input) INTEGER
!          The leading dimension of the array VS.  LDVS >= 1, and if
!          JOBVS = 'V', LDVS >= N.
!
!  RCONDE  (output) DOUBLE PRECISION
!          If SENSE = 'E' or 'B', RCONDE contains the reciprocal
!          condition number for the average of the selected eigenvalues.
!          Not referenced if SENSE = 'N' or 'V'.
!
!  RCONDV  (output) DOUBLE PRECISION
!          If SENSE = 'V' or 'B', RCONDV contains the reciprocal
!          condition number for the selected right invariant subspace.
!          Not referenced if SENSE = 'N' or 'E'.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,3*N).
!          Also, if SENSE = 'E' or 'V' or 'B',
!          LWORK >= N+2*SDIM*(N-SDIM), where SDIM is the number of
!          selected eigenvalues computed by this routine.  Note that
!          N+2*SDIM*(N-SDIM) <= N+N*N/2.
!          For good performance, LWORK must generally be larger.
!
!  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
!          Not referenced if SENSE = 'N' or 'E'.
!          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          LIWORK >= 1; if SENSE = 'V' or 'B', LIWORK >= SDIM*(N-SDIM).
!
!  BWORK   (workspace) LOGICAL array, dimension (N)
!          Not referenced if SORT = 'N'.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, and i is
!             <= N: the QR algorithm failed to compute all the
!                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
!                   contain those eigenvalues which have converged; if
!                   JOBVS = 'V', VS contains the transformation which
!                   reduces A to its partially converged Schur form.
!             = N+1: the eigenvalues could not be reordered because some
!                   eigenvalues were too close to separate (the problem
!                   is very ill-conditioned);
!             = N+2: after reordering, roundoff changed values of some
!                   complex eigenvalues so that leading eigenvalues in
!                   the Schur form no longer satisfy SELECT=.TRUE.  This
!                   could also be caused by underflow due to scaling.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            CURSL, LASTSL, LST2SL, SCALEA, WANTSB, WANTSE, &
                         WANTSN, WANTST, WANTSV, WANTVS
      INTEGER            HSWORK, I, I1, I2, IBAL, ICOND, IERR, IEVAL, &
                         IHI, ILO, INXT, IP, ITAU, IWRK, K, MAXB, &
                         MAXWRK, MINWRK
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SMLNUM
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY, &
                         DLASCL, DORGHR, DSWAP, DTRSEN, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      WANTVS = LSAME( JOBVS, 'V' )
      WANTST = LSAME( SORT, 'S' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
      if ( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) then
         INFO = -1
      else if ( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) then
         INFO = -2
      else if ( .NOT.( WANTSN .OR. WANTSE .OR. WANTSV .OR. WANTSB ) .OR. &
               ( .NOT.WANTST .AND. .NOT.WANTSN ) ) then
         INFO = -4
      else if ( N < 0 ) then
         INFO = -5
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -7
      else if ( LDVS < 1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) then
         INFO = -12
      end if
!
!     Compute workspace
!      (Note: Comments in the code beginning "RWorkspace:" describe the
!       minimal amount of real workspace needed at that point in the
!       code, as well as the preferred amount for good performance.
!       IWorkspace refers to integer workspace.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by DHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.
!       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
!       depends on SDIM, which is computed by the routine DTRSEN later
!       in the code.)
!
      MINWRK = 1
      if ( INFO == 0 .AND. LWORK.GE.1 ) then
         MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
         MINWRK = MAX( 1, 3*N )
         if ( .NOT.WANTVS ) then
            MAXB = MAX( ILAENV( 8, 'DHSEQR', 'SN', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'SN', N, 1, &
                N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, N+HSWORK, 1 )
         ELSE
            MAXWRK = MAX( MAXWRK, 2*N+( N-1 )* &
                     ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) )
            MAXB = MAX( ILAENV( 8, 'DHSEQR', 'SV', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'SV', N, 1, &
                N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, N+HSWORK, 1 )
         end if
         WORK( 1 ) = MAXWRK
      end if
      if ( LWORK < MINWRK ) then
         INFO = -16
      end if
      if ( LIWORK < 1 ) then
         INFO = -18
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEESX', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) then
         SDIM = 0
         RETURN
      end if
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM < SMLNUM ) then
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      else if ( ANRM.GT.BIGNUM ) then
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      end if
      if ( SCALEA ) &
         CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!     Permute the matrix to make it more nearly triangular
!     (RWorkspace: need N)
!
      IBAL = 1
      CALL DGEBAL( 'P', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
!
!     Reduce to upper Hessenberg form
!     (RWorkspace: need 3*N, prefer 2*N+N*NB)
!
      ITAU = N + IBAL
      IWRK = N + ITAU
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
                   LWORK-IWRK+1, IERR )
!
      if ( WANTVS ) then
!
!        Copy Householder vectors to VS
!
         CALL DLACPY( 'L', N, N, A, LDA, VS, LDVS )
!
!        Generate orthogonal matrix in VS
!        (RWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL DORGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
      end if
!
      SDIM = 0
!
!     Perform QR iteration, accumulating Schur vectors in VS if desired
!     (RWorkspace: need N+1, prefer N+HSWORK (see comments) )
!
      IWRK = ITAU
      CALL DHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, WR, WI, VS, LDVS, &
                   WORK( IWRK ), LWORK-IWRK+1, IEVAL )
      if ( IEVAL.GT.0 ) &
         INFO = IEVAL
!
!     Sort eigenvalues if desired
!
      if ( WANTST .AND. INFO == 0 ) then
         if ( SCALEA ) then
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WR, N, IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WI, N, IERR )
         end if
         DO I = 1, N
            BWORK( I ) = SELECT( WR( I ), WI( I ) )
         end do
!
!        Reorder eigenvalues, transform Schur vectors, and compute
!        reciprocal condition numbers
!        (RWorkspace: if SENSE is not 'N', need N+2*SDIM*(N-SDIM)
!                     otherwise, need N )
!        (IWorkspace: if SENSE is 'V' or 'B', need SDIM*(N-SDIM)
!                     otherwise, need 0 )
!
         CALL DTRSEN( SENSE, JOBVS, BWORK, N, A, LDA, VS, LDVS, WR, WI, &
                      SDIM, RCONDE, RCONDV, WORK( IWRK ), LWORK-IWRK+1, &
                      IWORK, LIWORK, ICOND )
         if ( .NOT.WANTSN ) &
            MAXWRK = MAX( MAXWRK, N+2*SDIM*( N-SDIM ) )
         if ( ICOND == -15 ) then
!
!           Not enough real workspace
!
            INFO = -16
         else if ( ICOND == -17 ) then
!
!           Not enough integer workspace
!
            INFO = -18
         else if ( ICOND.GT.0 ) then
!
!           DTRSEN failed to reorder or to restore standard Schur form
!
            INFO = ICOND + N
         end if
      end if
!
      if ( WANTVS ) then
!
!        Undo balancing
!        (RWorkspace: need N)
!
         CALL DGEBAK( 'P', 'R', N, ILO, IHI, WORK( IBAL ), N, VS, LDVS, &
                      IERR )
      end if
!
      if ( SCALEA ) then
!
!        Undo scaling for the Schur form of A
!
         CALL DLASCL( 'H', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL DCOPY( N, A, LDA+1, WR, 1 )
         if ( ( WANTSV .OR. WANTSB ) .AND. INFO == 0 ) then
            DUM( 1 ) = RCONDV
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
            RCONDV = DUM( 1 )
         end if
         if ( CSCALE == SMLNUM ) then
!
!           If scaling back towards underflow, adjust WI if an
!           offdiagonal element of a 2-by-2 block in the Schur form
!           underflows.
!
            if ( IEVAL.GT.0 ) then
               I1 = IEVAL + 1
               I2 = IHI - 1
               CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, &
                            IERR )
            else if ( WANTST ) then
               I1 = 1
               I2 = N - 1
            ELSE
               I1 = ILO
               I2 = IHI - 1
            end if
            INXT = I1 - 1
            DO 20 I = I1, I2
               if ( I < INXT ) &
                  GO TO 20
               if ( WI( I ) == ZERO ) then
                  INXT = I + 1
               ELSE
                  if ( A( I+1, I ) == ZERO ) then
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                  else if ( A( I+1, I ).NE.ZERO .AND. A( I, I+1 ) ==  &
                           ZERO ) then
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                     if ( I.GT.1 ) &
                        CALL DSWAP( I-1, A( 1, I ), 1, A( 1, I+1 ), 1 )
                     if ( N.GT.I+1 ) &
                        CALL DSWAP( N-I-1, A( I, I+2 ), LDA, &
                                    A( I+1, I+2 ), LDA )
                     CALL DSWAP( N, VS( 1, I ), 1, VS( 1, I+1 ), 1 )
                     A( I, I+1 ) = A( I+1, I )
                     A( I+1, I ) = ZERO
                  end if
                  INXT = I + 2
               end if
   20       CONTINUE
         end if
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-IEVAL, 1, &
                      WI( IEVAL+1 ), MAX( N-IEVAL, 1 ), IERR )
      end if
!
      if ( WANTST .AND. INFO == 0 ) then
!
!        Check if reordering successful
!
         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 30 I = 1, N
            CURSL = SELECT( WR( I ), WI( I ) )
            if ( WI( I ) == ZERO ) then
               if ( CURSL ) &
                  SDIM = SDIM + 1
               IP = 0
               if ( CURSL .AND. .NOT.LASTSL ) &
                  INFO = N + 2
            ELSE
               if ( IP == 1 ) then
!
!                 Last eigenvalue of conjugate pair
!
                  CURSL = CURSL .OR. LASTSL
                  LASTSL = CURSL
                  if ( CURSL ) &
                     SDIM = SDIM + 2
                  IP = -1
                  if ( CURSL .AND. .NOT.LST2SL ) &
                     INFO = N + 2
               ELSE
!
!                 First eigenvalue of conjugate pair
!
                  IP = 1
               end if
            end if
            LST2SL = LASTSL
            LASTSL = CURSL
   30    CONTINUE
      end if
!
      WORK( 1 ) = MAXWRK
      if ( WANTSV .OR. WANTSB ) then
         IWORK( 1 ) = SDIM*( N-SDIM )
      ELSE
         IWORK( 1 ) = 1
      end if
!
      RETURN
!
!     End of DGEESX
!
      END
