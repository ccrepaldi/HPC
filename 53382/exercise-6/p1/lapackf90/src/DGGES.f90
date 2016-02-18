      SUBROUTINE DGGES( JOBVSL, JOBVSR, SORT, DELCTG, N, A, LDA, B, LDB, &
                        SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, &
                        LDVSR, WORK, LWORK, BWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVSL, JOBVSR, SORT
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
!     ..
!     .. Array Arguments ..
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                         B( LDB, * ), BETA( * ), VSL( LDVSL, * ), &
                         VSR( LDVSR, * ), WORK( * )
!     ..
!     .. Function Arguments ..
      LOGICAL            DELCTG
      EXTERNAL           DELCTG
!     ..
!
!  Purpose
!  =======
!
!  DGGES computes for a pair of N-by-N real nonsymmetric matrices (A,B),
!  the generalized eigenvalues, the generalized real Schur form (S,T),
!  optionally, the left and/or right matrices of Schur vectors (VSL and
!  VSR). This gives the generalized Schur factorization
!
!           (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
!
!  Optionally, it also orders the eigenvalues so that a selected cluster
!  of eigenvalues appears in the leading diagonal blocks of the upper
!  quasi-triangular matrix S and the upper triangular matrix T.The
!  leading columns of VSL and VSR then form an orthonormal basis for the
!  corresponding left and right eigenspaces (deflating subspaces).
!
!  (If only the generalized eigenvalues are needed, use the driver
!  DGGEV instead, which is faster.)
!
!  A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
!  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
!  usually represented as the pair (alpha,beta), as there is a
!  reasonable interpretation for beta=0 or both being zero.
!
!  A pair of matrices (S,T) is in generalized real Schur form if T is
!  upper triangular with non-negative diagonal and S is block upper
!  triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
!  to real generalized eigenvalues, while 2-by-2 blocks of S will be
!  "standardized" by making the corresponding elements of T have the
!  form:
!          [  a  0  ]
!          [  0  b  ]
!
!  and the pair of corresponding 2-by-2 blocks in S and T will have a
!  complex conjugate pair of generalized eigenvalues.
!
!
!  Arguments
!  =========
!
!  JOBVSL  (input) CHARACTER*1
!          = 'N':  do not compute the left Schur vectors;
!          = 'V':  compute the left Schur vectors.
!
!  JOBVSR  (input) CHARACTER*1
!          = 'N':  do not compute the right Schur vectors;
!          = 'V':  compute the right Schur vectors.
!
!  SORT    (input) CHARACTER*1
!          Specifies whether or not to order the eigenvalues on the
!          diagonal of the generalized Schur form.
!          = 'N':  Eigenvalues are not ordered;
!          = 'S':  Eigenvalues are ordered (see DELZTG);
!
!  DELZTG  (input) LOGICAL FUNCTION of three DOUBLE PRECISION arguments
!          DELZTG must be declared EXTERNAL in the calling subroutine.
!          If SORT = 'N', DELZTG is not referenced.
!          If SORT = 'S', DELZTG is used to select eigenvalues to sort
!          to the top left of the Schur form.
!          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
!          DELZTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
!          one of a complex conjugate pair of eigenvalues is selected,
!          then both complex eigenvalues are selected.
!
!          Note that in the ill-conditioned case, a selected complex
!          eigenvalue may no longer satisfy DELZTG(ALPHAR(j),ALPHAI(j),
!          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2
!          in this case.
!
!  N       (input) INTEGER
!          The order of the matrices A, B, VSL, and VSR.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the first of the pair of matrices.
!          On exit, A has been overwritten by its generalized Schur
!          form S.
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
!          On entry, the second of the pair of matrices.
!          On exit, B has been overwritten by its generalized Schur
!          form T.
!
!  LDB     (input) INTEGER
!          The leading dimension of B.  LDB >= max(1,N).
!
!  SDIM    (output) INTEGER
!          If SORT = 'N', SDIM = 0.
!          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
!          for which DELZTG is true.  (Complex conjugate pairs for which
!          DELZTG is true for either eigenvalue count as 2.)
!
!  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
!  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
!  BETA    (output) DOUBLE PRECISION array, dimension (N)
!          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
!          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
!          and  BETA(j),j=1,...,N are the diagonals of the complex Schur
!          form (S,T) that would result if the 2-by-2 diagonal blocks of
!          the real Schur form of (A,B) were further reduced to
!          triangular form using 2-by-2 complex unitary transformations.
!          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
!          positive, then the j-th and (j+1)-st eigenvalues are a
!          complex conjugate pair, with ALPHAI(j+1) negative.
!
!          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
!          may easily over- or underflow, and BETA(j) may even be zero.
!          Thus, the user should avoid naively computing the ratio.
!          However, ALPHAR and ALPHAI will be always less than and
!          usually comparable with norm(A) in magnitude, and BETA always
!          less than and usually comparable with norm(B).
!
!  VSL     (output) DOUBLE PRECISION array, dimension (LDVSL,N)
!          If JOBVSL = 'V', VSL will contain the left Schur vectors.
!          Not referenced if JOBVSL = 'N'.
!
!  LDVSL   (input) INTEGER
!          The leading dimension of the matrix VSL. LDVSL >=1, and
!          if JOBVSL = 'V', LDVSL >= N.
!
!  VSR     (output) DOUBLE PRECISION array, dimension (LDVSR,N)
!          If JOBVSR = 'V', VSR will contain the right Schur vectors.
!          Not referenced if JOBVSR = 'N'.
!
!  LDVSR   (input) INTEGER
!          The leading dimension of the matrix VSR. LDVSR >= 1, and
!          if JOBVSR = 'V', LDVSR >= N.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= 8*N+16.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  BWORK   (workspace) LOGICAL array, dimension (N)
!          Not referenced if SORT = 'N'.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          = 1,...,N:
!                The QZ iteration failed.  (A,B) are not in Schur
!                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
!                be correct for j=INFO+1,...,N.
!          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
!                =N+2: after reordering, roundoff changed values of
!                      some complex eigenvalues so that leading
!                      eigenvalues in the Generalized Schur form no
!                      longer satisfy DELZTG=.TRUE.  This could also
!                      be caused due to scaling.
!                =N+3: reordering failed in DTGSEN.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, &
                         LQUERY, LST2SL, WANTST
      INTEGER            I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, &
                         ILO, IP, IRIGHT, IROWS, ITAU, IWRK, MAXWRK, &
                         MINWRK
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL, &
                         PVSR, SAFMAX, SAFMIN, SMLNUM
!     ..
!     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      DOUBLE PRECISION   Dif ( 2 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLABAD, &
                         DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGSEN, &
                         XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Decode the input arguments
!
      if ( LSAME( JOBVSL, 'N' ) ) then
         IJOBVL = 1
         ILVSL = .FALSE.
      else if ( LSAME( JOBVSL, 'V' ) ) then
         IJOBVL = 2
         ILVSL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVSL = .FALSE.
      end if
!
      if ( LSAME( JOBVSR, 'N' ) ) then
         IJOBVR = 1
         ILVSR = .FALSE.
      else if ( LSAME( JOBVSR, 'V' ) ) then
         IJOBVR = 2
         ILVSR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVSR = .FALSE.
      end if
!
      WANTST = LSAME( SORT, 'S' )
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( IJOBVL.LE.0 ) then
         INFO = -1
      else if ( IJOBVR.LE.0 ) then
         INFO = -2
      else if ( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) then
         INFO = -3
      else if ( N < 0 ) then
         INFO = -5
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -7
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -9
      else if ( LDVSL < 1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) then
         INFO = -15
      else if ( LDVSR < 1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) then
         INFO = -17
      end if
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.)
!
      MINWRK = 1
      if ( INFO == 0 .AND. ( LWORK.GE.1 .OR. LQUERY ) ) then
         MINWRK = 7*( N+1 ) + 16
         MAXWRK = 7*( N+1 ) + N*ILAENV( 1, 'DGEQRF', ' ', N, 1, N, 0 ) + &
                  16
         if ( ILVSL ) then
            MAXWRK = MAX( MAXWRK, 7*( N+1 )+N* &
                     ILAENV( 1, 'DORGQR', ' ', N, 1, N, -1 ) )
         end if
         WORK( 1 ) = MAXWRK
      end if
!
      if ( LWORK < MINWRK .AND. .NOT.LQUERY ) &
         INFO = -19
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGGES ', -INFO )
         RETURN
      else if ( LQUERY ) then
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
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      SMLNUM = SQRT( SAFMIN ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM < SMLNUM ) then
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      else if ( ANRM.GT.BIGNUM ) then
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      end if
      if ( ILASCL ) &
         CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = .FALSE.
      if ( BNRM.GT.ZERO .AND. BNRM < SMLNUM ) then
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      else if ( BNRM.GT.BIGNUM ) then
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      end if
      if ( ILBSCL ) &
         CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
!
!     Permute the matrix to make it more nearly triangular
!     (Workspace: need 6*N + 2*N space for storing balancing factors)
!
      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      CALL DGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), &
                   WORK( IRIGHT ), WORK( IWRK ), IERR )
!
!     Reduce B to triangular form (QR decomposition of B)
!     (Workspace: need N, prefer N*NB)
!
      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWRK
      IWRK = ITAU + IROWS
      CALL DGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), &
                   WORK( IWRK ), LWORK+1-IWRK, IERR )
!
!     Apply the orthogonal transformation to matrix A
!     (Workspace: need N, prefer N*NB)
!
      CALL DORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, &
                   WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), &
                   LWORK+1-IWRK, IERR )
!
!     Initialize VSL
!     (Workspace: need N, prefer N*NB)
!
      if ( ILVSL ) then
         CALL DLASET( 'Full', N, N, ZERO, ONE, VSL, LDVSL )
         CALL DLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, &
                      VSL( ILO+1, ILO ), LDVSL )
         CALL DORGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, &
                      WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      end if
!
!     Initialize VSR
!
      if ( ILVSR ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, VSR, LDVSR )
!
!     Reduce to generalized Hessenberg form
!     (Workspace: none needed)
!
      CALL DGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, &
                   LDVSL, VSR, LDVSR, IERR )
!
!     Perform QZ algorithm, computing Schur vectors if desired
!     (Workspace: need N)
!
      IWRK = ITAU
      CALL DHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, &
                   ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, &
                   WORK( IWRK ), LWORK+1-IWRK, IERR )
      if ( IERR.NE.0 ) then
         if ( IERR.GT.0 .AND. IERR.LE.N ) then
            INFO = IERR
         else if ( IERR.GT.N .AND. IERR.LE.2*N ) then
            INFO = IERR - N
         ELSE
            INFO = N + 1
         end if
         GO TO 50
      end if
!
!     Sort eigenvalues ALPHA/BETA if desired
!     (Workspace: need 4*N+16 )
!
      SDIM = 0
      if ( WANTST ) then
!
!        Undo scaling on eigenvalues before DELZTGing
!
         if ( ILASCL ) then
            CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, &
                         IERR )
            CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, &
                         IERR )
         end if
         if ( ILBSCL ) &
            CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
!
!        Select eigenvalues
!
         DO 10 I = 1, N
            BWORK( I ) = DELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
   10    CONTINUE
!
         CALL DTGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, &
                      ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, &
                      PVSR, DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, &
                      IERR )
         if ( IERR == 1 ) &
            INFO = N + 3
!
      end if
!
!     Apply back-permutation to VSL and VSR
!     (Workspace: none needed)
!
      if ( ILVSL ) &
         CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), &
                      WORK( IRIGHT ), N, VSL, LDVSL, IERR )
!
      if ( ILVSR ) &
         CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), &
                      WORK( IRIGHT ), N, VSR, LDVSR, IERR )
!
!     Check if unscaling would cause over/underflow, if so, rescale
!     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
!     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)
!
      if ( ILASCL ) then
         DO 20 I = 1, N
            if ( ALPHAI( I ).NE.ZERO ) then
               if ( ( ALPHAR( I ) / SAFMAX ).GT.( ANRMTO / ANRM ) .OR. &
                   ( SAFMIN / ALPHAR( I ) ).GT.( ANRM / ANRMTO ) ) then
                  WORK( 1 ) = ABS( A( I, I ) / ALPHAR( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               else if ( ( ALPHAI( I ) / SAFMAX ).GT. &
                        ( ANRMTO / ANRM ) .OR. &
                        ( SAFMIN / ALPHAI( I ) ).GT.( ANRM / ANRMTO ) ) &
                         THEN
                  WORK( 1 ) = ABS( A( I, I+1 ) / ALPHAI( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               end if
            end if
   20    CONTINUE
      end if
!
      if ( ILBSCL ) then
         DO 30 I = 1, N
            if ( ALPHAI( I ).NE.ZERO ) then
               if ( ( BETA( I ) / SAFMAX ).GT.( BNRMTO / BNRM ) .OR. &
                   ( SAFMIN / BETA( I ) ).GT.( BNRM / BNRMTO ) ) then
                  WORK( 1 ) = ABS( B( I, I ) / BETA( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               end if
            end if
   30    CONTINUE
      end if
!
!     Undo scaling
!
      if ( ILASCL ) then
         CALL DLASCL( 'H', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
      end if
!
      if ( ILBSCL ) then
         CALL DLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      end if
!
      if ( WANTST ) then
!
!        Check if reordering is correct
!
         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 40 I = 1, N
            CURSL = DELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
            if ( ALPHAI( I ) == ZERO ) then
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
   40    CONTINUE
!
      end if
!
   50 CONTINUE
!
      WORK( 1 ) = MAXWRK
!
      RETURN
!
!     End of DGGES
!
      END
