      SUBROUTINE DGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
                        BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                         B( LDB, * ), BETA( * ), VL( LDVL, * ), &
                         VR( LDVR, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  This routine is deprecated and has been replaced by routine DGGEV.
!
!  DGEGV computes for a pair of n-by-n real nonsymmetric matrices A and
!  B, the generalized eigenvalues (alphar +/- alphai*i, beta), and
!  optionally, the left and/or right generalized eigenvectors (VL and
!  VR).
!
!  A generalized eigenvalue for a pair of matrices (A,B) is, roughly
!  speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B
!  is singular.  It is usually represented as the pair (alpha,beta),
!  as there is a reasonable interpretation for beta=0, and even for
!  both being zero.  A good beginning reference is the book, "Matrix
!  Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)
!
!  A right generalized eigenvector corresponding to a generalized
!  eigenvalue  w  for a pair of matrices (A,B) is a vector  r  such
!  that  (A - w B) r = 0 .  A left generalized eigenvector is a vector
!  l such that l**H * (A - w B) = 0, where l**H is the
!  conjugate-transpose of l.
!
!  Note: this routine performs "full balancing" on A and B -- see
!  "Further Details", below.
!
!  Arguments
!  =========
!
!  JOBVL   (input) CHARACTER*1
!          = 'N':  do not compute the left generalized eigenvectors;
!          = 'V':  compute the left generalized eigenvectors.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N':  do not compute the right generalized eigenvectors;
!          = 'V':  compute the right generalized eigenvectors.
!
!  N       (input) INTEGER
!          The order of the matrices A, B, VL, and VR.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the first of the pair of matrices whose
!          generalized eigenvalues and (optionally) generalized
!          eigenvectors are to be computed.
!          On exit, the contents will have been destroyed.  (For a
!          description of the contents of A on exit, see "Further
!          Details", below.)
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
!          On entry, the second of the pair of matrices whose
!          generalized eigenvalues and (optionally) generalized
!          eigenvectors are to be computed.
!          On exit, the contents will have been destroyed.  (For a
!          description of the contents of B on exit, see "Further
!          Details", below.)
!
!  LDB     (input) INTEGER
!          The leading dimension of B.  LDB >= max(1,N).
!
!  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
!  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
!  BETA    (output) DOUBLE PRECISION array, dimension (N)
!          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
!          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
!          the j-th eigenvalue is real; if positive, then the j-th and
!          (j+1)-st eigenvalues are a complex conjugate pair, with
!          ALPHAI(j+1) negative.
!
!          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
!          may easily over- or underflow, and BETA(j) may even be zero.
!          Thus, the user should avoid naively computing the ratio
!          alpha/beta.  However, ALPHAR and ALPHAI will be always less
!          than and usually comparable with norm(A) in magnitude, and
!          BETA always less than and usually comparable with norm(B).
!
!  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
!          If JOBVL = 'V', the left generalized eigenvectors.  (See
!          "Purpose", above.)  Real eigenvectors take one column,
!          complex take two columns, the first for the real part and
!          the second for the imaginary part.  Complex eigenvectors
!          correspond to an eigenvalue with positive imaginary part.
!          Each eigenvector will be scaled so the largest component
!          will have abs(real part) + abs(imag. part) = 1, *except*
!          that for eigenvalues with alpha=beta=0, a zero vector will
!          be returned as the corresponding eigenvector.
!          Not referenced if JOBVL = 'N'.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the matrix VL. LDVL >= 1, and
!          if JOBVL = 'V', LDVL >= N.
!
!  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
!          If JOBVR = 'V', the right generalized eigenvectors.  (See
!          "Purpose", above.)  Real eigenvectors take one column,
!          complex take two columns, the first for the real part and
!          the second for the imaginary part.  Complex eigenvectors
!          correspond to an eigenvalue with positive imaginary part.
!          Each eigenvector will be scaled so the largest component
!          will have abs(real part) + abs(imag. part) = 1, *except*
!          that for eigenvalues with alpha=beta=0, a zero vector will
!          be returned as the corresponding eigenvector.
!          Not referenced if JOBVR = 'N'.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the matrix VR. LDVR >= 1, and
!          if JOBVR = 'V', LDVR >= N.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,8*N).
!          For good performance, LWORK must generally be larger.
!          To compute the optimal value of LWORK, call ILAENV to get
!          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:
!          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR;
!          The optimal LWORK is:
!              2*N + MAX( 6*N, N*(NB+1) ).
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          = 1,...,N:
!                The QZ iteration failed.  No eigenvectors have been
!                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
!                should be correct for j=INFO+1,...,N.
!          > N:  errors that usually indicate LAPACK problems:
!                =N+1: error return from DGGBAL
!                =N+2: error return from DGEQRF
!                =N+3: error return from DORMQR
!                =N+4: error return from DORGQR
!                =N+5: error return from DGGHRD
!                =N+6: error return from DHGEQZ (other than failed
!                                                iteration)
!                =N+7: error return from DTGEVC
!                =N+8: error return from DGGBAK (computing VL)
!                =N+9: error return from DGGBAK (computing VR)
!                =N+10: error return from DLASCL (various calls)
!
!  Further Details
!  ===============
!
!  Balancing
!  ---------
!
!  This driver calls DGGBAL to both permute and scale rows and columns
!  of A and B.  The permutations PL and PR are chosen so that PL*A*PR
!  and PL*B*R will be upper triangular except for the diagonal blocks
!  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as
!  possible.  The diagonal scaling matrices DL and DR are chosen so
!  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to
!  one (except for the elements that start out zero.)
!
!  After the eigenvalues and eigenvectors of the balanced matrices
!  have been computed, DGGBAK transforms the eigenvectors back to what
!  they would have been (in perfect arithmetic) if they had not been
!  balanced.
!
!  Contents of A and B on Exit
!  -------- -- - --- - -- ----
!
!  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or
!  both), then on exit the arrays A and B will contain the real Schur
!  form[*] of the "balanced" versions of A and B.  If no eigenvectors
!  are computed, then only the diagonal blocks will be correct.
!
!  [*] See DHGEQZ, DGEGS, or read the book "Matrix Computations",
!      by Golub & van Loan, pub. by Johns Hopkins U. Press.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ILIMIT, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, &
                         IN, IRIGHT, IROWS, ITAU, IWORK, JC, JR, LOPT, &
                         LWKMIN, LWKOPT, NB, NB1, NB2, NB3
      DOUBLE PRECISION   ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM, &
                         BNRM1, BNRM2, EPS, ONEPLS, SAFMAX, SAFMIN, &
                         SALFAI, SALFAR, SBETA, SCALE, TEMP
!     ..
!     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, &
                         DLASCL, DLASET, DORGQR, DORMQR, DTGEVC, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, MAX
!     ..
!     .. Executable Statements ..
!
!     Decode the input arguments
!
      if ( LSAME( JOBVL, 'N' ) ) then
         IJOBVL = 1
         ILVL = .FALSE.
      else if ( LSAME( JOBVL, 'V' ) ) then
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      end if
!
      if ( LSAME( JOBVR, 'N' ) ) then
         IJOBVR = 1
         ILVR = .FALSE.
      else if ( LSAME( JOBVR, 'V' ) ) then
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      end if
      ILV = ILVL .OR. ILVR
!
!     Test the input arguments
!
      LWKMIN = MAX( 8*N, 1 )
      LWKOPT = LWKMIN
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK == -1 )
      INFO = 0
      if ( IJOBVL.LE.0 ) then
         INFO = -1
      else if ( IJOBVR.LE.0 ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -5
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -7
      else if ( LDVL < 1 .OR. ( ILVL .AND. LDVL.LT.N ) ) then
         INFO = -12
      else if ( LDVR < 1 .OR. ( ILVR .AND. LDVR.LT.N ) ) then
         INFO = -14
      else if ( LWORK < LWKMIN .AND. .NOT.LQUERY ) then
         INFO = -16
      end if
!
      if ( INFO == 0 ) then
         NB1 = ILAENV( 1, 'DGEQRF', ' ', N, N, -1, -1 )
         NB2 = ILAENV( 1, 'DORMQR', ' ', N, N, N, -1 )
         NB3 = ILAENV( 1, 'DORGQR', ' ', N, N, N, -1 )
         NB = MAX( NB1, NB2, NB3 )
         LOPT = 2*N + MAX( 6*N, N*( NB+1 ) )
         WORK( 1 ) = LOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEGV ', -INFO )
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
!     Get machine constants
!
      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SAFMIN = DLAMCH( 'S' )
      SAFMIN = SAFMIN + SAFMIN
      SAFMAX = ONE / SAFMIN
      ONEPLS = ONE + ( 4*EPS )
!
!     Scale A
!
      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ANRM1 = ANRM
      ANRM2 = ONE
      if ( ANRM < ONE ) then
         if ( SAFMAX*ANRM < ONE ) then
            ANRM1 = SAFMIN
            ANRM2 = SAFMAX*ANRM
         end if
      end if
!
      if ( ANRM.GT.ZERO ) then
         CALL DLASCL( 'G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 10
            RETURN
         end if
      end if
!
!     Scale B
!
      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      BNRM1 = BNRM
      BNRM2 = ONE
      if ( BNRM < ONE ) then
         if ( SAFMAX*BNRM < ONE ) then
            BNRM1 = SAFMIN
            BNRM2 = SAFMAX*BNRM
         end if
      end if
!
      if ( BNRM.GT.ZERO ) then
         CALL DLASCL( 'G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 10
            RETURN
         end if
      end if
!
!     Permute the matrix to make it more nearly triangular
!     Workspace layout:  (8*N words -- "work" requires 6*N words)
!        left_permutation, right_permutation, work...
!
      ILEFT = 1
      IRIGHT = N + 1
      IWORK = IRIGHT + N
      CALL DGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), &
                   WORK( IRIGHT ), WORK( IWORK ), IINFO )
      if ( IINFO.NE.0 ) then
         INFO = N + 1
         GO TO 120
      end if
!
!     Reduce B to triangular form, and initialize VL and/or VR
!     Workspace layout:  ("work..." must have at least N words)
!        left_permutation, right_permutation, tau, work...
!
      IROWS = IHI + 1 - ILO
      if ( ILV ) then
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      end if
      ITAU = IWORK
      IWORK = ITAU + IROWS
      CALL DGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), &
                   WORK( IWORK ), LWORK+1-IWORK, IINFO )
      if ( IINFO.GE.0 ) &
         LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) then
         INFO = N + 2
         GO TO 120
      end if
!
      CALL DORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, &
                   WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), &
                   LWORK+1-IWORK, IINFO )
      if ( IINFO.GE.0 ) &
         LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) then
         INFO = N + 3
         GO TO 120
      end if
!
      if ( ILVL ) then
         CALL DLASET( 'Full', N, N, ZERO, ONE, VL, LDVL )
         CALL DLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, &
                      VL( ILO+1, ILO ), LDVL )
         CALL DORGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, &
                      WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, &
                      IINFO )
         if ( IINFO.GE.0 ) &
            LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         if ( IINFO.NE.0 ) then
            INFO = N + 4
            GO TO 120
         end if
      end if
!
      if ( ILVR ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )
!
!     Reduce to generalized Hessenberg form
!
      if ( ILV ) then
!
!        Eigenvectors requested -- work on whole matrix.
!
         CALL DGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, &
                      LDVL, VR, LDVR, IINFO )
      ELSE
         CALL DGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, &
                      B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IINFO )
      end if
      if ( IINFO.NE.0 ) then
         INFO = N + 5
         GO TO 120
      end if
!
!     Perform QZ algorithm
!     Workspace layout:  ("work..." must have at least 1 word)
!        left_permutation, right_permutation, work...
!
      IWORK = ITAU
      if ( ILV ) then
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      end if
      CALL DHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, &
                   ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, &
                   WORK( IWORK ), LWORK+1-IWORK, IINFO )
      if ( IINFO.GE.0 ) &
         LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) then
         if ( IINFO.GT.0 .AND. IINFO.LE.N ) then
            INFO = IINFO
         else if ( IINFO.GT.N .AND. IINFO.LE.2*N ) then
            INFO = IINFO - N
         ELSE
            INFO = N + 6
         end if
         GO TO 120
      end if
!
      if ( ILV ) then
!
!        Compute Eigenvectors  (DTGEVC requires 6*N words of workspace)
!
         if ( ILVL ) then
            if ( ILVR ) then
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            end if
         ELSE
            CHTEMP = 'R'
         end if
!
         CALL DTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, &
                      VR, LDVR, N, IN, WORK( IWORK ), IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 7
            GO TO 120
         end if
!
!        Undo balancing on VL and VR, rescale
!
         if ( ILVL ) then
            CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), &
                         WORK( IRIGHT ), N, VL, LDVL, IINFO )
            if ( IINFO.NE.0 ) then
               INFO = N + 8
               GO TO 120
            end if
            DO 50 JC = 1, N
               if ( ALPHAI( JC ) < ZERO ) &
                  GO TO 50
               TEMP = ZERO
               if ( ALPHAI( JC ) == ZERO ) then
                  DO JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
                  end do
               ELSE
                  DO 20 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ &
                            ABS( VL( JR, JC+1 ) ) )
   20             CONTINUE
               end if
               if ( TEMP < SAFMIN ) &
                  GO TO 50
               TEMP = ONE / TEMP
               if ( ALPHAI( JC ) == ZERO ) then
                  DO 30 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
   30             CONTINUE
               ELSE
                  DO 40 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                     VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
   40             CONTINUE
               end if
   50       CONTINUE
         end if
         if ( ILVR ) then
            CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), &
                         WORK( IRIGHT ), N, VR, LDVR, IINFO )
            if ( IINFO.NE.0 ) then
               INFO = N + 9
               GO TO 120
            end if
            DO 100 JC = 1, N
               if ( ALPHAI( JC ) < ZERO ) &
                  GO TO 100
               TEMP = ZERO
               if ( ALPHAI( JC ) == ZERO ) then
                  DO 60 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
   60             CONTINUE
               ELSE
                  DO 70 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ &
                            ABS( VR( JR, JC+1 ) ) )
   70             CONTINUE
               end if
               if ( TEMP < SAFMIN ) &
                  GO TO 100
               TEMP = ONE / TEMP
               if ( ALPHAI( JC ) == ZERO ) then
                  DO 80 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
   80             CONTINUE
               ELSE
                  DO 90 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                     VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
   90             CONTINUE
               end if
  100       CONTINUE
         end if
!
!        End of eigenvector calculation
!
      end if
!
!     Undo scaling in alpha, beta
!
!     Note: this does not give the alpha and beta for the unscaled
!     problem.
!
!     Un-scaling is limited to avoid underflow in alpha and beta
!     if they are significant.
!
      DO 110 JC = 1, N
         ABSAR = ABS( ALPHAR( JC ) )
         ABSAI = ABS( ALPHAI( JC ) )
         ABSB = ABS( BETA( JC ) )
         SALFAR = ANRM*ALPHAR( JC )
         SALFAI = ANRM*ALPHAI( JC )
         SBETA = BNRM*BETA( JC )
         ILIMIT = .FALSE.
         SCALE = ONE
!
!        Check for significant underflow in ALPHAI
!
         if ( ABS( SALFAI ) < SAFMIN .AND. ABSAI.GE. &
             MAX( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) then
            ILIMIT = .TRUE.
            SCALE = ( ONEPLS*SAFMIN / ANRM1 ) / &
                    MAX( ONEPLS*SAFMIN, ANRM2*ABSAI )
!
         else if ( SALFAI == ZERO ) then
!
!           If insignificant underflow in ALPHAI, then make the
!           conjugate eigenvalue real.
!
            if ( ALPHAI( JC ) < ZERO .AND. JC.GT.1 ) then
               ALPHAI( JC-1 ) = ZERO
            else if ( ALPHAI( JC ).GT.ZERO .AND. JC < N ) then
               ALPHAI( JC+1 ) = ZERO
            end if
         end if
!
!        Check for significant underflow in ALPHAR
!
         if ( ABS( SALFAR ) < SAFMIN .AND. ABSAR.GE. &
             MAX( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) then
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / ANRM1 ) / &
                    MAX( ONEPLS*SAFMIN, ANRM2*ABSAR ) )
         end if
!
!        Check for significant underflow in BETA
!
         if ( ABS( SBETA ) < SAFMIN .AND. ABSB.GE. &
             MAX( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) then
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / BNRM1 ) / &
                    MAX( ONEPLS*SAFMIN, BNRM2*ABSB ) )
         end if
!
!        Check for possible overflow when limiting scaling
!
         if ( ILIMIT ) then
            TEMP = ( SCALE*SAFMIN )*MAX( ABS( SALFAR ), ABS( SALFAI ), &
                   ABS( SBETA ) )
            if ( TEMP.GT.ONE ) &
               SCALE = SCALE / TEMP
            if ( SCALE < ONE ) &
               ILIMIT = .FALSE.
         end if
!
!        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.
!
         if ( ILIMIT ) then
            SALFAR = ( SCALE*ALPHAR( JC ) )*ANRM
            SALFAI = ( SCALE*ALPHAI( JC ) )*ANRM
            SBETA = ( SCALE*BETA( JC ) )*BNRM
         end if
         ALPHAR( JC ) = SALFAR
         ALPHAI( JC ) = SALFAI
         BETA( JC ) = SBETA
  110 CONTINUE
!
  120 CONTINUE
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of DGEGV
!
      END
