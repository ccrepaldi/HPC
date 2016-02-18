      SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, &
                         ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, &
                         IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, &
                         RCONDV, WORK, LWORK, IWORK, BWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          BALANC, JOBVL, JOBVR, SENSE
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   ABNRM, BBNRM
!     ..
!     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                         B( LDB, * ), BETA( * ), LSCALE( * ), &
                         RCONDE( * ), RCONDV( * ), RSCALE( * ), &
                         VL( LDVL, * ), VR( LDVR, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
!  the generalized eigenvalues, and optionally, the left and/or right
!  generalized eigenvectors.
!
!  Optionally also, it computes a balancing transformation to improve
!  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
!  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
!  the eigenvalues (RCONDE), and reciprocal condition numbers for the
!  right eigenvectors (RCONDV).
!
!  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
!  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
!  singular. It is usually represented as the pair (alpha,beta), as
!  there is a reasonable interpretation for beta=0, and even for both
!  being zero.
!
!  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
!  of (A,B) satisfies
!
!                   A * v(j) = lambda(j) * B * v(j) .
!
!  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
!  of (A,B) satisfies
!
!                   u(j)**H * A  = lambda(j) * u(j)**H * B.
!
!  where u(j)**H is the conjugate-transpose of u(j).
!
!
!  Arguments
!  =========
!
!  BALANC  (input) CHARACTER*1
!          Specifies the balance option to be performed.
!          = 'N':  do not diagonally scale or permute;
!          = 'P':  permute only;
!          = 'S':  scale only;
!          = 'B':  both permute and scale.
!          Computed reciprocal condition numbers will be for the
!          matrices after permuting and/or balancing. Permuting does
!          not change condition numbers (in exact arithmetic), but
!          balancing does.
!
!  JOBVL   (input) CHARACTER*1
!          = 'N':  do not compute the left generalized eigenvectors;
!          = 'V':  compute the left generalized eigenvectors.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N':  do not compute the right generalized eigenvectors;
!          = 'V':  compute the right generalized eigenvectors.
!
!  SENSE   (input) CHARACTER*1
!          Determines which reciprocal condition numbers are computed.
!          = 'N': none are computed;
!          = 'E': computed for eigenvalues only;
!          = 'V': computed for eigenvectors only;
!          = 'B': computed for eigenvalues and eigenvectors.
!
!  N       (input) INTEGER
!          The order of the matrices A, B, VL, and VR.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the matrix A in the pair (A,B).
!          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
!          or both, then A contains the first part of the real Schur
!          form of the "balanced" versions of the input A and B.
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
!          On entry, the matrix B in the pair (A,B).
!          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
!          or both, then B contains the second part of the real Schur
!          form of the "balanced" versions of the input A and B.
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
!          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less
!          than and usually comparable with norm(A) in magnitude, and
!          BETA always less than and usually comparable with norm(B).
!
!  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order as
!          their eigenvalues. If the j-th eigenvalue is real, then
!          u(j) = VL(:,j), the j-th column of VL. If the j-th and
!          (j+1)-th eigenvalues form a complex conjugate pair, then
!          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
!          Each eigenvector will be scaled so the largest component have
!          abs(real part) + abs(imag. part) = 1.
!          Not referenced if JOBVL = 'N'.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the matrix VL. LDVL >= 1, and
!          if JOBVL = 'V', LDVL >= N.
!
!  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
!          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!          after another in the columns of VR, in the same order as
!          their eigenvalues. If the j-th eigenvalue is real, then
!          v(j) = VR(:,j), the j-th column of VR. If the j-th and
!          (j+1)-th eigenvalues form a complex conjugate pair, then
!          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
!          Each eigenvector will be scaled so the largest component have
!          abs(real part) + abs(imag. part) = 1.
!          Not referenced if JOBVR = 'N'.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the matrix VR. LDVR >= 1, and
!          if JOBVR = 'V', LDVR >= N.
!
!  ILO,IHI (output) INTEGER
!          ILO and IHI are integer values such that on exit
!          A(i,j) = 0 and B(i,j) = 0 if i > j and
!          j = 1,...,ILO-1 or i = IHI+1,...,N.
!          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
!
!  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied
!          to the left side of A and B.  If PL(j) is the index of the
!          row interchanged with row j, and DL(j) is the scaling
!          factor applied to row j, then
!            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
!                      = DL(j)  for j = ILO,...,IHI
!                      = PL(j)  for j = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied
!          to the right side of A and B.  If PR(j) is the index of the
!          column interchanged with column j, and DR(j) is the scaling
!          factor applied to column j, then
!            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
!                      = DR(j)  for j = ILO,...,IHI
!                      = PR(j)  for j = IHI+1,...,N
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  ABNRM   (output) DOUBLE PRECISION
!          The one-norm of the balanced matrix A.
!
!  BBNRM   (output) DOUBLE PRECISION
!          The one-norm of the balanced matrix B.
!
!  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
!          If SENSE = 'E' or 'B', the reciprocal condition numbers of
!          the selected eigenvalues, stored in consecutive elements of
!          the array. For a complex conjugate pair of eigenvalues two
!          consecutive elements of RCONDE are set to the same value.
!          Thus RCONDE(j), RCONDV(j), and the j-th columns of VL and VR
!          all correspond to the same eigenpair (but not in general the
!          j-th eigenpair, unless all eigenpairs are selected).
!          If SENSE = 'V', RCONDE is not referenced.
!
!  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
!          If SENSE = 'V' or 'B', the estimated reciprocal condition
!          numbers of the selected eigenvectors, stored in consecutive
!          elements of the array. For a complex eigenvector two
!          consecutive elements of RCONDV are set to the same value. If
!          the eigenvalues cannot be reordered to compute RCONDV(j),
!          RCONDV(j) is set to 0; this can only occur when the true
!          value would be very small anyway.
!          If SENSE = 'E', RCONDV is not referenced.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,6*N).
!          If SENSE = 'E', LWORK >= 12*N.
!          If SENSE = 'V' or 'B', LWORK >= 2*N*N+12*N+16.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace) INTEGER array, dimension (N+6)
!          If SENSE = 'E', IWORK is not referenced.
!
!  BWORK   (workspace) LOGICAL array, dimension (N)
!          If SENSE = 'N', BWORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          = 1,...,N:
!                The QZ iteration failed.  No eigenvectors have been
!                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
!                should be correct for j=INFO+1,...,N.
!          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
!                =N+2: error return from DTGEVC.
!
!  Further Details
!  ===============
!
!  Balancing a matrix pair (A,B) includes, first, permuting rows and
!  columns to isolate eigenvalues, second, applying diagonal similarity
!  transformation to the rows and columns to make the rows and columns
!  as close in norm as possible. The computed reciprocal condition
!  numbers correspond to the balanced matrix. Permuting rows and columns
!  will not change the condition numbers (in exact arithmetic) but
!  diagonal scaling will.  For further explanation of balancing, see
!  section 4.11.1.2 of LAPACK Users' Guide.
!
!  An approximate error bound on the chordal distance between the i-th
!  computed generalized eigenvalue w and the corresponding exact
!  eigenvalue lambda is
!
!       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
!
!  An approximate error bound for the angle between the i-th computed
!  eigenvector VL(i) or VR(i) is given by
!
!       EPS * norm(ABNRM, BBNRM) / Dif (i).
!
!  For further explanation of the reciprocal condition numbers RCONDE
!  and RCONDV, see section 4.11 of LAPACK User's Guide.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY, PAIR, &
                         WANTSB, WANTSE, WANTSN, WANTSV
      CHARACTER          CHTEMP
      INTEGER            I, ICOLS, IERR, IJOBVL, IJOBVR, IN, IROWS, &
                         ITAU, IWRK, IWRK1, J, JC, JR, M, MAXWRK, &
                         MINWRK, MM
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, &
                         SMLNUM, TEMP
!     ..
!     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, &
                         DLASCL, DLASET, DORGQR, DORMQR, DTGEVC, DTGSNA, &
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
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( .NOT.( LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, &
          'S' ) .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' ) ) ) &
           THEN
         INFO = -1
      else if ( IJOBVL.LE.0 ) then
         INFO = -2
      else if ( IJOBVR.LE.0 ) then
         INFO = -3
      else if ( .NOT.( WANTSN .OR. WANTSE .OR. WANTSB .OR. WANTSV ) ) &
                THEN
         INFO = -4
      else if ( N < 0 ) then
         INFO = -5
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -7
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -9
      else if ( LDVL < 1 .OR. ( ILVL .AND. LDVL.LT.N ) ) then
         INFO = -14
      else if ( LDVR < 1 .OR. ( ILVR .AND. LDVR.LT.N ) ) then
         INFO = -16
      end if
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV. The workspace is
!       computed assuming ILO = 1 and IHI = N, the worst case.)
!
      MINWRK = 1
      if ( INFO == 0 .AND. ( LWORK.GE.1 .OR. LQUERY ) ) then
         MAXWRK = 5*N + N*ILAENV( 1, 'DGEQRF', ' ', N, 1, N, 0 )
         MINWRK = MAX( 1, 6*N )
         if ( WANTSE ) then
            MINWRK = MAX( 1, 12*N )
         else if ( WANTSV .OR. WANTSB ) then
            MINWRK = 2*N*N + 12*N + 16
            MAXWRK = MAX( MAXWRK, 2*N*N+12*N+16 )
         end if
         WORK( 1 ) = MAXWRK
      end if
!
      if ( LWORK < MINWRK .AND. .NOT.LQUERY ) then
         INFO = -26
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGGEVX', -INFO )
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
!     Permute and/or balance the matrix pair (A,B)
!     (Workspace: need 6*N)
!
      CALL DGGBAL( BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, &
                   WORK, IERR )
!
!     Compute ABNRM and BBNRM
!
      ABNRM = DLANGE( '1', N, N, A, LDA, WORK( 1 ) )
      if ( ILASCL ) then
         WORK( 1 ) = ABNRM
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, 1, 1, WORK( 1 ), 1, &
                      IERR )
         ABNRM = WORK( 1 )
      end if
!
      BBNRM = DLANGE( '1', N, N, B, LDB, WORK( 1 ) )
      if ( ILBSCL ) then
         WORK( 1 ) = BBNRM
         CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, 1, 1, WORK( 1 ), 1, &
                      IERR )
         BBNRM = WORK( 1 )
      end if
!
!     Reduce B to triangular form (QR decomposition of B)
!     (Workspace: need N, prefer N*NB )
!
      IROWS = IHI + 1 - ILO
      if ( ILV .OR. .NOT.WANTSN ) then
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      end if
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL DGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), &
                   WORK( IWRK ), LWORK+1-IWRK, IERR )
!
!     Apply the orthogonal transformation to A
!     (Workspace: need N, prefer N*NB)
!
      CALL DORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, &
                   WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), &
                   LWORK+1-IWRK, IERR )
!
!     Initialize VL and/or VR
!     (Workspace: need N, prefer N*NB)
!
      if ( ILVL ) then
         CALL DLASET( 'Full', N, N, ZERO, ONE, VL, LDVL )
         CALL DLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, &
                      VL( ILO+1, ILO ), LDVL )
         CALL DORGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, &
                      WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      end if
!
      if ( ILVR ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )
!
!     Reduce to generalized Hessenberg form
!     (Workspace: none needed)
!
      if ( ILV .OR. .NOT.WANTSN ) then
!
!        Eigenvectors requested -- work on whole matrix.
!
         CALL DGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, &
                      LDVL, VR, LDVR, IERR )
      ELSE
         CALL DGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, &
                      B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR )
      end if
!
!     Perform QZ algorithm (Compute eigenvalues, and optionally, the
!     Schur forms and Schur vectors)
!     (Workspace: need N)
!
      if ( ILV .OR. .NOT.WANTSN ) then
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      end if
!
      CALL DHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, &
                   ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, &
                   LWORK, IERR )
      if ( IERR.NE.0 ) then
         if ( IERR.GT.0 .AND. IERR.LE.N ) then
            INFO = IERR
         else if ( IERR.GT.N .AND. IERR.LE.2*N ) then
            INFO = IERR - N
         ELSE
            INFO = N + 1
         end if
         GO TO 130
      end if
!
!     Compute Eigenvectors and estimate condition numbers if desired
!     (Workspace: DTGEVC: need 6*N
!                 DTGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B',
!                         need N otherwise )
!
      if ( ILV .OR. .NOT.WANTSN ) then
         if ( ILV ) then
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
            CALL DTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, &
                         LDVL, VR, LDVR, N, IN, WORK, IERR )
            if ( IERR.NE.0 ) then
               INFO = N + 2
               GO TO 130
            end if
         end if
!
         if ( .NOT.WANTSN ) then
!
!           compute eigenvectors (DTGEVC) and estimate condition
!           numbers (DTGSNA). Note that the definition of the condition
!           number is not invariant under transformation (u,v) to
!           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
!           Schur form (S,T), Q and Z are orthogonal matrices. In order
!           to avoid using extra 2*N*N workspace, we have to recalculate
!           eigenvectors and estimate one condition numbers at a time.
!
            PAIR = .FALSE.
            DO 20 I = 1, N
!
               if ( PAIR ) then
                  PAIR = .FALSE.
                  GO TO 20
               end if
               MM = 1
               if ( I < N ) then
                  if ( A( I+1, I ).NE.ZERO ) then
                     PAIR = .TRUE.
                     MM = 2
                  end if
               end if
!
               DO 10 J = 1, N
                  BWORK( J ) = .FALSE.
   10          CONTINUE
               if ( MM == 1 ) then
                  BWORK( I ) = .TRUE.
               else if ( MM == 2 ) then
                  BWORK( I ) = .TRUE.
                  BWORK( I+1 ) = .TRUE.
               end if
!
               IWRK = MM*N + 1
               IWRK1 = IWRK + MM*N
!
!              Compute a pair of left and right eigenvectors.
!              (compute workspace: need up to 4*N + 6*N)
!
               if ( WANTSE .OR. WANTSB ) then
                  CALL DTGEVC( 'B', 'S', BWORK, N, A, LDA, B, LDB, &
                               WORK( 1 ), N, WORK( IWRK ), N, MM, M, &
                               WORK( IWRK1 ), IERR )
                  if ( IERR.NE.0 ) then
                     INFO = N + 2
                     GO TO 130
                  end if
               end if
!
               CALL DTGSNA( SENSE, 'S', BWORK, N, A, LDA, B, LDB, &
                            WORK( 1 ), N, WORK( IWRK ), N, RCONDE( I ), &
                            RCONDV( I ), MM, M, WORK( IWRK1 ), &
                            LWORK-IWRK1+1, IWORK, IERR )
!
   20       CONTINUE
         end if
      end if
!
!     Undo balancing on VL and VR and normalization
!     (Workspace: none needed)
!
      if ( ILVL ) then
         CALL DGGBAK( BALANC, 'L', N, ILO, IHI, LSCALE, RSCALE, N, VL, &
                      LDVL, IERR )
!
         DO 70 JC = 1, N
            if ( ALPHAI( JC ) < ZERO ) &
               GO TO 70
            TEMP = ZERO
            if ( ALPHAI( JC ) == ZERO ) then
               DO 30 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
   30          CONTINUE
            ELSE
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ &
                         ABS( VL( JR, JC+1 ) ) )
   40          CONTINUE
            end if
            if ( TEMP < SMLNUM ) &
               GO TO 70
            TEMP = ONE / TEMP
            if ( ALPHAI( JC ) == ZERO ) then
               DO 50 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   50          CONTINUE
            ELSE
               DO 60 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
                  VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
   60          CONTINUE
            end if
   70    CONTINUE
      end if
      if ( ILVR ) then
         CALL DGGBAK( BALANC, 'R', N, ILO, IHI, LSCALE, RSCALE, N, VR, &
                      LDVR, IERR )
         DO 120 JC = 1, N
            if ( ALPHAI( JC ) < ZERO ) &
               GO TO 120
            TEMP = ZERO
            if ( ALPHAI( JC ) == ZERO ) then
               DO 80 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
   80          CONTINUE
            ELSE
               DO 90 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ &
                         ABS( VR( JR, JC+1 ) ) )
   90          CONTINUE
            end if
            if ( TEMP < SMLNUM ) &
               GO TO 120
            TEMP = ONE / TEMP
            if ( ALPHAI( JC ) == ZERO ) then
               DO 100 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
  100          CONTINUE
            ELSE
               DO 110 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
                  VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
  110          CONTINUE
            end if
  120    CONTINUE
      end if
!
!     Undo scaling if necessary
!
      if ( ILASCL ) then
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
      end if
!
      if ( ILBSCL ) then
         CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      end if
!
  130 CONTINUE
      WORK( 1 ) = MAXWRK
!
      RETURN
!
!     End of DGGEVX
!
      END
