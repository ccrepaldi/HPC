      SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
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
!  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
!  the generalized eigenvalues, and optionally, the left and/or right
!  generalized eigenvectors.
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
!                   A * v(j) = lambda(j) * B * v(j).
!
!  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
!  of (A,B) satisfies
!
!                   u(j)**H * A  = lambda(j) * u(j)**H * B .
!
!  where u(j)**H is the conjugate-transpose of u(j).
!
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
!          On entry, the matrix A in the pair (A,B).
!          On exit, A has been overwritten.
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
!          On entry, the matrix B in the pair (A,B).
!          On exit, B has been overwritten.
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
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order as
!          their eigenvalues. If the j-th eigenvalue is real, then
!          u(j) = VL(:,j), the j-th column of VL. If the j-th and
!          (j+1)-th eigenvalues form a complex conjugate pair, then
!          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
!          Each eigenvector will be scaled so the largest component have
!          abs(real part)+abs(imag. part)=1.
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
!          abs(real part)+abs(imag. part)=1.
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
!          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
!                =N+2: error return from DTGEVC.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, &
                         IN, IRIGHT, IROWS, ITAU, IWRK, JC, JR, MAXWRK, &
                         MINWRK
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, &
                         SMLNUM, TEMP
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
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK == -1 )
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
         MAXWRK = 7*N + N*ILAENV( 1, 'DGEQRF', ' ', N, 1, N, 0 )
         MINWRK = MAX( 1, 8*N )
         WORK( 1 ) = MAXWRK
      end if
!
      if ( LWORK < MINWRK .AND. .NOT.LQUERY ) &
         INFO = -16
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGGEV ', -INFO )
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
!     Permute the matrices A, B to isolate eigenvalues if possible
!     (Workspace: need 6*N)
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
      if ( ILV ) then
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      end if
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
!     Initialize VL
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
!     Initialize VR
!
      if ( ILVR ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )
!
!     Reduce to generalized Hessenberg form
!     (Workspace: none needed)
!
      if ( ILV ) then
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
      IWRK = ITAU
      if ( ILV ) then
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      end if
      CALL DHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, &
                   ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, &
                   WORK( IWRK ), LWORK+1-IWRK, IERR )
      if ( IERR.NE.0 ) then
         if ( IERR.GT.0 .AND. IERR.LE.N ) then
            INFO = IERR
         else if ( IERR.GT.N .AND. IERR.LE.2*N ) then
            INFO = IERR - N
         ELSE
            INFO = N + 1
         end if
         GO TO 110
      end if
!
!     Compute Eigenvectors
!     (Workspace: need 6*N)
!
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
         CALL DTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, &
                      VR, LDVR, N, IN, WORK( IWRK ), IERR )
         if ( IERR.NE.0 ) then
            INFO = N + 2
            GO TO 110
         end if
!
!        Undo balancing on VL and VR and normalization
!        (Workspace: none needed)
!
         if ( ILVL ) then
            CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), &
                         WORK( IRIGHT ), N, VL, LDVL, IERR )
            DO 50 JC = 1, N
               if ( ALPHAI( JC ) < ZERO ) &
                  GO TO 50
               TEMP = ZERO
               if ( ALPHAI( JC ) == ZERO ) then
                  DO 10 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
   10             CONTINUE
               ELSE
                  DO 20 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ &
                            ABS( VL( JR, JC+1 ) ) )
   20             CONTINUE
               end if
               if ( TEMP < SMLNUM ) &
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
                         WORK( IRIGHT ), N, VR, LDVR, IERR )
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
               if ( TEMP < SMLNUM ) &
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
  110 CONTINUE
!
      WORK( 1 ) = MAXWRK
!
      RETURN
!
!     End of DGGEV
!
      END
