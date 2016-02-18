      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     December 8, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WI( * ), WORK( * ), WR( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
!  eigenvalues and, optionally, the left and/or right eigenvectors.
!
!  The right eigenvector v(j) of A satisfies
!                   A * v(j) = lambda(j) * v(j)
!  where lambda(j) is its eigenvalue.
!  The left eigenvector u(j) of A satisfies
!                u(j)**H * A = lambda(j) * u(j)**H
!  where u(j)**H denotes the conjugate transpose of u(j).
!
!  The computed eigenvectors are normalized to have Euclidean norm
!  equal to 1 and largest component real.
!
!  Arguments
!  =========
!
!  JOBVL   (input) CHARACTER*1
!          = 'N': left eigenvectors of A are not computed;
!          = 'V': left eigenvectors of A are computed.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N': right eigenvectors of A are not computed;
!          = 'V': right eigenvectors of A are computed.
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit, A has been overwritten.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  WR      (output) DOUBLE PRECISION array, dimension (N)
!  WI      (output) DOUBLE PRECISION array, dimension (N)
!          WR and WI contain the real and imaginary parts,
!          respectively, of the computed eigenvalues.  Complex
!          conjugate pairs of eigenvalues appear consecutively
!          with the eigenvalue having the positive imaginary part
!          first.
!
!  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order
!          as their eigenvalues.
!          If JOBVL = 'N', VL is not referenced.
!          If the j-th eigenvalue is real, then u(j) = VL(:,j),
!          the j-th column of VL.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
!          u(j+1) = VL(:,j) - i*VL(:,j+1).
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= 1; if
!          JOBVL = 'V', LDVL >= N.
!
!  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
!          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!          after another in the columns of VR, in the same order
!          as their eigenvalues.
!          If JOBVR = 'N', VR is not referenced.
!          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!          the j-th column of VR.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!          v(j+1) = VR(:,j) - i*VR(:,j+1).
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= 1; if
!          JOBVR = 'V', LDVR >= N.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,3*N), and
!          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
!          performance, LWORK must generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = i, the QR algorithm failed to compute all the
!                eigenvalues, and no eigenvectors have been computed;
!                elements i+1:N of WR and WI contain eigenvalues which
!                have converged.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K, &
                         MAXB, MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, &
                         SN
!     ..
!     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY, DLARTG, &
                         DLASCL, DORGHR, DROT, DSCAL, DTREVC, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2, &
                         DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK == -1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      if ( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) then
         INFO = -1
      else if ( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -5
      else if ( LDVL < 1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) then
         INFO = -9
      else if ( LDVR < 1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) then
         INFO = -11
      end if
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by DHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      MINWRK = 1
      if ( INFO == 0 .AND. ( LWORK.GE.1 .OR. LQUERY ) ) then
         MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
         if ( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) then
            MINWRK = MAX( 1, 3*N )
            MAXB = MAX( ILAENV( 8, 'DHSEQR', 'EN', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'EN', N, 1, &
                N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, N+1, N+HSWORK )
         ELSE
            MINWRK = MAX( 1, 4*N )
            MAXWRK = MAX( MAXWRK, 2*N+( N-1 )* &
                     ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) )
            MAXB = MAX( ILAENV( 8, 'DHSEQR', 'SV', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'SV', N, 1, &
                N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, N+1, N+HSWORK )
            MAXWRK = MAX( MAXWRK, 4*N )
         end if
         WORK( 1 ) = MAXWRK
      end if
      if ( LWORK < MINWRK .AND. .NOT.LQUERY ) then
         INFO = -13
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEEV ', -INFO )
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
!     Balance the matrix
!     (Workspace: need N)
!
      IBAL = 1
      CALL DGEBAL( 'B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 3*N, prefer 2*N+N*NB)
!
      ITAU = IBAL + N
      IWRK = ITAU + N
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
                   LWORK-IWRK+1, IERR )
!
      if ( WANTVL ) then
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         SIDE = 'L'
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
!
!        Generate orthogonal matrix in VL
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
         if ( WANTVR ) then
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            SIDE = 'B'
            CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         end if
!
      else if ( WANTVR ) then
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         SIDE = 'R'
         CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )
!
!        Generate orthogonal matrix in VR
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
      ELSE
!
!        Compute eigenvalues only
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
      end if
!
!     If INFO > 0 from DHSEQR, then quit
!
      if ( INFO.GT.0 ) &
         GO TO 50
!
      if ( WANTVL .OR. WANTVR ) then
!
!        Compute left and/or right eigenvectors
!        (Workspace: need 4*N)
!
         CALL DTREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                      N, NOUT, WORK( IWRK ), IERR )
      end if
!
      if ( WANTVL ) then
!
!        Undo balancing of left eigenvectors
!        (Workspace: need N)
!
         CALL DGEBAK( 'B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL, &
                      IERR )
!
!        Normalize left eigenvectors and make largest component real
!
         DO 20 I = 1, N
            if ( WI( I ) == ZERO ) then
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            else if ( WI( I ).GT.ZERO ) then
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ), &
                     DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO K = 1, N
                  WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2
               end do
               K = IDAMAX( N, WORK( IWRK ), 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            end if
   20    CONTINUE
      end if
!
      if ( WANTVR ) then
!
!        Undo balancing of right eigenvectors
!        (Workspace: need N)
!
         CALL DGEBAK( 'B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR, &
                      IERR )
!
!        Normalize right eigenvectors and make largest component real
!
         DO 40 I = 1, N
            if ( WI( I ) == ZERO ) then
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            else if ( WI( I ).GT.ZERO ) then
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ), &
                     DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = IDAMAX( N, WORK( IWRK ), 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            end if
   40    CONTINUE
      end if
!
!     Undo scaling if necessary
!
   50 CONTINUE
      if ( SCALEA ) then
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         if ( INFO.GT.0 ) then
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, &
                         IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, &
                         IERR )
         end if
      end if
!
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of DGEEV
!
      END
