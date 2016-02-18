      SUBROUTINE DGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR, &
                        ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, &
                        LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVSL, JOBVSR
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                         B( LDB, * ), BETA( * ), VSL( LDVSL, * ), &
                         VSR( LDVSR, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  This routine is deprecated and has been replaced by routine DGGES.
!
!  DGEGS computes for a pair of N-by-N real nonsymmetric matrices A, B:
!  the generalized eigenvalues (alphar +/- alphai*i, beta), the real
!  Schur form (A, B), and optionally left and/or right Schur vectors
!  (VSL and VSR).
!
!  (If only the generalized eigenvalues are needed, use the driver DGEGV
!  instead.)
!
!  A generalized eigenvalue for a pair of matrices (A,B) is, roughly
!  speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B
!  is singular.  It is usually represented as the pair (alpha,beta),
!  as there is a reasonable interpretation for beta=0, and even for
!  both being zero.  A good beginning reference is the book, "Matrix
!  Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)
!
!  The (generalized) Schur form of a pair of matrices is the result of
!  multiplying both matrices on the left by one orthogonal matrix and
!  both on the right by another orthogonal matrix, these two orthogonal
!  matrices being chosen so as to bring the pair of matrices into
!  (real) Schur form.
!
!  A pair of matrices A, B is in generalized real Schur form if B is
!  upper triangular with non-negative diagonal and A is block upper
!  triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
!  to real generalized eigenvalues, while 2-by-2 blocks of A will be
!  "standardized" by making the corresponding elements of B have the
!  form:
!          [  a  0  ]
!          [  0  b  ]
!
!  and the pair of corresponding 2-by-2 blocks in A and B will
!  have a complex conjugate pair of generalized eigenvalues.
!
!  The left and right Schur vectors are the columns of VSL and VSR,
!  respectively, where VSL and VSR are the orthogonal matrices
!  which reduce A and B to Schur form:
!
!  Schur form of (A,B) = ( (VSL)**T A (VSR), (VSL)**T B (VSR) )
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
!  N       (input) INTEGER
!          The order of the matrices A, B, VSL, and VSR.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the first of the pair of matrices whose generalized
!          eigenvalues and (optionally) Schur vectors are to be
!          computed.
!          On exit, the generalized Schur form of A.
!          Note: to avoid overflow, the Frobenius norm of the matrix
!          A should be less than the overflow threshold.
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
!          On entry, the second of the pair of matrices whose
!          generalized eigenvalues and (optionally) Schur vectors are
!          to be computed.
!          On exit, the generalized Schur form of B.
!          Note: to avoid overflow, the Frobenius norm of the matrix
!          B should be less than the overflow threshold.
!
!  LDB     (input) INTEGER
!          The leading dimension of B.  LDB >= max(1,N).
!
!  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
!  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
!  BETA    (output) DOUBLE PRECISION array, dimension (N)
!          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
!          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
!          j=1,...,N  and  BETA(j),j=1,...,N  are the diagonals of the
!          complex Schur form (A,B) that would result if the 2-by-2
!          diagonal blocks of the real Schur form of (A,B) were further
!          reduced to triangular form using 2-by-2 complex unitary
!          transformations.  If ALPHAI(j) is zero, then the j-th
!          eigenvalue is real; if positive, then the j-th and (j+1)-st
!          eigenvalues are a complex conjugate pair, with ALPHAI(j+1)
!          negative.
!
!          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
!          may easily over- or underflow, and BETA(j) may even be zero.
!          Thus, the user should avoid naively computing the ratio
!          alpha/beta.  However, ALPHAR and ALPHAI will be always less
!          than and usually comparable with norm(A) in magnitude, and
!          BETA always less than and usually comparable with norm(B).
!
!  VSL     (output) DOUBLE PRECISION array, dimension (LDVSL,N)
!          If JOBVSL = 'V', VSL will contain the left Schur vectors.
!          (See "Purpose", above.)
!          Not referenced if JOBVSL = 'N'.
!
!  LDVSL   (input) INTEGER
!          The leading dimension of the matrix VSL. LDVSL >=1, and
!          if JOBVSL = 'V', LDVSL >= N.
!
!  VSR     (output) DOUBLE PRECISION array, dimension (LDVSR,N)
!          If JOBVSR = 'V', VSR will contain the right Schur vectors.
!          (See "Purpose", above.)
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
!          The dimension of the array WORK.  LWORK >= max(1,4*N).
!          For good performance, LWORK must generally be larger.
!          To compute the optimal value of LWORK, call ILAENV to get
!          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:
!          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR
!          The optimal LWORK is  2*N + N*(NB+1).
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
!                The QZ iteration failed.  (A,B) are not in Schur
!                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
!                be correct for j=INFO+1,...,N.
!          > N:  errors that usually indicate LAPACK problems:
!                =N+1: error return from DGGBAL
!                =N+2: error return from DGEQRF
!                =N+3: error return from DORMQR
!                =N+4: error return from DORGQR
!                =N+5: error return from DGGHRD
!                =N+6: error return from DHGEQZ (other than failed
!                                                iteration)
!                =N+7: error return from DGGBAK (computing VSL)
!                =N+8: error return from DGGBAK (computing VSR)
!                =N+9: error return from DLASCL (various places)
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY
      INTEGER            ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, &
                         IRIGHT, IROWS, ITAU, IWORK, LOPT, LWKMIN, &
                         LWKOPT, NB, NB1, NB2, NB3
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, &
                         SAFMIN, SMLNUM
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, &
                         DLASCL, DLASET, DORGQR, DORMQR, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX
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
!     Test the input arguments
!
      LWKMIN = MAX( 4*N, 1 )
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
      else if ( LDVSL < 1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) then
         INFO = -12
      else if ( LDVSR < 1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) then
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
         LOPT = 2*N + N*( NB+1 )
         WORK( 1 ) = LOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEGS ', -INFO )
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
      SMLNUM = N*SAFMIN / EPS
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
!
      if ( ILASCL ) then
         CALL DLASCL( 'G', -1, -1, ANRM, ANRMTO, N, N, A, LDA, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 9
            RETURN
         end if
      end if
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
!
      if ( ILBSCL ) then
         CALL DLASCL( 'G', -1, -1, BNRM, BNRMTO, N, N, B, LDB, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 9
            RETURN
         end if
      end if
!
!     Permute the matrix to make it more nearly triangular
!     Workspace layout:  (2*N words -- "work..." not actually used)
!        left_permutation, right_permutation, work...
!
      ILEFT = 1
      IRIGHT = N + 1
      IWORK = IRIGHT + N
      CALL DGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), &
                   WORK( IRIGHT ), WORK( IWORK ), IINFO )
      if ( IINFO.NE.0 ) then
         INFO = N + 1
         GO TO 10
      end if
!
!     Reduce B to triangular form, and initialize VSL and/or VSR
!     Workspace layout:  ("work..." must have at least N words)
!        left_permutation, right_permutation, tau, work...
!
      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWORK
      IWORK = ITAU + IROWS
      CALL DGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), &
                   WORK( IWORK ), LWORK+1-IWORK, IINFO )
      if ( IINFO.GE.0 ) &
         LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) then
         INFO = N + 2
         GO TO 10
      end if
!
      CALL DORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, &
                   WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), &
                   LWORK+1-IWORK, IINFO )
      if ( IINFO.GE.0 ) &
         LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) then
         INFO = N + 3
         GO TO 10
      end if
!
      if ( ILVSL ) then
         CALL DLASET( 'Full', N, N, ZERO, ONE, VSL, LDVSL )
         CALL DLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, &
                      VSL( ILO+1, ILO ), LDVSL )
         CALL DORGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, &
                      WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, &
                      IINFO )
         if ( IINFO.GE.0 ) &
            LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         if ( IINFO.NE.0 ) then
            INFO = N + 4
            GO TO 10
         end if
      end if
!
      if ( ILVSR ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, VSR, LDVSR )
!
!     Reduce to generalized Hessenberg form
!
      CALL DGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, &
                   LDVSL, VSR, LDVSR, IINFO )
      if ( IINFO.NE.0 ) then
         INFO = N + 5
         GO TO 10
      end if
!
!     Perform QZ algorithm, computing Schur vectors if desired
!     Workspace layout:  ("work..." must have at least 1 word)
!        left_permutation, right_permutation, work...
!
      IWORK = ITAU
      CALL DHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, &
                   ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, &
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
         GO TO 10
      end if
!
!     Apply permutation to VSL and VSR
!
      if ( ILVSL ) then
         CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), &
                      WORK( IRIGHT ), N, VSL, LDVSL, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 7
            GO TO 10
         end if
      end if
      if ( ILVSR ) then
         CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), &
                      WORK( IRIGHT ), N, VSR, LDVSR, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 8
            GO TO 10
         end if
      end if
!
!     Undo scaling
!
      if ( ILASCL ) then
         CALL DLASCL( 'H', -1, -1, ANRMTO, ANRM, N, N, A, LDA, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 9
            RETURN
         end if
         CALL DLASCL( 'G', -1, -1, ANRMTO, ANRM, N, 1, ALPHAR, N, &
                      IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 9
            RETURN
         end if
         CALL DLASCL( 'G', -1, -1, ANRMTO, ANRM, N, 1, ALPHAI, N, &
                      IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 9
            RETURN
         end if
      end if
!
      if ( ILBSCL ) then
         CALL DLASCL( 'U', -1, -1, BNRMTO, BNRM, N, N, B, LDB, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 9
            RETURN
         end if
         CALL DLASCL( 'G', -1, -1, BNRMTO, BNRM, N, 1, BETA, N, IINFO )
         if ( IINFO.NE.0 ) then
            INFO = N + 9
            RETURN
         end if
      end if
!
   10 CONTINUE
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of DGEGS
!
      END
