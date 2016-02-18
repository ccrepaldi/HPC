      SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELSS computes the minimum norm solution to a real linear least
!  squares problem:
!
!  Minimize 2-norm(| b - A*x |).
!
!  using the singular value decomposition (SVD) of A. A is an M-by-N
!  matrix which may be rank-deficient.
!
!  Several right hand side vectors b and solution vectors x can be
!  handled in a single call; they are stored as the columns of the
!  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
!  X.
!
!  The effective rank of A is determined by treating as zero those
!  singular values which are less than RCOND times the largest singular
!  value.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A. N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X. NRHS >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the first min(m,n) rows of A are overwritten with
!          its right singular vectors, stored rowwise.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the M-by-NRHS right hand side matrix B.
!          On exit, B is overwritten by the N-by-NRHS solution
!          matrix X.  If m >= n and RANK = n, the residual
!          sum-of-squares for the solution in the i-th column is given
!          by the sum of squares of elements n+1:m in that column.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= max(1,max(M,N)).
!
!  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A in decreasing order.
!          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
!
!  RCOND   (input) DOUBLE PRECISION
!          RCOND is used to determine the effective rank of A.
!          Singular values S(i) <= RCOND*S(1) are treated as zero.
!          If RCOND < 0, machine precision is used instead.
!
!  RANK    (output) INTEGER
!          The effective rank of A, i.e., the number of singular values
!          which are greater than RCOND*S(1).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= 1, and also:
!          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
!          For good performance, LWORK should generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the algorithm for computing the SVD failed to converge;
!                if INFO = i, i off-diagonal elements of an intermediate
!                bidiagonal form did not converge to zero.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            BDSPAC, BL, CHUNK, I, IASCL, IBSCL, IE, IL, &
                         ITAU, ITAUP, ITAUQ, IWORK, LDWORK, MAXMN, &
                         MAXWRK, MINMN, MINWRK, MM, MNTHR
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   VDUM( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DBDSQR, DCOPY, DGEBRD, DGELQF, DGEMM, DGEMV, &
                         DGEQRF, DLABAD, DLACPY, DLASCL, DLASET, DORGBR, &
                         DORMBR, DORMLQ, DORMQR, DRSCL, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           ILAENV, DLAMCH, DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      MNTHR = ILAENV( 6, 'DGELSS', ' ', M, N, NRHS, -1 )
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -5
      else if ( LDB < MAX( 1, MAXMN ) ) then
         INFO = -7
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
         MAXWRK = 0
         MM = M
         if ( M.GE.N .AND. M.GE.MNTHR ) then
!
!           Path 1a - overdetermined, with many more rows than columns
!
            MM = N
            MAXWRK = MAX( MAXWRK, N+N*ILAENV( 1, 'DGEQRF', ' ', M, N, &
                     -1, -1 ) )
            MAXWRK = MAX( MAXWRK, N+NRHS* &
                     ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N, -1 ) )
         end if
         if ( M.GE.N ) then
!
!           Path 1 - overdetermined or exactly determined
!
!           Compute workspace needed for DBDSQR
!
            BDSPAC = MAX( 1, 5*N )
            MAXWRK = MAX( MAXWRK, 3*N+( MM+N )* &
                     ILAENV( 1, 'DGEBRD', ' ', MM, N, -1, -1 ) )
            MAXWRK = MAX( MAXWRK, 3*N+NRHS* &
                     ILAENV( 1, 'DORMBR', 'QLT', MM, NRHS, N, -1 ) )
            MAXWRK = MAX( MAXWRK, 3*N+( N-1 )* &
                     ILAENV( 1, 'DORGBR', 'P', N, N, N, -1 ) )
            MAXWRK = MAX( MAXWRK, BDSPAC )
            MAXWRK = MAX( MAXWRK, N*NRHS )
            MINWRK = MAX( 3*N+MM, 3*N+NRHS, BDSPAC )
            MAXWRK = MAX( MINWRK, MAXWRK )
         end if
         if ( N.GT.M ) then
!
!           Compute workspace needed for DBDSQR
!
            BDSPAC = MAX( 1, 5*M )
            MINWRK = MAX( 3*M+NRHS, 3*M+N, BDSPAC )
            if ( N.GE.MNTHR ) then
!
!              Path 2a - underdetermined, with many more columns
!              than rows
!
               MAXWRK = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
               MAXWRK = MAX( MAXWRK, M*M+4*M+2*M* &
                        ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS* &
                        ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )* &
                        ILAENV( 1, 'DORGBR', 'P', M, M, M, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+M+BDSPAC )
               if ( NRHS.GT.1 ) then
                  MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
               ELSE
                  MAXWRK = MAX( MAXWRK, M*M+2*M )
               end if
               MAXWRK = MAX( MAXWRK, M+NRHS* &
                        ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M, -1 ) )
            ELSE
!
!              Path 2 - underdetermined
!
               MAXWRK = 3*M + ( N+M )*ILAENV( 1, 'DGEBRD', ' ', M, N, &
                        -1, -1 )
               MAXWRK = MAX( MAXWRK, 3*M+NRHS* &
                        ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, 3*M+M* &
                        ILAENV( 1, 'DORGBR', 'P', M, N, M, -1 ) )
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MAXWRK = MAX( MAXWRK, N*NRHS )
            end if
         end if
         MAXWRK = MAX( MINWRK, MAXWRK )
         WORK( 1 ) = MAXWRK
      end if
!
      MINWRK = MAX( MINWRK, 1 )
      if ( LWORK < MINWRK .AND. .NOT.LQUERY ) &
         INFO = -12
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGELSS', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 ) then
         RANK = 0
         RETURN
      end if
!
!     Get machine parameters
!
      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM < SMLNUM ) then
!
!        Scale matrix norm up to SMLNUM
!
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      else if ( ANRM.GT.BIGNUM ) then
!
!        Scale matrix norm down to BIGNUM
!
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      else if ( ANRM == ZERO ) then
!
!        Matrix all zero. Return zero solution.
!
         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         CALL DLASET( 'F', MINMN, 1, ZERO, ZERO, S, 1 )
         RANK = 0
         GO TO 70
      end if
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO .AND. BNRM < SMLNUM ) then
!
!        Scale matrix norm up to SMLNUM
!
         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      else if ( BNRM.GT.BIGNUM ) then
!
!        Scale matrix norm down to BIGNUM
!
         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      end if
!
!     Overdetermined case
!
      if ( M.GE.N ) then
!
!        Path 1 - overdetermined or exactly determined
!
         MM = M
         if ( M.GE.MNTHR ) then
!
!           Path 1a - overdetermined, with many more rows than columns
!
            MM = N
            ITAU = 1
            IWORK = ITAU + N
!
!           Compute A=Q*R
!           (Workspace: need 2*N, prefer N+N*NB)
!
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                         LWORK-IWORK+1, INFO )
!
!           Multiply B by transpose(Q)
!           (Workspace: need N+NRHS, prefer N+NRHS*NB)
!
            CALL DORMQR( 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAU ), B, &
                         LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!
!           Zero out below R
!
            if ( N.GT.1 ) &
               CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
         end if
!
         IE = 1
         ITAUQ = IE + N
         ITAUP = ITAUQ + N
         IWORK = ITAUP + N
!
!        Bidiagonalize R in A
!        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
!
         CALL DGEBRD( MM, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                      WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                      INFO )
!
!        Multiply B by transpose of left bidiagonalizing vectors of R
!        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
!
         CALL DORMBR( 'Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK( ITAUQ ), &
                      B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!
!        Generate right bidiagonalizing vectors of R in A
!        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
!
         CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                      WORK( IWORK ), LWORK-IWORK+1, INFO )
         IWORK = IE + N
!
!        Perform bidiagonal QR iteration
!          multiply B by transpose of left singular vectors
!          compute right singular vectors in A
!        (Workspace: need BDSPAC)
!
         CALL DBDSQR( 'U', N, N, 0, NRHS, S, WORK( IE ), A, LDA, VDUM, &
                      1, B, LDB, WORK( IWORK ), INFO )
         if ( INFO.NE.0 ) &
            GO TO 70
!
!        Multiply B by reciprocals of singular values
!
         THR = MAX( RCOND*S( 1 ), SFMIN )
         if ( RCOND < ZERO ) &
            THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO I = 1, N
            if ( S( I ).GT.THR ) then
               CALL DRSCL( NRHS, S( I ), B( I, 1 ), LDB )
               RANK = RANK + 1
            ELSE
               CALL DLASET( 'F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            end if
         end do
!
!        Multiply B by right singular vectors
!        (Workspace: need N, prefer N*NRHS)
!
         if ( LWORK.GE.LDB*NRHS .AND. NRHS.GT.1 ) then
            CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, A, LDA, B, LDB, ZERO, &
                        WORK, LDB )
            CALL DLACPY( 'G', N, NRHS, WORK, LDB, B, LDB )
         else if ( NRHS.GT.1 ) then
            CHUNK = LWORK / N
            DO 20 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               CALL DGEMM( 'T', 'N', N, BL, N, ONE, A, LDA, B( 1, I ), &
                           LDB, ZERO, WORK, N )
               CALL DLACPY( 'G', N, BL, WORK, N, B( 1, I ), LDB )
   20       CONTINUE
         ELSE
            CALL DGEMV( 'T', N, N, ONE, A, LDA, B, 1, ZERO, WORK, 1 )
            CALL DCOPY( N, WORK, 1, B, 1 )
         end if
!
      else if ( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+ &
               MAX( M, 2*M-4, NRHS, N-3*M ) ) then
!
!        Path 2a - underdetermined, with many more columns than rows
!        and sufficient workspace for an efficient algorithm
!
         LDWORK = M
         if ( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ), &
             M*LDA+M+M*NRHS ) )LDWORK = LDA
         ITAU = 1
         IWORK = M + 1
!
!        Compute A=L*Q
!        (Workspace: need 2*M, prefer M+M*NB)
!
         CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                      LWORK-IWORK+1, INFO )
         IL = IWORK
!
!        Copy L to WORK(IL), zeroing out above it
!
         CALL DLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
         CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IL+LDWORK ), &
                      LDWORK )
         IE = IL + LDWORK*M
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M
!
!        Bidiagonalize L in WORK(IL)
!        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
!
         CALL DGEBRD( M, M, WORK( IL ), LDWORK, S, WORK( IE ), &
                      WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), &
                      LWORK-IWORK+1, INFO )
!
!        Multiply B by transpose of left bidiagonalizing vectors of L
!        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
!
         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, M, WORK( IL ), LDWORK, &
                      WORK( ITAUQ ), B, LDB, WORK( IWORK ), &
                      LWORK-IWORK+1, INFO )
!
!        Generate right bidiagonalizing vectors of R in WORK(IL)
!        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)
!
         CALL DORGBR( 'P', M, M, M, WORK( IL ), LDWORK, WORK( ITAUP ), &
                      WORK( IWORK ), LWORK-IWORK+1, INFO )
         IWORK = IE + M
!
!        Perform bidiagonal QR iteration,
!           computing right singular vectors of L in WORK(IL) and
!           multiplying B by transpose of left singular vectors
!        (Workspace: need M*M+M+BDSPAC)
!
         CALL DBDSQR( 'U', M, M, 0, NRHS, S, WORK( IE ), WORK( IL ), &
                      LDWORK, A, LDA, B, LDB, WORK( IWORK ), INFO )
         if ( INFO.NE.0 ) &
            GO TO 70
!
!        Multiply B by reciprocals of singular values
!
         THR = MAX( RCOND*S( 1 ), SFMIN )
         if ( RCOND < ZERO ) &
            THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 30 I = 1, M
            if ( S( I ).GT.THR ) then
               CALL DRSCL( NRHS, S( I ), B( I, 1 ), LDB )
               RANK = RANK + 1
            ELSE
               CALL DLASET( 'F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            end if
   30    CONTINUE
         IWORK = IE
!
!        Multiply B by right singular vectors of L in WORK(IL)
!        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)
!
         if ( LWORK.GE.LDB*NRHS+IWORK-1 .AND. NRHS.GT.1 ) then
            CALL DGEMM( 'T', 'N', M, NRHS, M, ONE, WORK( IL ), LDWORK, &
                        B, LDB, ZERO, WORK( IWORK ), LDB )
            CALL DLACPY( 'G', M, NRHS, WORK( IWORK ), LDB, B, LDB )
         else if ( NRHS.GT.1 ) then
            CHUNK = ( LWORK-IWORK+1 ) / M
            DO 40 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               CALL DGEMM( 'T', 'N', M, BL, M, ONE, WORK( IL ), LDWORK, &
                           B( 1, I ), LDB, ZERO, WORK( IWORK ), N )
               CALL DLACPY( 'G', M, BL, WORK( IWORK ), N, B( 1, I ), &
                            LDB )
   40       CONTINUE
         ELSE
            CALL DGEMV( 'T', M, M, ONE, WORK( IL ), LDWORK, B( 1, 1 ), &
                        1, ZERO, WORK( IWORK ), 1 )
            CALL DCOPY( M, WORK( IWORK ), 1, B( 1, 1 ), 1 )
         end if
!
!        Zero out below first M rows of B
!
         CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
         IWORK = ITAU + M
!
!        Multiply transpose(Q) by B
!        (Workspace: need M+NRHS, prefer M+NRHS*NB)
!
         CALL DORMLQ( 'L', 'T', N, NRHS, M, A, LDA, WORK( ITAU ), B, &
                      LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!
      ELSE
!
!        Path 2 - remaining underdetermined cases
!
         IE = 1
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M
!
!        Bidiagonalize A
!        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
!
         CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                      WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                      INFO )
!
!        Multiply B by transpose of left bidiagonalizing vectors
!        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
!
         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAUQ ), &
                      B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!
!        Generate right bidiagonalizing vectors in A
!        (Workspace: need 4*M, prefer 3*M+M*NB)
!
         CALL DORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ), &
                      WORK( IWORK ), LWORK-IWORK+1, INFO )
         IWORK = IE + M
!
!        Perform bidiagonal QR iteration,
!           computing right singular vectors of A in A and
!           multiplying B by transpose of left singular vectors
!        (Workspace: need BDSPAC)
!
         CALL DBDSQR( 'L', M, N, 0, NRHS, S, WORK( IE ), A, LDA, VDUM, &
                      1, B, LDB, WORK( IWORK ), INFO )
         if ( INFO.NE.0 ) &
            GO TO 70
!
!        Multiply B by reciprocals of singular values
!
         THR = MAX( RCOND*S( 1 ), SFMIN )
         if ( RCOND < ZERO ) &
            THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 50 I = 1, M
            if ( S( I ).GT.THR ) then
               CALL DRSCL( NRHS, S( I ), B( I, 1 ), LDB )
               RANK = RANK + 1
            ELSE
               CALL DLASET( 'F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            end if
   50    CONTINUE
!
!        Multiply B by right singular vectors of A
!        (Workspace: need N, prefer N*NRHS)
!
         if ( LWORK.GE.LDB*NRHS .AND. NRHS.GT.1 ) then
            CALL DGEMM( 'T', 'N', N, NRHS, M, ONE, A, LDA, B, LDB, ZERO, &
                        WORK, LDB )
            CALL DLACPY( 'F', N, NRHS, WORK, LDB, B, LDB )
         else if ( NRHS.GT.1 ) then
            CHUNK = LWORK / N
            DO 60 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               CALL DGEMM( 'T', 'N', N, BL, M, ONE, A, LDA, B( 1, I ), &
                           LDB, ZERO, WORK, N )
               CALL DLACPY( 'F', N, BL, WORK, N, B( 1, I ), LDB )
   60       CONTINUE
         ELSE
            CALL DGEMV( 'T', M, N, ONE, A, LDA, B, 1, ZERO, WORK, 1 )
            CALL DCOPY( N, WORK, 1, B, 1 )
         end if
      end if
!
!     Undo scaling
!
      if ( IASCL == 1 ) then
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, &
                      INFO )
      else if ( IASCL == 2 ) then
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, &
                      INFO )
      end if
      if ( IBSCL == 1 ) then
         CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      else if ( IBSCL == 2 ) then
         CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      end if
!
   70 CONTINUE
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of DGELSS
!
      END
