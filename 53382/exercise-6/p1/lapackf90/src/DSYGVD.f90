      SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
                         LWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
!  of a real generalized symmetric-definite eigenproblem, of the form
!  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
!  B are assumed to be symmetric and B is also positive definite.
!  If eigenvectors are desired, it uses a divide and conquer algorithm.
!
!  The divide and conquer algorithm makes very mild assumptions about
!  floating point arithmetic. It will work on machines with a guard
!  digit in add/subtract, or on those binary machines without guard
!  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!  Cray-2. It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Arguments
!  =========
!
!  ITYPE   (input) INTEGER
!          Specifies the problem type to be solved:
!          = 1:  A*x = (lambda)*B*x
!          = 2:  A*B*x = (lambda)*x
!          = 3:  B*A*x = (lambda)*x
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangles of A and B are stored;
!          = 'L':  Lower triangles of A and B are stored.
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          matrix Z of eigenvectors.  The eigenvectors are normalized
!          as follows:
!          if ITYPE = 1 or 2, Z**T*B*Z = I;
!          if ITYPE = 3, Z**T*inv(B)*Z = I.
!          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
!          or the lower triangle (if UPLO='L') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
!          On entry, the symmetric matrix B.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of B contains the
!          upper triangular part of the matrix B.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of B contains
!          the lower triangular part of the matrix B.
!
!          On exit, if INFO <= N, the part of B containing the matrix is
!          overwritten by the triangular factor U or L from the Cholesky
!          factorization B = U**T*U or B = L*L**T.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If N <= 1,               LWORK >= 1.
!          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1.
!          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
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
!          The dimension of the array IWORK.
!          If N <= 1,                LIWORK >= 1.
!          If JOBZ  = 'N' and N > 1, LIWORK >= 1.
!          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  DPOTRF or DSYEVD returned an error code:
!             <= N:  if INFO = i, DSYEVD failed to converge;
!                    i off-diagonal elements of an intermediate
!                    tridiagonal form did not converge to zero;
!             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
!                    minor of order i of B is not positive definite.
!                    The factorization of B could not be completed and
!                    no eigenvalues or eigenvectors were computed.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, WANTZ
      CHARACTER          TRANS
      INTEGER            LIOPT, LIWMIN, LOPT, LWMIN, NEIG
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DPOTRF, DSYEVD, DSYGST, DTRMM, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 .OR. LIWORK.EQ.-1 )
!
      INFO = 0
      if ( N.LE.1 ) then
         LIWMIN = 1
         LWMIN = 1
         LOPT = LWMIN
         LIOPT = LIWMIN
      ELSE
         if ( WANTZ ) then
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 6*N + 2*N**2
         ELSE
            LIWMIN = 1
            LWMIN = 2*N + 1
         end if
         LOPT = LWMIN
         LIOPT = LIWMIN
      end if
      if ( ITYPE < 0 .OR. ITYPE.GT.3 ) then
         INFO = -1
      else if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) then
         INFO = -2
      else if ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) then
         INFO = -3
      else if ( N < 0 ) then
         INFO = -4
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -6
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -8
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -11
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -13
      end if
!
      if ( INFO == 0 ) then
         WORK( 1 ) = LOPT
         IWORK( 1 ) = LIOPT
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSYGVD', -INFO )
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
!     Form a Cholesky factorization of B.
!
      CALL DPOTRF( UPLO, N, B, LDB, INFO )
      if ( INFO.NE.0 ) then
         INFO = N + INFO
         RETURN
      end if
!
!     Transform problem to standard eigenvalue problem and solve.
!
      CALL DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      CALL DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, &
                   INFO )
      LOPT = MAX( DBLE( LOPT ), DBLE( WORK( 1 ) ) )
      LIOPT = MAX( DBLE( LIOPT ), DBLE( IWORK( 1 ) ) )
!
      if ( WANTZ ) then
!
!        Backtransform eigenvectors to the original problem.
!
         NEIG = N
         if ( INFO.GT.0 ) &
            NEIG = INFO - 1
         if ( ITYPE == 1 .OR. ITYPE.EQ.2 ) then
!
!           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
!           backtransform eigenvectors: x = inv(L)'*y or inv(U)*y
!
            if ( UPPER ) then
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            end if
!
            CALL DTRSM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, &
                        B, LDB, A, LDA )
!
         else if ( ITYPE == 3 ) then
!
!           For B*A*x=(lambda)*x;
!           backtransform eigenvectors: x = L*y or U'*y
!
            if ( UPPER ) then
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            end if
!
            CALL DTRMM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, &
                        B, LDB, A, LDA )
         end if
      end if
!
      WORK( 1 ) = LOPT
      IWORK( 1 ) = LIOPT
!
      RETURN
!
!     End of DSYGVD
!
      END