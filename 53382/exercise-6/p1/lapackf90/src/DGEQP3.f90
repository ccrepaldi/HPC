      SUBROUTINE DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEQP3 computes a QR factorization with column pivoting of a
!  matrix A:  A*P = Q*R  using Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the upper triangle of the array contains the
!          min(M,N)-by-N upper trapezoidal matrix R; the elements below
!          the diagonal, together with the array TAU, represent the
!          orthogonal matrix Q as a product of min(M,N) elementary
!          reflectors.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  JPVT    (input/output) INTEGER array, dimension (N)
!          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
!          to the front of A*P (a leading column); if JPVT(J)=0,
!          the J-th column of A is a free column.
!          On exit, if JPVT(J)=K, then the J-th column of A*P was the
!          the K-th column of A.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= 3*N+1.
!          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
!          is the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real/complex scalar, and v is a real/complex vector
!  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
!  A(i+1:m,i), and tau in TAU(i).
!
!  Based on contributions by
!    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!    X. Sun, Computer Science Dept., Duke University, USA
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            INB, INBMIN, IXOVER
      PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB, &
                         NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQRF, DLAQP2, DLAQPS, DORMQR, DSWAP, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DNRM2
      EXTERNAL           ILAENV, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
!     ..
!     .. Executable Statements ..
!
      IWS = 3*N + 1
      MINMN = MIN( M, N )
!
!     Test input arguments
!     ====================
!
      INFO = 0
      NB = ILAENV( INB, 'DGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = 2*N + ( N+1 )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -4
      else if ( ( LWORK < IWS ) .AND. .NOT.LQUERY ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEQP3', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible.
!
      if ( MINMN == 0 ) then
         WORK( 1 ) = 1
         RETURN
      end if
!
!     Move initial columns up front.
!
      NFXD = 1
      DO 10 J = 1, N
         if ( JPVT( J ).NE.0 ) then
            if ( J.NE.NFXD ) then
               CALL DSWAP( M, A( 1, J ), 1, A( 1, NFXD ), 1 )
               JPVT( J ) = JPVT( NFXD )
               JPVT( NFXD ) = J
            ELSE
               JPVT( J ) = J
            end if
            NFXD = NFXD + 1
         ELSE
            JPVT( J ) = J
         end if
   10 CONTINUE
      NFXD = NFXD - 1
!
!     Factorize fixed columns
!     =======================
!
!     Compute the QR factorization of fixed columns and update
!     remaining columns.
!
      if ( NFXD.GT.0 ) then
         NA = MIN( M, NFXD )
!CC      CALL DGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         CALL DGEQRF( M, NA, A, LDA, TAU, WORK, LWORK, INFO )
         IWS = MAX( IWS, INT( WORK( 1 ) ) )
         if ( NA < N ) then
!CC         CALL DORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
!CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            CALL DORMQR( 'Left', 'Transpose', M, N-NA, NA, A, LDA, TAU, &
                         A( 1, NA+1 ), LDA, WORK, LWORK, INFO )
            IWS = MAX( IWS, INT( WORK( 1 ) ) )
         end if
      end if
!
!     Factorize free columns
!     ======================
!
      if ( NFXD < MINMN ) then
!
         SM = M - NFXD
         SN = N - NFXD
         SMINMN = MINMN - NFXD
!
!        Determine the block size.
!
         NB = ILAENV( INB, 'DGEQRF', ' ', SM, SN, -1, -1 )
         NBMIN = 2
         NX = 0
!
         if ( ( NB.GT.1 ) .AND. ( NB < SMINMN ) ) then
!
!           Determine when to cross over from blocked to unblocked code.
!
            NX = MAX( 0, ILAENV( IXOVER, 'DGEQRF', ' ', SM, SN, -1, &
                 -1 ) )
!
!
            if ( NX < SMINMN ) then
!
!              Determine if workspace is large enough for blocked code.
!
               MINWS = 2*SN + ( SN+1 )*NB
               IWS = MAX( IWS, MINWS )
               if ( LWORK < MINWS ) then
!
!                 Not enough workspace to use optimal NB: Reduce NB and
!                 determine the minimum value of NB.
!
                  NB = ( LWORK-2*SN ) / ( SN+1 )
                  NBMIN = MAX( 2, ILAENV( INBMIN, 'DGEQRF', ' ', SM, SN, &
                          -1, -1 ) )
!
!
               end if
            end if
         end if
!
!        Initialize partial column norms. The first N elements of work
!        store the exact column norms.
!
         DO 20 J = NFXD + 1, N
            WORK( J ) = DNRM2( SM, A( NFXD+1, J ), 1 )
            WORK( N+J ) = WORK( J )
   20    CONTINUE
!
         if ( ( NB.GE.NBMIN ) .AND. ( NB < SMINMN ) .AND. &
             ( NX < SMINMN ) ) then
!
!           Use blocked code initially.
!
            J = NFXD + 1
!
!           Compute factorization: while loop.
!
!
            TOPBMN = MINMN - NX
   30       CONTINUE
            if ( J.LE.TOPBMN ) then
               JB = MIN( NB, TOPBMN-J+1 )
!
!              Factorize JB columns among columns J:N.
!
               CALL DLAQPS( M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA, &
                            JPVT( J ), TAU( J ), WORK( J ), WORK( N+J ), &
                            WORK( 2*N+1 ), WORK( 2*N+JB+1 ), N-J+1 )
!
               J = J + FJB
               GO TO 30
            end if
         ELSE
            J = NFXD + 1
         end if
!
!        Use unblocked code to factor the last or only block.
!
!
         if ( J.LE.MINMN ) &
            CALL DLAQP2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ), &
                         TAU( J ), WORK( J ), WORK( N+J ), &
                         WORK( 2*N+1 ) )
!
      end if
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DGEQP3
!
      END
