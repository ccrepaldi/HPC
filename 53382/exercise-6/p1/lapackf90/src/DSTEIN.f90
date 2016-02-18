      SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
                         IWORK, IFAIL, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), &
                         IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEIN computes the eigenvectors of a real symmetric tridiagonal
!  matrix T corresponding to specified eigenvalues, using inverse
!  iteration.
!
!  The maximum number of iterations allowed for each eigenvector is
!  specified by an internal parameter MAXITS (currently set to 5).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N)
!          The (n-1) subdiagonal elements of the tridiagonal matrix
!          T, in elements 1 to N-1.  E(N) need not be set.
!
!  M       (input) INTEGER
!          The number of eigenvectors to be found.  0 <= M <= N.
!
!  W       (input) DOUBLE PRECISION array, dimension (N)
!          The first M elements of W contain the eigenvalues for
!          which eigenvectors are to be computed.  The eigenvalues
!          should be grouped by split-off block and ordered from
!          smallest to largest within the block.  ( The output array
!          W from DSTEBZ with ORDER = 'B' is expected here. )
!
!  IBLOCK  (input) INTEGER array, dimension (N)
!          The submatrix indices associated with the corresponding
!          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
!          the first submatrix from the top, =2 if W(i) belongs to
!          the second submatrix, etc.  ( The output array IBLOCK
!          from DSTEBZ is expected here. )
!
!  ISPLIT  (input) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to
!          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!          through ISPLIT( 2 ), etc.
!          ( The output array ISPLIT from DSTEBZ is expected here. )
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)
!          The computed eigenvectors.  The eigenvector associated
!          with the eigenvalue W(i) is stored in the i-th column of
!          Z.  Any vector which fails to converge is set to its current
!          iterate after MAXITS iterations.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  IFAIL   (output) INTEGER array, dimension (M)
!          On normal exit, all elements of IFAIL are zero.
!          If one or more eigenvectors fail to converge after
!          MAXITS iterations, then their indices are stored in
!          array IFAIL.
!
!  INFO    (output) INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, then i eigenvectors failed to converge
!               in MAXITS iterations.  Their indices are stored in
!               array IFAIL.
!
!  Internal Parameters
!  ===================
!
!  MAXITS  INTEGER, default = 5
!          The maximum number of iterations performed.
!
!  EXTRA   INTEGER, default = 2
!          The number of iterations performed after norm growth
!          criterion is satisfied, should be at least 1.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TEN, ODM3, ODM1
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 1.0D+1, &
                         ODM3 = 1.0D-3, ODM1 = 1.0D-1 )
      INTEGER            MAXITS, EXTRA
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
!     ..
!     .. Local Scalars ..
      INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1, &
                         INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1, &
                         JBLK, JMAX, NBLK, NRMCHK
      DOUBLE PRECISION   DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL, &
                         SCL, SEP, TOL, XJ, XJM, ZTR
!     ..
!     .. Local Arrays ..
      INTEGER            ISEED( 4 )
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH, DNRM2
      EXTERNAL           IDAMAX, DASUM, DDOT, DLAMCH, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DLAGTF, DLAGTS, DLARNV, DSCAL, &
                         XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE
!
      if ( N < 0 ) then
         INFO = -1
      else if ( M < 0 .OR. M.GT.N ) then
         INFO = -4
      else if ( LDZ < MAX( 1, N ) ) then
         INFO = -9
      ELSE
         DO 20 J = 2, M
            if ( IBLOCK( J ) < IBLOCK( J-1 ) ) then
               INFO = -6
               GO TO 30
            end if
            if ( IBLOCK( J ) == IBLOCK( J-1 ) .AND. W( J ) < W( J-1 ) ) &
                 THEN
               INFO = -5
               GO TO 30
            end if
   20    CONTINUE
   30    CONTINUE
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSTEIN', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 .OR. M.EQ.0 ) then
         RETURN
      else if ( N == 1 ) then
         Z( 1, 1 ) = ONE
         RETURN
      end if
!
!     Get machine constants.
!
      EPS = DLAMCH( 'Precision' )
!
!     Initialize seed for random number generator DLARNV.
!
      DO 40 I = 1, 4
         ISEED( I ) = 1
   40 CONTINUE
!
!     Initialize pointers.
!
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
!
!     Compute eigenvectors of matrix blocks.
!
      J1 = 1
      DO 160 NBLK = 1, IBLOCK( M )
!
!        Find starting and ending indices of block nblk.
!
         if ( NBLK == 1 ) then
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         end if
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         if ( BLKSIZ == 1 ) &
            GO TO 60
         GPIND = B1
!
!        Compute reorthogonalization criterion and stopping criterion.
!
         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+ &
                     ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM
!
         DTPCRT = SQRT( ODM1 / BLKSIZ )
!
!        Loop through eigenvalues of block nblk.
!
   60    CONTINUE
         JBLK = 0
         DO 150 J = J1, M
            if ( IBLOCK( J ).NE.NBLK ) then
               J1 = J
               GO TO 160
            end if
            JBLK = JBLK + 1
            XJ = W( J )
!
!           Skip all the work if the block size is one.
!
            if ( BLKSIZ == 1 ) then
               WORK( INDRV1+1 ) = ONE
               GO TO 120
            end if
!
!           If eigenvalues j and j-1 are too close, add a relatively
!           small perturbation.
!
            if ( JBLK.GT.1 ) then
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               if ( SEP < PERTOL ) &
                  XJ = XJM + PERTOL
            end if
!
            ITS = 0
            NRMCHK = 0
!
!           Get random starting vector.
!
            CALL DLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
!
!           Copy the matrix T so it won't be destroyed in factorization.
!
            CALL DCOPY( BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 )
!
!           Compute LU factors with partial pivoting  ( PT = LU )
!
            TOL = ZERO
            CALL DLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ), &
                         WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK, &
                         IINFO )
!
!           Update iteration count.
!
   70       CONTINUE
            ITS = ITS + 1
            if ( ITS.GT.MAXITS ) &
               GO TO 100
!
!           Normalize and scale the righthand side vector Pb.
!
            SCL = BLKSIZ*ONENRM*MAX( EPS, &
                  ABS( WORK( INDRV4+BLKSIZ ) ) ) / &
                  DASUM( BLKSIZ, WORK( INDRV1+1 ), 1 )
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
!
!           Solve the system LU = Pb.
!
            CALL DLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ), &
                         WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK, &
                         WORK( INDRV1+1 ), TOL, IINFO )
!
!           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
!           close enough.
!
            if ( JBLK == 1 ) &
               GO TO 90
            if ( ABS( XJ-XJM ).GT.ORTOL ) &
               GPIND = J
            if ( GPIND.NE.J ) then
               DO 80 I = GPIND, J - 1
                  ZTR = -DDOT( BLKSIZ, WORK( INDRV1+1 ), 1, Z( B1, I ), &
                        1 )
                  CALL DAXPY( BLKSIZ, ZTR, Z( B1, I ), 1, &
                              WORK( INDRV1+1 ), 1 )
   80          CONTINUE
            end if
!
!           Check the infinity norm of the iterate.
!
   90       CONTINUE
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            NRM = ABS( WORK( INDRV1+JMAX ) )
!
!           Continue for additional iterations after norm reaches
!           stopping criterion.
!
            if ( NRM < DTPCRT ) &
               GO TO 70
            NRMCHK = NRMCHK + 1
            if ( NRMCHK < EXTRA+1 ) &
               GO TO 70
!
            GO TO 110
!
!           If stopping criterion was not satisfied, update info and
!           store eigenvector number in array ifail.
!
  100       CONTINUE
            INFO = INFO + 1
            IFAIL( INFO ) = J
!
!           Accept iterate as jth eigenvector.
!
  110       CONTINUE
            SCL = ONE / DNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            if ( WORK( INDRV1+JMAX ) < ZERO ) &
               SCL = -SCL
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  120       CONTINUE
            DO 130 I = 1, N
               Z( I, J ) = ZERO
  130       CONTINUE
            DO 140 I = 1, BLKSIZ
               Z( B1+I-1, J ) = WORK( INDRV1+I )
  140       CONTINUE
!
!           Save the shift to check eigenvalue spacing at next
!           iteration.
!
            XJM = XJ
!
  150    CONTINUE
  160 CONTINUE
!
      RETURN
!
!     End of DSTEIN
!
      END
