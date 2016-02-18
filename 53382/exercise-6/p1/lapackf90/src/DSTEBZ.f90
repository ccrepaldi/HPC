      SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, &
                         M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, &
                         INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEBZ computes the eigenvalues of a symmetric tridiagonal
!  matrix T.  The user may ask for all eigenvalues, all eigenvalues
!  in the half-open interval (VL, VU], or the IL-th through IU-th
!  eigenvalues.
!
!  To avoid overflow, the matrix must be scaled so that its
!  largest element is no greater than overflow**(1/2) *
!  underflow**(1/4) in absolute value, and for greatest
!  accuracy, it should not be much smaller than that.
!
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!  Matrix", Report CS41, Computer Science Dept., Stanford
!  University, July 21, 1966.
!
!  Arguments
!  =========
!
!  RANGE   (input) CHARACTER
!          = 'A': ("All")   all eigenvalues will be found.
!          = 'V': ("Value") all eigenvalues in the half-open interval
!                           (VL, VU] will be found.
!          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
!                           entire matrix) will be found.
!
!  ORDER   (input) CHARACTER
!          = 'B': ("By Block") the eigenvalues will be grouped by
!                              split-off block (see IBLOCK, ISPLIT) and
!                              ordered from smallest to largest within
!                              the block.
!          = 'E': ("Entire matrix")
!                              the eigenvalues for the entire matrix
!                              will be ordered from smallest to
!                              largest.
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix T.  N >= 0.
!
!  VL      (input) DOUBLE PRECISION
!  VU      (input) DOUBLE PRECISION
!          If RANGE='V', the lower and upper bounds of the interval to
!          be searched for eigenvalues.  Eigenvalues less than or equal
!          to VL, or greater than VU, will not be returned.  VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!
!  IL      (input) INTEGER
!  IU      (input) INTEGER
!          If RANGE='I', the indices (in ascending order) of the
!          smallest and largest eigenvalues to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The absolute tolerance for the eigenvalues.  An eigenvalue
!          (or cluster) is considered to be located if it has been
!          determined to lie in an interval whose width is ABSTOL or
!          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
!          will be used, where |T| means the 1-norm of T.
!
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) off-diagonal elements of the tridiagonal matrix T.
!
!  M       (output) INTEGER
!          The actual number of eigenvalues found. 0 <= M <= N.
!          (See also the description of INFO=2,3.)
!
!  NSPLIT  (output) INTEGER
!          The number of diagonal blocks in the matrix T.
!          1 <= NSPLIT <= N.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          On exit, the first M elements of W will contain the
!          eigenvalues.  (DSTEBZ may use the remaining N-M elements as
!          workspace.)
!
!  IBLOCK  (output) INTEGER array, dimension (N)
!          At each row/column j where E(j) is zero or small, the
!          matrix T is considered to split into a block diagonal
!          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
!          block (from 1 to the number of blocks) the eigenvalue W(i)
!          belongs.  (DSTEBZ may use the remaining N-M elements as
!          workspace.)
!
!  ISPLIT  (output) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to ISPLIT(1),
!          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
!          etc., and the NSPLIT-th consists of rows/columns
!          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
!          (Only the first NSPLIT elements will actually be used, but
!          since the user cannot know a priori what value NSPLIT will
!          have, N words must be reserved for ISPLIT.)
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
!
!  IWORK   (workspace) INTEGER array, dimension (3*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  some or all of the eigenvalues failed to converge or
!                were not computed:
!                =1 or 3: Bisection failed to converge for some
!                        eigenvalues; these eigenvalues are flagged by a
!                        negative block number.  The effect is that the
!                        eigenvalues may not be as accurate as the
!                        absolute and relative tolerances.  This is
!                        generally caused by unexpectedly inaccurate
!                        arithmetic.
!                =2 or 3: RANGE='I' only: Not all of the eigenvalues
!                        IL:IU were found.
!                        Effect: M < IU+1-IL
!                        Cause:  non-monotonic arithmetic, causing the
!                                Sturm sequence to be non-monotonic.
!                        Cure:   recalculate, using RANGE='A', and pick
!                                out eigenvalues IL:IU.  In some cases,
!                                increasing the PARAMETER "FUDGE" may
!                                make things work.
!                = 4:    RANGE='I', and the Gershgorin interval
!                        initially used was too small.  No eigenvalues
!                        were computed.
!                        Probable cause: your machine has sloppy
!                                        floating-point arithmetic.
!                        Cure: Increase the PARAMETER "FUDGE",
!                              recompile, and try again.
!
!  Internal Parameters
!  ===================
!
!  RELFAC  DOUBLE PRECISION, default = 2.0e0
!          The relative tolerance.  An interval (a,b] lies within
!          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
!          where "ulp" is the machine precision (distance from 1 to
!          the next larger floating point number.)
!
!  FUDGE   DOUBLE PRECISION, default = 2
!          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
!          a value of 1 should work, but on machines with sloppy
!          arithmetic, this needs to be larger.  The default for
!          publicly released versions should be large enough to handle
!          the worst machine around.  Note that this has no effect
!          on accuracy of the solution.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         HALF = 1.0D0 / TWO )
      DOUBLE PRECISION   FUDGE, RELFAC
      PARAMETER          ( FUDGE = 2.0D0, RELFAC = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NCNVRG, TOOFEW
      INTEGER            IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, &
                         IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX, &
                         ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL, &
                         NWU
      DOUBLE PRECISION   ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN, &
                         TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
!     ..
!     .. Local Arrays ..
      INTEGER            IDUMMA( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, ILAENV, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAEBZ, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Decode RANGE
!
      if ( LSAME( RANGE, 'A' ) ) then
         IRANGE = 1
      else if ( LSAME( RANGE, 'V' ) ) then
         IRANGE = 2
      else if ( LSAME( RANGE, 'I' ) ) then
         IRANGE = 3
      ELSE
         IRANGE = 0
      end if
!
!     Decode ORDER
!
      if ( LSAME( ORDER, 'B' ) ) then
         IORDER = 2
      else if ( LSAME( ORDER, 'E' ) ) then
         IORDER = 1
      ELSE
         IORDER = 0
      end if
!
!     Check for Errors
!
      if ( IRANGE.LE.0 ) then
         INFO = -1
      else if ( IORDER.LE.0 ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( IRANGE == 2 ) then
         if ( VL.GE.VU ) &
            INFO = -5
      else if ( IRANGE == 3 .AND. ( IL < 1 .OR. IL.GT.MAX( 1, N ) ) ) &
                THEN
         INFO = -6
      else if ( IRANGE == 3 .AND. ( IU < MIN( N, IL ) .OR. IU.GT.N ) ) &
                THEN
         INFO = -7
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSTEBZ', -INFO )
         RETURN
      end if
!
!     Initialize error flags
!
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
!
!     Quick return if possible
!
      M = 0
      if ( N == 0 ) &
         RETURN
!
!     Simplifications:
!
      if ( IRANGE == 3 .AND. IL.EQ.1 .AND. IU.EQ.N ) &
         IRANGE = 1
!
!     Get machine constants
!     NB is the minimum vector length for vector bisection, or 0
!     if only scalar is to be done.
!
      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
      RTOLI = ULP*RELFAC
      NB = ILAENV( 1, 'DSTEBZ', ' ', N, -1, -1, -1 )
      if ( NB.LE.1 ) &
         NB = 0
!
!     Special Case when N=1
!
      if ( N == 1 ) then
         NSPLIT = 1
         ISPLIT( 1 ) = 1
         if ( IRANGE == 2 .AND. ( VL.GE.D( 1 ) .OR. VU < D( 1 ) ) ) then
            M = 0
         ELSE
            W( 1 ) = D( 1 )
            IBLOCK( 1 ) = 1
            M = 1
         end if
         RETURN
      end if
!
!     Compute Splitting Points
!
      NSPLIT = 1
      WORK( N ) = ZERO
      PIVMIN = ONE
!
!DIR$ NOVECTOR
      DO 10 J = 2, N
         TMP1 = E( J-1 )**2
         if ( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN.GT.TMP1 ) then
            ISPLIT( NSPLIT ) = J - 1
            NSPLIT = NSPLIT + 1
            WORK( J-1 ) = ZERO
         ELSE
            WORK( J-1 ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         end if
   10 CONTINUE
      ISPLIT( NSPLIT ) = N
      PIVMIN = PIVMIN*SAFEMN
!
!     Compute Interval and ATOLI
!
      if ( IRANGE == 3 ) then
!
!        RANGE='I': Compute the interval containing eigenvalues
!                   IL through IU.
!
!        Compute Gershgorin interval for entire (split) matrix
!        and use it as the initial interval
!
         GU = D( 1 )
         GL = D( 1 )
         TMP1 = ZERO
!
         DO 20 J = 1, N - 1
            TMP2 = SQRT( WORK( J ) )
            GU = MAX( GU, D( J )+TMP1+TMP2 )
            GL = MIN( GL, D( J )-TMP1-TMP2 )
            TMP1 = TMP2
   20    CONTINUE
!
         GU = MAX( GU, D( N )+TMP1 )
         GL = MIN( GL, D( N )-TMP1 )
         TNORM = MAX( ABS( GL ), ABS( GU ) )
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
!
!        Compute Iteration parameters
!
         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / &
                 LOG( TWO ) ) + 2
         if ( ABSTOL.LE.ZERO ) then
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         end if
!
         WORK( N+1 ) = GL
         WORK( N+2 ) = GL
         WORK( N+3 ) = GU
         WORK( N+4 ) = GU
         WORK( N+5 ) = GL
         WORK( N+6 ) = GU
         IWORK( 1 ) = -1
         IWORK( 2 ) = -1
         IWORK( 3 ) = N + 1
         IWORK( 4 ) = N + 1
         IWORK( 5 ) = IL - 1
         IWORK( 6 ) = IU
!
         CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E, &
                      WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT, &
                      IWORK, W, IBLOCK, IINFO )
!
         if ( IWORK( 6 ) == IU ) then
            WL = WORK( N+1 )
            WLU = WORK( N+3 )
            NWL = IWORK( 1 )
            WU = WORK( N+4 )
            WUL = WORK( N+2 )
            NWU = IWORK( 4 )
         ELSE
            WL = WORK( N+2 )
            WLU = WORK( N+4 )
            NWL = IWORK( 2 )
            WU = WORK( N+3 )
            WUL = WORK( N+1 )
            NWU = IWORK( 3 )
         end if
!
         if ( NWL < 0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) then
            INFO = 4
            RETURN
         end if
      ELSE
!
!        RANGE='A' or 'V' -- Set ATOLI
!
         TNORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ), &
                 ABS( D( N ) )+ABS( E( N-1 ) ) )
!
         DO 30 J = 2, N - 1
            TNORM = MAX( TNORM, ABS( D( J ) )+ABS( E( J-1 ) )+ &
                    ABS( E( J ) ) )
   30    CONTINUE
!
         if ( ABSTOL.LE.ZERO ) then
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         end if
!
         if ( IRANGE == 2 ) then
            WL = VL
            WU = VU
         ELSE
            WL = ZERO
            WU = ZERO
         end if
      end if
!
!     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
!     NWL accumulates the number of eigenvalues .le. WL,
!     NWU accumulates the number of eigenvalues .le. WU
!
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
!
      DO 70 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF
!
         if ( IN == 1 ) then
!
!           Special Case -- IN=1
!
            if ( IRANGE == 1 .OR. WL.GE.D( IBEGIN )-PIVMIN ) &
               NWL = NWL + 1
            if ( IRANGE == 1 .OR. WU.GE.D( IBEGIN )-PIVMIN ) &
               NWU = NWU + 1
            if ( IRANGE == 1 .OR. ( WL < D( IBEGIN )-PIVMIN .AND. WU.GE. &
                D( IBEGIN )-PIVMIN ) ) then
               M = M + 1
               W( M ) = D( IBEGIN )
               IBLOCK( M ) = JB
            end if
         ELSE
!
!           General Case -- IN > 1
!
!           Compute Gershgorin Interval
!           and use it as the initial interval
!
            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO
!
            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, D( J )+TMP1+TMP2 )
               GL = MIN( GL, D( J )-TMP1-TMP2 )
               TMP1 = TMP2
   40       CONTINUE
!
            GU = MAX( GU, D( IEND )+TMP1 )
            GL = MIN( GL, D( IEND )-TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
!
!           Compute ATOLI for the current submatrix
!
            if ( ABSTOL.LE.ZERO ) then
               ATOLI = ULP*MAX( ABS( GL ), ABS( GU ) )
            ELSE
               ATOLI = ABSTOL
            end if
!
            if ( IRANGE.GT.1 ) then
               if ( GU < WL ) then
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GO TO 70
               end if
               GL = MAX( GL, WL )
               GU = MIN( GU, WU )
               if ( GL.GE.GU ) &
                  GO TO 70
            end if
!
!           Set Up Initial Interval
!
            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, &
                         D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), &
                         IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM, &
                         IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
!
            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )
!
!           Compute Eigenvalues
!
            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) / &
                    LOG( TWO ) ) + 2
            CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, &
                         D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), &
                         IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT, &
                         IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
!
!           Copy Eigenvalues Into W and IBLOCK
!           Use -JB for block number for unconverged eigenvalues.
!
            DO 60 J = 1, IOUT
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )
!
!              Flag non-convergence.
!
               if ( J.GT.IOUT-IINFO ) then
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               end if
               DO 50 JE = IWORK( J ) + 1 + IWOFF, &
                       IWORK( J+IN ) + IWOFF
                  W( JE ) = TMP1
                  IBLOCK( JE ) = IB
   50          CONTINUE
   60       CONTINUE
!
            M = M + IM
         end if
   70 CONTINUE
!
!     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
!     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
!
      if ( IRANGE == 3 ) then
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
!
         if ( IDISCL.GT.0 .OR. IDISCU.GT.0 ) then
            DO 80 JE = 1, M
               if ( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) then
                  IDISCL = IDISCL - 1
               else if ( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) then
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               end if
   80       CONTINUE
            M = IM
         end if
         if ( IDISCL.GT.0 .OR. IDISCU.GT.0 ) then
!
!           Code to deal with effects of bad arithmetic:
!           Some low eigenvalues to be discarded are not in (WL,WLU],
!           or high eigenvalues to be discarded are not in (WUL,WU]
!           so just kill off the smallest IDISCL/largest IDISCU
!           eigenvalues, by simply finding the smallest/largest
!           eigenvalue(s).
!
!           (If N(w) is monotone non-decreasing, this should never
!               happen.)
!
            if ( IDISCL.GT.0 ) then
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     if ( IBLOCK( JE ).NE.0 .AND. &
                         ( W( JE ) < WKILL .OR. IW == 0 ) ) then
                        IW = JE
                        WKILL = W( JE )
                     end if
   90             CONTINUE
                  IBLOCK( IW ) = 0
  100          CONTINUE
            end if
            if ( IDISCU.GT.0 ) then
!
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     if ( IBLOCK( JE ).NE.0 .AND. &
                         ( W( JE ).GT.WKILL .OR. IW == 0 ) ) then
                        IW = JE
                        WKILL = W( JE )
                     end if
  110             CONTINUE
                  IBLOCK( IW ) = 0
  120          CONTINUE
            end if
            IM = 0
            DO 130 JE = 1, M
               if ( IBLOCK( JE ).NE.0 ) then
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               end if
  130       CONTINUE
            M = IM
         end if
         if ( IDISCL < 0 .OR. IDISCU.LT.0 ) then
            TOOFEW = .TRUE.
         end if
      end if
!
!     If ORDER='B', do nothing -- the eigenvalues are already sorted
!        by block.
!     If ORDER='E', sort the eigenvalues from smallest to largest
!
      if ( IORDER == 1 .AND. NSPLIT.GT.1 ) then
         DO 150 JE = 1, M - 1
            IE = 0
            TMP1 = W( JE )
            DO 140 J = JE + 1, M
               if ( W( J ) < TMP1 ) then
                  IE = J
                  TMP1 = W( J )
               end if
  140       CONTINUE
!
            if ( IE.NE.0 ) then
               ITMP1 = IBLOCK( IE )
               W( IE ) = W( JE )
               IBLOCK( IE ) = IBLOCK( JE )
               W( JE ) = TMP1
               IBLOCK( JE ) = ITMP1
            end if
  150    CONTINUE
      end if
!
      INFO = 0
      if ( NCNVRG ) &
         INFO = INFO + 1
      if ( TOOFEW ) &
         INFO = INFO + 2
      RETURN
!
!     End of DSTEBZ
!
      END
