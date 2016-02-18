      SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
                         INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEBAK forms the right or left eigenvectors of a real general matrix
!  by backward transformation on the computed eigenvectors of the
!  balanced matrix output by DGEBAL.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the type of backward transformation required:
!          = 'N', do nothing, return immediately;
!          = 'P', do backward transformation for permutation only;
!          = 'S', do backward transformation for scaling only;
!          = 'B', do backward transformations for both permutation and
!                 scaling.
!          JOB must be the same as the argument JOB supplied to DGEBAL.
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  V contains right eigenvectors;
!          = 'L':  V contains left eigenvectors.
!
!  N       (input) INTEGER
!          The number of rows of the matrix V.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          The integers ILO and IHI determined by DGEBAL.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  SCALE   (input) DOUBLE PRECISION array, dimension (N)
!          Details of the permutation and scaling factors, as returned
!          by DGEBAL.
!
!  M       (input) INTEGER
!          The number of columns of the matrix V.  M >= 0.
!
!  V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)
!          On entry, the matrix of right or left eigenvectors to be
!          transformed, as returned by DHSEIN or DTREVC.
!          On exit, V is overwritten by the transformed eigenvectors.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V. LDV >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
!
      INFO = 0
      if ( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) then
         INFO = -1
      else if ( .NOT.RIGHTV .AND. .NOT.LEFTV ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( ILO < 1 .OR. ILO.GT.MAX( 1, N ) ) then
         INFO = -4
      else if ( IHI < MIN( ILO, N ) .OR. IHI.GT.N ) then
         INFO = -5
      else if ( M < 0 ) then
         INFO = -7
      else if ( LDV < MAX( 1, N ) ) then
         INFO = -9
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEBAK', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
      if ( M == 0 ) &
         RETURN
      if ( LSAME( JOB, 'N' ) ) &
         RETURN
!
      if ( ILO == IHI ) &
         GO TO 30
!
!     Backward balance
!
      if ( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) then
!
         if ( RIGHTV ) then
            DO I = ILO, IHI
               S = SCALE( I )
               CALL DSCAL( M, S, V( I, 1 ), LDV )
            end do
         end if

         if ( LEFTV ) then
            DO 20 I = ILO, IHI
               S = ONE / SCALE( I )
               CALL DSCAL( M, S, V( I, 1 ), LDV )
   20       CONTINUE
         end if
!
      end if
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
   30 CONTINUE
      if ( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) then
         if ( RIGHTV ) then
            DO 40 II = 1, N
               I = II
               if ( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 40
               if ( I < ILO ) &
                  I = ILO - II
               K = SCALE( I )
               if ( K == I ) &
                  GO TO 40
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
         end if
!
         if ( LEFTV ) then
            DO 50 II = 1, N
               I = II
               if ( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 50
               if ( I < ILO ) &
                  I = ILO - II
               K = SCALE( I )
               if ( K == I ) &
                  GO TO 50
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   50       CONTINUE
         end if
      end if
!
      RETURN
!
!     End of DGEBAK
!
      END
