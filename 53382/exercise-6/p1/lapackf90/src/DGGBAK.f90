      SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, &
                         LDV, INFO )
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
      DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DGGBAK forms the right or left eigenvectors of a real generalized
!  eigenvalue problem A*x = lambda*B*x, by backward transformation on
!  the computed eigenvectors of the balanced pair of matrices output by
!  DGGBAL.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the type of backward transformation required:
!          = 'N':  do nothing, return immediately;
!          = 'P':  do backward transformation for permutation only;
!          = 'S':  do backward transformation for scaling only;
!          = 'B':  do backward transformations for both permutation and
!                  scaling.
!          JOB must be the same as the argument JOB supplied to DGGBAL.
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
!          The integers ILO and IHI determined by DGGBAL.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  LSCALE  (input) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and/or scaling factors applied
!          to the left side of A and B, as returned by DGGBAL.
!
!  RSCALE  (input) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and/or scaling factors applied
!          to the right side of A and B, as returned by DGGBAL.
!
!  M       (input) INTEGER
!          The number of columns of the matrix V.  M >= 0.
!
!  V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)
!          On entry, the matrix of right or left eigenvectors to be
!          transformed, as returned by DTGEVC.
!          On exit, V is overwritten by the transformed eigenvectors.
!
!  LDV     (input) INTEGER
!          The leading dimension of the matrix V. LDV >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  See R.C. Ward, Balancing the generalized eigenvalue problem,
!                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, K
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
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
      else if ( ILO < 1 ) then
         INFO = -4
      else if ( IHI < ILO .OR. IHI.GT.MAX( 1, N ) ) then
         INFO = -5
      else if ( M < 0 ) then
         INFO = -6
      else if ( LDV < MAX( 1, N ) ) then
         INFO = -10
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGGBAK', -INFO )
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
!        Backward transformation on right eigenvectors
!
         if ( RIGHTV ) then
            DO 10 I = ILO, IHI
               CALL DSCAL( M, RSCALE( I ), V( I, 1 ), LDV )
   10       CONTINUE
         end if
!
!        Backward transformation on left eigenvectors
!
         if ( LEFTV ) then
            DO 20 I = ILO, IHI
               CALL DSCAL( M, LSCALE( I ), V( I, 1 ), LDV )
   20       CONTINUE
         end if
      end if
!
!     Backward permutation
!
   30 CONTINUE
      if ( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) then
!
!        Backward permutation on right eigenvectors
!
         if ( RIGHTV ) then
            if ( ILO == 1 ) &
               GO TO 50
!
            DO 40 I = ILO - 1, 1, -1
               K = RSCALE( I )
               if ( K == I ) &
                  GO TO 40
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
!
   50       CONTINUE
            if ( IHI == N ) &
               GO TO 70
            DO 60 I = IHI + 1, N
               K = RSCALE( I )
               if ( K == I ) &
                  GO TO 60
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   60       CONTINUE
         end if
!
!        Backward permutation on left eigenvectors
!
   70    CONTINUE
         if ( LEFTV ) then
            if ( ILO == 1 ) &
               GO TO 90
            DO 80 I = ILO - 1, 1, -1
               K = LSCALE( I )
               if ( K == I ) &
                  GO TO 80
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   80       CONTINUE
!
   90       CONTINUE
            if ( IHI == N ) &
               GO TO 110
            DO 100 I = IHI + 1, N
               K = LSCALE( I )
               if ( K == I ) &
                  GO TO 100
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
  100       CONTINUE
         end if
      end if
!
  110 CONTINUE
!
      RETURN
!
!     End of DGGBAK
!
      END
