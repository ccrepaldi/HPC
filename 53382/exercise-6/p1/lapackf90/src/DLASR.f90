      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASR   performs the transformation
!
!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
!
!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
!
!  where A is an m by n real matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
!  and z = n when SIDE = 'R' or 'r' ):
!
!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
!
!     P = P( z - 1 )*...*P( 2 )*P( 1 ),
!
!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
!
!     P = P( 1 )*P( 2 )*...*P( z - 1 ),
!
!  where  P( k ) is a plane rotation matrix for the following planes:
!
!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
!        the plane ( k, k + 1 )
!
!     when  PIVOT = 'T' or 't'  ( Top pivot ),
!        the plane ( 1, k + 1 )
!
!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
!        the plane ( k, z )
!
!  c( k ) and s( k )  must contain the  cosine and sine that define the
!  matrix  P( k ).  The two by two plane rotation part of the matrix
!  P( k ), R( k ), is assumed to be of the form
!
!     R( k ) = (  c( k )  s( k ) ).
!              ( -s( k )  c( k ) )
!
!  This version vectorises across rows of the array A when SIDE = 'L'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'
!
!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C, S    (input) DOUBLE PRECISION arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R( k ) = (  c( k )  s( k ) ).
!                   ( -s( k )  c( k ) )
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      if ( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) then
         INFO = 1
      else if ( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, &
               'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) then
         INFO = 2
      else if ( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) &
                THEN
         INFO = 3
      else if ( M < 0 ) then
         INFO = 4
      else if ( N < 0 ) then
         INFO = 5
      else if ( LDA < MAX( 1, M ) ) then
         INFO = 9
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( ( M == 0 ) .OR. ( N.EQ.0 ) ) &
         RETURN
      if ( LSAME( SIDE, 'L' ) ) then
!
!        Form  P * A
!
         if ( LSAME( PIVOT, 'V' ) ) then
            if ( LSAME( DIRECT, 'F' ) ) then
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  end if
   20          CONTINUE
            else if ( LSAME( DIRECT, 'B' ) ) then
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  end if
   40          CONTINUE
            end if
         else if ( LSAME( PIVOT, 'T' ) ) then
            if ( LSAME( DIRECT, 'F' ) ) then
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  end if
   60          CONTINUE
            else if ( LSAME( DIRECT, 'B' ) ) then
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  end if
   80          CONTINUE
            end if
         else if ( LSAME( PIVOT, 'B' ) ) then
            if ( LSAME( DIRECT, 'F' ) ) then
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  end if
  100          CONTINUE
            else if ( LSAME( DIRECT, 'B' ) ) then
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  end if
  120          CONTINUE
            end if
         end if
      else if ( LSAME( SIDE, 'R' ) ) then
!
!        Form A * P'
!
         if ( LSAME( PIVOT, 'V' ) ) then
            if ( LSAME( DIRECT, 'F' ) ) then
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  end if
  140          CONTINUE
            else if ( LSAME( DIRECT, 'B' ) ) then
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  end if
  160          CONTINUE
            end if
         else if ( LSAME( PIVOT, 'T' ) ) then
            if ( LSAME( DIRECT, 'F' ) ) then
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  end if
  180          CONTINUE
            else if ( LSAME( DIRECT, 'B' ) ) then
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  end if
  200          CONTINUE
            end if
         else if ( LSAME( PIVOT, 'B' ) ) then
            if ( LSAME( DIRECT, 'F' ) ) then
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  end if
  220          CONTINUE
            else if ( LSAME( DIRECT, 'B' ) ) then
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) then
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  end if
  240          CONTINUE
            end if
         end if
      end if
!
      RETURN
!
!     End of DLASR
!
      END
