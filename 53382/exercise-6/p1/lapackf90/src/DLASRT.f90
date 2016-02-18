      SUBROUTINE DLASRT( ID, N, D, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
!     ..
!     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
      INFO = 0
      DIR = -1
      if ( LSAME( ID, 'D' ) ) then
         DIR = 0
      else if ( LSAME( ID, 'I' ) ) then
         DIR = 1
      end if
      if ( DIR == -1 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N.LE.1 ) &
         RETURN
!
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      if ( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) then
!
!        Do Insertion sort on D( START:ENDD )
!
         if ( DIR == 0 ) then
!
!           Sort into decreasing order
!
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  if ( D( J ).GT.D( J-1 ) ) then
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  end if
   20          CONTINUE
   30       CONTINUE
!
         ELSE
!
!           Sort into increasing order
!
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  if ( D( J ) < D( J-1 ) ) then
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  end if
   40          CONTINUE
   50       CONTINUE
!
         end if
!
      else if ( ENDD-START.GT.SELECT ) then
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         if ( D1 < D2 ) then
            if ( D3 < D1 ) then
               DMNMX = D1
            else if ( D3 < D2 ) then
               DMNMX = D3
            ELSE
               DMNMX = D2
            end if
         ELSE
            if ( D3 < D2 ) then
               DMNMX = D2
            else if ( D3 < D1 ) then
               DMNMX = D3
            ELSE
               DMNMX = D1
            end if
         end if
!
         if ( DIR == 0 ) then
!
!           Sort into decreasing order
!
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            if ( D( J ) < DMNMX ) &
               GO TO 70
   80       CONTINUE
            I = I + 1
            if ( D( I ).GT.DMNMX ) &
               GO TO 80
            if ( I < J ) then
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            end if
            if ( J-START.GT.ENDD-J-1 ) then
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            end if
         ELSE
!
!           Sort into increasing order
!
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            if ( D( J ).GT.DMNMX ) &
               GO TO 100
  110       CONTINUE
            I = I + 1
            if ( D( I ) < DMNMX ) &
               GO TO 110
            if ( I < J ) then
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            end if
            if ( J-START.GT.ENDD-J-1 ) then
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            end if
         end if
      end if
      if ( STKPNT.GT.0 ) &
         GO TO 10
      RETURN
!
!     End of DLASRT
!
      END
