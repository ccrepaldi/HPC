      SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
!     ..
!
!  Purpose
!  =======
!
!  DGTSV  solves the equation
!
!     A*X = B,
!
!  where A is an n by n tridiagonal matrix, by Gaussian elimination with
!  partial pivoting.
!
!  Note that the equation  A'*X = B  may be solved by interchanging the
!  order of the arguments DU and DL.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, DL must contain the (n-1) sub-diagonal elements of
!          A.
!
!          On exit, DL is overwritten by the (n-2) elements of the
!          second super-diagonal of the upper triangular matrix U from
!          the LU factorization of A, in DL(1), ..., DL(n-2).
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, D must contain the diagonal elements of A.
!
!          On exit, D is overwritten by the n diagonal elements of U.
!
!  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, DU must contain the (n-1) super-diagonal elements
!          of A.
!
!          On exit, DU is overwritten by the (n-1) elements of the first
!          super-diagonal of U.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N by NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N by NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
!               has not been computed.  The factorization has not been
!               completed unless i = N.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   FACT, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      if ( N < 0 ) then
         INFO = -1
      else if ( NRHS < 0 ) then
         INFO = -2
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -7
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGTSV ', -INFO )
         RETURN
      end if
!
      if ( N == 0 ) &
         RETURN
!
      if ( NRHS == 1 ) then
         DO 10 I = 1, N - 2
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) then
!
!              No row interchange required
!
               if ( D( I ).NE.ZERO ) then
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               end if
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            end if
   10    CONTINUE
         if ( N.GT.1 ) then
            I = N - 1
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) then
               if ( D( I ).NE.ZERO ) then
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               end if
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            end if
         end if
         if ( D( N ) == ZERO ) then
            INFO = N
            RETURN
         end if
      ELSE
         DO 40 I = 1, N - 2
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) then
!
!              No row interchange required
!
               if ( D( I ).NE.ZERO ) then
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 20 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   20             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               end if
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               DO 30 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   30          CONTINUE
            end if
   40    CONTINUE
         if ( N.GT.1 ) then
            I = N - 1
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) then
               if ( D( I ).NE.ZERO ) then
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 50 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   50             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               end if
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               DO 60 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   60          CONTINUE
            end if
         end if
         if ( D( N ) == ZERO ) then
            INFO = N
            RETURN
         end if
      end if
!
!     Back solve with the matrix U from the factorization.
!
      if ( NRHS.LE.2 ) then
         J = 1
   70    CONTINUE
         B( N, J ) = B( N, J ) / D( N )
         if ( N.GT.1 ) &
            B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 80 I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* &
                        B( I+2, J ) ) / D( I )
   80    CONTINUE
         if ( J < NRHS ) then
            J = J + 1
            GO TO 70
         end if
      ELSE
         DO 100 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            if ( N.GT.1 ) &
               B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / &
                             D( N-1 )
            DO 90 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* &
                           B( I+2, J ) ) / D( I )
   90       CONTINUE
  100    CONTINUE
      end if
!
      RETURN
!
!     End of DGTSV
!
      END
