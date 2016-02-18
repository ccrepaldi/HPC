      SUBROUTINE DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYR2  performs the symmetric rank 2 operation
!
!     A := alpha*x*y' + alpha*y*x' + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND. &
               .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      else if ( N < 0 )THEN
         INFO = 2
      else if ( INCX == 0 )THEN
         INFO = 5
      else if ( INCY == 0 )THEN
         INFO = 7
      else if ( LDA < MAX( 1, N ) )THEN
         INFO = 9
      end if
      if ( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR2 ', INFO )
         RETURN
      end if
!
!     Quick return if possible.
!
      if ( ( N == 0 ).OR.( ALPHA.EQ.ZERO ) ) &
         RETURN
!
!     Set up the start points in X and Y if the increments are not both
!     unity.
!
      if ( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         if ( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         end if
         if ( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         end if
         JX = KX
         JY = KY
      end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      if ( LSAME( UPLO, 'U' ) )THEN
!
!        Form  A  when A is stored in the upper triangle.
!
         if ( ( INCX == 1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               if ( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
               end if
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               if ( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1 &
                                           + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
               end if
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         end if
      ELSE
!
!        Form  A  when A is stored in the lower triangle.
!
         if ( ( INCX == 1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               if ( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               end if
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               if ( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = JX
                  IY    = JY
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1 &
                                           + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   70             CONTINUE
               end if
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         end if
      end if
!
      RETURN
!
!     End of DSYR2 .
!
      END
