      SUBROUTINE DLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, &
                         SCALE, CNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, KD, LDAB, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), CNORM( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLATBS solves one of the triangular systems
!
!     A *x = s*b  or  A'*x = s*b
!
!  with scaling to prevent overflow, where A is an upper or lower
!  triangular band matrix.  Here A' denotes the transpose of A, x and b
!  are n-element vectors, and s is a scaling factor, usually less than
!  or equal to 1, chosen so that the components of x will be less than
!  the overflow threshold.  If the unscaled problem will not cause
!  overflow, the Level 2 BLAS routine DTBSV is called.  If the matrix A
!  is singular (A(j,j) = 0 for some j), then s is set to 0 and a
!  non-trivial solution to A*x = 0 is returned.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  TRANS   (input) CHARACTER*1
!          Specifies the operation applied to A.
!          = 'N':  Solve A * x = s*b  (No transpose)
!          = 'T':  Solve A'* x = s*b  (Transpose)
!          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  NORMIN  (input) CHARACTER*1
!          Specifies whether CNORM has been set or not.
!          = 'Y':  CNORM contains the column norms on entry
!          = 'N':  CNORM is not set on entry.  On exit, the norms will
!                  be computed and stored in CNORM.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of subdiagonals or superdiagonals in the
!          triangular matrix A.  KD >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The upper or lower triangular band matrix A, stored in the
!          first KD+1 rows of the array. The j-th column of A is stored
!          in the j-th column of the array AB as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  X       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the right hand side b of the triangular system.
!          On exit, X is overwritten by the solution vector x.
!
!  SCALE   (output) DOUBLE PRECISION
!          The scaling factor s for the triangular system
!             A * x = s*b  or  A'* x = s*b.
!          If SCALE = 0, the matrix A is singular or badly scaled, and
!          the vector x is an exact or approximate solution to A*x = 0.
!
!  CNORM   (input or output) DOUBLE PRECISION array, dimension (N)
!
!          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
!          contains the norm of the off-diagonal part of the j-th column
!          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
!          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
!          must be greater than or equal to the 1-norm.
!
!          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
!          returns the 1-norm of the offdiagonal part of the j-th column
!          of A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -k, the k-th argument had an illegal value
!
!  Further Details
!  ======= =======
!
!  A rough bound on x is computed; if that is less than overflow, DTBSV
!  is called, otherwise, specific code is used which checks for possible
!  overflow or divide-by-zero at every operation.
!
!  A columnwise scheme is used for solving A*x = b.  The basic algorithm
!  if A is lower triangular is
!
!       x[1:n] := b[1:n]
!       for j = 1, ..., n
!            x(j) := x(j) / A(j,j)
!            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
!       end
!
!  Define bounds on the components of x after j iterations of the loop:
!     M(j) = bound on x[1:j]
!     G(j) = bound on x[j+1:n]
!  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
!
!  Then for iteration j+1 we have
!     M(j+1) <= G(j) / | A(j+1,j+1) |
!     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
!            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
!
!  where CNORM(j+1) is greater than or equal to the infinity-norm of
!  column j+1 of A, not counting the diagonal.  Hence
!
!     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
!                  1<=i<=j
!  and
!
!     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
!                                   1<=i< j
!
!  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTBSV if the
!  reciprocal of the largest M(j), j=1,..,n, is larger than
!  max(underflow, 1/overflow).
!
!  The bound on x(j) is also used to determine when a step in the
!  columnwise method can be performed without fear of overflow.  If
!  the computed bound is greater than a large constant, x is scaled to
!  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
!
!  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic
!  algorithm for A upper triangular is
!
!       for j = 1, ..., n
!            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
!       end
!
!  We simultaneously compute two bounds
!       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
!       M(j) = bound on x(i), 1<=i<=j
!
!  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
!  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
!  Then the bound on x(j) is
!
!       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
!
!            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
!                      1<=i<=j
!
!  and we can safely call DTBSV if 1/M(n) and 1/G(n) are both greater
!  than max(underflow, 1/overflow).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST, JLEN, MAIND
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS, &
                         TMAX, TSCAL, USCAL, XBND, XJ, XMAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DASUM, DDOT, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DTBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Test the input parameters.
!
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
               LSAME( TRANS, 'C' ) ) then
         INFO = -2
      else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) then
         INFO = -3
      else if ( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. &
               LSAME( NORMIN, 'N' ) ) then
         INFO = -4
      else if ( N < 0 ) then
         INFO = -5
      else if ( KD < 0 ) then
         INFO = -6
      else if ( LDAB < KD+1 ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DLATBS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Determine machine dependent parameters to control overflow.
!
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
!
      if ( LSAME( NORMIN, 'N' ) ) then
!
!        Compute the 1-norm of each column, not including the diagonal.
!
         if ( UPPER ) then
!
!           A is upper triangular.
!
            DO 10 J = 1, N
               JLEN = MIN( KD, J-1 )
               CNORM( J ) = DASUM( JLEN, AB( KD+1-JLEN, J ), 1 )
   10       CONTINUE
         ELSE
!
!           A is lower triangular.
!
            DO 20 J = 1, N
               JLEN = MIN( KD, N-J )
               if ( JLEN.GT.0 ) then
                  CNORM( J ) = DASUM( JLEN, AB( 2, J ), 1 )
               ELSE
                  CNORM( J ) = ZERO
               end if
   20       CONTINUE
         end if
      end if
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM.
!
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      if ( TMAX.LE.BIGNUM ) then
         TSCAL = ONE
      ELSE
         TSCAL = ONE / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      end if
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine DTBSV can be used.
!
      J = IDAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      if ( NOTRAN ) then
!
!        Compute the growth in A * x = b.
!
         if ( UPPER ) then
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = KD + 1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = 1
         end if
!
         if ( TSCAL.NE.ONE ) then
            GROW = ZERO
            GO TO 50
         end if
!
         if ( NOUNIT ) then
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, G(0) = max{x(i), i=1,...,n}.
!
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               if ( GROW.LE.SMLNUM ) &
                  GO TO 50
!
!              M(j) = G(j-1) / abs(A(j,j))
!
               TJJ = ABS( AB( MAIND, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               if ( TJJ+CNORM( J ).GE.SMLNUM ) then
!
!                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
!
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
!
!                 G(j) could overflow, set GROW to 0.
!
                  GROW = ZERO
               end if
   30       CONTINUE
            GROW = XBND
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               if ( GROW.LE.SMLNUM ) &
                  GO TO 50
!
!              G(j) = G(j-1)*( 1 + CNORM(j) )
!
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         end if
   50    CONTINUE
!
      ELSE
!
!        Compute the growth in A' * x = b.
!
         if ( UPPER ) then
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = KD + 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = 1
         end if
!
         if ( TSCAL.NE.ONE ) then
            GROW = ZERO
            GO TO 80
         end if
!
         if ( NOUNIT ) then
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, M(0) = max{x(i), i=1,...,n}.
!
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               if ( GROW.LE.SMLNUM ) &
                  GO TO 80
!
!              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
!
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
!
!              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
!
               TJJ = ABS( AB( MAIND, J ) )
               if ( XJ.GT.TJJ ) &
                  XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               if ( GROW.LE.SMLNUM ) &
                  GO TO 80
!
!              G(j) = ( 1 + CNORM(j) )*G(j-1)
!
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         end if
   80    CONTINUE
      end if
!
      if ( ( GROW*TSCAL ).GT.SMLNUM ) then
!
!        Use the Level 2 BLAS solve if the reciprocal of the bound on
!        elements of X is not too small.
!
         CALL DTBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, X, 1 )
      ELSE
!
!        Use a Level 1 BLAS solve, scaling intermediate results.
!
         if ( XMAX.GT.BIGNUM ) then
!
!           Scale X so that its components are less than or equal to
!           BIGNUM in absolute value.
!
            SCALE = BIGNUM / XMAX
            CALL DSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         end if
!
         if ( NOTRAN ) then
!
!           Solve A * x = b
!
            DO 110 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
!
               XJ = ABS( X( J ) )
               if ( NOUNIT ) then
                  TJJS = AB( MAIND, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  if ( TSCAL == ONE ) &
                     GO TO 100
               end if
               TJJ = ABS( TJJS )
               if ( TJJ.GT.SMLNUM ) then
!
!                    abs(A(j,j)) > SMLNUM:
!
                  if ( TJJ < ONE ) then
                     if ( XJ.GT.TJJ*BIGNUM ) then
!
!                          Scale x by 1/b(j).
!
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     end if
                  end if
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               else if ( TJJ.GT.ZERO ) then
!
!                    0 < abs(A(j,j)) <= SMLNUM:
!
                  if ( XJ.GT.TJJ*BIGNUM ) then
!
!                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
!                       to avoid overflow when dividing by A(j,j).
!
                     REC = ( TJJ*BIGNUM ) / XJ
                     if ( CNORM( J ).GT.ONE ) then
!
!                          Scale by 1/CNORM(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
                        REC = REC / CNORM( J )
                     end if
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  end if
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
!
!                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to A*x = 0.
!
                  DO 90 I = 1, N
                     X( I ) = ZERO
   90             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               end if
  100          CONTINUE
!
!              Scale x if necessary to avoid overflow when adding a
!              multiple of column j of A.
!
               if ( XJ.GT.ONE ) then
                  REC = ONE / XJ
                  if ( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) then
!
!                    Scale x by 1/(2*abs(x(j))).
!
                     REC = REC*HALF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  end if
               else if ( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) then
!
!                 Scale x by 1/2.
!
                  CALL DSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               end if
!
               if ( UPPER ) then
                  if ( J.GT.1 ) then
!
!                    Compute the update
!                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
!                                             x(j)* A(max(1,j-kd):j-1,j)
!
                     JLEN = MIN( KD, J-1 )
                     CALL DAXPY( JLEN, -X( J )*TSCAL, &
                                 AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 )
                     I = IDAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  end if
               else if ( J < N ) then
!
!                 Compute the update
!                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
!                                          x(j) * A(j+1:min(j+kd,n),j)
!
                  JLEN = MIN( KD, N-J )
                  if ( JLEN.GT.0 ) &
                     CALL DAXPY( JLEN, -X( J )*TSCAL, AB( 2, J ), 1, &
                                 X( J+1 ), 1 )
                  I = J + IDAMAX( N-J, X( J+1 ), 1 )
                  XMAX = ABS( X( I ) )
               end if
  110       CONTINUE
!
         ELSE
!
!           Solve A' * x = b
!
            DO 160 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               if ( CNORM( J ).GT.( BIGNUM-XJ )*REC ) then
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  REC = REC*HALF
                  if ( NOUNIT ) then
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  end if
                  TJJ = ABS( TJJS )
                  if ( TJJ.GT.ONE ) then
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  end if
                  if ( REC < ONE ) then
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  end if
               end if
!
               SUMJ = ZERO
               if ( USCAL == ONE ) then
!
!                 If the scaling needed for A in the dot product is 1,
!                 call DDOT to perform the dot product.
!
                  if ( UPPER ) then
                     JLEN = MIN( KD, J-1 )
                     SUMJ = DDOT( JLEN, AB( KD+1-JLEN, J ), 1, &
                            X( J-JLEN ), 1 )
                  ELSE
                     JLEN = MIN( KD, N-J )
                     if ( JLEN.GT.0 ) &
                        SUMJ = DDOT( JLEN, AB( 2, J ), 1, X( J+1 ), 1 )
                  end if
               ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
                  if ( UPPER ) then
                     JLEN = MIN( KD, J-1 )
                     DO 120 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( KD+I-JLEN, J )*USCAL )* &
                               X( J-JLEN-1+I )
  120                CONTINUE
                  ELSE
                     JLEN = MIN( KD, N-J )
                     DO 130 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( I+1, J )*USCAL )*X( J+I )
  130                CONTINUE
                  end if
               end if
!
               if ( USCAL == TSCAL ) then
!
!                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  if ( NOUNIT ) then
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     if ( TSCAL == ONE ) &
                        GO TO 150
                  end if
                  TJJ = ABS( TJJS )
                  if ( TJJ.GT.SMLNUM ) then
!
!                       abs(A(j,j)) > SMLNUM:
!
                     if ( TJJ < ONE ) then
                        if ( XJ.GT.TJJ*BIGNUM ) then
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = ONE / XJ
                           CALL DSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        end if
                     end if
                     X( J ) = X( J ) / TJJS
                  else if ( TJJ.GT.ZERO ) then
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     if ( XJ.GT.TJJ*BIGNUM ) then
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     end if
                     X( J ) = X( J ) / TJJS
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0, and compute a solution to A'*x = 0.
!
                     DO 140 I = 1, N
                        X( I ) = ZERO
  140                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  end if
  150             CONTINUE
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - sumj if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X( J ) = X( J ) / TJJS - SUMJ
               end if
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
         end if
         SCALE = SCALE / TSCAL
      end if
!
!     Scale the column norms by 1/TSCAL for return.
!
      if ( TSCAL.NE.ONE ) then
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      end if
!
      RETURN
!
!     End of DLATBS
!
      END
