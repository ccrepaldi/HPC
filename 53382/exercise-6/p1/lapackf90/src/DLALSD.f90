      SUBROUTINE DLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, &
                         RANK, WORK, IWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
      DOUBLE PRECISION   RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLALSD uses the singular value decomposition of A to solve the least
!  squares problem of finding X to minimize the Euclidean norm of each
!  column of A*X-B, where A is N-by-N upper bidiagonal, and X and B
!  are N-by-NRHS. The solution X overwrites B.
!
!  The singular values of A smaller than RCOND times the largest
!  singular value are treated as zero in solving the least squares
!  problem; in this case a minimum norm solution is returned.
!  The actual singular values are returned in D in ascending order.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Arguments
!  =========
!
!  UPLO   (input) CHARACTER*1
!         = 'U': D and E define an upper bidiagonal matrix.
!         = 'L': D and E define a  lower bidiagonal matrix.
!
!  SMLSIZ (input) INTEGER
!         The maximum size of the subproblems at the bottom of the
!         computation tree.
!
!  N      (input) INTEGER
!         The dimension of the  bidiagonal matrix.  N >= 0.
!
!  NRHS   (input) INTEGER
!         The number of columns of B. NRHS must be at least 1.
!
!  D      (input/output) DOUBLE PRECISION array, dimension (N)
!         On entry D contains the main diagonal of the bidiagonal
!         matrix. On exit, if INFO = 0, D contains its singular values.
!
!  E      (input) DOUBLE PRECISION array, dimension (N-1)
!         Contains the super-diagonal entries of the bidiagonal matrix.
!         On exit, E has been destroyed.
!
!  B      (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!         On input, B contains the right hand sides of the least
!         squares problem. On output, B contains the solution X.
!
!  LDB    (input) INTEGER
!         The leading dimension of B in the calling subprogram.
!         LDB must be at least max(1,N).
!
!  RCOND  (input) DOUBLE PRECISION
!         The singular values of A less than or equal to RCOND times
!         the largest singular value are treated as zero in solving
!         the least squares problem. If RCOND is negative,
!         machine precision is used instead.
!         For example, if diag(S)*X=B were the least squares problem,
!         where diag(S) is a diagonal matrix of singular values, the
!         solution would be X(i) = B(i) / S(i) if S(i) is greater than
!         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to
!         RCOND*max(S).
!
!  RANK   (output) INTEGER
!         The number of singular values of A greater than RCOND times
!         the largest singular value.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension at least
!         (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2),
!         where NLVL = max(0, INT(log_2 (N/(SMLSIZ+1))) + 1).
!
!  IWORK  (workspace) INTEGER array, dimension at least
!         (3*N*NLVL + 11*N)
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  The algorithm failed to compute an singular value while
!               working on the submatrix lying in rows and columns
!               INFO/(N+1) through MOD(INFO,N+1).
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!       California at Berkeley, USA
!     Osni Marques, LBNL/NERSC, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM, &
                         GIVPTR, I, ICMPQ1, ICMPQ2, IWK, J, K, NLVL, &
                         NM1, NSIZE, NSUB, NWORK, PERM, POLES, S, SIZEI, &
                         SMLSZP, SQRE, ST, ST1, U, VT, Z
      DOUBLE PRECISION   CS, EPS, ORGNRM, R, SN, TOL
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           IDAMAX, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLALSA, DLARTG, DLASCL, &
                         DLASDA, DLASDQ, DLASET, DLASRT, DROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      if ( N < 0 ) then
         INFO = -3
      else if ( NRHS < 1 ) then
         INFO = -4
      else if ( ( LDB < 1 ) .OR. ( LDB.LT.N ) ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DLALSD', -INFO )
         RETURN
      end if
!
      EPS = DLAMCH( 'Epsilon' )
!
!     Set up the tolerance.
!
      if ( ( RCOND.LE.ZERO ) .OR. ( RCOND.GE.ONE ) ) then
         RCOND = EPS
      end if
!
      RANK = 0
!
!     Quick return if possible.
!
      if ( N == 0 ) then
         RETURN
      else if ( N == 1 ) then
         if ( D( 1 ) == ZERO ) then
            CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, B, LDB )
         ELSE
            RANK = 1
            CALL DLASCL( 'G', 0, 0, D( 1 ), ONE, 1, NRHS, B, LDB, INFO )
            D( 1 ) = ABS( D( 1 ) )
         end if
         RETURN
      end if
!
!     Rotate the matrix if it is lower bidiagonal.
!
      if ( UPLO == 'L' ) then
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            if ( NRHS == 1 ) then
               CALL DROT( 1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN )
            ELSE
               WORK( I*2-1 ) = CS
               WORK( I*2 ) = SN
            end if
   10    CONTINUE
         if ( NRHS.GT.1 ) then
            DO 30 I = 1, NRHS
               DO 20 J = 1, N - 1
                  CS = WORK( J*2-1 )
                  SN = WORK( J*2 )
                  CALL DROT( 1, B( J, I ), 1, B( J+1, I ), 1, CS, SN )
   20          CONTINUE
   30       CONTINUE
         end if
      end if
!
!     Scale.
!
      NM1 = N - 1
      ORGNRM = DLANST( 'M', N, D, E )
      if ( ORGNRM == ZERO ) then
         CALL DLASET( 'A', N, NRHS, ZERO, ZERO, B, LDB )
         RETURN
      end if
!
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, INFO )
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
      if ( N.LE.SMLSIZ ) then
         NWORK = 1 + N*N
         CALL DLASET( 'A', N, N, ZERO, ONE, WORK, N )
         CALL DLASDQ( 'U', 0, N, N, 0, NRHS, D, E, WORK, N, WORK, N, B, &
                      LDB, WORK( NWORK ), INFO )
         if ( INFO.NE.0 ) then
            RETURN
         end if
         TOL = RCOND*ABS( D( IDAMAX( N, D, 1 ) ) )
         DO 40 I = 1, N
            if ( D( I ).LE.TOL ) then
               CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            ELSE
               CALL DLASCL( 'G', 0, 0, D( I ), ONE, 1, NRHS, B( I, 1 ), &
                            LDB, INFO )
               RANK = RANK + 1
            end if
   40    CONTINUE
         CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, &
                     WORK( NWORK ), N )
         CALL DLACPY( 'A', N, NRHS, WORK( NWORK ), N, B, LDB )
!
!        Unscale.
!
         CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
         CALL DLASRT( 'D', N, D, INFO )
         CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )
!
         RETURN
      end if
!
!     Book-keeping and setting up some constants.
!
      NLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
!
      SMLSZP = SMLSIZ + 1
!
      U = 1
      VT = 1 + SMLSIZ*N
      DIFL = VT + SMLSZP*N
      DIFR = DIFL + NLVL*N
      Z = DIFR + NLVL*N*2
      C = Z + NLVL*N
      S = C + N
      POLES = S + N
      GIVNUM = POLES + 2*NLVL*N
      BX = GIVNUM + 2*NLVL*N
      NWORK = BX + N*NRHS
!
      SIZEI = 1 + N
      K = SIZEI + N
      GIVPTR = K + N
      PERM = GIVPTR + N
      GIVCOL = PERM + NLVL*N
      IWK = GIVCOL + NLVL*N*2
!
      ST = 1
      SQRE = 0
      ICMPQ1 = 1
      ICMPQ2 = 0
      NSUB = 0
!
      DO 50 I = 1, N
         if ( ABS( D( I ) ) < EPS ) then
            D( I ) = SIGN( EPS, D( I ) )
         end if
   50 CONTINUE
!
      DO 60 I = 1, NM1
         if ( ( ABS( E( I ) ) < EPS ) .OR. ( I == NM1 ) ) then
            NSUB = NSUB + 1
            IWORK( NSUB ) = ST
!
!           Subproblem found. First determine its size and then
!           apply divide and conquer on it.
!
            if ( I < NM1 ) then
!
!              A subproblem with E(I) small for I < NM1.
!
               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            else if ( ABS( E( I ) ).GE.EPS ) then
!
!              A subproblem with E(NM1) not too small but I = NM1.
!
               NSIZE = N - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            ELSE
!
!              A subproblem with E(NM1) small. This implies an
!              1-by-1 subproblem at D(N), which is not solved
!              explicitly.
!
               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
               NSUB = NSUB + 1
               IWORK( NSUB ) = N
               IWORK( SIZEI+NSUB-1 ) = 1
               CALL DCOPY( NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N )
            end if
            ST1 = ST - 1
            if ( NSIZE == 1 ) then
!
!              This is a 1-by-1 subproblem and is not solved
!              explicitly.
!
               CALL DCOPY( NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N )
            else if ( NSIZE.LE.SMLSIZ ) then
!
!              This is a small subproblem and is solved by DLASDQ.
!
               CALL DLASET( 'A', NSIZE, NSIZE, ZERO, ONE, &
                            WORK( VT+ST1 ), N )
               CALL DLASDQ( 'U', 0, NSIZE, NSIZE, 0, NRHS, D( ST ), &
                            E( ST ), WORK( VT+ST1 ), N, WORK( NWORK ), &
                            N, B( ST, 1 ), LDB, WORK( NWORK ), INFO )
               if ( INFO.NE.0 ) then
                  RETURN
               end if
               CALL DLACPY( 'A', NSIZE, NRHS, B( ST, 1 ), LDB, &
                            WORK( BX+ST1 ), N )
            ELSE
!
!              A large problem. Solve it using divide and conquer.
!
               CALL DLASDA( ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ), &
                            E( ST ), WORK( U+ST1 ), N, WORK( VT+ST1 ), &
                            IWORK( K+ST1 ), WORK( DIFL+ST1 ), &
                            WORK( DIFR+ST1 ), WORK( Z+ST1 ), &
                            WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), &
                            IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), &
                            WORK( GIVNUM+ST1 ), WORK( C+ST1 ), &
                            WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), &
                            INFO )
               if ( INFO.NE.0 ) then
                  RETURN
               end if
               BXST = BX + ST1
               CALL DLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ), &
                            LDB, WORK( BXST ), N, WORK( U+ST1 ), N, &
                            WORK( VT+ST1 ), IWORK( K+ST1 ), &
                            WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), &
                            WORK( Z+ST1 ), WORK( POLES+ST1 ), &
                            IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, &
                            IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), &
                            WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), &
                            IWORK( IWK ), INFO )
               if ( INFO.NE.0 ) then
                  RETURN
               end if
            end if
            ST = I + 1
         end if
   60 CONTINUE
!
!     Apply the singular values and treat the tiny ones as zero.
!
      TOL = RCOND*ABS( D( IDAMAX( N, D, 1 ) ) )
!
      DO 70 I = 1, N
!
!        Some of the elements in D can be negative because 1-by-1
!        subproblems were not solved explicitly.
!
         if ( ABS( D( I ) ).LE.TOL ) then
            CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, WORK( BX+I-1 ), N )
         ELSE
            RANK = RANK + 1
            CALL DLASCL( 'G', 0, 0, D( I ), ONE, 1, NRHS, &
                         WORK( BX+I-1 ), N, INFO )
         end if
         D( I ) = ABS( D( I ) )
   70 CONTINUE
!
!     Now apply back the right singular vectors.
!
      ICMPQ2 = 1
      DO 80 I = 1, NSUB
         ST = IWORK( I )
         ST1 = ST - 1
         NSIZE = IWORK( SIZEI+I-1 )
         BXST = BX + ST1
         if ( NSIZE == 1 ) then
            CALL DCOPY( NRHS, WORK( BXST ), N, B( ST, 1 ), LDB )
         else if ( NSIZE.LE.SMLSIZ ) then
            CALL DGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE, &
                        WORK( VT+ST1 ), N, WORK( BXST ), N, ZERO, &
                        B( ST, 1 ), LDB )
         ELSE
            CALL DLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, WORK( BXST ), N, &
                         B( ST, 1 ), LDB, WORK( U+ST1 ), N, &
                         WORK( VT+ST1 ), IWORK( K+ST1 ), &
                         WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), &
                         WORK( Z+ST1 ), WORK( POLES+ST1 ), &
                         IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, &
                         IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), &
                         WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), &
                         IWORK( IWK ), INFO )
            if ( INFO.NE.0 ) then
               RETURN
            end if
         end if
   80 CONTINUE
!
!     Unscale and sort the singular values.
!
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
      CALL DLASRT( 'D', N, D, INFO )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )

      RETURN
      END
