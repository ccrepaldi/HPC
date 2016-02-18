      SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, &
                         WORK, IWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     December 1, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ, UPLO
      INTEGER            INFO, LDU, LDVT, N
!     ..
!     .. Array Arguments ..
      INTEGER            IQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), Q( * ), U( LDU, * ), &
                         VT( LDVT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DBDSDC computes the singular value decomposition (SVD) of a real
!  N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
!  using a divide and conquer method, where S is a diagonal matrix
!  with non-negative diagonal elements (the singular values of B), and
!  U and VT are orthogonal matrices of left and right singular vectors,
!  respectively. DBDSDC can be used to compute all singular values,
!  and optionally, singular vectors or singular vectors in compact form.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.  See DLASD3 for details.
!
!  The code currently call DLASDQ if singular values only are desired.
!  However, it can be slightly modified to compute singular values
!  using the divide and conquer method.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  B is upper bidiagonal.
!          = 'L':  B is lower bidiagonal.
!
!  COMPQ   (input) CHARACTER*1
!          Specifies whether singular vectors are to be computed
!          as follows:
!          = 'N':  Compute singular values only;
!          = 'P':  Compute singular values and compute singular
!                  vectors in compact form;
!          = 'I':  Compute singular values and singular vectors.
!
!  N       (input) INTEGER
!          The order of the matrix B.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the bidiagonal matrix B.
!          On exit, if INFO=0, the singular values of B.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the elements of E contain the offdiagonal
!          elements of the bidiagonal matrix whose SVD is desired.
!          On exit, E has been destroyed.
!
!  U       (output) DOUBLE PRECISION array, dimension (LDU,N)
!          If  COMPQ = 'I', then:
!             On exit, if INFO = 0, U contains the left singular vectors
!             of the bidiagonal matrix.
!          For other values of COMPQ, U is not referenced.
!
!  LDU     (input) INTEGER
!          The leading dimension of the array U.  LDU >= 1.
!          If singular vectors are desired, then LDU >= max( 1, N ).
!
!  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
!          If  COMPQ = 'I', then:
!             On exit, if INFO = 0, VT' contains the right singular
!             vectors of the bidiagonal matrix.
!          For other values of COMPQ, VT is not referenced.
!
!  LDVT    (input) INTEGER
!          The leading dimension of the array VT.  LDVT >= 1.
!          If singular vectors are desired, then LDVT >= max( 1, N ).
!
!  Q       (output) DOUBLE PRECISION array, dimension (LDQ)
!          If  COMPQ = 'P', then:
!             On exit, if INFO = 0, Q and IQ contain the left
!             and right singular vectors in a compact form,
!             requiring O(N log N) space instead of 2*N**2.
!             In particular, Q contains all the DOUBLE PRECISION data in
!             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
!             words of memory, where SMLSIZ is returned by ILAENV and
!             is equal to the maximum size of the subproblems at the
!             bottom of the computation tree (usually about 25).
!          For other values of COMPQ, Q is not referenced.
!
!  IQ      (output) INTEGER array, dimension (LDIQ)
!          If  COMPQ = 'P', then:
!             On exit, if INFO = 0, Q and IQ contain the left
!             and right singular vectors in a compact form,
!             requiring O(N log N) space instead of 2*N**2.
!             In particular, IQ contains all INTEGER data in
!             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
!             words of memory, where SMLSIZ is returned by ILAENV and
!             is equal to the maximum size of the subproblems at the
!             bottom of the computation tree (usually about 25).
!          For other values of COMPQ, IQ is not referenced.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
!          If COMPQ = 'N' then LWORK >= (4 * N).
!          If COMPQ = 'P' then LWORK >= (6 * N).
!          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).
!
!  IWORK   (workspace) INTEGER array, dimension (8*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  The algorithm failed to compute an singular value.
!                The update process of divide and conquer failed.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, IC, &
                         ICOMPQ, IERR, II, IS, IU, IUPLO, IVT, J, K, KK, &
                         MLVL, NM1, NSIZE, PERM, POLES, QSTART, SMLSIZ, &
                         SMLSZP, SQRE, START, WSTART, Z
      DOUBLE PRECISION   CS, EPS, ORGNRM, P, R, SN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLARTG, DLASCL, DLASD0, DLASDA, DLASDQ, &
                         DLASET, DLASR, DSWAP, XERBLA
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
      IUPLO = 0
      if ( LSAME( UPLO, 'U' ) ) &
         IUPLO = 1
      if ( LSAME( UPLO, 'L' ) ) &
         IUPLO = 2
      if ( LSAME( COMPQ, 'N' ) ) then
         ICOMPQ = 0
      else if ( LSAME( COMPQ, 'P' ) ) then
         ICOMPQ = 1
      else if ( LSAME( COMPQ, 'I' ) ) then
         ICOMPQ = 2
      ELSE
         ICOMPQ = -1
      end if
      if ( IUPLO == 0 ) then
         INFO = -1
      else if ( ICOMPQ < 0 ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( ( LDU < 1 ) .OR. ( ( ICOMPQ == 2 ) .AND. ( LDU.LT. &
               N ) ) ) then
         INFO = -7
      else if ( ( LDVT < 1 ) .OR. ( ( ICOMPQ == 2 ) .AND. ( LDVT.LT. &
               N ) ) ) then
         INFO = -9
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DBDSDC', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
      SMLSIZ = ILAENV( 9, 'DBDSDC', ' ', 0, 0, 0, 0 )
      if ( N == 1 ) then
         if ( ICOMPQ == 1 ) then
            Q( 1 ) = SIGN( ONE, D( 1 ) )
            Q( 1+SMLSIZ*N ) = ONE
         else if ( ICOMPQ == 2 ) then
            U( 1, 1 ) = SIGN( ONE, D( 1 ) )
            VT( 1, 1 ) = ONE
         end if
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      end if
      NM1 = N - 1
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
      WSTART = 1
      QSTART = 3
      if ( ICOMPQ == 1 ) then
         CALL DCOPY( N, D, 1, Q( 1 ), 1 )
         CALL DCOPY( N-1, E, 1, Q( N+1 ), 1 )
      end if
      if ( IUPLO == 2 ) then
         QSTART = 5
         WSTART = 2*N - 1
         DO I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            if ( ICOMPQ == 1 ) then
               Q( I+2*N ) = CS
               Q( I+3*N ) = SN
            else if ( ICOMPQ == 2 ) then
               WORK( I ) = CS
               WORK( NM1+I ) = -SN
            end if
         end do
      end if
!
!     If ICOMPQ = 0, use DLASDQ to compute the singular values.
!
      if ( ICOMPQ == 0 ) then
         CALL DLASDQ( 'U', 0, N, 0, 0, 0, D, E, VT, LDVT, U, LDU, U, &
                      LDU, WORK( WSTART ), INFO )
         GO TO 40
      end if
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
      if ( N.LE.SMLSIZ ) then
         if ( ICOMPQ == 2 ) then
            CALL DLASET( 'A', N, N, ZERO, ONE, U, LDU )
            CALL DLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
            CALL DLASDQ( 'U', 0, N, N, N, 0, D, E, VT, LDVT, U, LDU, U, &
                         LDU, WORK( WSTART ), INFO )
         else if ( ICOMPQ == 1 ) then
            IU = 1
            IVT = IU + N
            CALL DLASET( 'A', N, N, ZERO, ONE, Q( IU+( QSTART-1 )*N ), &
                         N )
            CALL DLASET( 'A', N, N, ZERO, ONE, Q( IVT+( QSTART-1 )*N ), &
                         N )
            CALL DLASDQ( 'U', 0, N, N, N, 0, D, E, &
                         Q( IVT+( QSTART-1 )*N ), N, &
                         Q( IU+( QSTART-1 )*N ), N, &
                         Q( IU+( QSTART-1 )*N ), N, WORK( WSTART ), &
                         INFO )
         end if
         GO TO 40
      end if
!
      if ( ICOMPQ == 2 ) then
         CALL DLASET( 'A', N, N, ZERO, ONE, U, LDU )
         CALL DLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
      end if
!
!     Scale.
!
      ORGNRM = DLANST( 'M', N, D, E )
      if ( ORGNRM == ZERO ) &
         RETURN
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, IERR )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, IERR )
!
      EPS = DLAMCH( 'Epsilon' )
!
      MLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
      SMLSZP = SMLSIZ + 1
!
      if ( ICOMPQ == 1 ) then
         IU = 1
         IVT = 1 + SMLSIZ
         DIFL = IVT + SMLSZP
         DIFR = DIFL + MLVL
         Z = DIFR + MLVL*2
         IC = Z + MLVL
         IS = IC + 1
         POLES = IS + 1
         GIVNUM = POLES + 2*MLVL
!
         K = 1
         GIVPTR = 2
         PERM = 3
         GIVCOL = PERM + MLVL
      end if
!
      DO 20 I = 1, N
         if ( ABS( D( I ) ) < EPS ) then
            D( I ) = SIGN( EPS, D( I ) )
         end if
   20 CONTINUE
!
      START = 1
      SQRE = 0
!
      DO 30 I = 1, NM1
         if ( ( ABS( E( I ) ) < EPS ) .OR. ( I == NM1 ) ) then
!
!        Subproblem found. First determine its size and then
!        apply divide and conquer on it.
!
            if ( I < NM1 ) then
!
!        A subproblem with E(I) small for I < NM1.
!
               NSIZE = I - START + 1
            else if ( ABS( E( I ) ).GE.EPS ) then
!
!        A subproblem with E(NM1) not too small but I = NM1.
!
               NSIZE = N - START + 1
            ELSE
!
!        A subproblem with E(NM1) small. This implies an
!        1-by-1 subproblem at D(N). Solve this 1-by-1 problem
!        first.
!
               NSIZE = I - START + 1
               if ( ICOMPQ == 2 ) then
                  U( N, N ) = SIGN( ONE, D( N ) )
                  VT( N, N ) = ONE
               else if ( ICOMPQ == 1 ) then
                  Q( N+( QSTART-1 )*N ) = SIGN( ONE, D( N ) )
                  Q( N+( SMLSIZ+QSTART-1 )*N ) = ONE
               end if
               D( N ) = ABS( D( N ) )
            end if
            if ( ICOMPQ == 2 ) then
               CALL DLASD0( NSIZE, SQRE, D( START ), E( START ), &
                            U( START, START ), LDU, VT( START, START ), &
                            LDVT, SMLSIZ, IWORK, WORK( WSTART ), INFO )
            ELSE
               CALL DLASDA( ICOMPQ, SMLSIZ, NSIZE, SQRE, D( START ), &
                            E( START ), Q( START+( IU+QSTART-2 )*N ), N, &
                            Q( START+( IVT+QSTART-2 )*N ), &
                            IQ( START+K*N ), Q( START+( DIFL+QSTART-2 )* &
                            N ), Q( START+( DIFR+QSTART-2 )*N ), &
                            Q( START+( Z+QSTART-2 )*N ), &
                            Q( START+( POLES+QSTART-2 )*N ), &
                            IQ( START+GIVPTR*N ), IQ( START+GIVCOL*N ), &
                            N, IQ( START+PERM*N ), &
                            Q( START+( GIVNUM+QSTART-2 )*N ), &
                            Q( START+( IC+QSTART-2 )*N ), &
                            Q( START+( IS+QSTART-2 )*N ), &
                            WORK( WSTART ), IWORK, INFO )
               if ( INFO.NE.0 ) then
                  RETURN
               end if
            end if
            START = I + 1
         end if
   30 CONTINUE
!
!     Unscale
!
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, IERR )
   40 CONTINUE
!
!     Use Selection Sort to minimize swaps of singular vectors
!
      DO 60 II = 2, N
         I = II - 1
         KK = I
         P = D( I )
         DO 50 J = II, N
            if ( D( J ).GT.P ) then
               KK = J
               P = D( J )
            end if
   50    CONTINUE
         if ( KK.NE.I ) then
            D( KK ) = D( I )
            D( I ) = P
            if ( ICOMPQ == 1 ) then
               IQ( I ) = KK
            else if ( ICOMPQ == 2 ) then
               CALL DSWAP( N, U( 1, I ), 1, U( 1, KK ), 1 )
               CALL DSWAP( N, VT( I, 1 ), LDVT, VT( KK, 1 ), LDVT )
            end if
         else if ( ICOMPQ == 1 ) then
            IQ( I ) = I
         end if
   60 CONTINUE
!
!     If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO
!
      if ( ICOMPQ == 1 ) then
         if ( IUPLO == 1 ) then
            IQ( N ) = 1
         ELSE
            IQ( N ) = 0
         end if
      end if
!
!     If B is lower bidiagonal, update U by those Givens rotations
!     which rotated B to be upper bidiagonal
!
      if ( ( IUPLO == 2 ) .AND. ( ICOMPQ.EQ.2 ) ) &
         CALL DLASR( 'L', 'V', 'B', N, N, WORK( 1 ), WORK( N ), U, LDU )
!
      RETURN
!
!     End of DBDSDC
!
      END
