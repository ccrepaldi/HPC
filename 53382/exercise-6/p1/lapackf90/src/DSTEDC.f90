      SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, &
                         LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEDC computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the divide and conquer method.
!  The eigenvectors of a full or band real symmetric matrix can also be
!  found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this
!  matrix to tridiagonal form.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.  See DLAED3 for details.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'I':  Compute eigenvectors of tridiagonal matrix also.
!          = 'V':  Compute eigenvectors of original dense symmetric
!                  matrix also.  On entry, Z contains the orthogonal
!                  matrix used to reduce the original matrix to
!                  tridiagonal form.
!
!  N       (input) INTEGER
!          The dimension of the symmetric tridiagonal matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the subdiagonal elements of the tridiagonal matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
!          On entry, if COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If  COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1.
!          If eigenvectors are desired, then LDZ >= max(1,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array,
!                                         dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
!          If COMPZ = 'V' and N > 1 then LWORK must be at least
!                         ( 1 + 3*N + 2*N*lg N + 3*N**2 ),
!                         where lg( N ) = smallest integer k such
!                         that 2**k >= N.
!          If COMPZ = 'I' and N > 1 then LWORK must be at least
!                         ( 1 + 4*N + N**2 ).
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
!          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
!          If COMPZ = 'V' and N > 1 then LIWORK must be at least
!                         ( 6 + 6*N + 5*N*lg N ).
!          If COMPZ = 'I' and N > 1 then LIWORK must be at least
!                         ( 3 + 5*N ).
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  The algorithm failed to compute an eigenvalue while
!                working on the submatrix lying in rows and columns
!                INFO/(N+1) through mod(INFO,N+1).
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Jeff Rutter, Computer Science Division, University of California
!     at Berkeley, USA
!  Modified by Francoise Tisseur, University of Tennessee.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            DTRTRW, END, I, ICOMPZ, II, J, K, LGN, LIWMIN, &
                         LWMIN, M, SMLSIZ, START, STOREZ
      DOUBLE PRECISION   EPS, ORGNRM, P, TINY
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLAED0, DLASCL, DLASET, DLASRT, &
                         DSTEQR, DSTERF, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX, MOD, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LQUERY = ( LWORK == -1 .OR. LIWORK.EQ.-1 )
!
      if ( LSAME( COMPZ, 'N' ) ) then
         ICOMPZ = 0
      else if ( LSAME( COMPZ, 'V' ) ) then
         ICOMPZ = 1
      else if ( LSAME( COMPZ, 'I' ) ) then
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      end if
      if ( N.LE.1 .OR. ICOMPZ.LE.0 ) then
         LIWMIN = 1
         LWMIN = 1
      ELSE
         LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) )
         if ( 2**LGN < N ) &
            LGN = LGN + 1
         if ( 2**LGN < N ) &
            LGN = LGN + 1
         if ( ICOMPZ == 1 ) then
            LWMIN = 1 + 3*N + 2*N*LGN + 3*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         else if ( ICOMPZ == 2 ) then
            LWMIN = 1 + 4*N + N**2
            LIWMIN = 3 + 5*N
         end if
      end if
      if ( ICOMPZ < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( ( LDZ < 1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, &
               N ) ) ) then
         INFO = -6
      else if ( LWORK < LWMIN .AND. .NOT.LQUERY ) then
         INFO = -8
      else if ( LIWORK < LIWMIN .AND. .NOT.LQUERY ) then
         INFO = -10
      end if
!
      if ( INFO == 0 ) then
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DSTEDC', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
      if ( N == 1 ) then
         if ( ICOMPZ.NE.0 ) &
            Z( 1, 1 ) = ONE
         RETURN
      end if
!
      SMLSIZ = ILAENV( 9, 'DSTEDC', ' ', 0, 0, 0, 0 )
!
!     If the following conditional clause is removed, then the routine
!     will use the Divide and Conquer routine to compute only the
!     eigenvalues, which requires (3N + 3N**2) real workspace and
!     (2 + 5N + 2N lg(N)) integer workspace.
!     Since on many architectures DSTERF is much faster than any other
!     algorithm for finding eigenvalues only, it is used here
!     as the default.
!
!     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
!
      if ( ICOMPZ == 0 ) then
         CALL DSTERF( N, D, E, INFO )
         RETURN
      end if
!
!     If N is smaller than the minimum divide size (SMLSIZ+1), then
!     solve the problem with another solver.
!
      if ( N.LE.SMLSIZ ) then
         if ( ICOMPZ == 0 ) then
            CALL DSTERF( N, D, E, INFO )
            RETURN
         else if ( ICOMPZ == 2 ) then
            CALL DSTEQR( 'I', N, D, E, Z, LDZ, WORK, INFO )
            RETURN
         ELSE
            CALL DSTEQR( 'V', N, D, E, Z, LDZ, WORK, INFO )
            RETURN
         end if
      end if
!
!     If COMPZ = 'V', the Z matrix must be stored elsewhere for later
!     use.
!
      if ( ICOMPZ == 1 ) then
         STOREZ = 1 + N*N
      ELSE
         STOREZ = 1
      end if
!
      if ( ICOMPZ == 2 ) then
         CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
      end if
!
!     Scale.
!
      ORGNRM = DLANST( 'M', N, D, E )
      if ( ORGNRM == ZERO ) &
         RETURN
!
      EPS = DLAMCH( 'Epsilon' )
!
      START = 1
!
!     while ( START <= N )
!
   10 CONTINUE
      if ( START.LE.N ) then
!
!     Let END be the position of the next subdiagonal entry such that
!     E( END ) <= TINY or END = N if no such subdiagonal exists.  The
!     matrix identified by the elements between START and END
!     constitutes an independent sub-problem.
!
         END = START
   20    CONTINUE
         if ( END < N ) then
            TINY = EPS*SQRT( ABS( D( END ) ) )*SQRT( ABS( D( END+1 ) ) )
            if ( ABS( E( END ) ).GT.TINY ) then
               END = END + 1
               GO TO 20
            end if
         end if
!
!        (Sub) Problem determined.  Compute its size and solve it.
!
         M = END - START + 1
         if ( M == 1 ) then
            START = END + 1
            GO TO 10
         end if
         if ( M.GT.SMLSIZ ) then
            INFO = SMLSIZ
!
!           Scale.
!
            ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, &
                         INFO )
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), &
                         M-1, INFO )
!
            if ( ICOMPZ == 1 ) then
               DTRTRW = 1
            ELSE
               DTRTRW = START
            end if
            CALL DLAED0( ICOMPZ, N, M, D( START ), E( START ), &
                         Z( DTRTRW, START ), LDZ, WORK( 1 ), N, &
                         WORK( STOREZ ), IWORK, INFO )
            if ( INFO.NE.0 ) then
               INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + &
                      MOD( INFO, ( M+1 ) ) + START - 1
               RETURN
            end if
!
!           Scale back.
!
            CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, &
                         INFO )
!
         ELSE
            if ( ICOMPZ == 1 ) then
!
!     Since QR won't update a Z matrix which is larger than the
!     length of D, we must solve the sub-problem in a workspace and
!     then multiply back into Z.
!
               CALL DSTEQR( 'I', M, D( START ), E( START ), WORK, M, &
                            WORK( M*M+1 ), INFO )
               CALL DLACPY( 'A', N, M, Z( 1, START ), LDZ, &
                            WORK( STOREZ ), N )
               CALL DGEMM( 'N', 'N', N, M, M, ONE, WORK( STOREZ ), LDZ, &
                           WORK, M, ZERO, Z( 1, START ), LDZ )
            else if ( ICOMPZ == 2 ) then
               CALL DSTEQR( 'I', M, D( START ), E( START ), &
                            Z( START, START ), LDZ, WORK, INFO )
            ELSE
               CALL DSTERF( M, D( START ), E( START ), INFO )
            end if
            if ( INFO.NE.0 ) then
               INFO = START*( N+1 ) + END
               RETURN
            end if
         end if
!
         START = END + 1
         GO TO 10
      end if
!
!     endwhile
!
!     If the problem split any number of times, then the eigenvalues
!     will not be properly ordered.  Here we permute the eigenvalues
!     (and the associated eigenvectors) into ascending order.
!
      if ( M.NE.N ) then
         if ( ICOMPZ == 0 ) then
!
!        Use Quick Sort
!
            CALL DLASRT( 'I', N, D, INFO )
!
         ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
            DO 40 II = 2, N
               I = II - 1
               K = I
               P = D( I )
               DO 30 J = II, N
                  if ( D( J ) < P ) then
                     K = J
                     P = D( J )
                  end if
   30          CONTINUE
               if ( K.NE.I ) then
                  D( K ) = D( I )
                  D( I ) = P
                  CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
               end if
   40       CONTINUE
         end if
      end if
!
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of DSTEDC
!
      END
