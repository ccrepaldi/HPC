      SUBROUTINE DLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK, &
                         WORK, INFO )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDU, LDVT, N, SMLSIZ, SQRE
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  Using a divide and conquer approach, DLASD0 computes the singular
!  value decomposition (SVD) of a real upper bidiagonal N-by-M
!  matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
!  The algorithm computes orthogonal matrices U and VT such that
!  B = U * S * VT. The singular values S are overwritten on D.
!
!  A related subroutine, DLASDA, computes only the singular values,
!  and optionally, the singular vectors in compact form.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         On entry, the row dimension of the upper bidiagonal matrix.
!         This is also the dimension of the main diagonal array D.
!
!  SQRE   (input) INTEGER
!         Specifies the column dimension of the bidiagonal matrix.
!         = 0: The bidiagonal matrix has column dimension M = N;
!         = 1: The bidiagonal matrix has column dimension M = N+1;
!
!  D      (input/output) DOUBLE PRECISION array, dimension (N)
!         On entry D contains the main diagonal of the bidiagonal
!         matrix.
!         On exit D, if INFO = 0, contains its singular values.
!
!  E      (input) DOUBLE PRECISION array, dimension (M-1)
!         Contains the subdiagonal entries of the bidiagonal matrix.
!         On exit, E has been destroyed.
!
!  U      (output) DOUBLE PRECISION array, dimension at least (LDQ, N)
!         On exit, U contains the left singular vectors.
!
!  LDU    (input) INTEGER
!         On entry, leading dimension of U.
!
!  VT     (output) DOUBLE PRECISION array, dimension at least (LDVT, M)
!         On exit, VT' contains the right singular vectors.
!
!  LDVT   (input) INTEGER
!         On entry, leading dimension of VT.
!
!  SMLSIZ (input) INTEGER
!         On entry, maximum size of the subproblems at the
!         bottom of the computation tree.
!
!  IWORK  INTEGER work array.
!         Dimension must be at least (8 * N)
!
!  WORK   DOUBLE PRECISION work array.
!         Dimension must be at least (3 * M**2 + 2 * M)
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an singular value did not converge
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
!     .. Local Scalars ..
      INTEGER            I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK, &
                         J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR, &
                         NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQREI
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASD1, DLASDQ, DLASDT, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      if ( N < 0 ) then
         INFO = -1
      else if ( ( SQRE < 0 ) .OR. ( SQRE.GT.1 ) ) then
         INFO = -2
      end if
!
      M = N + SQRE
!
      if ( LDU < N ) then
         INFO = -6
      else if ( LDVT < M ) then
         INFO = -8
      else if ( SMLSIZ < 3 ) then
         INFO = -9
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DLASD0', -INFO )
         RETURN
      end if
!
!     If the input matrix is too small, call DLASDQ to find the SVD.
!
      if ( N.LE.SMLSIZ ) then
         CALL DLASDQ( 'U', SQRE, N, M, N, 0, D, E, VT, LDVT, U, LDU, U, &
                      LDU, WORK, INFO )
         RETURN
      end if
!
!     Set up the computation tree.
!
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ = NDIMR + N
      IWK = IDXQ + N
      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), &
                   IWORK( NDIMR ), SMLSIZ )
!
!     For the nodes on bottom level of the tree, solve
!     their subproblems by DLASDQ.
!
      NDB1 = ( ND+1 ) / 2
      NCC = 0
      DO 30 I = NDB1, ND
!
!     IC : center row of each node
!     NL : number of rows of left  subproblem
!     NR : number of rows of right subproblem
!     NLF: starting row of the left   subproblem
!     NRF: starting row of the right  subproblem
!
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NLP1 = NL + 1
         NR = IWORK( NDIMR+I1 )
         NRP1 = NR + 1
         NLF = IC - NL
         NRF = IC + 1
         SQREI = 1
         CALL DLASDQ( 'U', SQREI, NL, NLP1, NL, NCC, D( NLF ), E( NLF ), &
                      VT( NLF, NLF ), LDVT, U( NLF, NLF ), LDU, &
                      U( NLF, NLF ), LDU, WORK, INFO )
         if ( INFO.NE.0 ) then
            RETURN
         end if
         ITEMP = IDXQ + NLF - 2
         DO 10 J = 1, NL
            IWORK( ITEMP+J ) = J
   10    CONTINUE
         if ( I == ND ) then
            SQREI = SQRE
         ELSE
            SQREI = 1
         end if
         NRP1 = NR + SQREI
         CALL DLASDQ( 'U', SQREI, NR, NRP1, NR, NCC, D( NRF ), E( NRF ), &
                      VT( NRF, NRF ), LDVT, U( NRF, NRF ), LDU, &
                      U( NRF, NRF ), LDU, WORK, INFO )
         if ( INFO.NE.0 ) then
            RETURN
         end if
         ITEMP = IDXQ + IC
         DO 20 J = 1, NR
            IWORK( ITEMP+J-1 ) = J
   20    CONTINUE
   30 CONTINUE
!
!     Now conquer each subproblem bottom-up.
!
      DO 50 LVL = NLVL, 1, -1
!
!        Find the first node LF and last node LL on the
!        current level LVL.
!
         if ( LVL == 1 ) then
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         end if
         DO 40 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            if ( ( SQRE == 0 ) .AND. ( I.EQ.LL ) ) then
               SQREI = SQRE
            ELSE
               SQREI = 1
            end if
            IDXQC = IDXQ + NLF - 1
            ALPHA = D( IC )
            BETA = E( IC )
            CALL DLASD1( NL, NR, SQREI, D( NLF ), ALPHA, BETA, &
                         U( NLF, NLF ), LDU, VT( NLF, NLF ), LDVT, &
                         IWORK( IDXQC ), IWORK( IWK ), WORK, INFO )
            if ( INFO.NE.0 ) then
               RETURN
            end if
   40    CONTINUE
   50 CONTINUE
!
      RETURN
!
!     End of DLASD0
!
      END
