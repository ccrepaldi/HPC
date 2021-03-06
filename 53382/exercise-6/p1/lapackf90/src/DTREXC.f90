      SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, &
                         INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DTREXC reorders the real Schur factorization of a real matrix
!  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
!  moved to row ILST.
!
!  The real Schur form T is reordered by an orthogonal similarity
!  transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
!  is updated by postmultiplying it with Z.
!
!  T must be in Schur canonical form (as returned by DHSEQR), that is,
!  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!  2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  COMPQ   (input) CHARACTER*1
!          = 'V':  update the matrix Q of Schur vectors;
!          = 'N':  do not update Q.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
!          On entry, the upper quasi-triangular matrix T, in Schur
!          Schur canonical form.
!          On exit, the reordered upper quasi-triangular matrix, again
!          in Schur canonical form.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
!          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!          orthogonal transformation matrix Z which reorders T.
!          If COMPQ = 'N', Q is not referenced.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q.  LDQ >= max(1,N).
!
!  IFST    (input/output) INTEGER
!  ILST    (input/output) INTEGER
!          Specify the reordering of the diagonal blocks of T.
!          The block with row index IFST is moved to row ILST, by a
!          sequence of transpositions between adjacent blocks.
!          On exit, if IFST pointed on entry to the second row of a
!          2-by-2 block, it is changed to point to the first row; ILST
!          always points to the first row of the block in its final
!          position (which may differ from its input value by +1 or -1).
!          1 <= IFST <= N; 1 <= ILST <= N.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          = 1:  two adjacent blocks were too close to swap (the problem
!                is very ill-conditioned); T may have been partially
!                reordered, and ILST points to the first row of the
!                current position of the block being moved.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            HERE, NBF, NBL, NBNEXT
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAEXC, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input arguments.
!
      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      if ( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDT < MAX( 1, N ) ) then
         INFO = -4
      else if ( LDQ < 1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) then
         INFO = -6
      else if ( IFST < 1 .OR. IFST.GT.N ) then
         INFO = -7
      else if ( ILST < 1 .OR. ILST.GT.N ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DTREXC', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N.LE.1 ) &
         RETURN
!
!     Determine the first row of specified block
!     and find out it is 1 by 1 or 2 by 2.
!
      if ( IFST.GT.1 ) then
         if ( T( IFST, IFST-1 ).NE.ZERO ) &
            IFST = IFST - 1
      end if
      NBF = 1
      if ( IFST < N ) then
         if ( T( IFST+1, IFST ).NE.ZERO ) &
            NBF = 2
      end if
!
!     Determine the first row of the final block
!     and find out it is 1 by 1 or 2 by 2.
!
      if ( ILST.GT.1 ) then
         if ( T( ILST, ILST-1 ).NE.ZERO ) &
            ILST = ILST - 1
      end if
      NBL = 1
      if ( ILST < N ) then
         if ( T( ILST+1, ILST ).NE.ZERO ) &
            NBL = 2
      end if
!
      if ( IFST == ILST ) &
         RETURN
!
      if ( IFST < ILST ) then
!
!        Update ILST
!
         if ( NBF == 2 .AND. NBL.EQ.1 ) &
            ILST = ILST - 1
         if ( NBF == 1 .AND. NBL.EQ.2 ) &
            ILST = ILST + 1
!
         HERE = IFST
!
   10    CONTINUE
!
!        Swap block with next one below
!
         if ( NBF == 1 .OR. NBF.EQ.2 ) then
!
!           Current block either 1 by 1 or 2 by 2
!
            NBNEXT = 1
            if ( HERE+NBF+1.LE.N ) then
               if ( T( HERE+NBF+1, HERE+NBF ).NE.ZERO ) &
                  NBNEXT = 2
            end if
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT, &
                         WORK, INFO )
            if ( INFO.NE.0 ) then
               ILST = HERE
               RETURN
            end if
            HERE = HERE + NBNEXT
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            if ( NBF == 2 ) then
               if ( T( HERE+1, HERE ) == ZERO ) &
                  NBF = 3
            end if
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            NBNEXT = 1
            if ( HERE+3.LE.N ) then
               if ( T( HERE+3, HERE+2 ).NE.ZERO ) &
                  NBNEXT = 2
            end if
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT, &
                         WORK, INFO )
            if ( INFO.NE.0 ) then
               ILST = HERE
               RETURN
            end if
            if ( NBNEXT == 1 ) then
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, &
                            WORK, INFO )
               HERE = HERE + 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               if ( T( HERE+2, HERE+1 ) == ZERO ) &
                  NBNEXT = 1
               if ( NBNEXT == 2 ) then
!
!                 2 by 2 Block did not split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, &
                               NBNEXT, WORK, INFO )
                  if ( INFO.NE.0 ) then
                     ILST = HERE
                     RETURN
                  end if
                  HERE = HERE + 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                               WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1, &
                               WORK, INFO )
                  HERE = HERE + 2
               end if
            end if
         end if
         if ( HERE < ILST ) &
            GO TO 10
!
      ELSE
!
         HERE = IFST
   20    CONTINUE
!
!        Swap block with next one above
!
         if ( NBF == 1 .OR. NBF.EQ.2 ) then
!
!           Current block either 1 by 1 or 2 by 2
!
            NBNEXT = 1
            if ( HERE.GE.3 ) then
               if ( T( HERE-1, HERE-2 ).NE.ZERO ) &
                  NBNEXT = 2
            end if
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                         NBF, WORK, INFO )
            if ( INFO.NE.0 ) then
               ILST = HERE
               RETURN
            end if
            HERE = HERE - NBNEXT
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            if ( NBF == 2 ) then
               if ( T( HERE+1, HERE ) == ZERO ) &
                  NBF = 3
            end if
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            NBNEXT = 1
            if ( HERE.GE.3 ) then
               if ( T( HERE-1, HERE-2 ).NE.ZERO ) &
                  NBNEXT = 2
            end if
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                         1, WORK, INFO )
            if ( INFO.NE.0 ) then
               ILST = HERE
               RETURN
            end if
            if ( NBNEXT == 1 ) then
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1, &
                            WORK, INFO )
               HERE = HERE - 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               if ( T( HERE, HERE-1 ) == ZERO ) &
                  NBNEXT = 1
               if ( NBNEXT == 2 ) then
!
!                 2 by 2 Block did not split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1, &
                               WORK, INFO )
                  if ( INFO.NE.0 ) then
                     ILST = HERE
                     RETURN
                  end if
                  HERE = HERE - 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                               WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1, &
                               WORK, INFO )
                  HERE = HERE - 2
               end if
            end if
         end if
         if ( HERE.GT.ILST ) &
            GO TO 20
      end if
      ILST = HERE
!
      RETURN
!
!     End of DTREXC
!
      END
