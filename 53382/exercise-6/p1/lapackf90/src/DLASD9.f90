      SUBROUTINE DLASD9( ICOMPQ, LDU, K, D, Z, VF, VL, DIFL, DIFR, &
                         DSIGMA, WORK, INFO )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
!     Courant Institute, NAG Ltd., and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, K, LDU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DIFL( * ), DIFR( LDU, * ), DSIGMA( * ), &
                         VF( * ), VL( * ), WORK( * ), Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASD9 finds the square roots of the roots of the secular equation,
!  as defined by the values in DSIGMA and Z.  It makes the
!  appropriate calls to DLASD4, and stores, for each  element in D,
!  the distance to its two nearest poles (elements in DSIGMA). It also
!  updates the arrays VF and VL, the first and last components of all
!  the right singular vectors of the original bidiagonal matrix.
!
!  DLASD9 is called from DLASD7.
!
!  Arguments
!  =========
!
!  ICOMPQ  (input) INTEGER
!          Specifies whether singular vectors are to be computed in
!          factored form in the calling routine:
!
!             ICOMPQ = 0             Compute singular values only.
!
!             ICOMPQ = 1             Compute singular vector matrices in
!                                    factored form also.
!  K       (input) INTEGER
!          The number of terms in the rational function to be solved by
!          DLASD4.  K >= 1.
!
!  D       (output) DOUBLE PRECISION array, dimension(K)
!          D(I) contains the updated singular values.
!
!  DSIGMA  (input) DOUBLE PRECISION array, dimension(K)
!          The first K elements of this array contain the old roots
!          of the deflated updating problem.  These are the poles
!          of the secular equation.
!
!  Z       (input) DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain the components
!          of the deflation-adjusted updating row vector.
!
!  VF      (input/output) DOUBLE PRECISION array, dimension(K)
!          On entry, VF contains  information passed through SBEDE8.f
!          On exit, VF contains the first K components of the first
!          components of all right singular vectors of the bidiagonal
!          matrix.
!
!  VL      (input/output) DOUBLE PRECISION array, dimension(K)
!          On entry, VL contains  information passed through SBEDE8.f
!          On exit, VL contains the first K components of the last
!          components of all right singular vectors of the bidiagonal
!          matrix.
!
!  DIFL    (output) DOUBLE PRECISION array, dimension (K).
!          On exit, DIFL(I) = D(I) - DSIGMA(I).
!
!  DIFR    (output) DOUBLE PRECISION array,
!                              dimension (LDU, 2) if ICOMPQ =1 and
!                              dimension (K) if ICOMPQ = 0.
!          On exit, DIFR(I, 1) = D(I) - DSIGMA(I+1), DIFR(K, 1) is not
!          defined and will not be referenced.
!
!          If ICOMPQ = 1, DIFR(1:K, 2) is an array containing the
!          normalizing factors for the right singular vector matrix.
!
!  WORK    (workspace) DOUBLE PRECISION array,
!                                 dimension at least (3 * K)
!          Workspace.
!
!  INFO    (output) INTEGER
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
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J
      DOUBLE PRECISION   DIFLJ, DIFRJ, DJ
!     DOUBLE PRECISION   DJP1
      DOUBLE PRECISION   DSIGJ, DSIGJP, RHO, TEMP
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DDOT, DLAMC3, DNRM2
      EXTERNAL           DDOT, DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLASCL, DLASD4, DLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      if ( ( ICOMPQ < 0 ) .OR. ( ICOMPQ.GT.1 ) ) then
         INFO = -1
      else if ( K < 1 ) then
         INFO = -3
      else if ( LDU < K ) then
         INFO = -2
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DLASD9', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( K == 1 ) then
         D( 1 ) = ABS( Z( 1 ) )
         DIFL( 1 ) = D( 1 )
         if ( ICOMPQ == 1 ) then
            DIFL( 2 ) = ONE
            DIFR( 1, 2 ) = ONE
         end if
         RETURN
      end if
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 10 I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   10 CONTINUE
!
!     Book keeping.
!
      IWK1 = 1
      IWK2 = IWK1 + K
      IWK3 = IWK2 + K
      IWK2I = IWK2 - 1
      IWK3I = IWK3 - 1
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Initialize WORK(IWK3).
!
      CALL DLASET( 'A', K, 1, ONE, ONE, WORK( IWK3 ), K )
!
!     Compute the updated singular values, the arrays DIFL, DIFR,
!     and the updated Z.
!
      DO 40 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), &
                      WORK( IWK2 ), INFO )
!
!        If the root finder fails, the computation is terminated.
!
         if ( INFO.NE.0 ) then
            RETURN
         end if
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )
         DIFR( J, 1 ) = -WORK( J+1 )
         DO 20 I = 1, J - 1
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* &
                              WORK( IWK2I+I ) / ( DSIGMA( I )- &
                              DSIGMA( J ) ) / ( DSIGMA( I )+ &
                              DSIGMA( J ) )
   20    CONTINUE
         DO 30 I = J + 1, K
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* &
                              WORK( IWK2I+I ) / ( DSIGMA( I )- &
                              DSIGMA( J ) ) / ( DSIGMA( I )+ &
                              DSIGMA( J ) )
   30    CONTINUE
   40 CONTINUE
!
!     Compute updated Z.
!
      DO 50 I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
   50 CONTINUE
!
!     Update VF and VL.
!
      DO 80 J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         if ( J < K ) then
            DIFRJ = -DIFR( J, 1 )
!           DJP1 = D( J+1 )
            DSIGJP = -DSIGMA( J+1 )
         end if
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         DO 60 I = 1, J - 1
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJ )-DIFLJ ) &
                         / ( DSIGMA( I )+DJ )
   60    CONTINUE
         DO I = J + 1, K
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJP )+DIFRJ ) &
                         / ( DSIGMA( I )+DJ )
         end do
         TEMP = DNRM2( K, WORK, 1 )
         WORK( IWK2I+J ) = DDOT( K, WORK, 1, VF, 1 ) / TEMP
         WORK( IWK3I+J ) = DDOT( K, WORK, 1, VL, 1 ) / TEMP
         if ( ICOMPQ == 1 ) then
            DIFR( J, 2 ) = TEMP
         end if
   80 CONTINUE

      CALL DCOPY( K, WORK( IWK2 ), 1, VF, 1 )
      CALL DCOPY( K, WORK( IWK3 ), 1, VL, 1 )

      RETURN
      END
