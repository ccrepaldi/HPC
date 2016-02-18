      SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
!     ..
!
!  Purpose
!  =======
!
!  DLASV2 computes the singular value decomposition of a 2-by-2
!  triangular matrix
!     [  F   G  ]
!     [  0   H  ].
!  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
!  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
!  right singular vectors for abs(SSMAX), giving the decomposition
!
!     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
!     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  G       (input) DOUBLE PRECISION
!          The (1,2) element of the 2-by-2 matrix.
!
!  H       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  SSMIN   (output) DOUBLE PRECISION
!          abs(SSMIN) is the smaller singular value.
!
!  SSMAX   (output) DOUBLE PRECISION
!          abs(SSMAX) is the larger singular value.
!
!  SNL     (output) DOUBLE PRECISION
!  CSL     (output) DOUBLE PRECISION
!          The vector (CSL, SNL) is a unit left singular vector for the
!          singular value abs(SSMAX).
!
!  SNR     (output) DOUBLE PRECISION
!  CSR     (output) DOUBLE PRECISION
!          The vector (CSR, SNR) is a unit right singular vector for the
!          singular value abs(SSMAX).
!
!  Further Details
!  ===============
!
!  Any input parameter may be aliased with any output parameter.
!
!  Barring over/underflow and assuming a guard digit in subtraction, all
!  output quantities are correct to within a few units in the last
!  place (ulps).
!
!  In IEEE arithmetic, the code works correctly if one matrix element is
!  infinite.
!
!  Overflow will not occur unless the largest singular value itself
!  overflows or is within a few ulps of overflow. (On machines with
!  partial overflow, like the Cray, overflow may occur if the largest
!  singular value is within a factor of 2 of overflow.)
!
!  Underflow is harmless if underflow is gradual. Otherwise, results
!  may correspond to a matrix modified by perturbations of size near
!  the underflow threshold.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   FOUR
      PARAMETER          ( FOUR = 4.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            GASMAL, SWAP
      INTEGER            PMAX
      DOUBLE PRECISION   A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M, &
                         MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Executable Statements ..
!
      FT = F
      FA = ABS( FT )
      HT = H
      HA = ABS( H )
!
!     PMAX points to the maximum absolute element of matrix
!       PMAX = 1 if F largest in absolute values
!       PMAX = 2 if G largest in absolute values
!       PMAX = 3 if H largest in absolute values
!
      PMAX = 1
      SWAP = ( HA.GT.FA )
      if ( SWAP ) then
         PMAX = 3
         TEMP = FT
         FT = HT
         HT = TEMP
         TEMP = FA
         FA = HA
         HA = TEMP
!
!        Now FA .ge. HA
!
      end if
      GT = G
      GA = ABS( GT )
      if ( GA == ZERO ) then
!
!        Diagonal matrix
!
         SSMIN = HA
         SSMAX = FA
         CLT = ONE
         CRT = ONE
         SLT = ZERO
         SRT = ZERO
      ELSE
         GASMAL = .TRUE.
         if ( GA.GT.FA ) then
            PMAX = 2
            if ( ( FA / GA ) < DLAMCH( 'EPS' ) ) then
!
!              Case of very large GA
!
               GASMAL = .FALSE.
               SSMAX = GA
               if ( HA.GT.ONE ) then
                  SSMIN = FA / ( GA / HA )
               ELSE
                  SSMIN = ( FA / GA )*HA
               end if
               CLT = ONE
               SLT = HT / GT
               SRT = ONE
               CRT = FT / GT
            end if
         end if
         if ( GASMAL ) then
!
!           Normal case
!
            D = FA - HA
            if ( D == FA ) then
!
!              Copes with infinite F or H
!
               L = ONE
            ELSE
               L = D / FA
            end if
!
!           Note that 0 .le. L .le. 1
!
            M = GT / FT
!
!           Note that abs(M) .le. 1/macheps
!
            T = TWO - L
!
!           Note that T .ge. 1
!
            MM = M*M
            TT = T*T
            S = SQRT( TT+MM )
!
!           Note that 1 .le. S .le. 1 + 1/macheps
!
            if ( L == ZERO ) then
               R = ABS( M )
            ELSE
               R = SQRT( L*L+MM )
            end if
!
!           Note that 0 .le. R .le. 1 + 1/macheps
!
            A = HALF*( S+R )
!
!           Note that 1 .le. A .le. 1 + abs(M)
!
            SSMIN = HA / A
            SSMAX = FA*A
            if ( MM == ZERO ) then
!
!              Note that M is very tiny
!
               if ( L == ZERO ) then
                  T = SIGN( TWO, FT )*SIGN( ONE, GT )
               ELSE
                  T = GT / SIGN( D, FT ) + M / T
               end if
            ELSE
               T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A )
            end if
            L = SQRT( T*T+FOUR )
            CRT = TWO / L
            SRT = T / L
            CLT = ( CRT+SRT*M ) / A
            SLT = ( HT / FT )*SRT / A
         end if
      end if
      if ( SWAP ) then
         CSL = SRT
         SNL = CRT
         CSR = SLT
         SNR = CLT
      ELSE
         CSL = CLT
         SNL = SLT
         CSR = CRT
         SNR = SRT
      end if
!
!     Correct signs of SSMAX and SSMIN
!
      if ( PMAX == 1 ) &
         TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )
      if ( PMAX == 2 ) &
         TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )
      if ( PMAX == 3 ) &
         TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H )
      SSMAX = SIGN( SSMAX, TSIGN )
      SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) )
      RETURN
!
!     End of DLASV2
!
      END
