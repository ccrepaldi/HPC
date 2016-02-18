      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
                       N4 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      if ( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000, &
              1100 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      if ( IZ == 90 .OR. IZ.EQ.122 ) then
!
!        ASCII character set
!
         if ( IC.GE.97 .AND. IC.LE.122 ) then
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if ( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         end if
!
      else if ( IZ == 233 .OR. IZ.EQ.169 ) then
!
!        EBCDIC character set
!
         if ( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) then
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if ( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) ) &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         end if
!
      else if ( IZ == 218 .OR. IZ.EQ.250 ) then
!
!        Prime machines:  ASCII+128
!
         if ( IC.GE.225 .AND. IC.LE.250 ) then
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if ( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         end if
      end if
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1 == 'S' .OR. C1.EQ.'D'
      CNAME = C1 == 'C' .OR. C1.EQ.'Z'
      if ( .NOT.( CNAME .OR. SNAME ) ) &
         RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      if ( C2 == 'GE' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         else if ( C3 == 'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3 == 'QLF' ) then
            if ( SNAME ) then
               NB = 32
            ELSE
               NB = 32
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NB = 32
            ELSE
               NB = 32
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NB = 32
            ELSE
               NB = 32
            end if
         else if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         end if
      else if ( C2 == 'PO' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         end if
      else if ( C2 == 'SY' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         else if ( SNAME .AND. C3 == 'TRD' ) then
            NB = 32
         else if ( SNAME .AND. C3 == 'GST' ) then
            NB = 64
         end if
      else if ( CNAME .AND. C2 == 'HE' ) then
         if ( C3 == 'TRF' ) then
            NB = 64
         else if ( C3 == 'TRD' ) then
            NB = 32
         else if ( C3 == 'GST' ) then
            NB = 64
         end if
      else if ( SNAME .AND. C2 == 'OR' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         end if
      else if ( CNAME .AND. C2 == 'UN' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         end if
      else if ( C2 == 'GB' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               if ( N4.LE.64 ) then
                  NB = 1
               ELSE
                  NB = 32
               end if
            ELSE
               if ( N4.LE.64 ) then
                  NB = 1
               ELSE
                  NB = 32
               end if
            end if
         end if
      else if ( C2 == 'PB' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               if ( N2.LE.64 ) then
                  NB = 1
               ELSE
                  NB = 32
               end if
            ELSE
               if ( N2.LE.64 ) then
                  NB = 1
               ELSE
                  NB = 32
               end if
            end if
         end if
      else if ( C2 == 'TR' ) then
         if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         end if
      else if ( C2 == 'LA' ) then
         if ( C3 == 'UUM' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         end if
      else if ( SNAME .AND. C2 == 'ST' ) then
         if ( C3 == 'EBZ' ) then
            NB = 1
         end if
      end if
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      if ( C2 == 'GE' ) then
         if ( C3 == 'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3 == 'QLF' ) then
            if ( SNAME ) then
               NBMIN = 2
            ELSE
               NBMIN = 2
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NBMIN = 2
            ELSE
               NBMIN = 2
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NBMIN = 2
            ELSE
               NBMIN = 2
            end if
         else if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NBMIN = 2
            ELSE
               NBMIN = 2
            end if
         end if
      else if ( C2 == 'SY' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NBMIN = 8
            ELSE
               NBMIN = 8
            end if
         else if ( SNAME .AND. C3 == 'TRD' ) then
            NBMIN = 2
         end if
      else if ( CNAME .AND. C2 == 'HE' ) then
         if ( C3 == 'TRD' ) then
            NBMIN = 2
         end if
      else if ( SNAME .AND. C2 == 'OR' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2
            end if
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2
            end if
         end if
      else if ( CNAME .AND. C2 == 'UN' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2
            end if
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2
            end if
         end if
      end if
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      if ( C2 == 'GE' ) then
         if ( C3 == 'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3 == 'QLF' ) then
            if ( SNAME ) then
               NX = 128
            ELSE
               NX = 128
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NX = 128
            ELSE
               NX = 128
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NX = 128
            ELSE
               NX = 128
            end if
         end if
      else if ( C2 == 'SY' ) then
         if ( SNAME .AND. C3 == 'TRD' ) then
            NX = 32
         end if
      else if ( CNAME .AND. C2 == 'HE' ) then
         if ( C3 == 'TRD' ) then
            NX = 32
         end if
      else if ( SNAME .AND. C2 == 'OR' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NX = 128
            end if
         end if
      else if ( CNAME .AND. C2 == 'UN' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NX = 128
            end if
         end if
      end if
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  900 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
 1000 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      if ( ILAENV == 1 ) then
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      end if
      RETURN
!
 1100 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      if ( ILAENV == 1 ) then
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      end if
      RETURN
!
!     End of ILAENV
!
      END
