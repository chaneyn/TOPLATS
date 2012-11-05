! ====================================================================
!
! 			function factln
!
! ====================================================================
!
! Returns the natural log of the i factorial.
!
! ====================================================================
!
!
!  Variable definitions:
!  
!    factln:     The return value (ln (i factorial))
!    i:          Argument for which factorial is calculated
!    jj:         Loop index for series of adding logs
!    x:          Real*8 value used to send to dlog function
! ====================================================================

      function factln(i)
      implicit none
      include "help/factln.h"

! ====================================================================
! Calculate the log of factorial by adding all logs up
! to the argument value.
! ====================================================================

      factln = 0.0

      do 100 jj=1,i

        x=jj
        factln = factln + dlog(x)

100   continue

      return

      end
