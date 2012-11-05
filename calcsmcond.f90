! ====================================================================
!
!                  subroutine calcsmcond
!
! ====================================================================
!
! Calculate the soil moisture conductance.
!
! ====================================================================

      subroutine calcsmcond(rzsm,tc,smcond,one,tw,zero)

      implicit none
      include "help/calcsmcond.h"

      if (rzsm.ge.tc) then

         smcond = one

      else if (rzsm.ge.tw) then

         smcond = (rzsm-tw)/(tc-tw)

      else

         smcond = zero

      endif

      if ( (smcond.ge.0.d0).and.(smcond.le.1.d0) ) then

         smcond=smcond

      else

         write (*,*) 'CALCSMCOND : smcond out of bounds ',smcond
         if (smcond.le.0.d0) smcond=zero
         if (smcond.ge.1.d0) smcond=one

      endif

      return

      end
