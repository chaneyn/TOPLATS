! ====================================================================
!
!                   subroutine calcvegcap
!
! ====================================================================
!
! Calculate vegetation capacity for the over story using linear
! interpolation for the canopy resistance.
!
! ====================================================================

      subroutine calcvegcap(smcond,zero,vegcap,epetd,resist,ravd)

      implicit none
      include "help/calcvegcap.h"

      if (smcond.gt.zero) then

         vegcap = epetd*(resist+ravd)/(resist/smcond + ravd)

      else

         vegcap = zero

      endif

      return

      end
