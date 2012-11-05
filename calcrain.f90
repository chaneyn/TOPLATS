! ====================================================================
!
!                subroutine calcrain
!
! ====================================================================
!
! Calculate the partitioning of precipitation into rain and snow.
! Do this based on air temperature.
!
! ====================================================================

      subroutine calcrain (tcel,snow,rain,precip_o,dt)

      implicit none
      include "help/calcrain.h"

      rain=0.d0
      snow=0.d0

      if (tcel.gt.(0.d0)) then

         snow=0.d0*dt
         rain=precip_o*dt

      endif

      if (tcel.le.(0.d0)) then

         rain=0.d0*dt
         snow=precip_o*dt

      endif

      return

      end
