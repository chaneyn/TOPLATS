! ====================================================================
!
!                   function clcdg
!
! ====================================================================
!
! Calculate the derivative of gravity drainage out of the surface
! and transmission zone with respect to soil moisture.
!
! ====================================================================

      function clcdg(sm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)

      implicit none
      include "help/clcdg.h"

! --------------------------------------------------------------------
! Calculate the derivative of the downward flux
! --------------------------------------------------------------------

      dgrz = xksrz * (2.d0 + 3.d0 * bcbeta)/(bcbeta* (thetas- thetar))

! ....................................................................
! Calculate the relative soil saturation.
! ....................................................................
    
      satrel = (sm- thetar)/(thetas - thetar)

      dgrz = dgrz * satrel**(2.d0+(2.d0/bcbeta))
   
      clcdg=dgrz

      return

      end
