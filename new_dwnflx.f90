! ====================================================================
!
!			subroutine dwnflx
!
! ====================================================================
!
! Subroutine dwnflx calculates the downward fluxes out of the root 
! and transmission zones.
!
! ====================================================================

      subroutine new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,&
       inc_frozen,thetar,thetas,&
       rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

      implicit none
      include "help/new_dwnflx.h"

! ====================================================================
! Calculate downward flux out of root and transmission zones using 
! ksat(theta) relations of Brooks and Corey.
! ====================================================================

      if (zrz.eq.zero) then

! ====================================================================
! If root zone is saturated by the water table then set drainage
! to zero.
! ====================================================================

         grz = zero

! ====================================================================
! If root zone not saturated then calculate the drainage.
! ====================================================================

      else

! ====================================================================
! If root zone not saturated then calculate the drainage.
! ====================================================================


! --------------------------------------------------------------------
! Calcualate relative saturation in the root zone.
! --------------------------------------------------------------------

         if (inc_frozen.eq.0) then

! ....................................................................
! Option 1 : Treat frozen soil water as liquid water.
! ....................................................................

            relsrz = (rzsm-thetar)/(thetas-thetar)

         else

! ....................................................................
! Option 2 : Treat frozen soil particles as solid soil.
! ....................................................................

            relsrz = (rzsm_u-thetar)/(thetas-thetar)

         endif

         if (relsrz.le.zero) relsrz=zero
         if (relsrz.ge.one) relsrz=one

! --------------------------------------------------------------------
! Calculate downward flux out of root zone.
! --------------------------------------------------------------------

         grz = xksrz*(relsrz**((two+three*bcbeta)/bcbeta))

      endif

! ====================================================================
! Now repeat calculation for transmission zone.
! ====================================================================

      if (ztz.eq.zero) then

! ====================================================================
! If transmission zone is saturated by the water table then set drainage
! to zero.
! ====================================================================

         gtz = zero

      else

! ====================================================================
! If transmission zone not saturated then calculate the drainage.
! ====================================================================


! --------------------------------------------------------------------
! Calcualate relative saturation in the transmission zone.
! --------------------------------------------------------------------

         if (inc_frozen.eq.0) then

! ....................................................................
! Option 1 : Treat frozen soil water as liquid water.
! ....................................................................

            relstz=(tzsm-thetar)/(thetas-thetar)

         else

! ....................................................................
! Option 2 : Treat frozen soil particles as solid soil.
! ....................................................................

            relstz=(tzsm_u-thetar)/(thetas-thetar)

         endif

         if (relstz.le.zero) relstz=zero
         if (relstz.ge.one) relstz=one

! --------------------------------------------------------------------
! Calculate downward flux out of transmission zone.
! --------------------------------------------------------------------

         gtz = xkstz*(relstz**((two+three*bcbeta)/bcbeta))

      endif

      return

      end
