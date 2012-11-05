! ====================================================================
!
!                   subroutine maxplevap
!
! ====================================================================
!
! Calculate the maximum plant evaporation.
!
! ====================================================================

      subroutine maxplevap(zrz,ztz,epetd,inc_frozen,srzrel,rzsm,thetas,&
       thetar,rzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xksrz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla)

      implicit none
      include "help/maxplevap.h"

! --------------------------------------------------------------------
! If the soil is saturated than the maximum flux of the water is not
! bounded by the plant/soil system.
! --------------------------------------------------------------------

      if (zrz.le.zero) then

         vegcap = epetd

      else

! --------------------------------------------------------------------
! If the soil is not saturated calculate the soil saturation.
! --------------------------------------------------------------------

         if (inc_frozen.eq.0) then

! ....................................................................
! No difference between ice particles and frozen soil water.
! ....................................................................

            srzrel = (rzsm-thetar)/(thetas-thetar)

         else

! ....................................................................
! Ice crystals are treated as mineral soil particles.
! ....................................................................

            srzrel = (rzsm_u-thetar)/(thetas-thetar)

         endif

! --------------------------------------------------------------------
! Double check whether relative saturation is between 1 and 0.
! --------------------------------------------------------------------

         if (srzrel.le.zero) srzrel=zero
         if (srzrel.ge.one) srzrel=one

! --------------------------------------------------------------------
! Calculate the soil water potential.
! --------------------------------------------------------------------

         psisoi=psic/(srzrel**(one/bcbeta))

! --------------------------------------------------------------------
! Calculate the saturated hydrauli! conductivity in the root zone.
! --------------------------------------------------------------------

         if (ikopt.eq.1) then

! ....................................................................
! Option 1 : No decline with depth.
! ....................................................................

            xksrz = xk0

         else

! ....................................................................
! Option 2 : Exponential decline with depth.
! ....................................................................

            xksrz = xk0*dexp(-ff*((zrz+ztz)/two))

         endif

! --------------------------------------------------------------------
! Calculate the unsaturated hydrauli! conductivity in the root zone.
! --------------------------------------------------------------------

         xkrz = xksrz*(srzrel**((two+three*bcbeta)/bcbeta))

! --------------------------------------------------------------------
! Calculate the soil resistance for evaporation.
! --------------------------------------------------------------------

         ressoi = one/(rtact*xkrz*rtdens)

! --------------------------------------------------------------------
! Calculate the maximum water vapor flux out of the plant.
! --------------------------------------------------------------------

         vegcap = (psisoi-psicri)/(ressoi+respla)

! --------------------------------------------------------------------
! Check whether the maximum water vapor flux out of the plant is
! positive.
! --------------------------------------------------------------------

         if (vegcap.lt.zero) vegcap = zero

      endif

      return

      end
