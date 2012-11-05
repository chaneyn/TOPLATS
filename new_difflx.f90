! ====================================================================
!
!			subroutine difflx
!
! ====================================================================
!
! Subroutine difflx calculates the diffusive fluxes out of the surface 
! and transmission zones (positive down).
!
! ====================================================================

      subroutine new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

      implicit none
      include "help/new_difflx.h"

      data tolsat / 0.001d0 /

! ====================================================================
! Calculate diffusive flux out of root and transmission zones using 
! D(theta) relations of Brooks and Corey.
! ====================================================================

! ====================================================================
! If root zone not saturated then calculate the drainage.
! ====================================================================

      if (inc_frozen.eq.0) then

! --------------------------------------------------------------------
! If frozen soil water is treated as if it were not frozen.
! --------------------------------------------------------------------

	 if (zrz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
! ....................................................................

	    difrz=diffuse(zrz,ztz,rzsm,(0.5*rzsm+0.5*tzsm),&
                          thetas,thetar,bcbeta,psic,xksrz,xkstz)

	 else

! ....................................................................
! The soil profile is saturated.
! ....................................................................

	    difrz=zero

	 endif

      else

! --------------------------------------------------------------------
! If frozen soil water is treated as solid soil particles.
! --------------------------------------------------------------------

         if (zrz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
! ....................................................................

	    difrz=diffuse(zrz,ztz,rzsm_u,(0.5*rzsm_u+0.5*tzsm_u),&
                          thetas,thetar,bcbeta,psic,xksrz,xkstz)

         else
! ....................................................................
! The soil profile is saturated.
! ....................................................................

	   difrz=zero

	 endif

      endif

! ====================================================================
! Now repeat calculation for transmission zone.
! ====================================================================

      if (inc_frozen.eq.0) then

! --------------------------------------------------------------------
! If frozen soil water is treated as if it were not frozen.
! --------------------------------------------------------------------

         if (ztz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
!
! Here we assume diffusion occurs over layer of thickness equal to
! transmission zone and again account for any nonlinearity in soil
! moisture profile by averaging gradients.
! ....................................................................

	    diftz=diffuse(ztz,ztz,tzsm,0.5*tzsm+0.5*thetas,&
                          thetas,thetar,bcbeta,psic,xkstz,xkstz)

         else

! ....................................................................
! The soil profile is saturated.
! ....................................................................

	    diftz=zero

         endif

      else

! --------------------------------------------------------------------
! If frozen soil water is treated as solid soil particles.
! --------------------------------------------------------------------

         if (ztz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
! ....................................................................

	    diftz=diffuse(ztz,ztz,tzsm_u,0.5*tzsm_u+0.5*(thetas),&
                          thetas,thetar,bcbeta,psic,xkstz,xkstz)

         else

! ....................................................................
! The soil profile is saturated.
! ....................................................................

            diftz=zero

         endif

      endif
 
      return

      end

! ====================================================================
!
!			subroutine diffuse
!
! ====================================================================
!
! Subroutine diffuse calculates the diffusive flux from layer 1 to layer 2
! in Richards' Equation using approach of Mahrt and Pan (1984).
!
! ====================================================================

      function diffuse(dz1,dz2,theta1,theta2,thetas,thetar,bcbeta,psic,&
                       xksat1,xksat2)

      implicit none
      include "help/diffuse.noup.h"

! ====================================================================
! Calculate diffusivity  using centered approx.
! ====================================================================

      theta = 0.5d0*(theta1 + theta2)

! ====================================================================
! Use harmoni! average since flow is perp. to layers.
! ====================================================================

      xksat = 1/(0.5/xksat1 + 0.5/xksat2)

      F1 = bcbeta * xksat * psic/(thetas-thetar) 

      satrel=(theta-thetar)/(thetas-thetar)

      if (satrel.lt.(0.d0)) satrel=0.d0

      DF = F1 * (satrel) ** (bcbeta + 2.)

! ====================================================================
! Calculate moisture gradient.
!
! dz1 and dz2 are layer thicknesses (not depths) over which diffusion
! is modeled.
! ====================================================================

      dz=0.5d0*(dz1+dz2)

      if (dz.gt.1.d-9) then

         grad= (theta1-theta2)/dz

      else

         grad=0.d0

      endif

! ====================================================================
! Calculate diffusive flux.
! ====================================================================

      difflx = DF * grad

      diffuse=difflx

      return

      end
