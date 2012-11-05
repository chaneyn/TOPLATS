! ====================================================================
!
!                   function calctc_j
!
! ====================================================================
!
! Calculate the thermal conductivity (thermc) of
! an unsaturated soil using Johansen's method
!
! INPUT:
!
! soil moisture (theta) (% vol),
! residual saturation (thetar) (% vol),
! porosity (thetas) (% vol),
! quartz content (quartz) (%),
! frozen soil flag (iffroz) (1=frozen,0=not).
! coarse soil flag (ifcoarse) (1=coarse,0=fine).
!
! ====================================================================

      function calctc_j(theta,thetar,thetas,quartz,iffrozen,ifcoarse)

      implicit none
      include "help/calctc_j.h"
      data soildens/2700.d0/,ko/2.d0/,kq/7.7d0/

! --------------------------------------------------------------------
! Double check the soil moisture inputs.
! --------------------------------------------------------------------

      if ( (theta.ge.0.d0).and.(theta.le.1.d0) ) then

         theta=theta

      else

         write (*,*) 'Soil moisture out of bounds in CALCTC_J'
         write (*,*) theta,thetar,thetas,quartz,iffrozen,ifcoarse
         stop

      endif

! --------------------------------------------------------------------
! Calculate relative saturation.
! --------------------------------------------------------------------

      satrel=(theta-thetar)/(thetas-thetar)

! --------------------------------------------------------------------
! Calculate bulk density .
! --------------------------------------------------------------------

      bulkdens = soildens*(1.d0 - thetas)

! --------------------------------------------------------------------
! Compute dry thermal conductivity (kdry).
! --------------------------------------------------------------------

      if (ifcoarse.eq.1) then

         kdry = 0.39d0*(thetas)**(-2.2d0)

      else

         kdry = 0.137d0*bulkdens + 64.7d0
         kdry = kdry/(soildens-0.947d0*bulkdens)

      endif
	
! --------------------------------------------------------------------
! Compute Kersten number (Ke).
! --------------------------------------------------------------------

      if (iffrozen.eq.1) then

         Ke = satrel

      else

         if (ifcoarse.eq.1) then

            if (satrel.gt.0.05) then

               Ke= 0.7d0*dlog10(satrel) + 1.d0

            else

               Ke = 0.d0

            endif

         else	

            if (satrel.gt.0.1) then

               Ke= dlog10(satrel) + 1.d0

            else

               Ke = 0.d0

            endif

         endif

      endif

! --------------------------------------------------------------------
! Compute solids thermal conductivity (ks).
! --------------------------------------------------------------------

      if (ifcoarse.eq.1.and.quartz.lt.0.20) then

         ks = kq**quartz * (ko+1.)**(1.d0-quartz)

      else

         ks = kq**quartz * ko**(1.d0-quartz)

      endif

! --------------------------------------------------------------------
! Compute saturated thermal conductivity (ksat).
! --------------------------------------------------------------------

      if (iffrozen.eq.1) then

         ksat = 2.2d0**thetas * ks**(1-thetas) * 0.269d0**theta	

      else

         ksat = 0.57d0**thetas * ks**(1-thetas)

      endif

! --------------------------------------------------------------------
! Compute thermal conductivity (thermc).
! --------------------------------------------------------------------

      thermc = (ksat-kdry)*Ke + kdry

      if ( (thermc.lt.100.d0).or.(thermc.gt.0.d0) ) then

         thermc=thermc

      else

         write (*,*) 'CALCTC_J : therm! unrealisti! ',thermc
         write (*,*) theta,thetar,thetas,quartz,iffrozen,ifcoarse
         stop

      endif

      calctc_j = thermc
	
      return

      end
