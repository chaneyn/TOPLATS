! ====================================================================
!
!               subroutine peteb_dv
!
! ====================================================================
!
! Solve the energy balance at potential rate for dry vegetation, over
! story.
!
! ====================================================================

      subroutine peteb_dv(thermc2,vpsat,heatcap,heatcap2,heatcapold,&
       rs_over,rain,snow,Swq,albd,emiss,thermc,f1par,f3vpd,f4temp,rescan,&
       ravd,rahd,tkd,tkmidd,rnetd,xled,epetd,hd,gd,dshd,&
       tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,&
       rld,toleb,maxnri,dt,i,PackWater,SurfWater,VaporMassFlux,TPack,&
       TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,&
       dens,za,zpd,z0h,RaSnow,appa,uzw,gact,alb_snow,row)

      implicit none
      include "SNOW.h"
      include "help/peteb_dv.h"

      if ( (Swq.le.(0.0002d0)).OR.(SNOW_RUN.eq.0) ) then

! ....................................................................
! If the total snow water equivalent on top of the over story layer
! is lower than 0.2 mm solve as if there were no snow.
! ....................................................................

         call nreb(1,albd,emiss,thermc,thermc2,heatcap,heatcap2,heatcapold,&
       f1par*f3vpd*f4temp*rescan+ravd,rahd,tkd,tkmidd,rnetd,xled,epetd,hd,gd,&
       dshd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rs_over,&
       rld,toleb,maxnri,dt,i)

      else

! ....................................................................
! If there is snow, solve the energy balance for the snow pack but
! do not change the snow pack variables yet, this call for the
! snow module is only done to estimate the potential evaporation.
! ....................................................................

         call inidum (dum1,PackWater,dum2,SurfWater,&
       dum3,Swq,dum4,VaporMassFlux,dum5,TPack,dum6,TSurf,&
       dum7,r_MeltEnergy,dum8,Outflow,dum9,xleact_snow,dum10,hact_snow,&
       dum11,rn_snow,dum12,dens)

         call calcsnowmelt(0,0,dt/3600.d0,za,zpd,z0h,RaSnow,roa,vppa,xlhv,&
       rs_over*(1.d0-alb_snow),rld,appa,rain,snow,tcel,vpsat-vppa,uzw,&
       dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11,1.d0,dum12,&
       gact)
!CVAL       0.d0)

! ....................................................................
! The sign convention of the snow model is heat fluxes postive
! downward, ours is positive upward.
! ....................................................................

         dum9=0.d0-dum9
         dum10=0.d0-dum10

! ....................................................................
! Assing the snow variables to the potential energy fluxes for the
! pixel.
! ....................................................................

         rnetd=dum11
         xled=dum9
         hd=dum10

! ....................................................................
! Calculate the snow pack temperature used to estimate the soil
! temperatures.
! ....................................................................

         tsnow=dum5+273.15d0
         if (dum3.lt.(0.005d0)) tsnow=dum6+273.15d0
         if (dum3.lt.(0.d0)) tsnow=tcel+273.15d0

! ....................................................................
! Calculate the soil temperatures under the snow pack.
! ....................................................................

         call nreb_snow(thermc,thermc2,heatcap,heatcap2,heatcapold,&
       tkd,tkmidd,tsnow,zdeep,Tdeepstep,zmid,dt,dum)

! ....................................................................
! Assing the snow variables to the potential energy fluxes for the
! pixel.
! ....................................................................

         epetd=dum9/(row*xlhv)
         gd=dum

      endif

      return

      end
