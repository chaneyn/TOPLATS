! ====================================================================
!
!               subroutine peteb_bs
!
! ====================================================================
!
! Solve the energy balance at potential rate for bare soil.
!
! ====================================================================

      subroutine peteb_bs(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       rain,snow,Swq,albd,emiss,ravd,rahd,tkd,tkmidd,rnetd,xled,epetd,hd,&
       gd,dshd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,rld,&
       toleb,maxnri,dt,i,tkw,tkmidw,rnetw,xlew,epetw,hw,gw,dshw,ravw,rahw,&
       PackWater,SurfWater,VaporMassFlux,TPack,TSurf,r_MeltEnergy,&
       Outflow,xleact_snow,hact_snow,rn_snow,dens,za,zpd,albw,&
       z0h,RaSnow,alb_snow,appa,vpsat,uzw,gact,row)

      implicit none
      include "SNOW.h"
      include "help/peteb_bs.h"

      if ( (Swq.le.(0.0002d0)).OR.(SNOW_RUN.eq.0) ) then

! --------------------------------------------------------------------
! If there is no snow solve the energy balance for bare soil at
! potential rate (no soil resistance assumed).
! --------------------------------------------------------------------

         call nreb(1,albd,emiss,thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       ravd,rahd,tkd,tkmidd,rnetd,xled,epetd,hd,gd,dshd,tcel,vppa,roa,psychr,&
       xlhv,zdeep,Tdeepstep,zmid,rsd,rld,toleb,maxnri,dt,i)
	
         call nreb(1,albw,emiss,thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       ravw,rahw,tkw,tkmidw,rnetw,xlew,epetw,hw,gw,dshw,tcel,vppa,roa,psychr,&
       xlhv,zdeep,Tdeepstep,zmid,rsd,rld,toleb,maxnri,dt,i)

      else

! --------------------------------------------------------------------
! If there is snow, solve the energy balance for the snow pack but
! do not change the snow pack variables yet, this call for the
! snow module is only done to estimate the potential evaporation.
! --------------------------------------------------------------------

         call inidum(dum1,PackWater,dum2,SurfWater,&
       dum3,Swq,dum4,VaporMassFlux,dum5,TPack,dum6,TSurf,&
       dum7,r_MeltEnergy,dum8,Outflow,dum9,xleact_snow,dum10,hact_snow,&
       dum11,rn_snow,dum12,dens)

         call calcsnowmelt(0,0,dt/3600.d0,za,zpd,z0h,RaSnow,roa,vppa,xlhv,&
       rsd*(1.d0-alb_snow),rld,appa,rain,snow,tcel,vpsat-vppa,uzw,&
       dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11,1.d0,dum12,&
       gact)
!CVAL       0.d0)

! ....................................................................
! The sign convention of the snow model is heat fluxes postive
! downward, ours is positive upward.
! ....................................................................

         dum9=0.d0-dum9
         dum10=0.d0-dum10
         epetd=dum9/(row*xlhv)
         epetw=dum9/(row*xlhv)

! ....................................................................
! Calculate the snow pack temperature used to estimate the soil
! temperatures.
! ....................................................................

         tsnow=dum5+273.15d0
         if (dum3.lt.(0.005d0)) tsnow=dum6+273.15d0
         if (dum3.lt.(0.d0)) tsnow=tcel+273.15d0

! --------------------------------------------------------------------
! Calculate the soil temperatures under the snow pack.
! --------------------------------------------------------------------

         call nreb_snow(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       tkd,tkmidd,tsnow,zdeep,Tdeepstep,zmid,dt,dum)

! ....................................................................
! Assing the snow variables to the potential energy fluxes for the
! pixel.
! ....................................................................

         xled=dum9
         hd=dum10
         rnetd=dum11
         gd=dum

         call nreb_snow(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       tkw,tkmidw,tsnow,zdeep,Tdeepstep,zmid,dt,dum)

! ....................................................................
! Assing the snow variables to the potential energy fluxes for the
! pixel.
! ....................................................................

         gw=dum
         xlew=dum9
         hw=dum10
         rnetw=dum11

      endif

      return

      end
