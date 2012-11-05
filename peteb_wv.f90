! ====================================================================
!
!               subroutine peteb_wv
!
! ====================================================================
!
! Solve the energy balance at potential rate for wet vegetation, over
! story.
!
! ====================================================================

      subroutine peteb_wv(thermc2,heatcap,heatcap2,heatcapold,&
       rs_over,rain,snow,Swq,albw,emiss,&
       thermc,ravw,rahw,tkw,tkmidw,rnetw,xlew,epetw,hw,&
       gw,dshw,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,&
       zmid,rld,toleb,maxnri,dt,i,iopstab,tkact,zww,za,uzw,zpd,&
       z0m,tkel,press,rib,z0h,PackWater,SurfWater,&
       VaporMassFlux,TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,&
       hact_snow,rn_snow,dens,RaSnow,alb_snow,appa,vpsat,gact,row,ipix)

      implicit none
      include "SNOW.h"
      include "help/peteb_wv.h"
      integer :: ipix

      if ( (Swq.le.(0.0002d0)).OR.(SNOW_RUN.eq.0) ) then

! ....................................................................
! If the total snow water equivalent on top of the over story layer
! is lower than 0.2 mm solve as if there were no snow.
! ....................................................................
         do iter=1,2

        !if (ipix.eq.500)print*,"peteb_wv.f90",albw,emiss,thermc,thermc2,heatcap,heatcap2,heatcapold,&
       !ravw,rahw,tkw,tkmidw,rnetw,xlew,epetw,hw,gw,dshw,tcel,vppa,roa,psychr,&
       !xlhv,zdeep,Tdeepstep,zmid,rs_over,rld,toleb,maxnri,dt,i
            call nreb(1,albw,emiss,thermc,thermc2,heatcap,heatcap2,heatcapold,&
       ravw,rahw,tkw,tkmidw,rnetw,xlew,epetw,hw,gw,dshw,tcel,vppa,roa,psychr,&
       xlhv,zdeep,Tdeepstep,zmid,rs_over,rld,toleb,maxnri,dt,i)
        !if (ipix .eq.500)print*,"peteb_wv.f90",albw,emiss,thermc,thermc2,heatcap,heatcap2,heatcapold,&
       !ravw,rahw,tkw,tkmidw,rnetw,xlew,epetw,hw,gw,dshw,tcel,vppa,roa,psychr,&
       !xlhv,zdeep,Tdeepstep,zmid,rs_over,rld,toleb,maxnri,dt,i

! ....................................................................
! Check for large temperature differences in wet vegetation.
! If large, need to recompute richardson number and call nreb again.
! ....................................................................

            if (iopstab.eq.1.and.i.gt.1) then

               tktmp = tkact
               tkact = tkw
               call stabcor(zww,za,uzw,zpd,z0m,tkel,press,tkact,vppa,rib)
               rahw=calcra(uzw,zww,za,zpd,z0m,z0h,rib)
               ravw=rahw
               tkact = tktmp

            endif

         enddo

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

         call calcsnowmelt(0,0,dt/3600.d0,za,zpd,z0h,&
       RaSnow,roa,vppa,xlhv,rs_over*(1.d0-alb_snow),rld,appa,rain,&
       snow,tcel,vpsat-vppa,uzw,dum1,dum2,dum3,dum4,dum5,&
       dum6,dum7,dum8,dum9,dum10,dum11,1.d0,dum12,gact)
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

         rnetw=dum11
         xlew=dum9
         hw=dum10
         epetw=dum9/(row*xlhv)

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

         call nreb_snow(thermc,thermc2,heatcap,heatcap2,&
       heatcapold,tkw,tkmidw,tsnow,zdeep,Tdeepstep,zmid,dt,dum)
         gw=dum

      endif

      return

      end
