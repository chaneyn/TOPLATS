! ====================================================================
!
!                      subroutine peteb_2_wv
!
! ====================================================================
!
! Solve the energy balance for over story (wet vegetation) under the
! assumption that the radiation balances of both the over and under
! story are linked with each other through the long wave radiation
! terms.
!
! ====================================================================

      subroutine peteb_2_wv(rlwdn,rswdn,rlwup,ccc,&
      thermc2,heatcap,heatcap2,heatcapold,tkmidw_tmp,tkw_tmp,r_diff,tt,&
      rain,snow,Swq,tkw,tkmidw,f1,f2,f3,thermc,dt,i,&
      emiss,ravw,rahw,rnetw,xlew,epetw,hw,gw,dshw,tcel,vppa,roa,&
      psychr,xlhv,zdeep,rld,Tdeepstep,zmid,toleb,maxnri,iopstab,tkact,&
      zww,za,uzw,zpd,z0m,press,tkel,rib,z0h,PackWater,&
      SurfWater,VaporMassFlux,TPack,TSurf,r_MeltEnergy,&
      Outflow,xleact_snow,hact_snow,rn_snow,dens,&
      alb_snow,RaSnow,appa,vpsat,gact,row)

      implicit none
      include "SNOW.h"
      include "help/peteb_2_wv.h"

      if ( (Swq.le.(0.0002d0)).OR.(SNOW_RUN.eq.0) ) then

! ....................................................................
! If the total snow water equivalent on top of the over story layer
! is lower than 0.2 mm solve as if there were no snow.
! ....................................................................

        tkw_tmp=tkw
        tkmidw_tmp=tkmidw

        do iter=1,2

           call nreb(1,f1+f2-f3,2.d0*emiss,thermc,thermc2,heatcap,heatcap2,&
      heatcapold,ravw,rahw,tkw_tmp,tkmidw_tmp,rnetw,xlew,epetw,hw,gw,dshw,&
      tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rswdn,rlwup+rld,toleb,&
      maxnri,dt,i)

! ....................................................................
! Check for large temperature differences in wet vegetation.
! If large, need to recompute richardson number and call nreb again.
! ....................................................................

           if (iopstab.eq.1.and.i.gt.1) then

             tktmp = tkact
             tkact = tkw_tmp
             call stabcor(zww,za,uzw,zpd,z0m,tkel,press,tkact,vppa,rib)
             rahw=calcra(uzw,zww,za,zpd,z0m,z0h,rib)
             ravw=rahw
             tkact = tktmp

           endif

        enddo

! ....................................................................
! Calculate the difference between the input over story skin
! temperature and the new estimated over story skin temperature.
! ....................................................................

        r_diff=tkw_tmp-tt

      else

! ....................................................................
! If there is snow, solve the energy balance for the snow pack but
! do not change the snow pack variables yet, this call for the
! snow module is only done to estimate the potential evaporation.
! ....................................................................

        tkmidw_tmp=tkmidw

        call inidum(dum1,PackWater,dum2,SurfWater,&
      dum3,Swq,dum4,VaporMassFlux,dum5,TPack,dum6,TSurf,&
      dum7,r_MeltEnergy,dum8,Outflow,dum9,xleact_snow,dum10,hact_snow,&
      dum11,rn_snow,dum12,dens)

        albsnow=f1+f2-f3
        albsnow=alb_snow
        call calcsnowmelt(0,0,dt/3600.d0,za,zpd,z0h,&
      RaSnow,roa,vppa,xlhv,rswdn*(1.d0-albsnow),rld+rlwup,appa,rain,&
      snow,tcel,vpsat-vppa,uzw,dum1,dum2,dum3,dum4,dum5,&
      dum6,dum7,dum8,dum9,dum10,dum11,2.d0,dum12,gact)
!CVAL                        0.d0)

! ....................................................................
! The sign convention of the snow model is heat fluxes postive
! downward, ours is positive upward.
! ....................................................................

        dum9=0.d0-dum9
        dum10=0.d0-dum10

! ....................................................................
! Calculate the difference between the input over story skin
! temperature and the new estimated over story skin temperature.
! ....................................................................

        r_diff=dum6+273.15-tt

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
      heatcapold,tkw,tkmidw_tmp,tsnow,zdeep,Tdeepstep,zmid,dt,dum)

! ....................................................................
! Assing the snow variables to the potential energy fluxes for the
! pixel.
! ....................................................................

        rnetw=dum11
        gw=dum
        xlew=dum9
        hw=dum10
        epetw=dum9/(row*xlhv)

      endif

      return

      end
