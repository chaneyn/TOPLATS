! ====================================================================
!
!                         subroutine peteb_2_dv
!
! ====================================================================
!
! Solve the energy balance for over story (dry vegetation) under the
! assumption that the radiation balances of both the over and under
! story are linked with each other through the long wave radiation
! terms.
!
! ====================================================================

      subroutine peteb_2_dv(rlwdn,rswdn,rlwup,ccc,&
       thermc2,heatcap,heatcap2,heatcapold,tkmidd_tmp,r_diff,tt,&
       rain,snow,Swq,tkmidd,f1,f2,f3,emiss,thermc,&
       f1par,f3vpd,f4temp,rescan,ravd,rahd,tkd,rnetd,xled,&
       epetd,hd,gd,dshd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,&
       zmid,rld,toleb,maxnri,dt,i,PackWater,SurfWater,&
       VaporMassFlux,TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow,&
       rn_snow,dens,alb_snow,za,zpd,z0h,RaSnow,appa,vpsat,uzw,gact,row)

      implicit none
      include "SNOW.h"
      include "help/peteb_2_dv.h"

      if ( (Swq.le.(0.0002d0)).OR.(SNOW_RUN.eq.1) ) then

! ....................................................................
! If the total snow water equivalent on top of the over story layer
! is lower than 0.2 mm solve as if there were no snow.
! ....................................................................

         tkmidd_tmp=tkmidd

         call nreb(1,f1+f2-f3,2.d0*emiss,thermc,thermc2,heatcap,heatcap2,&
       heatcapold,f1par*f3vpd*f4temp*rescan+ravd,rahd,tkd,tkmidd_tmp,rnetd,&
       xled,epetd,hd,gd,dshd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,&
       rswdn,rld+rlwup,toleb,maxnri,dt,i)

! ....................................................................
! Calculate the difference between the input over story skin
! temperature and the new estimated over story skin temperature.
! ....................................................................

         r_diff=tkd-tt

      else

! ....................................................................
! If there is snow, solve the energy balance for the snow pack but
! do not change the snow pack variables yet, this call for the
! snow module is only done to estimate the potential evaporation.
! ....................................................................

         tkmidd_tmp=tkmidd

! ....................................................................
! Assign the snow pack variables to dummy variables.
! ....................................................................

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
!CVAL                           0.d0)

! ....................................................................
! The sign convention of the snow model is heat fluxes postive
! downward, ours is positive upward.
! ....................................................................

         dum9=0.d0-dum9
         dum10=0.d0-dum10

! ....................................................................
! Calculate the snow pack temperature used to estimate the soil
! temperatures.
! ....................................................................

         tsnow=dum5+273.15d0
         if (dum3.lt.(0.005d0)) tsnow=dum6+273.15d0
         if (dum3.lt.(0.d0)) tsnow=tcel+273.15d0

! .................................................................
! Calculate the soil temperatures under the snow pack.
! ....................................................................

         call nreb_snow(thermc,thermc2,heatcap,heatcap2,heatcapold,&
       tkd,tkmidd_tmp,tsnow,zdeep,Tdeepstep,zmid,dt,dum)

! ....................................................................
! Assing the snow variables to the potential energy fluxes for the
! pixel.
! ....................................................................

         gd=dum
         rnetd=dum11
         xled=dum9
         hd=dum10
         epetd=dum9/(row*xlhv)

! ....................................................................
! Calculate the difference between the input over story skin
! temperature and the new estimated over story skin temperature.
! ....................................................................

         r_diff=dum6+273.15-tt

      endif

      return

      end
