! ====================================================================
!
!                    subroutine land_os
!
! ====================================================================
!
! Solve the energy balance for the over story assuming the radiation
! balances of under and over story are not linked with each other.
!
! ====================================================================

      subroutine land_os(rain,snow,thermc2,heatcap,heatcap2,heatcapold,&
       tkactd,tkmidactd,canclos,ievcon,xlhv,row,ivgtyp,xleactd,evtact,&
       bsdew,ioppet,iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,&
       tzsm,smold,rzsmold,tzsmold,iopthermc,thermc1,thetar,&
       thetas,psic,bcbeta,quartz,ifcoarse,heatcap1,rocpsoil,cph2o,&
       roa,cp,roi,thermc,rzdthetaudtemp,iopgveg,iopthermc_v,tcbeta,&
       xlai,tkact,i_2l,f1,f2,f3,emiss,rescan,ravd,rahd,rnactd,&
       hactd,gactd,dshactd,tcel,vppa,psychr,zdeep,Tdeepstep,&
       rsd,r_lup,rld,toleb,maxnri,dt,i,albd,r_sdn,rnetpn,&
       gbspen,rnetd,xled,hd,gd,dshd,tkd,tkmidd,rnact,xleact,hact,&
       gact,dshact,rnetw,xlew,hw,gw,dshw,tkw,tkmidw,dc,fw,tdiff,inc_frozen,&
       ipix,initer,PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,dens,precip_o,za,&
       zpd,z0h,RaSnow,appa,vpsat,uzw,rn_snow,alb_snow)

      implicit none
      include "SNOW.h"
      include "help/land_os.h"

      if ( ( (Swq.gt.(0.0002d0)).or.&
               ((tcel.le.(0.d0)).and.(dt*precip_o.gt.(0.0002d0)) )).AND.&
           (SNOW_RUN.eq.1) ) then

! --------------------------------------------------------------------
! If the snow water equivalent on the overstory is higher than 0.2
! millimeter or if there is rainfall when freezing air temperature
! (in these cases assumed to be snowfall) than solve the energy
! balance with the snow model.
! --------------------------------------------------------------------

! ....................................................................
! Calculate the partitioning of precipitation into rain and snow.
! This is purely bases on air temperature.
! ....................................................................

         call calcrain (tcel,snow,rain,precip_o,dt)

! ....................................................................
! Solve the energy balance for the snow on top of the over story.
! ....................................................................

         call calcsnowmelt(0,0,dt/3600.d0,za,zpd,z0h,RaSnow,roa,vppa,xlhv,&
       r_sdn*(1.d0-alb_snow),rld,appa,rain,snow,tcel,vpsat-vppa,uzw,&
       PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,1.d0,dens,gact)
!CVAL                           rn_snow,1.d0,dens,0.d0)

! ....................................................................
! The sign convention of the snow model is heat fluxes postive
! downward, ours is positive upward.
! ....................................................................

         xleact_snow=0.d0-xleact_snow
         hact_snow=0.d0-hact_snow

! ....................................................................
! Calculate the snow pack temperature used to estimate the soil
! temperatures.
! ....................................................................

         tsnow=TPack+273.15d0
         if (Swq.lt.(0.005d0)) tsnow=TSurf+273.15d0
         if (Swq.lt.(0.d0)) tsnow=tcel+273.15d0

! ....................................................................
! Calculate the ground heat flux and the soil temperature under
! the snow pack for the over story.
! ....................................................................

         call nreb_snow(thermc,thermc2,heatcap,heatcap2,heatcapold,&
       tkactd,tkmidactd,tsnow,zdeep,Tdeepstep,zmid,dt,gactd)

! ....................................................................
! Assing the snow variables to the actual energy fluxes for the pixel.
! ....................................................................

         rnact=rn_snow
         xleact=xleact_snow
         hact=hact_snow
         gact=gactd
         dshact=0.d0
         tkact=tkactd
         tkmid=tkmidactd

!         write (390,*) gact,DeltaColdContent
!         write (391,*) i,gact
!         write (392,*) i,DeltaColdContent

      else

! --------------------------------------------------------------------
! If there is no snow fall or no snow pack present solve the
! energy balance for the over story vegetation and set all snow
! variables to zero (this is only a double check).
! --------------------------------------------------------------------

         call engact(canclos,ievcon,xlhv,row,ivgtyp,xleactd,evtact,&
       bsdew,ioppet,iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,&
       tzsm,smold,rzsmold,tzsmold,iopthermc,thermc1,thermc2,thetar,&
       thetas,psic,bcbeta,quartz,ifcoarse,heatcap1,heatcap2,&
       heatcapold,rocpsoil,cph2o,roa,cp,roi,thermc,heatcap,rzdthetaudtemp,&
       iopgveg,iopthermc_v,tcbeta,xlai,tkact,tkactd,&
       tkmidactd,i_2l,f1,f2,f3,emiss,rescan,ravd,rahd,rnactd,&
       hactd,gactd,dshactd,tcel,vppa,psychr,zdeep,Tdeepstep,&
       rsd,r_lup,rld,toleb,maxnri,dt,i,albd,r_sdn,rnetpn,&
       gbspen,rnetd,xled,hd,gd,dshd,tkd,tkmidd,rnact,xleact,hact,&
       gact,dshact,rnetw,xlew,hw,gw,&
       dshw,tkw,tkmidw,dc,fw,tdiff,inc_frozen,ipix,initer)

! ....................................................................
! Set all snow variables to zero.
! ....................................................................

         call zero_snowvar(PackWater,SurfWater,Swq,VaporMassFlux,&
       TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow,dens)

      endif

      return

      end
