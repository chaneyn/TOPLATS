! ====================================================================
!
!                    subroutine land_us
!
! ====================================================================
!
! Solve the energy balance for the under story assuming the radiation
! balances of under and over story are not linked with each other.
!
! ====================================================================

      subroutine land_us(rain,snow,thermc1,thermc2,heatcap_moss,heatcap,&
       heatcap1,heatcap2,heatcapold,tkactd_us,tkmidactd_us,&
       tskinactd_moss,tkactd_moss,tkmidactd_moss,canclos,ievcon_us,&
       rnact_us,xleact_us,hact_us,gact_us,dshact_us,tkact_us,&
       tkmid_us,rnactd_us,rnetw_us,xleactd_us,xlew_us,hactd_us,hw_us,&
       gactd_us,gw_us,dshactd_us,dshw_us,tkw_us,tkmidw_us,xlai_us,&
       dc_us,fw_us,trlup,ipix,xlhv_ic,row,evtact_us,iffroz_us,tkmid,zmid,&
       zrzmax,smtmp,rzsm,tzsm,smold,rzsmold,tzsmold,iopthermc,thetar,thetas,&
       psic,bcbeta,quartz,ifcoarse,rocpsoil,cph2o,roa,cp,roi,thermc,inc_frozen,&
       rzdthetaudtemp,iopgveg,thermc_us,iopthermc_v,tcbeta_us,xlai,f3,albd_us,&
       emiss_us,ravd_us,rahd_us,rescan_us,tcel_ic,vppa_ic,roa_ic,psychr_ic,&
       zdeep,Tdeepstep,r_sdn,r_ldn,toleb,maxnri,dt,i,&
       rld,rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us,tkmidd_us,initer,&
       ievcon_moss,xleactd_moss,bsdew_moss,evtact_moss,thermc_moss,&
       r_mossm,tskinact_moss,tkact_moss,tkmid_moss,hactd_moss,gactd_moss,&
       dshactd_moss,rav_moss,rah_moss,r_moss_depth,alb_moss,&
       rnactd_moss,emiss_moss,eact_moss,rnet_pot_moss,xle_p_moss,h_p_moss,&
       g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,zmoss,&
       thetas_moss,rnact_moss,xleact_moss,hact_moss,gact_moss,dshact_moss,gact,&
       Swq_us,precip_u,za,zpd,z0h,RaSnow,alb_snow,appa,&
       vpsat_ic,uzw,PackWater_us,SurfWater_us,VaporMassFlux_us,&
       TPack_us,TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us,hact_snow_us,&
       rn_snow_us,dens_us,heatcap_us,tkel_ic,eps,ds_p_moss,i_und,i_moss,i_2l)

      implicit none
      include "SNOW.h"
      integer ievcon_us,ipix,iffroz_us,iopthermc,ifcoarse,inc_frozen,i
      integer iopgveg,iopthermc_v,initer,ievcon_moss,i_und,i_moss,i_2l
      integer maxnri
      real*8 rain,snow,thermc1,thermc2,heatcap_moss,heatcap,heatcap1
      real*8 heatcap2,heatcapold,tkactd_us,tkmidactd_us,tskinactd_moss
      real*8 tkactd_moss,tkmidactd_moss,canclos,rnact_us,xleact_us
      real*8 hact_us,gact_us,dshact_us,tkact_us,tkmid_us,rnactd_us
      real*8 rnetw_us,xleactd_us,xlew_us,hactd_us,hw_us,gactd_us,gw_us
      real*8 dshactd_us,dshw_us,tkw_us,tkmidw_us,xlai_us,dc_us,fw_us
      real*8 trlup,xlhv_ic,row,evtact_us,tkmid,zmid,zrzmax,smtmp
      real*8 rzsm,tzsm,smold,rzsmold,tzsmold,thetar,thetas,psic
      real*8 bcbeta,quartz,rocpsoil,cph2o,roa,cp,roi,thermc
      real*8 rzdthetaudtemp,thermc_us,tcbeta_us,xlai,f3,albd_us
      real*8 emiss_us,ravd_us,rahd_us,rescan_us,tcel_ic,vppa_ic
      real*8 roa_ic,psychr_ic,zdeep,Tdeepstep,r_sdn,r_ldn,toleb
      real*8 dt,rld,rnetd_us,xled_us,hd_us,gd_us,dshd_us
      real*8 tkd_us,tkmidd_us,xleactd_moss,bsdew_moss,evtact_moss
      real*8 thermc_moss,tskinact_moss,tkact_moss,tkmid_moss
      real*8 hactd_moss,gactd_moss,dshactd_moss,rav_moss,rah_moss
      real*8 r_moss_depth,alb_moss,rnactd_moss,emiss_moss,eact_moss
      real*8 rnet_pot_moss,xle_p_moss,h_p_moss,g_p_moss,tk_p_moss
      real*8 tkmid_p_moss,tskin_p_moss,zmoss,r_mossm,thetas_moss
      real*8 rnact_moss,xleact_moss,hact_moss,gact_moss,dshact_moss
      real*8 gact,Swq_us,precip_u,za,zpd,z0h,RaSnow
      real*8 alb_snow,appa,vpsat_ic,uzw,PackWater_us,SurfWater_us
      real*8 VaporMassFlux_us,TPack_us,TSurf_us,r_MeltEnergy_us
      real*8 Outflow_us,xleact_snow_us,hact_snow_us,rn_snow_us
      real*8 dens_us,heatcap_us,tkel_ic,eps,ds_p_moss
      real*8 tsnow,told0,told1,told2,hold0,hold1,hold2
      real*8 tcold0,tcold1,tcold2,dum
      real*8 zero,one,two,three,four,five,six
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/

      if ( ( (Swq_us.gt.(0.0002d0)).or.&
               ((tcel_ic.le.(0.d0)).and.&
               (dt*precip_u.gt.(0.0002d0)) ) ).AND.&
           (SNOW_RUN.eq.1) )then

! --------------------------------------------------------------------
! If the snow water equivalent for the under story is higher than 0.2
! millimeter or if there is precipitation under freezing air
! temperature (assumed to be snow) solve the snow energy balance for
! the under story.
! --------------------------------------------------------------------

! ....................................................................
! Calculate the partitioning of precipitation into rain and snow.
! This is purely bases on air temperature.
! ....................................................................

         call calcrain(tcel_ic,snow,rain,precip_u,dt)

! ....................................................................
! Solve the energy balance for the snow on top of the under story.
! ....................................................................

         call calcsnowmelt(0,0,dt/3600.d0,za,zpd,z0h,&
       RaSnow,roa_ic,vppa_ic,xlhv_ic,r_sdn*(1.d0-alb_snow),rld,appa,rain,&
       snow,tcel_ic,vpsat_ic-vppa_ic,uzw,PackWater_us,SurfWater_us,&
       Swq_us,VaporMassFlux_us,TPack_us,TSurf_us,r_MeltEnergy_us,&
       Outflow_us,xleact_snow_us,hact_snow_us,rn_snow_us,1.d0,dens_us,gact)

! ....................................................................
! The sign convention of the snow model is heat fluxes postive
! downward, ours is positive upward.
! ....................................................................

         xleact_snow_us=0.d0-xleact_snow_us
         hact_snow_us=0.d0-hact_snow_us

! ....................................................................
! Calculate the snow pack temperature used to estimate the soil
! temperatures.
! ....................................................................

         tsnow=TPack_us+273.15d0
         if (Swq_us.lt.(0.005d0)) tsnow=TSurf_us+273.15d0
         if (Swq_us.lt.(0.d0)) tsnow=tcel_ic+273.15d0

      else

! ....................................................................
! Set all snow variables for the under story to zero.
! ....................................................................

         call zero_snowvar(PackWater_us,SurfWater_us,Swq_us,VaporMassFlux_us,&
       TPack_us,TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us,hact_snow_us,&
       dens_us)

      endif

      return

      end
