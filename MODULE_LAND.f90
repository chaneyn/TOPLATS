MODULE MODULE_LAND

  USE VARIABLES

  USE MODULE_SNOW

  implicit none

  contains

! ====================================================================
!
!                       subroutine land
!
! ====================================================================
!
! Subroutine land calculates the land surface water balance.
!
! ====================================================================
    subroutine land(newstorm,ipix,i,dt,inc_frozen,i_2l,&

! Factor to multiply the regional parameters with

       mul_fac,&

! General vegetation parameters

       canclos,extinct,i_und,i_moss,ivgtyp,f_moss,f_und,&

! Snow pack variables

       PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,PackWater_us,&
       SurfWater_us,Swq_us,VaporMassFlux_us,TPack_us,&
       TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us,&
       hact_snow_us,rn_snow_us,dens,dens_us,&

! Albedos of the over story, under story,&
! and moss layer

       albd_us,alb_moss,alb_snow,albd,&

! Meteorological data

       rsd,rld,tcel,vppa,psychr,xlhv,tkel,zww,za,uzw,press,&
       pptms,appa,vpsat,tcel_ic,vppa_ic,psychr_ic,&
       xlhv_ic,tkel_ic,vpsat_ic,precip_o,precip_u,&

! Temperature variables

       tkmid,tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,tkmidpet,tkmidpet_us,tkmidpet_moss,tsoilold,Tdeepstep,&

! Energy fluxes

       dshact,rnetpn,gbspen,epetd,evtact,ievcon,bsdew,gact,rnact,xleact,&
       hact,epetd_us,dshact_moss,xle_act_moss,rnetd,xled,hd,&
       gd,dshd,tkd,tkmidd,rnetw,xlew,hw,gw,dshw,tkw,tkmidw,ievcon_us,rnact_us,&
       xleact_us,hact_us,gact_us,dshact_us,rnetw_us,xlew_us,hw_us,gw_us,&
       dshw_us,tkw_us,tkmidw_us,evtact_us,rnetd_us,xled_us,&
       hd_us,gd_us,dshd_us,tkd_us,tkmidd_us,ievcon_moss,bsdew_moss,&
       evtact_moss,rnet_pot_moss,xle_p_moss,h_p_moss,g_p_moss,tk_p_moss,&
       tkmid_p_moss,tskin_p_moss,eact_moss,rnact_moss,&
       xleact_moss,hact_moss,gact_moss,ds_p_moss,&

! Soil parameters

       thetar,thetas,psic,bcbeta,quartz,ifcoarse,rocpsoil,tcbeta,&
       tcbeta_us,bulk_dens,a_ice,b_ice,xk0,bcgamm,srespar1,srespar2,&
       srespar3,zdeep,zmid,zrzmax,&

! Moss parameters

       r_moss_depth,thetas_moss,srespar1_moss,srespar2_moss,srespar3_moss,&
       eps,emiss_moss,zpd_moss,rib_moss,z0m_moss,z0h_moss,epet_moss,&
       a_ice_moss,b_ice_moss,bulk_dens_moss,&

! Vegetation parameters

       xlai,xlai_us,emiss,zpd,zpd_us,z0m,z0h,f1par,f3vpd,f4temp,f1par_us,&
       f3vpd_us,f4temp_us,rescan,tc,tw,tc_us,tw_us,rescan_us,rtact,&
       rtdens,psicri,respla,f1,f2,f3,emiss_us,&

! Constants

       row,cph2o,roa,cp,roi,toleb,maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,rib,RaSnow,&

! Water balance variables

       rzsm,tzsm,rzsm1,tzsm1,rzsm_u,tzsm_u,rzsm1_u,tzsm1_u,rzsm_f,tzsm_f,&
       rzsm1_f,tzsm1_f,r_mossm,r_mossm1,r_mossm_f,r_mossm1_f,r_mossm_u,&
       r_mossm1_u,zrz,ztz,r_mossmold,smold,rzsmold,tzsmold,rzdthetaudtemp,&
       rzdthetaidt,tzdthetaidt,zw,zbar,zmoss,capflx,difrz,diftz,grz,gtz,pnet,&
       cuminf,sorp,cc,deltrz,xinact,satxr,xinfxr,runtot,irntyp,sesq,&
       corr,idifind,dc,fw,dc_us,fw_us,wcip1,par,dewrun,dsrz,rzrhs,&
       dstz,tzrhs,&

! Storm parameters

       istmst,intstm,istmst_moss,intstm_moss,intstp,istorm,&

! Topmodel parameters

       ff,atanb,xlamda,&


! Regional saturation parameters

       fwcat,fwreg,pr3sat,perrg2,pr2sat,pr2uns,perrg1,pr1sat,&
       pr1rzs,pr1tzs,pr1uns,persxr,perixr,persac,peruac,perusc,&

! Different option paramters

       iopthermc,iopgveg,iopthermc_v,iopsmini,ikopt,irestype,ioppet,&
       iopveg)

      implicit none
      include "SNOW.h" !Switch to use snow module when it is finished
      include "help/land.h"!Remove when variables are changed
      real*8 gold

      data tolinf/1.0d-09/
      initer=2

! ====================================================================
! Initialize the rain and snowfall.
! ====================================================================

      rain=0.d0
      snow=0.d0
      f1=0.d0
      f2=0.d0
      f3=0.d0

      gold=gact

! ====================================================================
! Calculate the incoming solar radiation for the under story and over
! story layers.
! ====================================================================

      call calc_rs(canclos,extinct,i_und,i_moss,Swq_us,&
                   albd_us,alb_moss,alb_snow,rsd,rs_over,rs_under)

! ====================================================================
! Initialize soil moisture for the calculation of the thermodynami!
! parameters, as a centered difference.
! ====================================================================

      call sm_cen_dif(iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm,smold,&
                      rzsmold,tzsmold)

! ====================================================================
! Calculate the soil thermal parameters.
! ====================================================================

      call soiltherm(iopthermc,thermc1,thermc2,rzsm,smtmp,&
       thetar,thetas,psic,bcbeta,tkmid,quartz,ifcoarse,&
       heatcap1,heatcap2,heatcapold,rocpsoil,row,cph2o,roa,cp,roi,&
       smold,thermc,heatcap,inc_frozen,rzdthetaudtemp)

! ====================================================================
! Modify the thermal parameters for soils under vegetation.
! ====================================================================

      if (ivgtyp.ne.0) then
         call soiladapt(iopgveg,thermc,iopthermc_v,tcbeta,&
                        xlai,thermc1,heatcap,heatcap1,zero)

      endif

      heatcap_us=heatcap
      thermc_us=thermc
      !call soiladapt(iopgveg,thermc_us,iopthermc_v,tcbeta_us,&
      !               xlai_us,thermc1,heatcap,heatcap1,zero)

! ====================================================================
! Initialize actual temperatures.
! ====================================================================

      tkactd = tkact
      tkmidactd = tkmid
      tkactd_us = tkact_us
      tkmidactd_us = tkmid_us
      tskinactd_moss = tskinact_moss
      tkactd_moss = tkact_moss
      tkmidactd_moss = tkmid_moss

! ====================================================================
! Determine the frozen and unfrozen part of the soil water if
! the representation of frozen soil processes is requested.
! ====================================================================

      rzdthetaidt=0.d0
      tzdthetaidt=0.d0

      if (inc_frozen.eq.1) then

         call ice_change(ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,bulk_dens,&
       a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,thetas,thetar,&
       rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold)

      endif

      rzsm_f=rzsm1_f
      tzsm_f=tzsm1_f

! ====================================================================
! Update local water table depth, root and transmission zone soil
! moisture.
! ====================================================================

      call states(zw0,inc_frozen,i_moss,0.5d0*(tskinact_moss+tskinact_moss),&
       r_mossm_u,r_mossm_f,r_mossm,zw,zbar,ff,atanb,xlamda,psic,&
       zrz,ztz,rzsm1,tzsm1,thetas,zrzmax,iopsmini,thetar,bcbeta,rzsm1_u,&
       tzsm1_u,rzsm1_f,tzsm1_f,tsoilold,bulk_dens,a_ice,b_ice,&
       row,rzsmold,tzsmold,r_mossmold,rzsm,tzsm,r_mossm1,zmoss,r_moss_depth,&
       thetas_moss,rzsm_u,rzsm_f,tzsm_u,tzsm_f,r_mossm1_u,&
       r_mossm1_f,i,a_ice_moss,b_ice_moss,bulk_dens_moss)

! ====================================================================
! Calculate the infiltration.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : No moss layer on the surface, infiltration capacity
!            depends on cumulative infiltration, rainfall, ...
! --------------------------------------------------------------------

         call infilt(pnet,i_moss,i_und,PackWater_us,SurfWater_us,Swq_us,&
       Outflow_us,dt,PackWater,SurfWater,Swq,Outflow,&
       istmst,cuminf,inc_frozen,rzsmst,rzsm,rzsm_u,thetas,thetar,&
       tolinf,sorp,xk0,psic,bcgamm,bcbeta,deltrz,cc,zw,xinact,satxr,xinfxr,&
       intstm,xinfcp,runtot,irntyp)

! ====================================================================
! Calculate actual rate of evaporation.
! ====================================================================

      if (ivgtyp.eq.0) then

! --------------------------------------------------------------------
! Option 1 : Bare soil surface.
! --------------------------------------------------------------------

         if ( (Swq.le.(0.d0)).or.(SNOW_RUN.eq.0) ) then

! --------------------------------------------------------------------
! In case of absence of a snow pack, solve the energy balance for
! the bare soil.
! --------------------------------------------------------------------
            call ebsres(inc_frozen,irestype,rsoil,rzsm,srespar1,tkact,&
       srespar2,rzsm_u,srespar3,ravd,iffroz,thetas,tkmid,&
       zmid,zrzmax,smtmp,tzsm,smold,rzsmold,tzsmold,&
       iopthermc,thermc1,thermc2,thetar,heatcapold,psic,bcbeta,&
       quartz,heatcap1,ifcoarse,heatcap2,rocpsoil,row,cph2o,roa,cp,&
       roi,thermc,heatcap,rzdthetaudtemp,dshact,albd,emiss,rahd,ebscap,&
       tcel,vppa,psychr,xlhv,zdeep,Tdeepstep,rsd,rld,toleb,maxnri,dt,i,tkel,&
       zww,za,uzw,zpd,z0m,press,rib,rnetpn,gbspen,epetd,evtact,ievcon,&
       bsdew,z0h,ioppet)

         else

! --------------------------------------------------------------------
! If there is snow, the evaporation is determined by the solution of
! the snow model.
! --------------------------------------------------------------------

            ebscap=epetd
            evtact = ebscap
            ievcon = 3

            call calcrain (tcel,snow,rain,pptms,dt)

            call calcsnowmelt(0,0,dt/3600.d0,za,zpd,z0h,RaSnow,roa,vppa,xlhv,&
            rsd*(1.d0-alb_snow),rld,appa,rain,snow,tcel,vpsat-vppa,uzw,&
            PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,r_MeltEnergy,&
            Outflow,xleact_snow,hact_snow,rn_snow,1.d0,dens,gold)

            xleact_snow=0.d0-xleact_snow
            hact_snow=0.d0-hact_snow
            tsnow=TPack+273.15d0

            if (Swq.lt.(0.005d0)) tsnow=TSurf+273.15d0
            if (Swq.lt.(0.d0)) tsnow=tcel+273.15d0

            call nreb_snow(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       tkactd,tkmidactd,tsnow,zdeep,Tdeepstep,zmid,dt,dum)

            rnact=rn_snow
            xleact=xleact_snow
            hact=hact_snow
            gact=gactd
            dshact=0.d0
            tkact=tkactd
            tkmid=tkmidactd

         endif

      else

! --------------------------------------------------------------------
! Option 2: Vegetated surface.
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Calculate the transpiration from dry canopy for vegetated pixels.
! --------------------------------------------------------------------

         call transv(epetd,epetd_us,i_und,iopveg,f1par,f3vpd,f4temp,&
       rescan,inc_frozen,ivgtyp,rzsm,rzsm_u,tc,tw,smcond,tzsm,tzsm_u,&
       tc_us,tw_us,smcond_us,f1par_us,f3vpd_us,f4temp_us,rescan_us,&
       vegcap,ravd,vegcap_us,ravd_us,zrz,srzrel,thetas,thetar,psisoi,&
       psic,bcbeta,ikopt,xksrz,xk0,ff,ressoi,rtact,rtdens,psicri,&
       respla,xkrz,ztz,stzrel,xkstz,xktz,Swq,evtact,ievcon,Swq_us,evtact_us,&
       ievcon_us,bsdew,i,ipix)

      endif

! ====================================================================
! balance transmission and surface root zones
! ====================================================================

      if (inc_frozen.eq.0) then

         rzsm_test=rzsm
         tzsm_test=tzsm
         rzsm_u_test=rzsm
         tzsm_u_test=tzsm
         thetas_add=thetas

      endif

      if (inc_frozen.eq.1) then

         rzsm_test=rzsm_u
         tzsm_test=tzsm_u
         rzsm_u_test=rzsm_u
         tzsm_u_test=tzsm_u
         thetas_add=thetas-rzsm_f

      endif

      call tz_and_rzbal(i,newstorm,inc_frozen,ikopt,ivgtyp,&
       dt,rzsm_test,tzsm_test,rzsm_u_test,tzsm_u_test,rzsm1,tzsm1,&
       zrz,ztz,zw0,zrzmax,&
       evtact,evtact_us,bsdew,dewrun,grz,gtz,diftz,difrz,&
       satxr,runtot,xinact,cuminf,&
       ff,thetar,thetas_add,bcbeta,xk0,psic,&
       Swq,Swq_us,&
       dc,i_und,i_moss,fw,dc_us,fw_us,evrz_moss,f_und,dstz,dsrz,&
       tzrhs,rzrhs)

      if (inc_frozen.eq.1) then

         rzsm_u=rzsm1
         tzsm_u=tzsm1
         rzsm1_u=rzsm1
         tzsm1_u=tzsm1
         rzsm1=rzsm1_u+rzsm_f
         tzsm1=tzsm1_u+tzsm_f

      endif

! ====================================================================
! Calculate the percent of land surface in various surface states.
! ====================================================================

      call sursat(fwcat,fwreg,zw,psic,fw,mul_fac,pr3sat,zrzmax,perrg2,rzsm1,&
       thetas,pr2sat,pr2uns,perrg1,tzsm1,pr1sat,pr1rzs,pr1tzs,pr1uns,satxr,&
       persxr,xinfxr,perixr,ievcon,persac,peruac,perusc,i,ipix)

! ====================================================================
! Calculate the actual surface energy fluxes - if PET values
! are used then skip this.
! ====================================================================

      if (ioppet.ne.2) then

! --------------------------------------------------------------------
! Calculate the incoming solar radiation for the overstory.
! --------------------------------------------------------------------

         r_sdn=rs_over
         dshactd=zero
! --------------------------------------------------------------------
! Solve the over story energy balance.
! --------------------------------------------------------------------

         call land_os(rain,snow,thermc2,heatcap,heatcap2,heatcapold,&
       tkactd,tkmidactd,canclos,ievcon,xlhv,row,ivgtyp,xleactd,evtact,&
       bsdew,ioppet,iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,&
       tzsm,smold,rzsmold,tzsmold,iopthermc,thermc1,thetar,thetas,psic,&
       bcbeta,quartz,ifcoarse,heatcap1,rocpsoil,cph2o,roa,cp,roi,thermc,&
       rzdthetaudtemp,iopgveg,iopthermc_v,tcbeta,xlai,tkact,&
       i_2l,f1,f2,f3,emiss,rescan,ravd,rahd,rnactd,&
       hactd,gactd,dshactd,tcel,vppa,psychr,zdeep,Tdeepstep,&
       rsd,r_lup,rld,toleb,maxnri,dt,i,albd,r_sdn,rnetpn,&
       gbspen,rnetd,xled,hd,gd,dshd,tkd,tkmidd,rnact,xleact,hact,&
       gact,dshact,rnetw,xlew,hw,gw,dshw,tkw,tkmidw,dc,fw,tdiff,inc_frozen,&
       ipix,initer,PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,dens,precip_o,za,&
       zpd,z0h,RaSnow,appa,vpsat,uzw,rn_snow,alb_snow)

! --------------------------------------------------------------------
! Calculate the incoming solar radiation for the under story and
! solve the under story energy balance.
! --------------------------------------------------------------------

         r_sdn=rs_under

         call land_us(rain,snow,thermc1,thermc2,heatcap_moss,heatcap,&
       heatcap1,heatcap2,heatcapold,tkactd_us,tkmidactd_us,&
       tskinactd_moss,tkactd_moss,tkmidactd_moss,canclos,ievcon_us,&
       rnact_us,xleact_us,hact_us,gact_us,dshact_us,tkact_us,&
       tkmid_us,rnactd_us,rnetw_us,xleactd_us,xlew_us,hactd_us,hw_us,&
       gactd_us,gw_us,dshactd_us,dshw_us,tkw_us,tkmidw_us,xlai_us,&
       dc_us,fw_us,trlup,ipix,xlhv_ic,row,evtact_us,iffroz_us,tkmid,zmid,&
       zrzmax,smtmp,rzsm,tzsm,smold,rzsmold,tzsmold,iopthermc,&
       thetar,thetas,psic,bcbeta,quartz,ifcoarse,rocpsoil,cph2o,roa,cp,roi,&
       thermc,inc_frozen,rzdthetaudtemp,iopgveg,thermc_us,iopthermc_v,tcbeta_us,&
       xlai,f3,albd_us,emiss_us,ravd_us,rahd_us,rescan_us,tcel_ic,vppa_ic,&
       roa_ic,psychr_ic,zdeep,Tdeepstep,r_sdn,r_ldn,toleb,maxnri,dt,i,&
       rld,rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us,tkmidd_us,initer,&
       ievcon_moss,xleactd_moss,bsdew_moss,evtact_moss,thermc_moss,&
       r_mossm,tskinact_moss,tkact_moss,tkmid_moss,hactd_moss,gactd_moss,&
       dshactd_moss,rav_moss,rah_moss,r_moss_depth,alb_moss,&
       rnactd_moss,emiss_moss,eact_moss,rnet_pot_moss,xle_p_moss,h_p_moss,&
       g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,zmoss,&
       thetas_moss,rnact_moss,xleact_moss,hact_moss,&
       gact_moss,dshact_moss,gold,Swq_us,precip_u,za,zpd_us,&
       z0h,RaSnow,alb_snow,appa,vpsat_ic,uzw,PackWater_us,&
       SurfWater_us,VaporMassFlux_us,TPack_us,TSurf_us,&
       r_MeltEnergy_us,Outflow_us,xleact_snow_us,hact_snow_us,&
       rn_snow_us,dens_us,heatcap_us,tkel_ic,eps,ds_p_moss,i_und,i_moss,i_2l)

      endif

! ====================================================================
! Print some results out.
! ====================================================================

      if (dens_us.gt.(0.d0)) then

         p1=PackWater_us*1000.d0
         p2=SurfWater_us*1000.d0
         p3=Swq_us*1000.d0
         r_dens=(p3/(p1+p2+p3))*dens_us+((p1+p2)/(p1+p2+p3))*1.d0

      endif

      if (dens.gt.(0.d0)) then

         p1=PackWater*1000.d0
         p2=SurfWater*1000.d0
         p3=Swq*1000.d0
         r_dens=(p3/(p1+p2+p3))*dens+((p1+p2)/(p1+p2+p3))*1.d0

      endif

5432  format (1i5," ",f11.6," ",2(f5.1," "),f7.3," ",2(f7.2," "),2(f6.3," "),f5.1,f6.1)

! ====================================================================
! Check whether the snow water equivalent on top of the overstory
! does not exceed its water holding capacity.
! ====================================================================

      if (ivgtyp.gt.0) then

            if (Swq.ge.wcip1) wcip1=Swq

      endif

      return
    end subroutine land

! ====================================================================
!
!                       subroutine calc_rs
!
! ====================================================================
!
! Calculate the incoming solar radiation for the under and over story
! under the assumption of only one reflection.
!
! ====================================================================

    subroutine calc_rs(canclos,extinct,i_und,i_moss,Swq_us,&
                         albd_us,alb_moss,alb_snow,rsd,rs_over,rs_under)

      implicit none
      include "help/calc_rs.h"

! ====================================================================
! Calculate the incoming solar radiation for the under story and over
! story layers.
! ====================================================================

      refus=0.d0
      ccc=canclos
      thr=extinct

! --------------------------------------------------------------------
! Determine what albedo of the understory is : moss, snow or normal
! vegetation.
! --------------------------------------------------------------------

      if (i_und.gt.0) refus=albd_us
      if (i_moss.gt.0) refus=alb_moss
      if (Swq_us.gt.(0.d0)) refus=alb_snow

! --------------------------------------------------------------------
! Calculate the incoming radiation under the assumption of only
! one reflection.
! --------------------------------------------------------------------

      if ( (i_und.gt.0).or.(i_moss.gt.0) ) then

         rs_over=(1.d0+refus*thr)*rsd
         rs_under=rsd*(thr*ccc+1.d0-ccc)

      endif

      if ( (i_und.eq.0).and.(i_moss.eq.0) ) then

         rs_over=rsd
         rs_under=rsd

      endif

      return
    end subroutine calc_rs

! ====================================================================
!
!            subroutine sm_cen_dif
!
! ====================================================================
!
! Initialize soil moisture for the calculation of the thermodynami!
! parameters, as a centered difference.
!
! ====================================================================

    subroutine sm_cen_dif(iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm,smold,&
                            rzsmold,tzsmold)

      implicit none
      include "help/sm_cen_dif.h"

      iffroz=0
      if (tkmid.lt.273.15) iffroz=1

      if (zmid.ge.zrzmax) then

         smtmp=0.5*rzsm + 0.5*tzsm
         smold=0.5*rzsmold + 0.5*tzsmold

      else

         smtmp=rzsm
         smold=rzsmold

      endif

      return

    end subroutine sm_cen_dif

! ====================================================================
!
!                     subroutine soiltherm
!
! ====================================================================
!
! Calculate the soil thermal parameters.
!
! ====================================================================

    subroutine soiltherm(iopthermc,thermc1,thermc2,rzsm,smtmp,&
       thetar,thetas,psic,bcbeta,tkmid,quartz,ifcoarse,&
       heatcap1,heatcap2,heatcapold,rocpsoil,row,cph2o,roa,cp,roi,&
       smold,thermc,heatcap,inc_frozen,rzdthetaudtemp)

      implicit none
      include "help/soiltherm.h"

! ====================================================================
! Calculate the termal conductivity.
! ====================================================================

      if (iopthermc.eq.1) then

! --------------------------------------------------------------------&
! McCumber-Pielke method
! --------------------------------------------------------------------&

         thermc1 = calctc_m(rzsm,thetar,thetas,psic,bcbeta)
         thermc2 = calctc_m(smtmp,thetar,thetas,psic,bcbeta)

      else

! --------------------------------------------------------------------&
! Johansen's method
! --------------------------------------------------------------------&

         iffroz=0
         if (tkmid.lt.273.15) iffroz=1
         thermc1 = calctc_j(rzsm,thetar,thetas,quartz,iffroz,ifcoarse)
         thermc2 = calctc_j(smtmp,thetar,thetas,quartz,iffroz,ifcoarse)

      endif

! ====================================================================
! Calculate the heat capacity.
! ====================================================================

      heatcap1 = calchc(rzsm,thetas,rocpsoil,row,cph2o,roa,cp,tkmid,roi,&
                        inc_frozen,rzdthetaudtemp)

! --------------------------------------------------------------------&
! For the heat capacity of the transmission zone : tzsmpet will
! lead to overestation of the heat capacity, use average
! --------------------------------------------------------------------&

      heatcap2 = calchc(smtmp,thetas,rocpsoil,row,cph2o,roa,cp,tkmid,roi,&
                        inc_frozen,rzdthetaudtemp)
      heatcapold = calchc(smold,thetas,rocpsoil,row,cph2o,roa,cp,tkmid,roi,&
                          inc_frozen,rzdthetaudtemp)
      thermc=thermc1
      heatcap=heatcap1

      return

    end subroutine soiltherm

! ====================================================================
!
!             subroutine soiladapt
!
! ====================================================================
!
! Adapt the thermal parameters for vegetated surfaces.
!
! ====================================================================

    subroutine soiladapt(iopgveg,thermc,iopthermc_v,tcbeta,&
                           xlai,thermc1,heatcap,heatcap1,zero)

      implicit none
      include "help/soiladapt.h"

! ====================================================================
! Check the thermal conductivity inputs.
! ====================================================================

      if ( (thermc.lt.100.d0).and.(thermc.ge.0.d0) )then

         xlai=xlai

      else

         write (*,*) 'Input therm! unrealisti! in SOILADAPT ',thermc
         write (*,*) iopgveg,thermc,iopthermc_v,tcbeta,&
                     xlai,thermc1,heatcap,heatcap1,zero
         stop

      endif

! ====================================================================
! Modify the thermal parameters for soils under vegetation.
! ====================================================================

      if (iopgveg.eq.0) then

! --------------------------------------------------------------------&
! Option 1 : Assume no ground heat flux under vegetation.
! --------------------------------------------------------------------&

         thermc = zero

      else

! --------------------------------------------------------------------&
! Option 2 : Assume ground heat flux under vegetation, and an
!            exponential decay in thermal conducivity of the soil
!            under vegetation (Choudhury et al., 1987) with LAI.
! -------------------------------------------------------------------&

         if (iopthermc_v.ne.1) then

            tau = dexp(-tcbeta * xlai)
            thermc = tau * thermc1
            heatcap = heatcap1

         else

            thermc = dexp(-tcbeta*xlai)*7.0
            heatcap = 0.d0

         endif

      endif

! ====================================================================
! Check the thermal conductivity outputs.
! ====================================================================

      if ( (thermc.lt.100.d0).and.(thermc.gt.0.d0) )then

         thermc=thermc

      else

         write (*,*) 'Corrected therm! unrealisti! in SOILADAPT ',thermc
         write (*,*) iopgveg,thermc,iopthermc_v,tcbeta,&
                     xlai,thermc1,heatcap,heatcap1,zero
         stop

      endif

      return

    end subroutine soiladapt

! ====================================================================
!
!                       subroutine states
!
! ====================================================================
!
! Subroutine to update water table depths and root and transmission
! zone soil moisture to the conditions at the end of the previous
! time step.
!
! ====================================================================

    subroutine states(zw0,inc_frozen,i_moss,tkmid_moss,&
       r_mossm_u,r_mossm_f,r_mossm,zw,zbar,ff,atanb,xlamda,psic,&
       zrz,ztz,rzsm1,tzsm1,thetas,zrzmax,iopsmini,thetar,bcbeta,rzsm1_u,&
       tzsm1_u,rzsm1_f,tzsm1_f,tsoilold,bulk_dens,a_ice,b_ice,&
       row,rzsmold,tzsmold,r_mossmold,rzsm,tzsm,r_mossm1,zmoss,r_moss_depth,&
       thetas_moss,rzsm_u,rzsm_f,tzsm_u,tzsm_f,r_mossm1_u,&
       r_mossm1_f,i,a_ice_moss,b_ice_moss,bulk_dens_moss)

      implicit none
      include "help/states.h"
    
! ====================================================================
! Calculate the frozen moss content.
! ====================================================================

      if (inc_frozen.eq.1) then

         if (i_moss.eq.1) then

            if (tkmid_moss.lt.(273.15d0)) then

               r_mossm_u=0.d0
               r_mossm_u=1000.d0*bulk_dens_moss*a_ice_moss*(273.15d0-tkmid_moss)**b_ice_moss/row

               if (r_mossm_u.ge.r_mossm) then

                  r_mossm_u=r_mossm

               endif

               r_mossm_f=r_mossm-r_mossm_u

            endif

            if (tkmid_moss.ge.(273.15d0)) then

               r_mossm_f=0.d0
               r_mossm_u=r_mossm

            endif

         endif

      endif

! ====================================================================
! Update local water table depth.
! ====================================================================

      zw0 = zw - psic
      zw = zbar-(one/ff)*(atanb-xlamda)
      
!cw! minimum size for root zone and transmission zone is 1 cm

      if (zw-psic.lt.(zrzmax+0.01).and.zw-psic.gt.zrzmax) then

         zw=zrzmax + psic - .00001

      endif
         
      if (zw-psic.lt.(0.01).and.zw-psic.gt.zero) then

         zw=psic

      endif
!cw!
! ====================================================================
! First time step call initsm to get spatially variable initial
! conditions based on brooks-corey or user specified values.
! ====================================================================

      if (i.eq.1) then

         call initsm(zw,psic,zrz,ztz,rzsm1,tzsm1,thetas,&
       zrzmax,iopsmini,thetar,bcbeta,rzsm1_u,tzsm1_u,rzsm1_f,tzsm1_f,&
       inc_frozen,tsoilold,bulk_dens,a_ice,b_ice,row)

         rzsm_u=rzsm1_u
         tzsm_u=tzsm1_u
         rzsm_f=rzsm1_f
         tzsm_f=tzsm1_f

      endif

      if (inc_frozen.eq.0) then

! ====================================================================
! Update the values for the old and new timestep for in the option
! that the frozen soil water is treated as liquid water.
! ====================================================================

         rzsmold = rzsm
         tzsmold = tzsm
         r_mossmold = r_mossm

! --------------------------------------------------------------------&
! Update root and transmission zone depths and soil moisture.
! Also update upper and lower soil layer soil moistures.
! --------------------------------------------------------------------&

         if ((zw-psic).le.zero) then

! ....................................................................
! For saturated areas.
! ....................................................................

            zrz = zero
            rzsm = thetas
            ztz = zero
            tzsm = thetas


         else if ((zw-psic).lt.zrzmax) then

! ....................................................................
! For not saturated areas where the ground water level reaches
! the root zone.
! ....................................................................

            zrz = zw-psic
            rzsm = rzsm1
            ztz = zero
            tzsm = thetas

         else

! ....................................................................
! For not saturated areas where the ground water level is in the
! transmission zone.
! ....................................................................

            zrz = zrzmax
            ztz = zw-psic-zrz
            rzsm = rzsm1
            tzsm = tzsm1

         endif

     
! ....................................................................
! Update the values for moss moisture content of the old and new
! timestep.
! ....................................................................

         r_mossm1=r_mossm

         if (thetas_moss.gt.0.d0) then

            zmoss=r_moss_depth*r_mossm/thetas_moss

         else

            zmoss=0.d0

         endif

      else

! ====================================================================
! Update the values for the old and new timestep for in the option
! that the frozen soil water is treated as solid soil particles.
! ====================================================================

         rzsmold = rzsm
         tzsmold = tzsm
         r_mossmold = r_mossm

! --------------------------------------------------------------------&
! Update root and transmission zone depths and soil moisture.
! Also update upper and lower soil layer soil moistures.
! --------------------------------------------------------------------&

         if ((zw-psic).le.zero) then

! ....................................................................
! For saturated areas.
! ....................................................................

            zrz = zero
            rzsm_u = thetas-rzsm_f
            ztz = zero
            tzsm_u = thetas-tzsm_f
            rzsm=rzsm_u+rzsm_f
            tzsm=tzsm_u+tzsm_f

         else if ((zw-psic).lt.zrzmax) then

! ....................................................................
! For not saturated areas where the ground water level reaches
! the root zone.
! ....................................................................

            zrz = zw-psic
            rzsm_u = rzsm1_u
            rzsm=rzsm_u+rzsm_f
            ztz = zero
            tzsm_u = thetas-tzsm_f
            tzsm=tzsm_u+tzsm_f

         else

! ....................................................................
! For not saturated areas where the ground water level is in the
! transmission zone.
! ....................................................................

            zrz = zrzmax
            ztz = zw-psic-zrz
            rzsm_u = rzsm1_u
            tzsm_u = tzsm1_u
            rzsm=rzsm_u+rzsm_f
            tzsm=tzsm_u+tzsm_f

         endif

      
  
! ....................................................................
! Update the values for moss moisture content of the old and new
! timestep.
! ....................................................................

         r_mossm1_u=r_mossm_u
         r_mossm=r_mossm_u+r_mossm_f
!CVAL         zmoss=r_moss_depth*r_mossm/(thetas_moss-r_mossm_f)

         if (thetas_moss.gt.0.d0) then

            zmoss=r_moss_depth*r_mossm/(thetas_moss)

         else

            zmoss=0.d0

         endif

         r_mossm1_f=r_mossm_f

      endif

!CVAL      write (129,*)
!real(r_mossm_f/r_mossm),real(rzsm_f/rzsm),real(tzsm_f/tzsm)
!CVAL      write (130,*) real(zw-psic)

      return

    end subroutine states

! ====================================================================
!
!                       subroutine infilt
!
! ====================================================================
!
! Subroutine infilt calculates the actual rate of infiltration 
! and the rate of infiltration excess runoff.
!
! ====================================================================   

    subroutine infilt(pnet,i_moss,i_und,PackWater_us,SurfWater_us,Swq_us,&
       Outflow_us,dt,PackWater,SurfWater,Swq,Outflow,&
       istmst,cuminf,inc_frozen,rzsmst,rzsm,rzsm_u,thetas,thetar,&
       tolinf,sorp,xk0,psic,bcgamm,bcbeta,deltrz,cc,&
       zw,xinact,satxr,xinfxr,intstm,xinfcp,runtot,irntyp)

      implicit none
      include "help/infilt.h"

! ====================================================================
! Calculate the precipitation input to the soil column.
! ====================================================================

      precipitation=pnet

      if ( (i_moss+i_und) .gt.0) then

         if ((PackWater_us+SurfWater_us+Swq_us).gt.(0.d0))then

! --------------------------------------------------------------------
! In case of snow on top of the under story : the precipitation is
! the liquid water outflow from the snow pack.
! --------------------------------------------------------------------

            precipitation=Outflow_us/dt

         endif

      else

         if ((PackWater+SurfWater+Swq).gt.(0.d0))then

! --------------------------------------------------------------------
! In case of snow on top of the over story : the precipitation is
! the liquid water outflow from the snow pack.
! --------------------------------------------------------------------

            precipitation=Outflow/dt

         endif

      endif

      if(istmst.eq.1)then

! ====================================================================
! If this is the first step of storm event: reset cumulative    
! infiltration, initial root zone soil moisture, sorptivity
! and dimensionless gravity parameter, cc.
! ====================================================================

         if (inc_frozen.eq.0) then

! ....................................................................
! Option 1 : Treat frozen water as liquid water.
! ....................................................................

            call reset_inf_pars(cuminf,zero,rzsmst,rzsm,thetas,tolinf,&
       sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)

         else

! ....................................................................
! Option 2 : Treat frozen water as soil particles.
! ....................................................................

            call reset_inf_pars(cuminf,zero,rzsmst,rzsm_u,thetas,tolinf,&
       sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)

         endif

      endif

! ====================================================================
! If this is the first step of the storm event all the infiltration
! parameters are calculated now, if this is not the first step
! the parameters are calculated from the previous timestep.
! ====================================================================

      if ((zw-psic).le.zero) then

! ====================================================================
! If surface is saturated then set infiltration to zero.  Also    
! calclulate saturation excess runoff and set infiltration
! excess to zero.
! ====================================================================

         xinact = zero
         satxr = precipitation
         xinfxr = zero

      else

! ====================================================================
! If surface is unsaturated then calculate infiltration.  Set
! saturation excess runoff to zero and calculate infiltration
! excess runoff.
! ====================================================================

         if (intstm.eq.1) then

! --------------------------------------------------------------------
! In interstorm period let bare soil infiltration rate
! equal zero, and let vegetated soil infiltration rate
! equal pnetms to allow for surface wetting due to throughfall
! of condensation.
! --------------------------------------------------------------------

            xinact = precipitation

         else if(cuminf.le.zero)then

! --------------------------------------------------------------------
! If first time step with infiltration then all precipitation is
! infiltrated to avoid division by zero.
! --------------------------------------------------------------------

            xinact = precipitation

         else

! --------------------------------------------------------------------
! Calculate infitration capacity as a function of cumulative      
! infiltration for all other time steps of storm.
! --------------------------------------------------------------------

            xinfcp = cc*xk0*(one+(one/(((one+((four*cc*xk0*cuminf)/&
                                     (sorp**two)))**0.5d0)-one)))

! --------------------------------------------------------------------
! Take the actual infiltration rate as the minimum of the 
! precipitation rate or the infiltration capacity.
! --------------------------------------------------------------------

            if (precipitation.gt.xinfcp) then

               xinact = xinfcp

            else

               xinact = precipitation

            endif

         endif

! --------------------------------------------------------------------
! Calculate infiltration excess runoff; set saturation 
! excess runoff to zero.
! --------------------------------------------------------------------

         xinfxr = precipitation - xinact
         satxr = zero

      endif

! ====================================================================
! Set the value of flag used in output image for type of
! runoff and find total runoff.
! ====================================================================

      runtot = precipitation - xinact
      irntyp = 0

      if (xinfxr.gt.zero) irntyp=1
      if (satxr.gt.zero) irntyp=2

      return

    end subroutine infilt

! ====================================================================
!
!                       subroutine ebsres
!
! ====================================================================
!
! Subroutine ebares calculates the actual rate of evaporation 
! from bare soils using a soil resistance parameterization.
!
! ====================================================================

    subroutine ebsres(inc_frozen,irestype,rsoil,rzsm,srespar1,tkact,&
       srespar2,rzsm_u,srespar3,ravd,iffroz,thetas,tkmid,&
       zmid,zrzmax,smtmp,tzsm,smold,rzsmold,tzsmold,&
       iopthermc,thermc1,thermc2,thetar,heatcapold,psic,bcbeta,&
       quartz,heatcap1,ifcoarse,heatcap2,rocpsoil,row,cph2o,roa,cp,&
       roi,thermc,heatcap,rzdthetaudtemp,dshact,albd,emiss,rahd,ebscap,&
       tcel,vppa,psychr,xlhv,zdeep,Tdeepstep,rsd,rld,toleb,maxnri,dt,i,tkel,&
       zww,za,uzw,zpd,z0m,press,rib,rnetpn,gbspen,epetd,evtact,ievcon,&
       bsdew,z0h,ioppet)

      implicit none
      include "help/ebsres.h"

      data TOLSTAB/0.1/
      data MAXITER/10/

! ====================================================================
! Calculate the bare soil resistance to evaporation.
! ====================================================================

      if (inc_frozen.eq.0) then

! --------------------------------------------------------------------
! No distinction between ice and unfrozen water is made.
! --------------------------------------------------------------------

         call calcrsoil(irestype,rsoil,srespar1,&
                        srespar2,srespar3,thetas,rzsm,tkact)

      else

! --------------------------------------------------------------------
! Treat the ice as mineral soil in the calculation of soil
! resistance.
! --------------------------------------------------------------------

         call calcrsoil(irestype,rsoil,srespar1,&
                        srespar2,srespar3,thetas,rzsm_u,tkact)

      endif

! --------------------------------------------------------------------
! Total resistance = soil resistance plus aerodynami! resistance.
! --------------------------------------------------------------------

      raveff = rsoil + ravd

! ====================================================================
! Initialize soil moisture for the calculation of the thermodynamic
! parameters, as a centered difference.
! ====================================================================

      call sm_cen_dif(iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm,smold,&
                      rzsmold,tzsmold)

! ====================================================================
! Calculate the soil thermal parameters.
! ====================================================================

      call soiltherm(iopthermc,thermc1,thermc2,rzsm,smtmp,&
       thetar,thetas,psic,bcbeta,tkmid,quartz,ifcoarse,&
       heatcap1,heatcap2,heatcapold,rocpsoil,row,cph2o,roa,cp,roi,&
       smold,thermc,heatcap,inc_frozen,rzdthetaudtemp)

! ====================================================================
! Initialize temperatures for the solution of the energy balance
! equations.
! ====================================================================

      dumtk = tkact
      dumtkmid = tkmid
      dumds = dshact

! ====================================================================
! If energy balance method is used then call subroutine to 
! do newton-raphson iterations to find the energy balance.
! ====================================================================

      if (ioppet.eq.0) then

         ttemp = dumtk
         tmidtemp = dumtkmid

! --------------------------------------------------------------------
! Solve the energy balance at potential evaporation rate.
! --------------------------------------------------------------------

         call nreb(1,albd,emiss,thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       raveff,rahd,dumtk,dumtkmid,dumrn,dumxle,ebscap,dumh,dumg,dumds,tcel,&
       vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,rld,toleb,maxnri,dt,i)

! --------------------------------------------------------------------
! Check for large change in skin temperature which affects
! magnitude of stability correction.
! --------------------------------------------------------------------

         iter = 0

200      if ((abs(ttemp-dumtk)).gt.TOLSTAB.and.iter.lt.MAXITER) then

            iter = iter + 1

! --------------------------------------------------------------------
! As long as there is no convergence in skin temperature of the
! bare soil, keep on iterating for the skin temperature.
! --------------------------------------------------------------------

            ttemp = dumtk
            tacttemp = tkact
            tkact = dumtk

! --------------------------------------------------------------------
! Recalculate the stability correction.
! --------------------------------------------------------------------

            call stabcor(zww,za,uzw,zpd,z0m,tkel,press,tkact,vppa,rib)
            tkact = tacttemp
            dumtkmid = tmidtemp

! --------------------------------------------------------------------
! Recalculate the aerodynami! resistances.
! --------------------------------------------------------------------

            rahd = calcra(uzw,zww,za,zpd,z0m,z0h,rib)

            ravd = rahd
            raveff = rsoil + ravd

! --------------------------------------------------------------------
! Resolve the energy balance at potential evaporation rate.
! --------------------------------------------------------------------

            call nreb(1,albd,emiss,thermc1,thermc2,heatcap1,heatcap2,&
       heatcapold,raveff,rahd,dumtk,dumtkmid,dumrn,dumxle,ebscap,dumh,dumg,&
       dumds,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,rld,toleb,&
       maxnri,dt,i)

! --------------------------------------------------------------------
! Check again whether convergence has been achieved.
! --------------------------------------------------------------------

            goto 200

         endif

! ====================================================================
! If unreasonable skin temperatures print the variables out and
! abort the program.
! ====================================================================

         if (dumtk.lt.100.d0.or.dumtk.gt.400.d0) then

            print*,"ebsres:ERROR:skin temp outside reas. bounds Tskin="&
                ,dumtk
            print*,'thermc1 =             ',thermc1
            print*,'thermc2 =             ',thermc2
            print*,'heatcap1 =             ',heatcap1
            print*,'heatcap2 =             ',heatcap2
            print*,'dt =                  ',dt
            print*,'tmidk =             ',dumtkmid
            print*,'rav =             ',raveff
            print*,'roa =             ',roa
            stop

         endif

      else

! ====================================================================
! If penman-montieth is used then just recalculate evaporation
! using rsoil like canopy resistance.
! ====================================================================

         vpsat = 611.d0*dexp((17.27d0*tcel)/(237.3d0+tcel))
         vpdef = vpsat-vppa
         dvpsdt = 4098.d0*vpsat/((237.3d0+tcel)**two)
         pstar = psychr*(one+rsoil/ravd)
         ebsmf = ((dvpsdt*(rnetpn-gbspen))+&
                  ((cp*roa*vpdef)/ravd)) /&
                 ((dvpsdt+pstar)*xlhv)
         ebscap = ebsmf/row

      endif

! ====================================================================
! Make sure evaporative capacity is not negative.
! ====================================================================

      if (ebscap.lt.zero) ebscap=zero

! ====================================================================
! Set actual evaporation to minimum of potential evaportation from
! the bare soil and the soil controlled evaporation rate.
! Also add the land area to the soil controlled or atmospheric
! controlled percentage.
! ====================================================================

      if (ebscap.le.epetd) then

         evtact = ebscap
         ievcon = 1

      else

         evtact = epetd
         ievcon = 2

      endif

! ====================================================================
! If actual evaporation is negative then condensation occuring.
! ====================================================================

      if (evtact.lt.zero) then

         bsdew = -evtact
         evtact = zero

      else

         bsdew = zero

      endif

      return

    end subroutine ebsres

! ====================================================================
!
!                subroutine calcrain
!
! ====================================================================
!
! Calculate the partitioning of precipitation into rain and snow.
! Do this based on air temperature.
!
! ====================================================================

    subroutine calcrain (tcel,snow,rain,precip_o,dt)

      implicit none
      include "help/calcrain.h"

      rain=0.d0
      snow=0.d0

      if (tcel.gt.(0.d0)) then

         snow=0.d0*dt
         rain=precip_o*dt

      endif

      if (tcel.le.(0.d0)) then

         rain=0.d0*dt
         snow=precip_o*dt

      endif

      return

    end subroutine calcrain

! ====================================================================
!
!                       subroutine transv
!
! ====================================================================
!
! Subroutine to find actual transpiration from dry canopy.
!
! ====================================================================

    subroutine transv(epetd,epetd_us,i_und,iopveg,f1par,f3vpd,f4temp,&
       rescan,inc_frozen,ivgtyp,rzsm,rzsm_u,tc,tw,smcond,tzsm,tzsm_u,&
       tc_us,tw_us,smcond_us,f1par_us,f3vpd_us,f4temp_us,rescan_us,&
       vegcap,ravd,vegcap_us,ravd_us,zrz,srzrel,thetas,thetar,psisoi,&
       psic,bcbeta,ikopt,xksrz,xk0,ff,ressoi,rtact,rtdens,psicri,&
       respla,xkrz,ztz,stzrel,xkstz,xktz,Swq,evtact,ievcon,Swq_us,evtact_us,&
       ievcon_us,bsdew,i,ipix)

      implicit none
      include "help/transv.h"
      integer :: i,ipix

! ====================================================================
! Set potential transpiration to zero if negative for both under s.or.&
! and over story.
! ====================================================================

      if (epetd.le.zero) epetd=zero
      if ( (i_und.gt.0) .and. (epetd_us.le.zero) ) epetd_us=zero

! ====================================================================
! Calculate the maximum water vapor flux out of the plants possible
! by using a soil moisture conductance parameterization.
! ====================================================================

      if (iopveg.eq.0) then

! --------------------------------------------------------------------&
! Set the multiplier for the linear interpolation between 
! critical and wilting soil moisture condition based on 
! soil moisture in either the upper or lower layer.
! --------------------------------------------------------------------&

         resist = f1par*f3vpd*f4temp*rescan

         if (inc_frozen.eq.0) then

! --------------------------------------------------------------------&
! Assuming no difference between ice particles and unfrozen water.
! --------------------------------------------------------------------&

            if (ivgtyp.eq.1) then

! ....................................................................
! Vegetation roots are in the upper layer.  Linear increase with
! soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(rzsm,tc,smcond,one,tw,zero)

            else

! ....................................................................
! Vegetation roots are in the lower layer.  Linear increase with
! soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(tzsm,tc,smcond,one,tw,zero)

            endif

         else

! --------------------------------------------------------------------&
! Treating the frozen soil water as soil particles.
! --------------------------------------------------------------------&

            if (ivgtyp.eq.1) then

! ....................................................................
! Vegetation roots are in the upper layer.  Linear increase with
! liquid soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(rzsm_u,tc,smcond,one,tw,zero)

            else

! ....................................................................
! Vegetation roots are in the lower layer.  Linear increase with
! liquid soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(tzsm_u,tc,smcond,one,tw,zero)

            endif

         endif

! ====================================================================
! Calculate the soil moisture effect on stomatal resistance for
! the under story.  The roots of the under story are always assumed
! to be in the upper soil layer.
! ====================================================================

         if (inc_frozen.eq.0) then

! --------------------------------------------------------------------&
! Assuming no difference between ice particles and unfrozen water.
! --------------------------------------------------------------------&

            if (i_und.gt.0) then

! ....................................................................
! Linear increase with liquid soil moisture between wilting and critical
! soil moisture, 1 if the soil moisture is higher than the wilting soil
! moisture and 0 if the soil moisture is lower than the critical soil
! moisture.
! ....................................................................

               call calcsmcond(rzsm,tc_us,smcond_us,one,tw_us,zero)

               resist_us=f1par_us*f3vpd_us*f4temp_us*rescan_us

           endif

        else

! --------------------------------------------------------------------&
! Treating the frozen soil water as soil particles.
! --------------------------------------------------------------------&

           if (i_und.gt.0) then

! ....................................................................
! Linear increase with liquid soil moisture between wilting and critical
! soil moisture, 1 if the soil moisture is higher than the wilting soil
! moisture and 0 if the soil moisture is lower than the critical soil
! moisture.
! ....................................................................

               call calcsmcond(rzsm_u,tc_us,smcond_us,one,tw_us,zero)

               resist_us=f1par_us*f3vpd_us*f4temp_us*rescan_us

            endif

         endif

! ====================================================================
! Calculate vegetation capacity for the over story using linear
! interpolation for the canopy resistance.
! ====================================================================

         call calcvegcap(smcond,zero,vegcap,epetd,resist,ravd,smcond)

! ====================================================================
! Calculate vegetation capacity for the under story using linear
! interpolation for the canopy resistance.
! ====================================================================

         if (i_und.gt.0) then

            call calcvegcap(smcond_us,zero,vegcap_us,&
                            epetd_us,resist_us,ravd_us,smcond_us)

         endif

! ====================================================================
! Calculate the maximum water vapor flux out of the plants possible
! by using a root resistivity parameterization.
! ====================================================================

      else

! ====================================================================
! Calculate the maximum plant evaporation for the under story.
! ====================================================================

         if (i_und.gt.0) then

            call maxplevap(zrz,0.d0,epetd,inc_frozen,srzrel,rzsm,thetas,&
       thetar,rzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xksrz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap_us,psicri,respla,xkrz)

         endif

! ====================================================================
! Calculate the maximum plant evaporation for the over story for
! vegetation with its roots in the upper soil layer.
! ====================================================================

         if (ivgtyp.eq.1) then

            call maxplevap(zrz,0.d0,epetd,inc_frozen,srzrel,rzsm,thetas,&
       thetar,rzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xksrz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla,xkrz)

! ====================================================================
! Calculate the maximum plant evaporation for the over story for
! vegetation with its roots in the upper soil layer.
! ====================================================================

         else

            call maxplevap(ztz,zrz,epetd,inc_frozen,stzrel,tzsm,thetas,&
       thetar,tzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xkstz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla,xktz)

         endif

      endif

! ====================================================================
! Set actual transpiration to the minimum of potential 
! transpiration or vegetation capacity.
! ====================================================================

      call acttrans(Swq,vegcap,epetd,evtact,ievcon,zrz)

      if (i_und.gt.0) then

        call acttrans(Swq_us,vegcap_us,epetd_us,evtact_us,ievcon_us,zrz)

      endif 
! ====================================================================
! Since condensation occurs in canopy, set bare soil condensation
! to zero.
! ====================================================================

      bsdew = zero

      return

    end subroutine transv

! ====================================================================
!
!                       subroutine tz_and_rzbal
!
! ====================================================================
!
! Calculate the root zone and transmission zone water balances.
!
! ====================================================================
!
! Created by Wade Crow 3/1/1999
! 
! ====================================================================
     subroutine tz_and_rzbal(i,newstorm,inc_frozen,ikopt,ivgtyp,dt,&

! Surface zone and transmission zone soil moistures 

       rzsm,tzsm,rzsm_u,tzsm_u,rzsm1,tzsm1,&

! Layer geometry

       zrz,ztz,zw0,zrzmax,&

! Water fluxes

       evtact,evtact_us,bsdew,dewrun,grz,gtz,diftz,difrz,&
       satxr,runtot,xinact,cuminf,&

! Soil parameters

       ff,thetar,thetas,bcbeta,xk0,psic,&

! Snow
       Swq,Swq_us,&

! Understory/moss

       dc,i_und,i_moss,fw,dc_us,fw_us,evrz_moss,f_und,&

! Storage changes

       dstz,dsrz,&

! Summations of fluxes

       tzrhs,rzrhs)

      implicit none
      include "help/tz_and_rzbal.h"
      integer :: test_flag
      
 
      rzsm0 = rzsm
      tzsm0 = tzsm
      rzsm1old = rzsm
      rzsm1new = rzsm
      tzsm1new = tzsm
      tzsm1old = tzsm
      evtran_rz = zero
      evtran_tz = zero
      cor_flx_tz = zero
      cor_flx_rz = zero
      ddifrzdth1 = zero
      ddifrzdth2 = zero
      dgrzdth1 = zero
      dgrzdth2 = zero
      dgtzdth1 = zero
      dgtzdth2 = zero
      ddiftzdth1 = zero
      ddiftzdth2 = zero
      dewrun = zero
      iter = 1
      max_iter = 10
      num_exp_iter = 60

! ====================================================================
! Calculate change of hydrauli! conductivity with depth
! ====================================================================

      if(ikopt.eq.1)then

        xksrz=xk0
        xkstz=xk0

      else

        xksrz=xk0*dexp(-ff*(zrz/two))
        xkstz=xk0*dexp(-ff*(zrz+(ztz/two)))

      end if

! ====================================================================
! Calculate the second derivative of soil moisture with 
! respect to time.
! ====================================================================

      call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u) 

      call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

      call clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,ikopt,&
       xksrz,xkstz,ff,zrz,ztz,bcbeta,thetas,thetar,psic,tzsm)

      dgrzdth1= clcdg(rzsm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)
      dgtzdth2= clcdg(tzsm,ikopt,xkstz,ff,ztz,bcbeta,thetas,thetar)

      if (zrz.gt.zero) then
        
        d2rzsmdt2 = (-dgrzdth1 - ddifrzdth1)*(-grz - difrz + xinact - evtran_rz)*(dt*dt)/(zrz*zrz)

      else

        d2rzsmdt2 = zero

      endif

      if (ztz.gt.zero) then

        d2tzsmdt2 = (ddifrzdth2 -dgtzdth2- ddiftzdth2)*(grz + difrz - gtz - diftz - evtran_tz)*(dt*dt)/(ztz*ztz)

      else

        d2tzsmdt2 = zero

      endif

! ====================================================================
! Modify the numerical approach based on linearity of problem.
! ====================================================================
  
      if (abs(d2rzsmdt2).lt..15d0.and.abs(d2tzsmdt2).lt..15d0) then
 
        explicit_weight = .5d0
        implicit_weight = .5d0
        non_linear_flag = 0
        tol = .001d0

      else if (abs(d2rzsmdt2).gt.5.or.abs(d2tzsmdt2).gt.5) then
        explicit_weight = zero
        implicit_weight = one
        tol = .0001d0
        non_linear_flag = 1

      else

        explicit_weight = zero
        implicit_weight = one
        tol = .0001d0
        non_linear_flag = 0
      endif

! ====================================================================
! Turn on and off numerical testing procedure
        test_flag = 0
! ====================================================================

! ====================================================================
! Identify created and destroyed surface and transmission zones
! - define correction fluxes to account for created water
! don't change the order of these.
! ====================================================================

! ====================================================================
! New surface zone created.
! ====================================================================


      if (zw0.le.zero.and.zrz.gt.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_rz = -zrz*(thetas-thetar) 

      endif

! ====================================================================
! New transmission zone created
! ====================================================================

      if (zw0.le.zrzmax.and.ztz.gt.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_tz = -ztz*(thetas-thetar)

      endif

! ====================================================================
! Transmission zone destroyed
! ====================================================================

      if (zw0.gt.zrzmax.and.ztz.le.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_rz = (zw0-zrzmax)*(thetas-thetar)

      endif

! ====================================================================
! Surface zone or surface and transmission zone destroyed.
! ====================================================================

      if (zw0.gt.zero.and.zrz.le.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_rz = zw0*(thetas-thetar)

      endif


!================================================
! shut off correction fluxes
!=================================================

      cor_flx_rz = 0
      cor_flx_tz = 0

      
! ====================================================================
! Calculate evaporation and transpiration from surface and transmission
! zone.
! ====================================================================

      call clc_evrz(evtran_rz,Swq,Swq_us,ivgtyp,evtact,dc,i_und,&
       i_moss,fw,evtact_us,dc_us,fw_us,evrz_moss,dummy,f_und)

      if (ivgtyp.eq.2) then

        evtran_tz=evtact*dc*(1-fw)

      else

        evtran_tz  = 0

      endif

! ====================================================================
! Decide if dew goes into surface surface.
! ====================================================================

      if (ivgtyp.eq.0) then 

        dewrz = bsdew

      else

        dewrz = zero

      endif

! ====================================================================
! Case I - water table is at surface
! ====================================================================
 
      if (zrz.eq.zero) then

        case_flag = 1

        rzsm1 = thetas
        tzsm1 = thetas
        difrz = zero
        diftz = zero
        grz = zero
        gtz = zero
        evtran_rz = zero
        evtran_tz = zero
        dewrun = dewrz
        dewrz = zero 
        satxr = satxr + dewrun
        runtot = runtot + dewrun 
        
      endif


! ====================================================================
! Case II - water table is in transmission zone
! ====================================================================

      if (ztz.gt.zero.and.non_linear_flag.eq.0) then

        case_flag = 2.1

500     dgrzdth1= clcdg(rzsm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)
        dgtzdth2= clcdg(tzsm,ikopt,xkstz,ff,ztz,bcbeta,thetas,thetar)
        dgrzdth2 = zero
        dgtzdth1 = zero 

        call clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,ikopt,xksrz,xkstz,ff,&
       zrz,ztz,bcbeta,thetas,thetar,psic,tzsm)

        ddiftzdth1 = zero 

        call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)
          
        call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

        srzflx = dt * (dewrz + xinact - difrz - grz - evtran_rz)
        stzflx = dt * (difrz - diftz  + grz - gtz - evtran_tz)
     
        F11 = rzsm0-rzsm1new+srzflx/zrz
        F22 = tzsm0-tzsm1new+stzflx/ztz

        dF1dtheta1 = dt/zrz*(-ddifrzdth1 - dgrzdth1) - one
        dF1dtheta2 = dt/zrz*(-ddifrzdth2 - dgrzdth2) 
        dF2dtheta1 = dt/ztz*(ddifrzdth1 + dgrzdth1 - dgtzdth1 - ddiftzdth1)
        dF2dtheta2 = dt/ztz*(ddifrzdth2 + dgrzdth2 - dgtzdth2 - ddiftzdth2) - one
        
! --------------------------------------------------------------------&
! Solve 2 x 2 matrix */&
!
! dF1dtheta1(rzsm)   dF1dtheta2(tzsm)           del_rzsm
! F1(rzsm)
! dF2dtheta1(rzsm)   dF2dtheta2(tzsm)           del_tzsm
! F2(tzsm)
! --------------------------------------------------------------------&

        del_rzsm = one/(dF2dtheta1-(dF1dtheta1*dF2dtheta2/dF1dtheta2))*(F22 -(dF2dtheta2*F11/dF1dtheta2))
        del_tzsm = (F11 - dF1dtheta1*del_rzsm)/dF1dtheta2
        rzsm1new = -(del_rzsm) + rzsm1old
        tzsm1new = -(del_tzsm) + tzsm1old

        if (rzsm1new.gt.thetas) then

          rzsm1new=thetas

        endif

        if (rzsm1new.lt.thetar) then

          rzsm1new=thetar
!         write(*,*) rzsm1new

        endif

        if (tzsm1new.gt.thetas) then

          tzsm1new=thetas

        endif

        if (tzsm1new.lt.thetar) then

          tzsm1new=thetar

        endif

        if ((((rzsm1new.lt.rzsm1old-tol).or.(rzsm1new.gt.rzsm1old+tol)).and.(iter.le.max_iter)).or.&
       (((tzsm1new.lt.tzsm1old-tol).or.(tzsm1new.gt.tzsm1old+tol)).and.(iter.le.max_iter))) then

          rzsm=implicit_weight*rzsm1new + explicit_weight*rzsm0
          tzsm=implicit_weight*tzsm1new + explicit_weight*tzsm0
          iter=iter+1
          rzsm1old = rzsm1new
          tzsm1old = tzsm1new
          goto 500

        endif

        rzsm1 = rzsm0 + dt/zrz*(dewrz -difrz - grz + xinact - evtran_rz)
        tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

      endif

! ====================================================================
! Case III - water table is in surface zone.
! ====================================================================

      if (zrz.gt.zero.and.zrz.lt.zrzmax.and.non_linear_flag.eq.0) then

        case_flag = 3.1

        diftz = zero
        gtz = zero
        evtran_tz  = zero
        tzsm = thetas

600     dgrzdth1= clcdg(rzsm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)

        call clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,ikopt,xksrz,xkstz,ff,zrz,ztz,bcbeta,thetas,&
       thetar,psic,tzsm)

        call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

        call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

        srzflx = dt * (dewrz + xinact - difrz - grz - evtran_rz)
        F11 = rzsm0-rzsm1new+srzflx/zrz    
        dF1dtheta1 = dt/zrz*(-ddifrzdth1 - dgrzdth1) - one
        rzsm1new = rzsm1old - F11/dF1dtheta1       

        if (rzsm1new.gt.thetas) then

          rzsm1new=thetas

        endif

        if (rzsm1new.lt.thetar) then

          rzsm1new=thetar

        endif

        if ((((rzsm1new.lt.rzsm1old-tol).or.(rzsm1new.gt.rzsm1old+tol)).and.(iter.le.max_iter))) then

          rzsm=implicit_weight*rzsm1new + explicit_weight*rzsm0
          iter=iter+1
          rzsm1old = rzsm1new
          goto 600

        endif

        rzsm1 = rzsm0 + dt/zrz*(dewrz -difrz - grz + xinact - evtran_rz)
        tzsm1 = thetas

      endif

! ====================================================================
! Back-up Numerics.
! ====================================================================
    
      if (tzsm1.gt.thetas.or.tzsm1.lt.thetar.or.rzsm1.lt.thetar.or.rzsm1.gt.thetas.or.&
       iter.gt.max_iter.or.non_linear_flag.eq.1.or.test_flag.eq.1) then
        
        grz_sum = zero
        difrz_sum = zero
        srzflx = zero
        gtz_sum = zero
        diftz_sum = zero
        stzflx = zero
        tzsm = tzsm0
        rzsm = rzsm0

! --------------------------------------------------------------------&
! Case II - water table is in tranmission zone
! --------------------------------------------------------------------&

        if (ztz.gt.zero) then

          case_flag = 2.2

          do 700 num=1,num_exp_iter

            call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

            call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

            dtsrzflx = dt * (dewrz + xinact - difrz  - grz - evtran_rz)
            srzflx = srzflx + dtsrzflx*(num_exp_iter**(-one))
            difrz_sum = difrz_sum + difrz*(num_exp_iter**(-one))
            grz_sum = grz_sum + grz*(num_exp_iter**(-one))
            rzsm = rzsm0 + srzflx/zrz
            dtstzflx = dt * (difrz  + grz - diftz - gtz - evtran_tz)
            stzflx = stzflx + dtstzflx*(num_exp_iter**(-one))
            diftz_sum = diftz_sum + diftz*(num_exp_iter**(-one))
            gtz_sum = gtz_sum + gtz*(num_exp_iter**(-one))            
            tzsm = tzsm0 + stzflx/ztz

            if (rzsm.gt.thetas) then

              rzsm=thetas

            endif

            if (rzsm.lt.thetar) then

              rzsm=thetar

            endif

            if (tzsm.gt.thetas) then

              tzsm=thetas

            endif
            if (tzsm.lt.thetar) then

              tzsm=thetar

            endif

 700      continue
          difrz =  difrz_sum
          grz =  grz_sum
          diftz = diftz_sum
          gtz = gtz_sum
          rzsm1 = rzsm0 + dt/zrz*(dewrz - difrz - grz + xinact - evtran_rz)
          tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

        endif
     
! --------------------------------------------------------------------&
! Case III - water table is in surface zone.
! --------------------------------------------------------------------&

        if (zrz.gt.zero.and.zrz.lt.zrzmax) then

          case_flag = 3.2

          diftz = zero
          gtz = zero
          evtran_tz  = zero 
          tzsm = thetas

         do 800 num=1,num_exp_iter

            call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

            call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

            dtsrzflx = dt * (dewrz + xinact - difrz  - grz - evtran_rz)
            srzflx = srzflx + dtsrzflx*(num_exp_iter**(-one))
            difrz_sum = difrz_sum + difrz*(num_exp_iter**(-one))
            grz_sum = grz_sum + grz*(num_exp_iter**(-one))
            rzsm = rzsm0 + srzflx/zrz

            if (rzsm.gt.thetas) then

              rzsm=thetas

            endif

            if (rzsm.lt.thetar) then

              rzsm=thetar

            endif

 800      continue

          difrz =  difrz_sum
          grz =  grz_sum
          diftz = zero
          gtz = zero
          rzsm1 = rzsm0 + dt/zrz*(dewrz - difrz - grz + xinact - evtran_rz)
          tzsm1 = thetas

        endif

      endif
 
! ====================================================================
! Balance checks and corrections.
! ====================================================================

      if (rzsm1.gt.thetas) then
      
        grz = grz + (rzsm1-thetas)*zrz/dt

        if (ztz.gt.zero) then

          tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

        endif

        rzsm1 = thetas

      endif

      if (tzsm1.gt.thetas) then

        gtz = gtz + (tzsm1-thetas)*ztz/dt
        tzsm1=thetas

      endif
   
      if (rzsm1.lt.thetar) then

        difrz = difrz + (rzsm1-thetar)*zrz/dt

        if (ztz.gt.zero) then

           tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

        endif

        rzsm1=thetar

      endif

      if (tzsm1.lt.thetar) then

        diftz = diftz + (tzsm1-thetar)*ztz/dt
        tzsm1=thetar

      endif

! ====================================================================
! Calculate storage changes (used to check water balance in unit 95).
! ====================================================================

      dstz = ztz*(tzsm1-tzsm0)
      dsrz = zrz*(rzsm1-rzsm0)

! ====================================================================
! Sum fluxes (used to check water balance in unit 95).
! ====================================================================

      if (ztz.gt.zero) then

        rzrhs = dt*(xinact-difrz + dewrz-evtran_rz-grz)
        tzrhs = dt*(grz-gtz-evtran_tz+ difrz - diftz )

      else
      
        rzrhs = dt*(xinact-difrz + dewrz-evtran_rz-grz)
        tzrhs = dt*(-gtz-evtran_tz-diftz )
 
      endif

! ======================================================================
! Corrections for "created" and "destroyed" layers
! ======================================================================

      difrz = difrz + cor_flx_rz/dt
      diftz = diftz + cor_flx_tz/dt

! ====================================================================
! Update cumulative infiltration.
! ====================================================================

      cuminf = cuminf + xinact * dt
 
      rzsm = rzsm0
      tzsm = tzsm0

      return
    end subroutine tz_and_rzbal

! ====================================================================
!
!                       subroutine sursat
!
! ====================================================================
!
! Define land surface saturation states for the region.
!
! ====================================================================

    subroutine sursat(fwcat,fwreg,zw,psic,fw,mul_fac,pr3sat,zrzmax,perrg2,&
       rzsm1,thetas,pr2sat,pr2uns,perrg1,tzsm1,pr1sat,pr1rzs,pr1tzs,pr1uns,&
       satxr,persxr,xinfxr,perixr,ievcon,persac,peruac,perusc,i,ipix)

      implicit none
      include "help/sursat.h"
      integer :: i,ipix

      data tolsat / 0.0001d0 /

! ====================================================================
! Calculate regional fraction of wet canopy.
! ====================================================================

      !fwcat = fwcat + fw*mul_fac
      fwcat = fw*mul_fac
      !fwreg = fwreg + fw*mul_fac
      fwreg = fw*mul_fac

! ====================================================================
! Define land surface saturation states for the region:
!
! Region 3:  Water Table at surface.
! Region 2:  Water Table in root zone.
! Region 1:  Water Table below root zone.
! ====================================================================

      if ((zw-psic).le.zero) then

! --------------------------------------------------------------------&
! First find areas in region 3.
! --------------------------------------------------------------------&

         !pr3sat = pr3sat + one*mul_fac
         pr3sat = one*mul_fac

      else if (((zw-psic).lt.zrzmax).and.((zw-psic).gt.zero)) then

! --------------------------------------------------------------------&
! For all pixels not in area 3 : first see if the water table is
! in the root zone and if the root zone is not saturated (region 2).
! --------------------------------------------------------------------&

         !perrg2 = perrg2 + one*mul_fac
         perrg2 = one*mul_fac

         if (rzsm1.ge.thetas-tolsat) then

            !pr2sat = pr2sat + one*mul_fac
            pr2sat = one*mul_fac

         else

            !pr2uns = pr2uns + one*mul_fac
            pr2uns = one*mul_fac

         endif

      else

! --------------------------------------------------------------------&
! If a pixel is not in in region 3 or 2 it has to be in region 1.
! Ssplit into four possibilities for root and transmission zone
! saturation:
! --------------------------------------------------------------------&

         !perrg1 = perrg1 + one
         perrg1 = one

         if ((rzsm1.ge.thetas-tolsat).and.(tzsm1.ge.thetas-tolsat)) then

! ....................................................................
! 1) Boot root and transmsission zone are saturated.
! ....................................................................

            !pr1sat = pr1sat + one*mul_fac
            pr1sat = one*mul_fac

         else if ((rzsm1.ge.thetas-tolsat).and.(tzsm1.lt.thetas-tolsat)) then

! ....................................................................
! 2) Root zone is saturated and transmsission zone is not
!    saturated.
! ....................................................................

            !pr1rzs = pr1rzs + one*mul_fac
            pr1rzs = one*mul_fac

         else if ((rzsm1.lt.thetas-tolsat).and.(tzsm1.ge.thetas-tolsat)) then 

! ....................................................................
! 3) Root zone is not saturated and transmsission zone is
!    saturated.
! ....................................................................

            !pr1tzs = pr1tzs + one*mul_fac
            pr1tzs = one*mul_fac

         else

! ....................................................................
! 4) Both root and transmsission zone are not saturated.
! ....................................................................

            !pr1uns = pr1uns + one*mul_fac
            pr1uns = one*mul_fac

         endif

      endif

! ====================================================================
! Determine fractions of land surface contribtuting saturation
! or infiltration excess runoff.
! ====================================================================

      if (satxr.gt.zero) then

         !persxr = persxr + one*mul_fac
         persxr = one*mul_fac

      else if (xinfxr.gt.zero) then

         !perixr = perixr + one*mul_fac
         perixr = one*mul_fac

      endif

! ====================================================================
! Determine areal fractions of bare soil evaporation
! controls - check for atmospheri! contolled (saturated),&
! atmospheri! contolled (unsaturated) and soil controlled.
! ====================================================================

      if (ievcon.eq.3) then

         !persac = persac + one*mul_fac
         persac = one*mul_fac

      else if (ievcon.eq.2) then

         !peruac = peruac + one*mul_fac
         peruac = one*mul_fac

      else 

         !perusc = perusc + one*mul_fac
         perusc = one*mul_fac

      endif

      return

    end subroutine sursat

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

      if (((Swq.gt.(0.0002d0)).or.&
           ((tcel.le.(0.d0)).and.(dt*precip_o.gt.(0.0002d0)))).AND.&
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

    end subroutine land_os

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

    end subroutine land_us




END MODULE MODULE_LAND
