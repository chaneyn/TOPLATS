! ====================================================================
!
!			subroutine land
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
      include "SNOW.h"
      include "help/land.h"
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

      end
