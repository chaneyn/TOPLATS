MODULE MODULE_LAND

  USE MODULE_VARIABLES
  
  USE MODULE_SHARED
  
  USE MODULE_FROZEN_SOIL

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
    subroutine land(ipix,i,&

! Meteorological data

       tcel,vppa,psychr,xlhv,tkel,&
       appa,vpsat,&

! Energy fluxes

       epetd,bsdew,&
       rnetd,&
       tkd,tkmidd,&

! Moss parameters

       r_moss_depth,thetas_moss,srespar1_moss,srespar2_moss,srespar3_moss,&
       eps,emiss_moss,zpd_moss,rib_moss,z0m_moss,z0h_moss,epet_moss,&
       a_ice_moss,b_ice_moss,bulk_dens_moss,&

! Vegetation parameters

       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,f1,f2,f3,&

! Constants

       roa,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,rib,RaSnow,&

! Water balance variables

       rzsm_u,tzsm_u,r_mossmold,deltrz,dc_us,fw_us,dewrun,&
       CELL_VARS,GRID_MET,GRID_VEG,GRID_VARS,GRID_SOIL,CAT,GLOBAL)

      implicit none
      include "help/land.h"!Remove when variables are changed

      type (GRID_MET_template) :: GRID_MET
      type (GRID_VEG_template) :: GRID_VEG
      type (GRID_VARS_template) :: GRID_VARS
      type (GRID_SOIL_template) :: GRID_SOIL
      type (CATCHMENT_template) :: CAT
      type (GLOBAL_template) :: GLOBAL
      type (CELL_VARS_template) :: CELL_VARS
      real*8 gold

      data tolinf/1.0d-09/
      initer=2
! ====================================================================
! Temporarily reassgin variables from old to new format
    !Meteorological data
      !tcel
      !vppa
      !psychr
      !xlhv
      !tkel
      !appa
      !vpsat
      
    !Energy fluxes and states
      !epetd
      !bsdew
      !rnetd
      !xled
      !hd
      !gd
      !dshd
      !tkd
      !tkmidd
      !rnetw
      !xlew
      !hw
      !gw
      !dshw
      !tkw
      !tkmidw
      
    !Vegetation parameters
     
      !f1par
      !f3vpd
      !f4temp
      !f1par_us
      !f3vpd_us
      !f4temp_us
      !f1
      !f2
      !f3
    
    !Constants
      !roa
      !roa_ic

    !Energy balance variables
      !ravd
      !rahd
      !ravd_us
      !rahd_us
      !rav_moss
      !rah_moss
      !rib
      !RaSnow
      
    !Water balance variables
      !rzsm_u
      !tzsm_u
      !r_mossold
      !deltrz 
      !dc_us
      !fw_us
      !dewrun
                  
    !Different option parameters

! ====================================================================
! Initialize the rain and snowfall.
! ====================================================================

      rain=0.d0
      snow=0.d0
      f1=0.d0
      f2=0.d0
      f3=0.d0

      gold=GRID_VARS%gact

! ====================================================================
! Calculate the incoming solar radiation for the under story and over
! story layers.
! ====================================================================

      call calc_rs(GRID_VEG,GRID_VEG%i_und,GRID_VEG%i_moss,GRID_VARS%Swq_us,&
                   GRID_VEG%alb_moss,GRID_VARS%alb_snow,GRID_MET%rsd,rs_over,rs_under)

! ====================================================================
! Initialize soil moisture for the calculation of the thermodynami!
! parameters, as a centered difference.
! ====================================================================

      call sm_cen_dif(iffroz,GRID_VARS%tkmid,GRID_SOIL%zmid,GLOBAL%zrzmax,smtmp,GRID_VARS%rzsm,GRID_VARS%tzsm,GRID_VARS%smold,&
                      GRID_VARS%rzsmold,GRID_VARS%tzsmold)

! ====================================================================
! Calculate the soil thermal parameters.
! ====================================================================

      call soiltherm(GLOBAL%iopthermc,thermc1,thermc2,GRID_VARS%rzsm,smtmp,&
       GRID_SOIL%thetar,GRID_SOIL%thetas,GRID_SOIL%psic,GRID_SOIL%bcbeta,GRID_VARS%tkmid,&
       GRID_SOIL%quartz,GRID_SOIL%ifcoarse,&
       heatcap1,heatcap2,heatcapold,GRID_SOIL%rocpsoil,GRID_VARS%row,GRID_VARS%cph2o,roa,GRID_VARS%cp,GRID_VARS%roi,&
       GRID_VARS%smold,thermc,heatcap,GLOBAL%inc_frozen,GRID_VARS%rzdthetaudtemp)

! ====================================================================
! Modify the thermal parameters for soils under vegetation.
! ====================================================================

      if (GRID_VEG%ivgtyp.ne.0) then
         call soiladapt(GLOBAL%iopgveg,thermc,GLOBAL%iopthermc_v,GRID_VEG%tcbeta,&
                        GRID_VEG%xlai,thermc1,heatcap,heatcap1,zero)

      endif

      heatcap_us=heatcap
      thermc_us=thermc
      !call soiladapt(iopgveg,thermc_us,iopthermc_v,tcbeta_us,&
      !               xlai_us,thermc1,heatcap,heatcap1,zero)

! ====================================================================
! Initialize actual temperatures.
! ====================================================================

      tkactd = GRID_VARS%tkact
      tkmidactd = GRID_VARS%tkmid
      tkactd_us = CELL_VARS%tkact_us
      tkmidactd_us = CELL_VARS%tkmid_us
      tskinactd_moss = CELL_VARS%tskinact_moss
      tkactd_moss = CELL_VARS%tkact_moss
      tkmidactd_moss = CELL_VARS%tkmid_moss

! ====================================================================
! Determine the frozen and unfrozen part of the soil water if
! the representation of frozen soil processes is requested.
! ====================================================================

      GRID_VARS%rzdthetaidt=0.d0
      GRID_VARS%tzdthetaidt=0.d0

      if (GLOBAL%inc_frozen.eq.1) then

         call ice_change(ipix,GRID_VARS%rzdthetaidt,GRID_VARS%tzdthetaidt,GRID_VEG%f_moss,GRID_VEG%f_und,&
       GRID_VARS%tkmidpet,CELL_VARS%tkmidpet_us,CELL_VARS%tkmidpet_moss,GRID_VARS%rzsm1_f,GRID_VARS%tzsm1_f,GRID_SOIL%bulk_dens,&
       GRID_SOIL%a_ice,GRID_SOIL%b_ice,GRID_VARS%row,GRID_VARS%roi,GRID_VARS%rzsm1_u,GRID_VARS%rzsm1,&
       GRID_VARS%tzsm1_u,GRID_VARS%tzsm1,GRID_SOIL%thetas,GRID_SOIL%thetar,&
       GRID_VARS%rzdthetaudtemp,GLOBAL%dt,GRID_VARS%rzsm_f,GRID_VARS%tzsm_f,CELL_VARS%tsoilold)

      endif

      GRID_VARS%rzsm_f=GRID_VARS%rzsm1_f
      GRID_VARS%tzsm_f=GRID_VARS%tzsm1_f

! ====================================================================
! Update local water table depth, root and transmission zone soil
! moisture.
! ====================================================================

      call states(zw0,GLOBAL%inc_frozen,GRID_VEG%i_moss,0.5d0*(CELL_VARS%tskinact_moss+CELL_VARS%tskinact_moss),&
       GRID_VARS%r_mossm_u,GRID_VARS%r_mossm_f,GRID_VARS%r_mossm,GRID_VARS%zw,CAT%zbar,&
       CAT%ff,GRID_VARS%atanb,CAT%xlamda,GRID_SOIL%psic,&
       GRID_VARS%zrz,GRID_VARS%ztz,GRID_VARS%rzsm1,GRID_VARS%tzsm1,GRID_SOIL%thetas,GLOBAL%zrzmax,GLOBAL%iopsmini,&
       GRID_SOIL%thetar,GRID_SOIL%bcbeta,GRID_VARS%rzsm1_u,&
       GRID_VARS%tzsm1_u,GRID_VARS%rzsm1_f,GRID_VARS%tzsm1_f,CELL_VARS%tsoilold,&
       GRID_SOIL%bulk_dens,GRID_SOIL%a_ice,GRID_SOIL%b_ice,&
       GRID_VARS%row,GRID_VARS%rzsmold,GRID_VARS%tzsmold,r_mossmold,GRID_VARS%rzsm,GRID_VARS%tzsm,&
       GRID_VARS%r_mossm1,GRID_VARS%zmoss,r_moss_depth,&
       thetas_moss,rzsm_u,GRID_VARS%rzsm_f,tzsm_u,GRID_VARS%tzsm_f,GRID_VARS%r_mossm1_u,&
       GRID_VARS%r_mossm1_f,i,a_ice_moss,b_ice_moss,bulk_dens_moss)

! ====================================================================
! Calculate the infiltration.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : No moss layer on the surface, infiltration capacity
!            depends on cumulative infiltration, rainfall, ...
! --------------------------------------------------------------------

         call infilt(GRID_VARS%pnet,GRID_VEG%i_moss,GRID_VEG%i_und,GRID_VARS%PackWater_us,GRID_VARS%SurfWater_us,GRID_VARS%Swq_us,&
       GRID_VARS%Outflow_us,GLOBAL%dt,GRID_VARS%PackWater,GRID_VARS%SurfWater,GRID_VARS%Swq,GRID_VARS%Outflow,&
       GRID_VARS%istmst,GRID_VARS%cuminf,GLOBAL%inc_frozen,rzsmst,GRID_VARS%rzsm,rzsm_u,GRID_SOIL%thetas,GRID_SOIL%thetar,&
       tolinf,GRID_VARS%sorp,GRID_SOIL%xk0,GRID_SOIL%psic,GRID_SOIL%bcgamm,GRID_SOIL%bcbeta,deltrz,GRID_VARS%cc,&
       GRID_VARS%zw,GRID_VARS%xinact,GRID_VARS%satxr,GRID_VARS%xinfxr,&
       GRID_VARS%intstm,xinfcp,GRID_VARS%runtot,GRID_VARS%irntyp)

! ====================================================================
! Calculate actual rate of evaporation.
! ====================================================================

      if (GRID_VEG%ivgtyp.eq.0) then

! --------------------------------------------------------------------
! Option 1 : Bare soil surface.
! --------------------------------------------------------------------

         if ( (GRID_VARS%Swq.le.(0.d0))) then

! --------------------------------------------------------------------
! In case of absence of a snow pack, solve the energy balance for
! the bare soil.
! --------------------------------------------------------------------
            call ebsres(GLOBAL%inc_frozen,GLOBAL%irestype,rsoil,GRID_VARS%rzsm,GRID_SOIL%srespar1,GRID_VARS%tkact,&
       GRID_SOIL%srespar2,rzsm_u,GRID_SOIL%srespar3,ravd,iffroz,GRID_SOIL%thetas,GRID_VARS%tkmid,&
       GRID_SOIL%zmid,GLOBAL%zrzmax,smtmp,GRID_VARS%tzsm,GRID_VARS%smold,GRID_VARS%rzsmold,GRID_VARS%tzsmold,&
       GLOBAL%iopthermc,thermc1,thermc2,GRID_SOIL%thetar,heatcapold,GRID_SOIL%psic,GRID_SOIL%bcbeta,&
       GRID_SOIL%quartz,heatcap1,GRID_SOIL%ifcoarse,heatcap2,GRID_SOIL%rocpsoil,GRID_VARS%row,&
       GRID_VARS%cph2o,roa,GRID_VARS%cp,GRID_VARS%roi,thermc,heatcap,GRID_VARS%rzdthetaudtemp,&
       GRID_VARS%dshact,GRID_VEG%albd,GRID_VEG%emiss,rahd,ebscap,tcel,vppa,psychr,xlhv,&
       GRID_SOIL%zdeep,GRID_SOIL%Tdeepstep,GRID_MET%rsd,GRID_MET%rld,GLOBAL%toleb,GLOBAL%maxnri,GLOBAL%dt,i,tkel,&
       GRID_VEG%zww,GRID_VEG%za,GRID_MET%uzw,GRID_VEG%zpd,GRID_VEG%z0m,GRID_MET%press,rib,GRID_VARS%rnetpn,&
       GRID_VARS%gbspen,epetd,GRID_VARS%evtact,GRID_VARS%ievcon,&
       bsdew,GRID_VEG%z0h,GLOBAL%ioppet)

         else

! --------------------------------------------------------------------
! If there is snow, the evaporation is determined by the solution of
! the snow model.
! --------------------------------------------------------------------

            ebscap=epetd
            GRID_VARS%evtact = ebscap
            GRID_VARS%ievcon = 3

            call calcrain (tcel,snow,rain,GRID_MET%pptms,GLOBAL%dt)

            call calcsnowmelt(0,0,GLOBAL%dt/3600.d0,GRID_VEG%za,GRID_VEG%zpd,GRID_VEG%z0h,RaSnow,roa,vppa,xlhv,&
            GRID_MET%rsd*(1.d0-GRID_VARS%alb_snow),GRID_MET%rld,appa,rain,snow,tcel,vpsat-vppa,GRID_MET%uzw,&
            GRID_VARS%PackWater,GRID_VARS%SurfWater,GRID_VARS%Swq,GRID_VARS%VaporMassFlux,&
            GRID_VARS%TPack,GRID_VARS%TSurf,GRID_VARS%r_MeltEnergy,&
            GRID_VARS%Outflow,GRID_VARS%xleact_snow,GRID_VARS%hact_snow,GRID_VARS%rn_snow,1.d0,GRID_VARS%dens,gold)

            GRID_VARS%xleact_snow=0.d0-GRID_VARS%xleact_snow
            GRID_VARS%hact_snow=0.d0-GRID_VARS%hact_snow
            tsnow=GRID_VARS%TPack+273.15d0

            if (GRID_VARS%Swq.lt.(0.005d0)) tsnow=GRID_VARS%TSurf+273.15d0
            if (GRID_VARS%Swq.lt.(0.d0)) tsnow=tcel+273.15d0

            call nreb_snow(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       tkactd,tkmidactd,tsnow,GRID_SOIL%zdeep,GRID_SOIL%Tdeepstep,GRID_SOIL%zmid,GLOBAL%dt,dum)

            GRID_VARS%rnact=GRID_VARS%rn_snow
            GRID_VARS%xleact=GRID_VARS%xleact_snow
            GRID_VARS%hact=GRID_VARS%hact_snow
            GRID_VARS%gact=gactd
            GRID_VARS%dshact=0.d0
            GRID_VARS%tkact=tkactd
            GRID_VARS%tkmid=tkmidactd

         endif

      else

! --------------------------------------------------------------------
! Option 2: Vegetated surface.
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Calculate the transpiration from dry canopy for vegetated pixels.
! --------------------------------------------------------------------

         call transv(epetd,CELL_VARS%epetd_us,GRID_VEG%i_und,GLOBAL%iopveg,f1par,f3vpd,f4temp,&
       GRID_VEG%rescan,GLOBAL%inc_frozen,GRID_VEG%ivgtyp,GRID_VARS%rzsm,rzsm_u,GRID_VEG%tc,GRID_VEG%tw,&
       smcond,GRID_VARS%tzsm,tzsm_u,&
       GRID_VEG%tc_us,GRID_VEG%tw_us,smcond_us,f1par_us,f3vpd_us,f4temp_us,GRID_VEG%rescan_us,&
       vegcap,ravd,vegcap_us,ravd_us,GRID_VARS%zrz,srzrel,GRID_SOIL%thetas,GRID_SOIL%thetar,psisoi,&
       GRID_SOIL%psic,GRID_SOIL%bcbeta,GLOBAL%ikopt,xksrz,GRID_SOIL%xk0,CAT%ff,ressoi,GRID_VEG%rtact,&
       GRID_VEG%rtdens,GRID_VEG%psicri,&
       GRID_VEG%respla,xkrz,GRID_VARS%ztz,stzrel,xkstz,xktz,GRID_VARS%Swq,GRID_VARS%evtact,&
       GRID_VARS%ievcon,GRID_VARS%Swq_us,CELL_VARS%evtact_us,CELL_VARS%ievcon_us,bsdew,i,ipix)

      endif

! ====================================================================
! balance transmission and surface root zones
! ====================================================================

      if (GLOBAL%inc_frozen.eq.0) then

         rzsm_test=GRID_VARS%rzsm
         tzsm_test=GRID_VARS%tzsm
         rzsm_u_test=GRID_VARS%rzsm
         tzsm_u_test=GRID_VARS%tzsm
         thetas_add=GRID_SOIL%thetas

      endif

      if (GLOBAL%inc_frozen.eq.1) then

         rzsm_test=rzsm_u
         tzsm_test=tzsm_u
         rzsm_u_test=rzsm_u
         tzsm_u_test=tzsm_u
         thetas_add=GRID_SOIL%thetas-GRID_VARS%rzsm_f

      endif

      call tz_and_rzbal(i,GLOBAL%newstorm,GLOBAL%inc_frozen,GLOBAL%ikopt,GRID_VEG%ivgtyp,&
       GLOBAL%dt,rzsm_test,tzsm_test,rzsm_u_test,tzsm_u_test,GRID_VARS%rzsm1,GRID_VARS%tzsm1,&
       GRID_VARS%zrz,GRID_VARS%ztz,zw0,GLOBAL%zrzmax,&
       GRID_VARS%evtact,CELL_VARS%evtact_us,bsdew,dewrun,GRID_VARS%grz,GRID_VARS%gtz,GRID_VARS%diftz,GRID_VARS%difrz,&
       GRID_VARS%satxr,GRID_VARS%runtot,GRID_VARS%xinact,GRID_VARS%cuminf,&
       CAT%ff,GRID_SOIL%thetar,thetas_add,GRID_SOIL%bcbeta,GRID_SOIL%xk0,GRID_SOIL%psic,&
       GRID_VARS%Swq,GRID_VARS%Swq_us,&
       GRID_VARS%dc,GRID_VEG%i_und,GRID_VEG%i_moss,GRID_VARS%fw,dc_us,fw_us,evrz_moss,GRID_VEG%f_und,GRID_VARS%dstz,GRID_VARS%dsrz,&
       GRID_VARS%tzrhs,GRID_VARS%rzrhs)

      if (GLOBAL%inc_frozen.eq.1) then

         rzsm_u=GRID_VARS%rzsm1
         tzsm_u=GRID_VARS%tzsm1
         GRID_VARS%rzsm1_u=GRID_VARS%rzsm1
         GRID_VARS%tzsm1_u=GRID_VARS%tzsm1
         GRID_VARS%rzsm1=GRID_VARS%rzsm1_u+GRID_VARS%rzsm_f
         GRID_VARS%tzsm1=GRID_VARS%tzsm1_u+GRID_VARS%tzsm_f

      endif

! ====================================================================
! Calculate the actual surface energy fluxes - if PET values
! are used then skip this.
! ====================================================================

      if (GLOBAL%ioppet.ne.2) then

! --------------------------------------------------------------------
! Calculate the incoming solar radiation for the overstory.
! --------------------------------------------------------------------

         r_sdn=rs_over
         dshactd=zero
! --------------------------------------------------------------------
! Solve the over story energy balance.
! --------------------------------------------------------------------

         call land_os(rain,snow,thermc2,heatcap,heatcap2,heatcapold,&
       tkactd,tkmidactd,GRID_VEG%canclos,GRID_VARS%ievcon,xlhv,GRID_VARS%row,GRID_VEG%ivgtyp,xleactd,GRID_VARS%evtact,&
       bsdew,GLOBAL%ioppet,iffroz,GRID_VARS%tkmid,GRID_SOIL%zmid,GLOBAL%zrzmax,smtmp,GRID_VARS%rzsm,&
       GRID_VARS%tzsm,GRID_VARS%smold,GRID_VARS%rzsmold,GRID_VARS%tzsmold,GLOBAL%iopthermc,thermc1,&
       GRID_SOIL%thetar,GRID_SOIL%thetas,GRID_SOIL%psic,&
       GRID_SOIL%bcbeta,GRID_SOIL%quartz,GRID_SOIL%ifcoarse,heatcap1,GRID_SOIL%rocpsoil,GRID_VARS%cph2o,roa,&
       GRID_VARS%cp,GRID_VARS%roi,thermc,&
       GRID_VARS%rzdthetaudtemp,GLOBAL%iopgveg,GLOBAL%iopthermc_v,GRID_VEG%tcbeta,GRID_VEG%xlai,GRID_VARS%tkact,&
       GLOBAL%i_2l,f1,f2,f3,GRID_VEG%emiss,GRID_VEG%rescan,ravd,rahd,rnactd,&
       hactd,gactd,dshactd,tcel,vppa,psychr,GRID_SOIL%zdeep,GRID_SOIL%Tdeepstep,&
       GRID_MET%rsd,r_lup,GRID_MET%rld,GLOBAL%toleb,GLOBAL%maxnri,GLOBAL%dt,i,GRID_VEG%albd,r_sdn,GRID_VARS%rnetpn,&
       GRID_VARS%gbspen,rnetd,GRID_VEG%xled,GRID_VEG%hd,GRID_VEG%gd,GRID_VEG%dshd,tkd,tkmidd,&
       GRID_VARS%rnact,GRID_VARS%xleact,GRID_VARS%hact,&
       GRID_VARS%gact,GRID_VARS%dshact,GRID_VEG%rnetw,GRID_VEG%xlew,GRID_VEG%hw,GRID_VEG%gw,&
       GRID_VEG%dshw,GRID_VEG%tkw,GRID_VEG%tkmidw,GRID_VARS%dc,GRID_VARS%fw,tdiff,GLOBAL%inc_frozen,&
       ipix,initer,GRID_VARS%PackWater,GRID_VARS%SurfWater,GRID_VARS%Swq,GRID_VARS%VaporMassFlux,GRID_VARS%TPack,&
       GRID_VARS%TSurf,GRID_VARS%r_MeltEnergy,GRID_VARS%Outflow,GRID_VARS%xleact_snow,GRID_VARS%hact_snow,&
       GRID_VARS%dens,GRID_VARS%precip_o,GRID_VEG%za,&
       GRID_VEG%zpd,GRID_VEG%z0h,RaSnow,appa,vpsat,GRID_MET%uzw,GRID_VARS%rn_snow,GRID_VARS%alb_snow)

! --------------------------------------------------------------------
! Calculate the incoming solar radiation for the under story and
! solve the under story energy balance.
! --------------------------------------------------------------------

         r_sdn=rs_under

         call land_us(rain,snow,thermc1,thermc2,heatcap_moss,heatcap,&
       heatcap1,heatcap2,heatcapold,tkactd_us,tkmidactd_us,&
       tskinactd_moss,tkactd_moss,tkmidactd_moss,GRID_VEG%canclos,CELL_VARS%ievcon_us,&
       CELL_VARS%rnact_us,CELL_VARS%xleact_us,CELL_VARS%hact_us,CELL_VARS%gact_us,CELL_VARS%dshact_us,CELL_VARS%tkact_us,&
       CELL_VARS%tkmid_us,rnactd_us,CELL_VARS%rnetw_us,xleactd_us,CELL_VARS%xlew_us,hactd_us,CELL_VARS%hw_us,&
       gactd_us,CELL_VARS%gw_us,dshactd_us,CELL_VARS%dshw_us,CELL_VARS%tkw_us,CELL_VARS%tkmidw_us,GRID_VEG%xlai_us,&
       dc_us,fw_us,trlup,ipix,CELL_VARS%xlhv_ic,GRID_VARS%row,CELL_VARS%evtact_us,iffroz_us,GRID_VARS%tkmid,GRID_SOIL%zmid,&
       GLOBAL%zrzmax,smtmp,GRID_VARS%rzsm,GRID_VARS%tzsm,GRID_VARS%smold,GRID_VARS%rzsmold,GRID_VARS%tzsmold,GLOBAL%iopthermc,&
       GRID_SOIL%thetar,GRID_SOIL%thetas,GRID_SOIL%psic,GRID_SOIL%bcbeta,&
       GRID_SOIL%quartz,GRID_SOIL%ifcoarse,GRID_SOIL%rocpsoil,GRID_VARS%cph2o,roa,GRID_VARS%cp,GRID_VARS%roi,&
       thermc,GLOBAL%inc_frozen,GRID_VARS%rzdthetaudtemp,GLOBAL%iopgveg,thermc_us,GLOBAL%iopthermc_v,GRID_VEG%tcbeta_us,&
       GRID_VEG%xlai,f3,GRID_VEG%albd_us,GRID_VEG%emiss_us,ravd_us,rahd_us,GRID_VEG%rescan_us,CELL_VARS%tcel_ic,CELL_VARS%vppa_ic,&
       roa_ic,CELL_VARS%psychr_ic,GRID_SOIL%zdeep,GRID_SOIL%Tdeepstep,r_sdn,r_ldn,GLOBAL%toleb,GLOBAL%maxnri,GLOBAL%dt,i,&
       GRID_MET%rld,CELL_VARS%rnetd_us,CELL_VARS%xled_us,CELL_VARS%hd_us,CELL_VARS%gd_us,CELL_VARS%dshd_us,&
       CELL_VARS%tkd_us,CELL_VARS%tkmidd_us,initer,&
       CELL_VARS%ievcon_moss,xleactd_moss,CELL_VARS%bsdew_moss,CELL_VARS%evtact_moss,thermc_moss,&
       GRID_VARS%r_mossm,CELL_VARS%tskinact_moss,CELL_VARS%tkact_moss,CELL_VARS%tkmid_moss,hactd_moss,gactd_moss,&
       dshactd_moss,rav_moss,rah_moss,r_moss_depth,GRID_VEG%alb_moss,&
       rnactd_moss,emiss_moss,CELL_VARS%eact_moss,CELL_VARS%rnet_pot_moss,CELL_VARS%xle_p_moss,CELL_VARS%h_p_moss,&
       CELL_VARS%g_p_moss,CELL_VARS%tk_p_moss,CELL_VARS%tkmid_p_moss,CELL_VARS%tskin_p_moss,GRID_VARS%zmoss,&
       thetas_moss,CELL_VARS%rnact_moss,CELL_VARS%xleact_moss,CELL_VARS%hact_moss,&
       CELL_VARS%gact_moss,CELL_VARS%dshact_moss,gold,GRID_VARS%Swq_us,GRID_VARS%precip_u,GRID_VEG%za,GRID_VEG%zpd_us,&
       GRID_VEG%z0h,RaSnow,GRID_VARS%alb_snow,appa,CELL_VARS%vpsat_ic,GRID_MET%uzw,GRID_VARS%PackWater_us,&
       GRID_VARS%SurfWater_us,GRID_VARS%VaporMassFlux_us,GRID_VARS%TPack_us,GRID_VARS%TSurf_us,&
       GRID_VARS%r_MeltEnergy_us,GRID_VARS%Outflow_us,GRID_VARS%xleact_snow_us,GRID_VARS%hact_snow_us,&
       GRID_VARS%rn_snow_us,GRID_VARS%dens_us,heatcap_us,CELL_VARS%tkel_ic,eps,CELL_VARS%ds_p_moss,&
       GRID_VEG%i_und,GRID_VEG%i_moss,GLOBAL%i_2l)

      endif

! ====================================================================
! Print some results out.
! ====================================================================

      if (GRID_VARS%dens_us.gt.(0.d0)) then

         p1=GRID_VARS%PackWater_us*1000.d0
         p2=GRID_VARS%SurfWater_us*1000.d0
         p3=GRID_VARS%Swq_us*1000.d0
         r_dens=(p3/(p1+p2+p3))*GRID_VARS%dens_us+((p1+p2)/(p1+p2+p3))*1.d0

      endif

      if (GRID_VARS%dens.gt.(0.d0)) then

         p1=GRID_VARS%PackWater*1000.d0
         p2=GRID_VARS%SurfWater*1000.d0
         p3=GRID_VARS%Swq*1000.d0
         r_dens=(p3/(p1+p2+p3))*GRID_VARS%dens+((p1+p2)/(p1+p2+p3))*1.d0

      endif

5432  format (1i5," ",f11.6," ",2(f5.1," "),f7.3," ",2(f7.2," "),2(f6.3," "),f5.1,f6.1)

! ====================================================================
! Check whether the snow water equivalent on top of the overstory
! does not exceed its water holding capacity.
! ====================================================================

      if (GRID_VEG%ivgtyp.gt.0) then

            if (GRID_VARS%Swq.ge.GRID_VARS%wcip1) GRID_VARS%wcip1=GRID_VARS%Swq

      endif

      return
    end subroutine land

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

      endif

!CVAL      write (129,*)
!real(r_mossm_f/r_mossm),real(rzsm_f/rzsm),real(tzsm_f/tzsm)
!CVAL      write (130,*) real(zw-psic)

      return

    end subroutine states

! ====================================================================
!
!			subroutine initsm
!
! ====================================================================
!
! Subroutine to assign spatially variable initial conditions in
! the root and transmission zones for the first time step in the
! storm.
!
! ====================================================================

      subroutine initsm(zw,psic,zrz,ztz,rzsm1,tzsm1,thetas,&
       zrzmax,iopsmini,thetar,bcbeta,rzsm1_u,tzsm1_u,rzsm1_f,tzsm1_f,&
       inc_frozen,tsoilold,bulk_dens,a_ice,b_ice,row)

      implicit none
      include "help/initsm.h"

! ====================================================================
! Update root and transmission zone depths and soil moisture.
! ====================================================================

      if ((zw-psic).le.zero) then

! --------------------------------------------------------------------
! If the root zone is saturated (Region 3).
! --------------------------------------------------------------------

         zrz = zero
         ztz = zero
         rzsm1 = thetas
         tzsm1 = thetas

      else if ((zw-psic).lt.zrzmax) then

! --------------------------------------------------------------------
! If the transmission zone is saturated and root zone is
! unsaturated (Region 2).
! --------------------------------------------------------------------

         zrz = zw-psic

         if (iopsmini.eq.0) then

            rzsm1 = thetar+(thetas-thetar)*((psic/zw)**bcbeta)

         endif

         ztz = zero
         tzsm1 = thetas

      else

! --------------------------------------------------------------------
! If the transmission and root zone are both unsaturated (Region 1).
! --------------------------------------------------------------------

         zrz = zrzmax
         ztz = zw-psic-zrz
 
         if (iopsmini.eq.0) then

            rzsm1 = thetar+(thetas-thetar)*((psic/(zw-0.5*zrz))**bcbeta)
            tzsm1 = thetar+(thetas-thetar)*((psic/(0.5*ztz+psic))** bcbeta)

         endif

      endif

! ====================================================================
! Calculate the frozen water and liquid water fractions in the
! soil if requested.
! ====================================================================

      rzsm1_u=0.d0
      tzsm1_u=0.d0
      rzsm1_f=0.d0
      tzsm1_f=0.d0

      if (inc_frozen.eq.1) then

         if (tsoilold.gt.(273.15d0)) then

! --------------------------------------------------------------------
! If the soil is not frozen all soil water is liquid.
! --------------------------------------------------------------------

            rzsm1_u=rzsm1
            tzsm1_u=tzsm1
            rzsm1_f=0.d0
            tzsm1_f=0.d0

         else

! --------------------------------------------------------------------
! In case of frozen soil calculate the unfrozen soil water.
! --------------------------------------------------------------------

            ttt=273.15d0-tsoilold
            rzsm1_u=1000.d0*bulk_dens*a_ice*ttt**b_ice/row
            tzsm1_u=1000.d0*bulk_dens*a_ice*ttt**b_ice/row

            if (rzsm1_u.gt.rzsm1) then

! ....................................................................
! Check whether the unfrozen soil water in the root zone
! is not higher than the root zone soil moisture.  Adapt
! if necessary.
! ....................................................................

               rzsm1_u=rzsm1
               rzsm1_f=0.d0

            else

! ....................................................................
! Frozen water content = total water content minus unfrozen water
! content.
! ....................................................................

               rzsm1_f=rzsm1-rzsm1_u

            endif

! ....................................................................
! Check whether the unfrozen soil water in the transmission zone
! is not higher than the transmission zone soil moisture.  Adapt
! if necessary.
! ....................................................................

            if (tzsm1_u.gt.tzsm1) then

               tzsm1_u=tzsm1
               tzsm1_f=0.d0

            else

! ....................................................................
! Frozen water content = total water content minus unfrozen water
! content.
! ....................................................................

               tzsm1_f=tzsm1-tzsm1_u

            endif

         endif

! --------------------------------------------------------------------
! Frozen water content = total water content minus unfrozen water
! content.  (This, in some cases, has already been calculated
! here above but doint this calculation again dos not affect the
! results.
! --------------------------------------------------------------------

         rzsm1_u=rzsm1-rzsm1_f
         tzsm1_u=tzsm1-tzsm1_f

         if ( (thetas-rzsm1_f).le.(thetar+0.d0)) then

! ....................................................................
! Check whether the frozen water content in the root zone is not
! higher than the saturated water content.
! ....................................................................

            rtdif=0.0d0+thetar-(thetas-rzsm1_f)
            rzsm1_f=rzsm1_f-rtdif
            rzsm1_u=rzsm1_u+rtdif

         endif

         if ( (thetas-tzsm1_f).le.(thetar)+0.d0) then

! ....................................................................
! Check whether the frozen water content in the transmission zone is not
! higher than the saturated water content.
! ....................................................................

            rtdif=0.0d0+thetar-(thetas-tzsm1_f)
            tzsm1_f=tzsm1_f-rtdif
            tzsm1_u=tzsm1_u+rtdif

         endif

         if (rzsm1_u.le.(thetar+0.d0)) then

! ....................................................................
! Check whether the unfrozen water content in the root zone is not
! lower than the residual water content.
! ....................................................................

            rtdif=thetar+0.01d0-rzsm1_u
            rzsm1_u=rzsm1_u+rtdif
            rzsm1_f=rzsm1_f-rtdif

         endif

         if (tzsm1_u.le.(thetar+0.d0)) then

! ....................................................................
! Check whether the unfrozen water content in the transmission zone is
! not lower than the residual water content.
! ....................................................................

            rtdif=thetar+0.01d0-tzsm1_u
            tzsm1_u=tzsm1_u+rtdif
            tzsm1_f=tzsm1_f-rtdif

         endif

      endif

      return

      end subroutine initsm

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

      if ( (i_und) .gt.0) then

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
!                       subroutine reset_inf_pars
!
! ====================================================================
!
! Reset cumulative infiltration, initial root zone soil moisture,&
! sorptivity and dimensionless gravity parameter, cc.
!
! ====================================================================

      subroutine reset_inf_pars(cuminf,zero,rzsmst,rzsm,thetas,tolinf,&
       sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)

      implicit none
      include "help/reset_inf_pars.h"

! ====================================================================
! Calculate the root zone soil moisture
! ====================================================================

      cuminf = zero

      rzsmst = rzsm

      if (rzsmst.ge.thetas)  then

         rzsmst=thetas-tolinf

      endif

! ====================================================================
! Calculate sorptivity and gravity parameters at the first step
! of the storm event.
! ====================================================================

      sorp = (((two*xk0*((thetas-rzsmst)**two)*psic)&
                   /(thetas-thetar))*((one/(bcgamm+0.5d0*bcbeta-one))+&
                    ((thetas-thetar)/ (thetas-rzsmst))))**0.5d0

      deltrz = rzsmst-thetar
      if(deltrz.le.zero) deltrz=zero
      cc = 0.5d0*(one+ ((deltrz/(thetas-thetar))**(bcgamm/bcbeta)))

      return

      end subroutine reset_inf_pars


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

         call calcvegcap(smcond,zero,vegcap,epetd,resist,ravd)

! ====================================================================
! Calculate vegetation capacity for the under story using linear
! interpolation for the canopy resistance.
! ====================================================================

         if (i_und.gt.0) then

            call calcvegcap(smcond_us,zero,vegcap_us,&
                            epetd_us,resist_us,ravd_us)

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
       two,three,ressoi,rtact,rtdens,vegcap_us,psicri,respla)

         endif

! ====================================================================
! Calculate the maximum plant evaporation for the over story for
! vegetation with its roots in the upper soil layer.
! ====================================================================

         if (ivgtyp.eq.1) then

            call maxplevap(zrz,0.d0,epetd,inc_frozen,srzrel,rzsm,thetas,&
       thetar,rzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xksrz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla)

! ====================================================================
! Calculate the maximum plant evaporation for the over story for
! vegetation with its roots in the upper soil layer.
! ====================================================================

         else

            call maxplevap(ztz,zrz,epetd,inc_frozen,stzrel,tzsm,thetas,&
       thetar,tzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xkstz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla)

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
!                   subroutine maxplevap
!
! ====================================================================
!
! Calculate the maximum plant evaporation.
!
! ====================================================================

      subroutine maxplevap(zrz,ztz,epetd,inc_frozen,srzrel,rzsm,thetas,&
       thetar,rzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xksrz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla)

      implicit none
      include "help/maxplevap.h"

! --------------------------------------------------------------------
! If the soil is saturated than the maximum flux of the water is not
! bounded by the plant/soil system.
! --------------------------------------------------------------------

      if (zrz.le.zero) then

         vegcap = epetd

      else

! --------------------------------------------------------------------
! If the soil is not saturated calculate the soil saturation.
! --------------------------------------------------------------------

         if (inc_frozen.eq.0) then

! ....................................................................
! No difference between ice particles and frozen soil water.
! ....................................................................

            srzrel = (rzsm-thetar)/(thetas-thetar)

         else

! ....................................................................
! Ice crystals are treated as mineral soil particles.
! ....................................................................

            srzrel = (rzsm_u-thetar)/(thetas-thetar)

         endif

! --------------------------------------------------------------------
! Double check whether relative saturation is between 1 and 0.
! --------------------------------------------------------------------

         if (srzrel.le.zero) srzrel=zero
         if (srzrel.ge.one) srzrel=one

! --------------------------------------------------------------------
! Calculate the soil water potential.
! --------------------------------------------------------------------

         psisoi=psic/(srzrel**(one/bcbeta))

! --------------------------------------------------------------------
! Calculate the saturated hydrauli! conductivity in the root zone.
! --------------------------------------------------------------------

         if (ikopt.eq.1) then

! ....................................................................
! Option 1 : No decline with depth.
! ....................................................................

            xksrz = xk0

         else

! ....................................................................
! Option 2 : Exponential decline with depth.
! ....................................................................

            xksrz = xk0*dexp(-ff*((zrz+ztz)/two))

         endif

! --------------------------------------------------------------------
! Calculate the unsaturated hydrauli! conductivity in the root zone.
! --------------------------------------------------------------------

         xkrz = xksrz*(srzrel**((two+three*bcbeta)/bcbeta))

! --------------------------------------------------------------------
! Calculate the soil resistance for evaporation.
! --------------------------------------------------------------------

         ressoi = one/(rtact*xkrz*rtdens)

! --------------------------------------------------------------------
! Calculate the maximum water vapor flux out of the plant.
! --------------------------------------------------------------------

         vegcap = (psisoi-psicri)/(ressoi+respla)

! --------------------------------------------------------------------
! Check whether the maximum water vapor flux out of the plant is
! positive.
! --------------------------------------------------------------------

         if (vegcap.lt.zero) vegcap = zero

      endif

      return

      end subroutine maxplevap

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
!			subroutine difflx
!
! ====================================================================
!
! Subroutine difflx calculates the diffusive fluxes out of the surface 
! and transmission zones (positive down).
!
! ====================================================================

      subroutine new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

      implicit none
      include "help/new_difflx.h"

      data tolsat / 0.001d0 /

! ====================================================================
! Calculate diffusive flux out of root and transmission zones using 
! D(theta) relations of Brooks and Corey.
! ====================================================================

! ====================================================================
! If root zone not saturated then calculate the drainage.
! ====================================================================

      if (inc_frozen.eq.0) then

! --------------------------------------------------------------------
! If frozen soil water is treated as if it were not frozen.
! --------------------------------------------------------------------

	 if (zrz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
! ....................................................................

	    difrz=diffuse(zrz,ztz,rzsm,(0.5*rzsm+0.5*tzsm),&
                          thetas,thetar,bcbeta,psic,xksrz,xkstz)

	 else

! ....................................................................
! The soil profile is saturated.
! ....................................................................

	    difrz=zero

	 endif

      else

! --------------------------------------------------------------------
! If frozen soil water is treated as solid soil particles.
! --------------------------------------------------------------------

         if (zrz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
! ....................................................................

	    difrz=diffuse(zrz,ztz,rzsm_u,(0.5*rzsm_u+0.5*tzsm_u),&
                          thetas,thetar,bcbeta,psic,xksrz,xkstz)

         else
! ....................................................................
! The soil profile is saturated.
! ....................................................................

	   difrz=zero

	 endif

      endif

! ====================================================================
! Now repeat calculation for transmission zone.
! ====================================================================

      if (inc_frozen.eq.0) then

! --------------------------------------------------------------------
! If frozen soil water is treated as if it were not frozen.
! --------------------------------------------------------------------

         if (ztz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
!
! Here we assume diffusion occurs over layer of thickness equal to
! transmission zone and again account for any nonlinearity in soil
! moisture profile by averaging gradients.
! ....................................................................

	    diftz=diffuse(ztz,ztz,tzsm,0.5*tzsm+0.5*thetas,&
                          thetas,thetar,bcbeta,psic,xkstz,xkstz)

         else

! ....................................................................
! The soil profile is saturated.
! ....................................................................

	    diftz=zero

         endif

      else

! --------------------------------------------------------------------
! If frozen soil water is treated as solid soil particles.
! --------------------------------------------------------------------

         if (ztz.gt.zero) then

! ....................................................................
! The soil profile is not saturated.
! ....................................................................

	    diftz=diffuse(ztz,ztz,tzsm_u,0.5*tzsm_u+0.5*(thetas),&
                          thetas,thetar,bcbeta,psic,xkstz,xkstz)

         else

! ....................................................................
! The soil profile is saturated.
! ....................................................................

            diftz=zero

         endif

      endif
 
      return

      end subroutine new_difflx

! ====================================================================
!
!			subroutine diffuse
!
! ====================================================================
!
! Subroutine diffuse calculates the diffusive flux from layer 1 to layer 2
! in Richards' Equation using approach of Mahrt and Pan (1984).
!
! ====================================================================

      function diffuse(dz1,dz2,theta1,theta2,thetas,thetar,bcbeta,psic,&
                       xksat1,xksat2)

      implicit none
      include "help/diffuse.noup.h"

! ====================================================================
! Calculate diffusivity  using centered approx.
! ====================================================================

      theta = 0.5d0*(theta1 + theta2)

! ====================================================================
! Use harmoni! average since flow is perp. to layers.
! ====================================================================

      xksat = 1/(0.5/xksat1 + 0.5/xksat2)

      F1 = bcbeta * xksat * psic/(thetas-thetar) 

      satrel=(theta-thetar)/(thetas-thetar)

      if (satrel.lt.(0.d0)) satrel=0.d0

      DF = F1 * (satrel) ** (bcbeta + 2.)

! ====================================================================
! Calculate moisture gradient.
!
! dz1 and dz2 are layer thicknesses (not depths) over which diffusion
! is modeled.
! ====================================================================

      dz=0.5d0*(dz1+dz2)

      if (dz.gt.1.d-9) then

         grad= (theta1-theta2)/dz

      else

         grad=0.d0

      endif

! ====================================================================
! Calculate diffusive flux.
! ====================================================================

      difflx = DF * grad

      diffuse=difflx

      return

      end function diffuse

! ====================================================================
!
!			subroutine dwnflx
!
! ====================================================================
!
! Subroutine dwnflx calculates the downward fluxes out of the root 
! and transmission zones.
!
! ====================================================================

      subroutine new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,&
       inc_frozen,thetar,thetas,&
       rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

      implicit none
      include "help/new_dwnflx.h"

! ====================================================================
! Calculate downward flux out of root and transmission zones using 
! ksat(theta) relations of Brooks and Corey.
! ====================================================================

      if (zrz.eq.zero) then

! ====================================================================
! If root zone is saturated by the water table then set drainage
! to zero.
! ====================================================================

         grz = zero

! ====================================================================
! If root zone not saturated then calculate the drainage.
! ====================================================================

      else

! ====================================================================
! If root zone not saturated then calculate the drainage.
! ====================================================================


! --------------------------------------------------------------------
! Calcualate relative saturation in the root zone.
! --------------------------------------------------------------------

         if (inc_frozen.eq.0) then

! ....................................................................
! Option 1 : Treat frozen soil water as liquid water.
! ....................................................................

            relsrz = (rzsm-thetar)/(thetas-thetar)

         else

! ....................................................................
! Option 2 : Treat frozen soil particles as solid soil.
! ....................................................................

            relsrz = (rzsm_u-thetar)/(thetas-thetar)

         endif

         if (relsrz.le.zero) relsrz=zero
         if (relsrz.ge.one) relsrz=one

! --------------------------------------------------------------------
! Calculate downward flux out of root zone.
! --------------------------------------------------------------------

         grz = xksrz*(relsrz**((two+three*bcbeta)/bcbeta))

      endif

! ====================================================================
! Now repeat calculation for transmission zone.
! ====================================================================

      if (ztz.eq.zero) then

! ====================================================================
! If transmission zone is saturated by the water table then set drainage
! to zero.
! ====================================================================

         gtz = zero

      else

! ====================================================================
! If transmission zone not saturated then calculate the drainage.
! ====================================================================


! --------------------------------------------------------------------
! Calcualate relative saturation in the transmission zone.
! --------------------------------------------------------------------

         if (inc_frozen.eq.0) then

! ....................................................................
! Option 1 : Treat frozen soil water as liquid water.
! ....................................................................

            relstz=(tzsm-thetar)/(thetas-thetar)

         else

! ....................................................................
! Option 2 : Treat frozen soil particles as solid soil.
! ....................................................................

            relstz=(tzsm_u-thetar)/(thetas-thetar)

         endif

         if (relstz.le.zero) relstz=zero
         if (relstz.ge.one) relstz=one

! --------------------------------------------------------------------
! Calculate downward flux out of transmission zone.
! --------------------------------------------------------------------

         gtz = xkstz*(relstz**((two+three*bcbeta)/bcbeta))

      endif

      return

      end subroutine new_dwnflx

! ====================================================================
!
!                  subroutine clc_evrz
!
! ====================================================================
!
! Calculate the evaporation/condensation into the root zone.
!
! ====================================================================

      subroutine clc_evrz(evrz,Swq,Swq_us,ivgtyp,evtact,dc,i_und,&
       i_moss,fw,evtact_us,dc_us,fw_us,evrz_moss,dummy,f_und)

      implicit none
      include "help/clc_evrz.h"

      evrz=zero

      if ((Swq+Swq_us).le.(0.d0)) then

! ....................................................................
! If the vegetation is bare soil and there is no snow layer present
! all the evaporative demand comes from the bare soil and under story
! or moss is not represented.
! ....................................................................

         if (ivgtyp.eq.0) evrz=evtact*dc

         if (ivgtyp.eq.1) then

            if (i_und.eq.0)then

! ....................................................................
! In case of vegetation with lower roots all the evaporative demand
! for the soil comes from the over story if there is no under story
! represented and if there is no snow.
! ....................................................................

               evrz=evtact*dc*(1.d0-fw)
               dummy=evrz

            else

               if (i_und.gt.0) then

! ....................................................................
! If there is under story and no snow part of the evporative demand
! comes from the over story and part from the under story if under
! story is represented.
! ....................................................................

                  evrz=(evtact*dc*(1.d0-fw))+ &
                       f_und*(evtact_us*dc_us*(1.d0-fw_us))

                  dummy=evrz

               endif

               if (i_moss.gt.0) then

! ....................................................................
! If there is a moss layer and no snow then all the evaporative demand
! for the soil comes from the over story layer.  This is justified
! by the fact that roots under a moss layer  re sufficiently
! vertically distributed to make this assumption hold.
! ....................................................................

                  evrz=(evtact*dc*(1.d0-fw))
                  dummy=evrz

               endif

            endif

         endif

         if (ivgtyp.eq.2) then

            if ( (i_und.eq.0).and.&
                 (i_moss.eq.0) ) then

! ....................................................................
! In case of lower layer roots and no under story
! than there is no evaporative demand for the root zone.
! ....................................................................

               evrz=zero
               dummy=evrz

            else

               if (i_und.gt.0) then

! ....................................................................
! In case of lower layer roots and an under story layer than
! the evaporative demand for the root zone comes entirely from
! the under story.
! ....................................................................

                  evrz=zero+&
                       f_und*(evtact_us*dc_us*(1.d0-fw_us))
                  dummy=evrz

               endif

               if (i_moss.gt.0) then

! ....................................................................
! In case of lower layer roots and a moss layer
! than there is no evaporative demand for the root zone.
! ....................................................................

                  evrz=zero
                  dummy=evrz

               endif

            endif

         endif

      else

! ....................................................................
! In case of snow on top of the under story all the evaporative demand
! comes from the over story.
! ....................................................................

         evrz=zero

         if (Swq.le.0.d0) then

            if (ivgtyp.eq.1) then

               evrz=evrz+evtact*dc*(1.d0-fw)*(1.-f_und)

            endif

            if (ivgtyp.eq.0) evrz=evtact*dc

         endif
! ....................................................................
! In case of snow on top of the over story all the evaporative demand
! comes from the under story.
! ....................................................................

         if (Swq_us.le.0.d0) then

            if (i_und.gt.0) then

               evrz=evrz+evtact_us*dc_us*(1.d0-fw_us)*f_und

            endif

         endif

         dummy=evrz

      endif

! ....................................................................
! Water for evaporation of the over story is for 50 % extracted
! from soil water and 50 % from water in the moss layer.
! ....................................................................

      evrz_moss=0.d0

      return

      end subroutine clc_evrz

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
!       include "SNOW.h"
      include "help/land_os.h"

      if (((Swq.gt.(0.0002d0)).or.&
           ((tcel.le.(0.d0)).and.(dt*precip_o.gt.(0.0002d0))))) then

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
!       include "SNOW.h"
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

      if ( ( (Swq_us.gt.(0.0002d0)).or.&
               ((tcel_ic.le.(0.d0)).and.&
               (dt*precip_u.gt.(0.0002d0)) ) ))then

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


! ====================================================================
!
!		subroutine engact
!
! ====================================================================
!
! Subroutine to calculate the actual surface energy fluxes.
!
! ====================================================================

      subroutine engact(canclos,ievcon,xlhv,row,ivgtyp,xleactd,evtact,&
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

      implicit none
      include "help/engact.h"

      ccc=canclos

! ====================================================================
! If evaporation is land surface controled then recalculate
! actual fluxes.
! ====================================================================

      if (ievcon.eq.1) then

! --------------------------------------------------------------------
! Calculate actual latent heat flux from dry canopy or bare soil.
! --------------------------------------------------------------------

         if (ivgtyp.eq.0) then

            xleactd = xlhv * row * (evtact-bsdew)

         else

            xleactd = xlhv * row * evtact

         endif

! --------------------------------------------------------------------
! Solve energy balance again if this option is specified.
! --------------------------------------------------------------------

         if (ioppet.eq.0) then

            call sm_cen_dif(iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm,smold,&
                            rzsmold,tzsmold)

! ====================================================================
! Calculate the soil thermal parameters.
! ====================================================================

            call soiltherm(iopthermc,thermc1,thermc2,rzsm,smtmp,&
       thetar,thetas,psic,bcbeta,tkmid,quartz,ifcoarse,&
       heatcap1,heatcap2,heatcapold,rocpsoil,row,cph2o,roa,cp,roi,&
       smold,thermc,heatcap,inc_frozen,rzdthetaudtemp)

! --------------------------------------------------------------------
! Correct the thermal parameters for soils under vegetation.
! --------------------------------------------------------------------

            if (ivgtyp.ne.0) then

               call soiladapt(iopgveg,thermc,iopthermc_v,tcbeta,&
                              xlai,thermc1,heatcap,heatcap1,zero)

            endif

! ====================================================================
! Initialize actual temperatures.
! ====================================================================

            tkactd = tkact
            tkmidactd = tkmid

! ====================================================================
! Solve the energy balance for actual temperatures and fluxes.
! ====================================================================

            tdum1=tkactd
            tdum2=tkmidactd

            if (i_2l.eq.1) then

! --------------------------------------------------------------------
! Solve assuming that the radiation budgets for over story and under
! story are related with each other through the long wave radiation
! term.
! --------------------------------------------------------------------

               call nreb(0,f1+f2-f3,2.*emiss,thermc,thermc2,heatcap,heatcap2,&
       heatcapold,rescan+ravd,rahd,tdum1,tdum2,rnactd,xleactd,evtact,hactd,&
       gactd,dshactd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,&
       r_lup+rld,toleb,maxnri,dt,i)

            else

! --------------------------------------------------------------------
! Solve  assuming that the radiation budget for over and under story
! are independent from each other, which means that the incoming long
! wave radiation for both layers are equal.
! --------------------------------------------------------------------

               call nreb(0,albd,emiss,thermc,thermc2,heatcap,heatcap2,&
       heatcapold,rescan+ravd,rahd,tdum1,tdum2,rnactd,xleactd,evtact,hactd,&
       gactd,dshactd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,r_sdn,&
       rld,toleb,maxnri,dt,i)

            endif

! ====================================================================
! For penmen-montieth case just solve for sensible heat flux
! since Rnet, G and latent heat are known.
! ====================================================================

         else if(ioppet.eq.1) then

! --------------------------------------------------------------------
! Solve for sensible heat.
! --------------------------------------------------------------------

            hactd = rnetpn - xleactd - gbspen

! --------------------------------------------------------------------
! Set values for Rnet and Ground Heat Flux.
! --------------------------------------------------------------------

            rnactd = rnetpn
            gactd = gbspen

         endif

! ====================================================================
! If the evapotranspiration is atmosphere controled then set
! the fluxes equal to potential fluxes.
! ====================================================================

      else

         rnactd = rnetd
         xleactd = xled
         hactd = hd
         gactd = gd
         dshactd = dshd
         tdum1 = tkd
         tdum2 = tkmidd

      endif

! ====================================================================
! Compute average actual surface energy fluxes including any
! wet canopy fluxes.
! ====================================================================

      if (initer.eq.2) then

         tkactd=tdum1
         tkmidactd=tdum2

         rnact = rnactd*dc*(1-fw) + rnetw*(one-dc*(1-fw))
         xleact = xleactd*dc*(1-fw) + xlew*(one-dc*(1-fw))
         hact = hactd*dc*(1-fw) + hw*(one-dc*(1-fw))
         gact = gactd*dc*(1-fw) + gw*(one-dc*(1-fw))
         dshact = dshactd*dc*(1-fw) + dshw*(one-dc*(1-fw))
         tkact = tkactd*dc*(1-fw) + tkw*(one-dc*(1-fw))
         tkmid = tkmidactd*dc*(1-fw) + tkmidw*(one-dc*(1-fw))

      endif

      tdiff=tdum1*dc*(1-fw) + tkw*(one-dc*(1-fw))

! ====================================================================
! In case of unreasonable soil temperatures, which will lead to
! unreasonable ground heat flux values, abort the program and
! print the most important variables out.
! ====================================================================

      if(gactd.gt.15000) then

         write (*,*) 'ENGACT error '
         write (*,*) 'tkactd =',tkactd
         write (*,*) 'tkw =',tkw
         write (*,*) 'tkmidactd =',tkmidactd
         write (*,*) 'tkmidw =',tkmidw
         write (*,*) 'gactd =',gactd
         write (*,*) 'gw =',gw
         write (*,*) 'gact    =',gact
         write (*,*) 'dc    =',dc
         write (*,*) 'fw    =',fw 
         write (*,*) 'gactd    =',gactd
         write (*,*) 'gw    =',gw
         write (*,*) 'dc    =',dc
         write (*,*) 'fw    =',fw
         write (*,*) 'gact = ',gact
         write (*,*) 'therm! =',thermc
         write (*,*) 'thermc1 =',thermc1
         write (*,*) 'thermc2 =',thermc2
         write (*,*) 'heatcap1 =',heatcap1
         write (*,*) 'heatcap2 =',heatcap2
         write (*,*) 'heatcapold =',heatcapold
         write (*,*) 'smtmp = ',smtmp
         write (*,*) 'smold = ',smold
         write (*,*) 'ipix = ',ipix
         write (*,*) 'rzsm = ',rzsm
         write (*,*) 'tzsm = ',tzsm
         write (*,*) 'thetar = ',thetar
         write (*,*) 'thetas = ',thetas
         write (*,*) 'quartz = ',quartz
         write (*,*) 'iffroz = ',iffroz
         write (*,*) 'ifcoarse = ',ifcoarse
         write (*,*) 'rocpsoil = ',rocpsoil
         write (*,*) 'row = ',row
         write (*,*) 'cph2o = ',cph2o
         write (*,*) 'roa = ',roa
         write (*,*) 'cp = ',cp
         write (*,*) 'zmid = ',zmid
         write (*,*) 'zrzmax = ',zrzmax
         write (*,*) 'rzsmold = ',rzsmold
         write (*,*) 'tzsmold = ',tzsmold

         stop

      endif

      return

      end subroutine engact


! ====================================================================
!
!            subroutine calcrsoil
!
! ====================================================================
!
! Calculate the bare soil resistance to evaporation.
!
! ====================================================================

      subroutine calcrsoil(irestype,rsoil,&
                           srespar1,srespar2,srespar3,&
                           thetas,rzsm,tkact)

      implicit none
      include "help/calcrsoil.h"

      if (irestype.eq.2) then

! ....................................................................
! Formulation of Sun (1982).
! ....................................................................

         rsoil = (srespar1*(((thetas)/(rzsm))**srespar2)) + srespar3

      endif

      if (irestype.eq.3) then

! ....................................................................
! Formulation of Kondo (1990).
! ....................................................................

         D0 = 0.229d-04
         Datm = D0 * (tkact/273.16)**srespar3
         rsoil = (srespar1* ((thetas - rzsm) **srespar2)) / Datm

      endif

      if (irestype.eq.4) then

! ....................................................................
! Formulation of Camillo and Gurney (1986).
! ....................................................................

         rsoil = srespar1*(thetas-rzsm) - srespar2

         if (rsoil.lt.(0.d0)) rsoil=0.d0

      endif

      if (irestype.eq.5) then

! ....................................................................
! Formulation of Passerat (1986).
! ....................................................................

         rsoil = srespar1* dexp(srespar2*rzsm/ srespar3)

      endif

      return

      end subroutine calcrsoil

! ====================================================================
!
!                  subroutine calcsmcond
!
! ====================================================================
!
! Calculate the soil moisture conductance.
!
! ====================================================================

      subroutine calcsmcond(rzsm,tc,smcond,one,tw,zero)

      implicit none
      include "help/calcsmcond.h"

      if (rzsm.ge.tc) then

         smcond = one

      else if (rzsm.ge.tw) then

         smcond = (rzsm-tw)/(tc-tw)

      else

         smcond = zero

      endif

      if ( (smcond.ge.0.d0).and.(smcond.le.1.d0) ) then

         smcond=smcond

      else

         write (*,*) 'CALCSMCOND : smcond out of bounds ',smcond
         if (smcond.le.0.d0) smcond=zero
         if (smcond.ge.1.d0) smcond=one

      endif

      return

      end subroutine calcsmcond

! ====================================================================
!
!                   subroutine calcvegcap
!
! ====================================================================
!
! Calculate vegetation capacity for the over story using linear
! interpolation for the canopy resistance.
!
! ====================================================================

      subroutine calcvegcap(smcond,zero,vegcap,epetd,resist,ravd)

      implicit none
      include "help/calcvegcap.h"

      if (smcond.gt.zero) then

         vegcap = epetd*(resist+ravd)/(resist/smcond + ravd)

      else

         vegcap = zero

      endif

      return

      end subroutine calcvegcap

!====================================================================
!
!                   subroutine acttrans
!
! ====================================================================
!
! Set actual transpiration to the minimum of potential
! transpiration or vegetation capacity.
!
! ====================================================================

      subroutine acttrans(Swq,vegcap,epet,evtact,ievcon,zrz)

      implicit none
      include "help/acttrans.h"

      if (Swq.le.(0.d0)) then

! --------------------------------------------------------------------
! If no snow present present on top of the over story the transpiration
! is determined by the plant.
! --------------------------------------------------------------------

         if (vegcap.lt.epet) then

            evtact = vegcap
            ievcon = 1

         else

            evtact = epet

            if (zrz.gt.(0.d0)) then

               ievcon = 2

            else

               ievcon = 3

            endif

         endif

      else

! --------------------------------------------------------------------
! If there is snow the evaporation is determined by the flux out
! of the snow pack.
! --------------------------------------------------------------------

         evtact = epet
         ievcon = 3

      endif

      return

      end subroutine acttrans

! ====================================================================
!
!            subroutine zero_snowvar
!
! ====================================================================
!
! If there is no snow all the snow pack variablesa are set to zero.
!
! ====================================================================

      subroutine zero_snowvar(PackWater,SurfWater,Swq,VaporMassFlux,&
       TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow,dens)

      implicit none
      include "help/zero_snowvar.h"

      PackWater=0.d0
      SurfWater=0.d0
      Swq=0.d0
      VaporMassFlux=0.d0
      TPack=0.d0
      TSurf=0.d0
      r_MeltEnergy=0.d0
      Outflow=0.d0
      xleact_snow=0.d0
      hact_snow=0.d0
      dens=0.d0

      return

      end subroutine zero_snowvar

! ====================================================================
!
!                         function clcddif
!
! ====================================================================
!
! Calculate the derivative of the diffusive flux out of the surface zone
! and transmission zone with respect to soil moisture.
!
! ====================================================================

      subroutine clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,&
         ikopt,xksrz,xkstz,ff,&
         zrz,ztz,bcbeta,thetas,thetar,psic,tzsm)
        
      implicit none
      include "help/clcddif.h"

! ====================================================================
! Calculate the derivatives of the diffusive flux with respect to
! soil moisture
! ====================================================================
    
      if(zrz.eq.0.and.ztz.eq.0)then
       ddifrzdth1=0.0
       ddifrzdth2=0.0
      else     
        ddifrzdth1 = (.75*bcbeta*(2+bcbeta)*psic*(.5*rzsm - .5*tzsm)&
               *((.75*rzsm + .25*tzsm - thetar)/(-thetar + thetas))**(1+bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas)**2) 

        ddifrzdth1 = ddifrzdth1 + (.5*bcbeta*psic*((.75*rzsm - thetar + &
                .25*tzsm)/(- thetar + thetas))**(2+ bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas))

        ddifrzdth2 = (.25*bcbeta*(2+bcbeta)*psic*(.5*rzsm - .5*tzsm)&
               *((.75*rzsm + .25*tzsm - thetar)/(-thetar + thetas))**(1+bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas)**2)

        ddifrzdth2 = ddifrzdth2 - (.5*bcbeta*psic*&
               ((.75*rzsm - thetar + .25*tzsm)/(-thetar + thetas))**(2+ bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas))

      endif

      if(ztz.eq.0.)then
        ddiftzdth2 = 0.0
      else
        ddiftzdth2 = (.75*bcbeta*(2+bcbeta)*xkstz*psic*(.5*tzsm - .5*thetas)&
               *((.75*tzsm + .25*thetas - thetar)/(-thetar + thetas))**(1+bcbeta))&
               /(ztz*(-thetar + thetas)**2)  
      
        ddiftzdth2 = ddiftzdth2 + (.5*bcbeta*xkstz*psic*&
               ((.75*tzsm - thetar + .25*thetas)/(-thetar + thetas))**(2 + bcbeta))&
               /(ztz*(-thetar + thetas))
      endif
 
      return

      end subroutine clcddif

! ====================================================================
!
!                   function clcdg
!
! ====================================================================
!
! Calculate the derivative of gravity drainage out of the surface
! and transmission zone with respect to soil moisture.
!
! ====================================================================

      function clcdg(sm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)

      implicit none
      include "help/clcdg.h"

! --------------------------------------------------------------------
! Calculate the derivative of the downward flux
! --------------------------------------------------------------------

      dgrz = xksrz * (2.d0 + 3.d0 * bcbeta)/(bcbeta* (thetas- thetar))

! ....................................................................
! Calculate the relative soil saturation.
! ....................................................................
    
      satrel = (sm- thetar)/(thetas - thetar)

      dgrz = dgrz * satrel**(2.d0+(2.d0/bcbeta))
   
      clcdg=dgrz

      return

      end function clcdg

END MODULE MODULE_LAND
