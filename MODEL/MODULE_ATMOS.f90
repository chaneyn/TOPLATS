MODULE MODULE_ATMOS

USE MODULE_VARIABLES

USE MODULE_SHARED

USE MODULE_SNOW

implicit none

contains

! ====================================================================
!
!		subroutine atmos
!
! ====================================================================
!
! Subroutine to read in and pass meteorological data (e.g.
! rainfall) and calculate the potential evaporation for bare
! soil and the unstressed (potential) transpiration for vegetation.
!
! ====================================================================
!
! Note: Sign Convention: All radiative fluxes directed toward the 
!        surface are positive (e.g. net radiation).  All non-radiative 
!        fluxes directed away from surface are positive (e.g. latent,&
!        sensible and soil heat fluxes)
!
! ====================================================================



  subroutine atmos(ipix,i,GRID_VEG,CELL_VARS,&

! Meteorological data

       GRID_MET,tcel,vppa,psychr,xlhv,tkel,uzw,&
       appa,vpsat,&
       twet,&
       qv,ra,&

! Temperature variables

       GRID_VARS,tkmid,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,&
       
! Energy fluxes and states

       epetd,epetd_us,rnetd,&
       tkd,tkmidd,&
       tsoilold,&

       GRID_SOIL,&

! Vegetation parameters

       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,&
       f1,f2,f3,&
       f3vpdpar_us,&
       f4temppar_us,&

! Constants

       roa,roa_ic,&

! Energy balance variables

       ravd,rahd,rah_moss,&
       RaSnow,rib_us,ravw,ravw_us,rahw,rahw_us,&

       GLOBAL)

    implicit none
   ! include "help/atmos.h" !take this out when variables are fixed

    integer ipix,i

    real*8 tcel,vppa,psychr,xlhv,tkel
    real*8 uzw,appa,vpsat
    real*8 twet
    real*8 qv,ra,tkmid,tkmid_us,tkact_us
    real*8 tskinact_moss,tkact_moss,tkmid_moss
    real*8 epetd,epetd_us,rnetd
    real*8 tkd,tkmidd
    real*8 tsoilold
    real*8 f1par
    real*8 f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us
    real*8 f1,f2,f3
    real*8 f3vpdpar_us,f4temppar_us
    real*8 roa,roa_ic,ravd,rahd
    real*8 rah_moss,RaSnow,rib_us,ravw,ravw_us
    real*8 rahw,rahw_us
    real*8 zero,one,two,three,four,five,six,rrr,rrrr,vpdef

    data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
          3.d0,4.d0,5.d0,6.d0/

    type (GRID_MET_template) :: GRID_MET
    type (GRID_VEG_template) :: GRID_VEG
    type (GRID_VARS_template) :: GRID_VARS
    type (GRID_SOIL_template) :: GRID_SOIL
    type (GLOBAL_template) :: GLOBAL
    type (CELL_VARS_template) :: CELL_VARS

! Temporarily changing over variables from old to new format
uzw = GRID_MET%uzw
tkmid = GRID_VARS%tkmid
!rnetd, tkmidd, tkd, tcel, vppa, roa, are problems


! ====================================================================
! Define the albedo for the snow layer.
! ====================================================================

    GRID_VARS%alb_snow = 0.75

! ====================================================================
! Calculate the temperature of the deep soil layer.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : Assuming a constant temperature.
! --------------------------------------------------------------------

      if ( (GRID_SOIL%amp.eq.(0.d0)).and.(GRID_SOIL%phase.eq.(0.d0)).and.&
      (GRID_SOIL%shift.eq.(0.d0)) ) then

         GRID_SOIL%Tdeepstep=GRID_SOIL%tdeep

      else

! --------------------------------------------------------------------
! Option 2 : Assuming a cosine wave form.
! --------------------------------------------------------------------

         rrr=real(i)

         GRID_SOIL%Tdeepstep=GRID_SOIL%tdeep + GRID_SOIL%amp*cos ( rrr*GRID_SOIL%phase - GRID_SOIL%shift )

      endif

! ====================================================================
! Calculate the temperature in the canopy given the temperature
! above the canopy.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : Assuming there is no temperature difference.
! --------------------------------------------------------------------

      if ( (GRID_VEG%Tslope1.eq.(0.d0)).and.(GRID_VEG%Tint1.eq.(0.d0)).and.&
           (GRID_VEG%Tslope2.eq.(0.d0)).and.(GRID_VEG%Tint2.eq.(0.d0)).and.&
           (GRID_VEG%Tsep.eq.(0.d0)) ) then

         GRID_VARS%Tincan=GRID_MET%tdry

      else

! --------------------------------------------------------------------
! Option 2 : Assuming the temperature difference depends on the
! total incoming radiation.
! --------------------------------------------------------------------

         rrrr=GRID_MET%rld+GRID_MET%rsd

         if (rrrr.ge.GRID_VEG%Tsep) then

            GRID_VARS%Tincan=GRID_MET%tdry+GRID_VEG%Tint2+GRID_VEG%Tslope2*rrrr

         endif

         if (rrrr.lt.GRID_VEG%Tsep) then

            GRID_VARS%Tincan=GRID_MET%tdry+GRID_VEG%Tint1+GRID_VEG%Tslope1*rrrr

         endif

      endif

! ====================================================================
! Calculate the dew point temperature in the canopy given the dew point
! temperature above the canopy.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : Assuming there is no temperature difference.
! --------------------------------------------------------------------

      if ( (GRID_VEG%Twslope1.eq.(0.d0)).and.(GRID_VEG%Twint1.eq.(0.d0)).and.&
           (GRId_VEG%Twslope2.eq.(0.d0)).and.(GRID_VEG%Twint2.eq.(0.d0)).and.&
           (GRID_VEG%Twsep.eq.(0.d0)) ) then

         CELL_VARS%twet_ic=twet

         !if (GRID_VEG%r_moss_depth.lt.0.d0) stop

      else

! --------------------------------------------------------------------
! Option 2 : Assuming the temperature difference depends on the
! total incoming radiation.
! --------------------------------------------------------------------

         rrrr=GRID_MET%rld+GRID_MET%rsd

         if (rrrr.ge.GRID_VEG%Twsep) then

            CELL_VARS%twet_ic=twet+GRId_VEG%Twint2+GRID_VEG%Twslope2*rrrr

         endif

         if (rrrr.lt.GRID_VEG%Twsep) then

            CELL_VARS%twet_ic=twet+GRID_VEG%Twint1+GRID_VEG%Twslope1*rrrr

         endif

      endif

! ====================================================================
! Initialize temperature variables.
! ====================================================================

      tcel=GRID_MET%tdry
      tkel=tcel+273.15d0
      CELL_VARS%tcel_ic=GRID_VARS%Tincan
      CELL_VARS%tkel_ic=CELL_VARS%tcel_ic+273.15d0

! ====================================================================
! If first time step, use air temperature to initialize mid soil temp.
! ====================================================================

      if (i.eq.1) then

         call inittk(GRID_SOIL,GRID_VEG,GRID_VARS,GRID_SOIL%tdeep,&
       GRID_SOIL%tmid0,GRID_VEG%tmid0_moss,tkmid,&
       tkmid_us,tkmid_moss,tkel,&
       GRID_VEG%tk0moss,GRID_VARS%tkact,tkact_us,tkact_moss,&
       tskinact_moss,GRID_VARS%dshact,&
       CELL_VARS%dshact_us,CELL_VARS%dshact_moss,GRID_VARS%tkpet,GRID_VARS%tkmidpet,&
       CELL_VARS%tkmidpet_us,CELL_VARS%tkmidpet_moss,&
       GRID_VARS%dspet,CELL_VARS%dspet_us,CELL_VARS%dspet_moss,GRID_VARS%TSurf,GRID_VARS%TPack,&
       GRID_VARS%TSurf_us,GRID_VARS%TPack_us)

      endif

      tsoilold=GRID_VARS%tkmidpet

! ====================================================================
! Vapor pressure variables -- use different method depending
! if input includes wet bulb temperature or relative humidity
!
! If running as a two layer method READ THE HUMIDITY IN AS
! DEW POINT TEMPERATURES !
! ====================================================================

      appa=100.d0*GRID_MET%press
      vpsat=611.d0*dexp((17.27d0*tcel)/(237.3d0+tcel))
      CELL_VARS%vpsat_ic=611.d0*dexp((17.27d0*CELL_VARS%tcel_ic)/(237.3d0+CELL_VARS%tcel_ic))

      if (GLOBAL%iopwv.eq.0) then

         vppa=611.0d0*dexp((17.27d0*(twet))/(237.3d0+(twet)))
         GRID_MET%rh=100.*vppa/vpsat
         CELL_VARS%vppa_ic=611.0d0*dexp((17.27d0*(CELL_VARS%twet_ic))/&
         (237.3d0+(CELL_VARS%twet_ic)))
         GRID_VARS%rh_ic=100.*CELL_VARS%vppa_ic/CELL_VARS%vpsat_ic

      else

         vppa=0.01*GRID_MET%rh*vpsat
         CELL_VARS%vppa_ic=0.01*GRID_VARS%rh_ic*CELL_VARS%vpsat_ic

      endif

      qv=0.622d0*(vppa/appa)
      CELL_VARS%qv_ic=0.622d0*(CELL_VARS%vppa_ic/appa)

! ====================================================================
! Calculate wind if two components are input -- check to make sure
! wind is positive number.
! ====================================================================

      if (uzw.lt.(0.)) then

         write(*,*) 'uzw is negative - time step ',i,' pixel ',ipix

      endif

! ====================================================================
! Calculate thermodynamic values for air and water.
! ====================================================================

      ra=287.d0*(one+0.608d0*qv)
      roa=appa/(ra*tkel)
      xlhv =2.501d6-2370.d0*tcel
      psychr=(GRID_VARS%cp*appa)/(0.622d0*xlhv)

      CELL_VARS%ra_ic=287.d0*(one+0.608d0*CELL_VARS%qv_ic)
      roa_ic=appa/(CELL_VARS%ra_ic*CELL_VARS%tkel_ic)
      CELL_VARS%xlhv_ic=2.501d6-2370.d0*CELL_VARS%tcel_ic
      CELL_VARS%psychr_ic=(GRID_VARS%cp*appa)/(0.622d0*CELL_VARS%xlhv_ic)



! ====================================================================
! Now, if requested, read in the Richardson Number for stability
! correction for aerodynamic resistance.  If no stability
! correction or first time step then set the Richardson number to zero.
!
! You can only use the stability correction if you are solving
! for a skin temperature i.e. GLOBAL%ioppet = 0
! 
! Do this for overstory, understory and moss.
! ====================================================================

      if (GLOBAL%iopstab.eq.1.and.i.gt.1.and.GLOBAL%ioppet.eq.0) then 
         
         call stabcor(GRID_VEG%zww,GRID_VEG%za,uzw,GRID_VEG%zpd,GRId_VEG%z0m,&
       tkel,GRID_MET%press,&
       GRID_VARS%tkact,vppa,GRID_VARS%rib)
         
      else

         GRID_VARS%rib= zero

      endif


! ====================================================================
! Calculate aerodynami! resistances to heat and mass transfer
! (assumed equal)--including stability correction if rib!=0.
!
! Again, do this for over story, under story and moss.
! Also do this for snow.
! ====================================================================

      rahd=calcra(uzw,GRID_VEG%zww,GRID_VEG%za,GRID_VEG%zpd,GRID_VEG%z0m,&
      GRID_VEG%z0h,GRID_VARS%rib)
      rahw=calcra(uzw,GRID_VEG%zww,GRID_VEG%za,GRId_VEG%zpd,GRID_VEG%z0m,&
      GRID_VEG%z0h,GRID_VARS%rib)

      ravd=rahd
      ravw=rahw

      RaSnow=calcra(uzw,GRID_VEG%zww,GRID_VEG%za,GRID_VEG%zpd,0.005d0,0.0005d0,1.d0)

! ====================================================================
! Choose option to calculate potentials with Penman and Penman-Monteith,&
! or by solving the nonlinear energy balance equations.
!
! TAKE CARE !
!
! If there is a moss layer OR if the program is run as a two-layer
! program, either the detailed or simplified model, SOLVE THE
! NONLINEAR ENERGY BALANCE EQUATIONS.  The program has not been set
! up for Penman or Penman-Monteith in these cases yet !
! ====================================================================
      GRID_MET%uzw = uzw

      if(GLOBAL%ioppet.eq.0)then

            call peteb(ipix,i,GLOBAL%dt,GLOBAL%inc_frozen,&

! General vegetation parameters

       GRID_VEG,GRID_VEG%i_und,GRID_VEG%i_moss,GRID_VEG%ivgtyp,&

! Snow pack variables

       GRID_VARS%PackWater,GRID_VARS%SurfWater,GRID_VARS%Swq,GRID_VARS%VaporMassFlux,GRID_VARS%TPack,&
       GRID_VARS%TSurf,&
       GRID_VARS%r_MeltEnergy,GRID_VARS%Outflow,GRID_VARS%xleact_snow,GRID_VARS%hact_snow,&
       GRID_VARS%rn_snow,GRID_VARS%PackWater_us,&
       GRID_VARS%SurfWater_us,GRID_VARS%Swq_us,GRID_VARS%VaporMassFlux_us,GRID_VARS%TPack_us,&
       GRID_VARS%TSurf_us,&
       GRID_VARS%r_MeltEnergy_us,GRID_VARS%Outflow_us,GRID_VARS%xleact_snow_us,GRID_VARS%hact_snow_us,&
       GRID_VARS%rn_snow_us,GRID_VARS%dens,GRID_VARS%dens_us,&

! Albedos of the over story, under story,&
! and moss layer

       GRID_VEG%alb_moss,GRID_VARS%alb_snow,GRID_VEG%albd,GRID_VEG%albw,GRID_VEG%albw_us,&

! Meteorological data

       GRID_MET,GRID_MET%rsd,GRID_MET%rld,tcel,vppa,psychr,xlhv,tkel,GRID_VEG%zww,GRID_VEG%za,uzw,GRID_MET%press,&
       appa,vpsat,CELL_VARS%tcel_ic,CELL_VARS%vppa_ic,CELL_VARS%psychr_ic,CELL_VARS%xlhv_ic,&
       CELL_VARS%tkel_ic,CELL_VARS%vpsat_ic,&

! Temperature variables

       GRID_VARS,tkmid,GRID_VARS%tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,GRID_SOIL%Tdeepstep,&

! Energy fluxes and states

       GRID_VARS%dshact,epetd,GRID_VARS%gact,epetd_us,CELL_VARS%dshact_moss,CELL_VARS%xle_act_moss,&
       rnetd,GRID_VEG%xled,GRID_VEG%hd,&
       GRID_VEG%gd,GRID_VEG%dshd,tkd,tkmidd,GRID_VEG%rnetw,GRID_VEG%xlew,GRID_VEG%hw,GRID_VEG%gw,&
       GRID_VEG%dshw,GRID_VEG%tkw,&
       GRID_VEG%tkmidw,CELL_VARS%tskinactd_moss,CELL_VARS%tkactd_moss,CELL_VARS%tkmidactd_moss,&
       CELL_VARS%ds_p_moss,GRID_VARS%epetw,&
       CELL_VARS%dshact_us,CELL_VARS%rnetw_us,CELL_VARS%xlew_us,CELL_VARS%hw_us,CELL_VARS%gw_us,&
       CELL_VARS%dshw_us,CELL_VARS%tkw_us,CELL_VARS%tkmidw_us,CELL_VARS%epetw_us,&
       CELL_VARS%rnetd_us,CELL_VARS%xled_us,CELL_VARS%hd_us,CELL_VARS%gd_us,&
       CELL_VARS%dshd_us,CELL_VARS%tkd_us,&
       CELL_VARS%tkmidd_us,CELL_VARS%rnet_pot_moss,CELL_VARS%xle_p_moss,&
       CELL_VARS%h_p_moss,CELL_VARS%g_p_moss,CELL_VARS%tk_p_moss,CELL_VARS%tkmid_p_moss,&
       CELL_VARS%tskin_p_moss,CELL_VARS%eact_moss,&

! Soil parameters

       GRID_SOIL,GRID_SOIL%thetar,GRID_SOIL%thetas,GRID_SOIL%psic,GRID_SOIL%bcbeta,GRID_SOIL%quartz,&
       GRID_SOIL%ifcoarse,GRID_SOIL%rocpsoil,GRID_VEG%tcbeta,&
       GRID_VEG%tcbeta_us,GRID_SOIL%zdeep,GRID_SOIL%zmid,GLOBAL%zrzmax,&

! Moss parameters

       GRID_VEG%r_moss_depth,GRID_VEG%eps,GRID_VEG%emiss_moss,GRID_VEG%zpd_moss,CELL_VARS%rib_moss,&
       GRID_VEG%z0m_moss,GRID_VEG%z0h_moss,CELL_VARS%epet_moss,&

! Vegetation parameters

       GRID_VEG%xlai,GRID_VEG%xlai_us,GRID_VEG%emiss,GRID_VEG%zpd,GRID_VEG%zpd_us,GRID_VEG%z0m,&
       GRID_VEG%z0h,GRID_VEG%z0m_us,GRID_VEG%z0h_us,&
       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,GRID_VEG%rescan,&
       GRID_VEG%rescan_us,f1,f2,f3,GRID_VEG%emiss_us,GRID_VEG%rsmin,GRID_VEG%rsmax,&
       GRID_VEG%rsmin_us,GRID_VEG%rsmax_us,GRID_VEG%Rpl,GRID_VEG%Rpl_us,GRID_VEG%f3vpdpar,f3vpdpar_us,&
       GRID_VEG%trefk,GRID_VEG%f4temppar,&
       GRID_VEG%trefk_us,f4temppar_us,&

! Constants

       GRID_VARS%row,GRID_VARS%cph2o,roa,GRID_VARS%cp,GRID_VARS%roi,GLOBAL%toleb,&
       GLOBAL%maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,CELL_VARS%ravd_us,CELL_VARS%rahd_us,CELL_VARS%rav_moss,rah_moss,&
       GRID_VARS%rib,RaSnow,rib_us,ravw,ravw_us,rahw,rahw_us,&

! Water balance variables

       GRID_VARS%rzsm,GRID_VARS%tzsm,GRID_VARS%rzsm1,GRID_VARS%tzsm1,GRID_VARS%r_mossm,&
       GRID_VARS%zrz,GRID_VARS%smold,GRID_VARS%rzdthetaudtemp,GLOBAL%smpet0,&

! Different option paramters

       GLOBAL%iopthermc,GLOBAL%iopgveg,GLOBAL%iopthermc_v,GLOBAL%iopstab,GLOBAL%iopsmini,GLOBAL)

      else if(GLOBAL%ioppet.eq.1)then

        call petpen(GRID_VEG,GRID_MET,GRID_VARS,tcel,vpsat,vpdef,f1par,GRID_VEG%albd,&
       GRID_VEG%xlai,GRID_MET%rsd,GRID_VEG%rsmin,GRID_VEG%rsmax,GRID_VEG%Rpl,&
       tkel,vppa,f3vpd,GRID_VEG%f3vpdpar,f4temp,GRID_VEG%trefk,&
       GRID_VEG%f4temppar,GRID_VARS%rnetpn,GRID_VARS%gbspen,rnetd,GRID_VEG%rnetw,GRID_VEG%gd,GRID_VEG%gw,&
       GRID_VEG%rescan,ravd,xlhv,&
       GRID_VARS%row,epetd,GRID_VARS%epetw,ravw,psychr,GRID_VEG%xled,GRID_VEG%xlew,GRID_VEG%hd,&
       GRID_VEG%hw,GRID_VARS%cp,roa)
 
      endif

      GRID_MET%uzw = uzw
      GRID_VARS%tkmid = tkmid

      return

  end subroutine atmos

    

! ====================================================================
!
!			subroutine inittk
!
! ====================================================================
!
! Subroutine to assign initial soil temperatures
! and set heat storage variables for the first time step in the
! simulation.
!
! ====================================================================

  subroutine inittk(GRID_SOIL,GRID_VEG,GRID_VARS,tdeep,tmid0,tmid0_moss,tkmid,&
       tkmid_us,tkmid_moss,tkel,tk0moss,tkact,tkact_us,&
       tkact_moss,tskinact_moss,dshact,dshact_us,dshact_moss,tkpet,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,dspet,dspet_us,dspet_moss,&
       Tsurf,Tpack,Tsurf_us,Tpack_us)

      implicit none
      include "help/inittk.h" 
      
      type (GRID_SOIL_template) :: GRID_SOIL
      type (GRID_VEG_template) :: GRID_VEG
      type (GRID_VARS_template) :: GRID_VARS

!Temporary storing of variables
tdeep = GRID_SOIL%tdeep
tmid0 = GRID_SOIL%tmid0
tmid0_moss = GRID_VEG%tmid0_moss
tkmid = GRID_VARS%tkmid

!tkmid_us, tkel,tkact, tkact_us, tkact_moss, tskinact_moss, 
!dshact_us (i), dshact_moss (i), tkmidpet_us (i), tkmidpet_moss (i)
! dspet_us (i), dspet_moss (i), TSurf, TPack, TSurf_us, TPack_us
! is not in struct

tk0moss = GRID_VEG%tk0moss
dshact = GRID_VARS%dshact
tkpet = GRID_VARS%tkpet
tkmidpet = GRID_VARS%tkmidpet
dspet = GRID_VARS%dspet




! ====================================================================
! Initialize average intermediate soil temperatures.
! ====================================================================

      ttemp=tdeep
      tkmidtmp=tmid0
      tkmidtmp_us=tmid0
      tkmidtmp_moss=tmid0_moss

! ====================================================================
! Initialize distributed intermediate soil temperatures.
! ====================================================================

      tkmid = tkmidtmp
      tkmid_us = tkmidtmp_us
      tkmid_moss = tkmidtmp_moss
      tkact = tkel
      tkact_us = tkel
      tkact_moss = tk0moss
      tskinact_moss = tkel

! ====================================================================
! Initialize storages.
! ====================================================================

      dshact= zero
      dshact_us= zero
      dshact_moss= zero

! ====================================================================
! Initialize skin temperatures.
! ====================================================================

      tkpet=ttemp
      tkmidpet = tkmidtmp
      tkmidpet_us = tkmidtmp_us
      tkmidpet_moss = tkmidtmp_moss

! ====================================================================
! Initialize storages.
! ====================================================================

      dspet=zero
      dspet_us=zero
      dspet_moss=zero

! ====================================================================
! Initialize snow temperatures.
! ====================================================================

      Tsurf=tkel-273.15d0
      Tpack=tkel-273.15d0
      Tsurf=-10.d0
      Tpack=-10.d0
      Tsurf_us=tkel-273.15d0
      Tpack_us=tkel-273.15d0

      GRID_SOIL%tdeep = tdeep
      GRID_SOIL%tmid0 = tmid0
      GRID_VEG%tmid0_moss = tmid0_moss
      GRID_VARS%tkmid = tkmid

      GRID_VEG%tk0moss = tk0moss
      GRID_VARS%dshact = dshact
      GRID_VARS%tkpet = tkpet
      GRID_VARS%tkmidpet = tkmidpet
      GRID_VARS%dspet = dspet

      return

  end subroutine inittk
  
! ====================================================================
!
!		subroutine petpen
!
! ====================================================================
!
! Subroutine to calculate the potential evapotranspiration
! using penman's equation for bare soil and penman-monteith
! for vegetated areas.
!
! ====================================================================
!
! NOTE: sign convention: all radiative fluxes directed toward the 
!         surface are positive (e.g. net radiation).  All 
!         non-radiative fluxes directed away from surface are 
!         positive (e.g. latent, sensible and soil heat fluxes)
!
! ====================================================================

  subroutine petpen(GRID_VEG,GRID_MET,GRID_VARS,tcel,vpsat,vpdef,f1par,&
       albd,xlai,rsd,rsmin,rsmax,Rpl,tkel,vppa,f3vpd,f3vpdpar,f4temp,trefk,&
       f4temppar,rnetpn,gbspen,rnetd,rnetw,gd,gw,rescan,ravd,xlhv,&
       row,epetd,epetw,ravw,psychr,xled,xlew,hd,hw,cp,roa)

      implicit none
      include "help/petpen.h"

      type (GRID_MET_template) :: GRID_MET
      type (GRID_VEG_template) :: GRID_VEG
      type (GRID_VARS_template) :: GRID_VARS

      albd = GRID_VEG%albd
      xlai = GRID_VEG%xlai
      rsd = GRID_MET%rsd
      rsmin = GRID_VEG%rsmin
      rsmax = GRID_VEG%rsmax
      Rpl = GRID_VEG%Rpl
      trefk = GRID_VEG%trefk
      rnetpn = GRID_VARS%rnetpn
      gbspen = GRID_VARS%gbspen
      rescan = GRID_VEG%rescan
      row = GRID_VARS%row
      cp = GRID_VARS%cp



! ====================================================================
! Calculate the atmospheri! vapor pressure deficit and slope of
! vapor pressure-temperature curve.
! ====================================================================

      vpsat = 611.d0 * dexp((17.27d0*tcel)/(237.3d0+tcel))
      vpdef = vpsat - vppa
      dvpsdt = 4098.d0 * vpsat/((237.3d0+tcel)**two)

! ====================================================================
! Calculate the environmental controls on stomatal
! ====================================================================

! --------------------------------------------------------------------
! Photosynthetically active radiation (PAR).
! --------------------------------------------------------------------
     

      f1par = 1

! --------------------------------------------------------------------
! Vapor pressure deficit.
! --------------------------------------------------------------------
      f3vpd = clcf3vpd(tkel,vppa,f3vpdpar)

! --------------------------------------------------------------------
! Air temperature.
! --------------------------------------------------------------------

      f4temp = clcf4temp(tkel,trefk,f4temppar)

! ====================================================================
! Assign inputs to each net radiation variable and ground heat flux.
! ====================================================================

      rnetd = rnetpn
      rnetw = rnetpn
      gd = gbspen
      gw = gbspen

! ====================================================================
! Calculate potential evapotranspiration rate from dry canopy or
! bare soil.
! ====================================================================

      pstar = psychr * (one+f1par*f3vpd*f4temp*rescan/ravd)

      if (1.eq.0) then

! --------------------------------------------------------------------
! This is the original formulation without environmental control
! on stomatal resistance that we do not want to reach.
! --------------------------------------------------------------------
        epetdmf = ((dvpsdt*(rnetd-gd)) +&
                  ((cp*roa*vpdef)/ravd)) /&
                  ((dvpsdt+pstar)*xlhv)

      endif

      epetdmf = ((dvpsdt*(rnetd-gd)) +&
                 ((cp*roa*vpdef)/&
                  (f1par*f3vpd*f4temp*rescan+ravd))) /&
                ((dvpsdt+pstar)*xlhv)


      epetd=epetdmf/row

! ====================================================================
! Calculate potential evaporation rate from wet canopy.
! ====================================================================
     
      epetwmf = ((dvpsdt*(rnetw-gw)) +&
                  ((cp*roa*vpdef)/ravw)) /&
                 ((dvpsdt+psychr)*xlhv)
      epetw = epetwmf/row

! ====================================================================
! Calculate potential latent heat fluxes.
! ====================================================================

      xled = epetdmf*xlhv
      xlew = epetwmf*xlhv
   
! ====================================================================
! Calculate sensible heat fluxes using potential evapotranspiration
! rate.
! ====================================================================

      hd = rnetd-xled-gd
      hw = rnetw-xlew-gw
      
      GRID_VEG%albd = albd
      GRID_VEG%xlai = xlai
      GRID_MET%rsd = rsd
      GRID_VEG%rsmin = rsmin
      GRID_VEG%rsmax = rsmax
      GRID_VEG%Rpl = Rpl
      GRID_VEG%trefk = trefk
      GRID_VARS%rnetpn = rnetpn
      GRID_VARS%gbspen = gbspen
      GRID_VEG%rescan = rescan
      GRID_VARS%row = row
      GRID_VARS%cp = cp    

      return

      end subroutine petpen

! ====================================================================
!
!		subroutine peteb
!
! ====================================================================
!
! Subroutine to calculate the potential evapotranspiration
! using the combination method (i.e. solving the
! set of nonlinear energy balance equations iteratively).
! three energy balances are maintained: one for bare soil, one
! for wet vegetation and one for dry vegetation.
!
! This is the subroutine that solves the energy balance under the
! assumption that the incoming long wave radiation for both under
! and over story is equal.
!
! ====================================================================
!
! NOTE: sign convention: all radiative fluxes directed toward the 
!         surface are positive (e.g. net radiation).  All 
!         non-radiative fluxes directed away from surface are 
!         positive (e.g. latent, sensible and soil heat fluxes)
!
! ====================================================================

      subroutine peteb(ipix,i,dt,inc_frozen,&

! General vegetation parameters

       GRID_VEG,i_und,i_moss,ivgtyp,&

! Snow pack variables

       PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,PackWater_us,&
       SurfWater_us,Swq_us,VaporMassFlux_us,TPack_us,TSurf_us,r_MeltEnergy_us,&
       Outflow_us,xleact_snow_us,hact_snow_us,rn_snow_us,dens,dens_us,&

! Albedos of the over story, under story,&
! and moss layer

       alb_moss,alb_snow,albd,albw,albw_us,&

! Meteorological data

       GRID_MET,rsd,rld,tcel,vppa,psychr,xlhv,tkel,zww,za,uzw,press,&
       appa,vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,&

! Temperature variables

       GRID_VARS,tkmid,tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,Tdeepstep,&

! Energy fluxes and states

       dshact,epetd,gact,epetd_us,dshact_moss,xle_act_moss,rnetd,xled,hd,&
       gd,dshd,tkd,tkmidd,rnetw,xlew,hw,gw,dshw,tkw,&
       tkmidw,tskinactd_moss,tkactd_moss,tkmidactd_moss,ds_p_moss,epetw,&
       dshact_us,rnetw_us,xlew_us,hw_us,gw_us,&
       dshw_us,tkw_us,tkmidw_us,epetw_us,&
       rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us,&
       tkmidd_us,rnet_pot_moss,xle_p_moss,&
       h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,eact_moss,&

! Soil parameters

       GRID_SOIL,thetar,thetas,psic,bcbeta,quartz,ifcoarse,&
       rocpsoil,tcbeta,tcbeta_us,zdeep,zmid,zrzmax,&

! Moss parameters

       r_moss_depth,eps,emiss_moss,zpd_moss,rib_moss,&
       z0m_moss,z0h_moss,epet_moss,&

! Vegetation parameters

       xlai,xlai_us,emiss,zpd,zpd_us,z0m,z0h,z0m_us,z0h_us,&
       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,rescan,rescan_us,&
       f1,f2,f3,emiss_us,rsmin,rsmax,rsmin_us,rsmax_us,Rpl,Rpl_us,f3vpdpar,&
       f3vpdpar_us,trefk,f4temppar,trefk_us,f4temppar_us,&

! Constants

       row,cph2o,roa,cp,roi,toleb,maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,&
       rib,RaSnow,rib_us,ravw,ravw_us,rahw,rahw_us,&

! Water balance variables

       rzsm,tzsm,rzsm1,tzsm1,r_mossm,zrz,smold,rzdthetaudtemp,smpet0,&

! Different option paramters

       iopthermc,iopgveg,iopthermc_v,iopstab,iopsmini,GLOBAL)

      implicit none
      include "help/peteb.h"


      type (GRID_MET_template),intent(inout) :: GRID_MET
      type (GRID_VEG_template),intent(inout) :: GRID_VEG
      type (GRID_VARS_template),intent(inout) :: GRID_VARS
      type (GRID_SOIL_template),intent(inout) :: GRID_SOIL
      type (GLOBAL_template),intent(inout) :: GLOBAL


!General Vegetation parameters
!canclos = GRID_VEG%canclos
!extinct = GRID_VEG%extinct
i_und = GRID_VEG%i_und
i_moss = GRID_VEG%i_moss
ivgtyp = GRID_VEG%ivgtyp

!Snow Pack variables
PackWater = GRID_VARS%PackWater
SurfWater = GRID_VARS%SurfWater
Swq = GRID_VARS%Swq
VaporMassFlux = GRID_VARS%VaporMassFlux
r_MeltEnergy = GRID_VARS%r_MeltEnergy
Outflow = GRID_VARS%Outflow
PackWater_us = GRID_VARS%PackWater_us
SurfWater_us = GRID_VARS%SurfWater_us
Swq_us = GRID_VARS%Swq_us
VaporMassFlux_us = GRID_VARS%VaporMassFlux_us
r_MeltEnergy_us = GRID_VARS%r_MeltEnergy_us
Outflow_us = GRID_VARS%Outflow_us

!Albedos of the over story, under story, and moss layer
albd_us = GRID_VEG%albd_us
alb_moss = GRID_VEG%alb_moss
albd = GRID_VEG%albd
albw = GRID_VEG%albw
albw_us = GRID_VEG%albw_us
alb_snow = GRID_VARS%alb_snow

!Meteorological data
rsd = GRID_MET%rsd
rld = GRID_MET%rld
zww = GRID_VEG%zww
za = GRID_VEG%za
uzw = GRID_MET%uzw
press = GRID_MET%press

!Temperature variables
tkmid = GRID_VARS%tkmid
tkact = GRID_VARS%tkact

!Energy fluxes and states
dshact = GRID_VARS%dshact
gact = GRID_VARS%gact
rnetd = GRID_VEG%rnetd
xled = GRID_VEG%xled
hd = GRID_VEG%hd
gd = GRID_VEG%gd
dshd = GRID_VEG%dshd
tkd = GRID_VEG%tkd
tkmidd = GRID_VEG%tkmidd
rnetw = GRID_VEG%rnetw
xlew = GRID_VEG%xlew
hw = GRID_VEG%hw
gw = GRID_VEG%gw
dshw = GRID_VEG%dshw
tkw = GRID_VEG%tkw
tkmidw = GRID_VEG%tkmidw
epetw = GRID_VARS%epetw

!Soil Parameters
thetar = GRID_SOIL%thetar
thetas = GRID_SOIL%thetas
psic = GRID_SOIL%psic
bcbeta = GRID_SOIL%bcbeta
quartz = GRID_SOIL%quartz
ifcoarse = GRID_SOIL%ifcoarse 
rocpsoil = GRID_SOIL%rocpsoil
tcbeta = GRID_VEG%tcbeta
tcbeta_us = GRID_VEG%tcbeta_us
zdeep = GRID_SOIL%zdeep
zmid = GRID_SOIL%zmid
zrzmax = GLOBAL%zrzmax

!Moss Parameters
r_moss_depth = GRID_VEG%r_moss_depth
eps = GRID_VEG%eps
emiss_moss = GRID_VEG%emiss_moss
zpd_moss = GRID_VEG%zpd_moss
z0m_moss = GRID_VEG%z0m_moss
z0h_moss = GRID_VEG%z0h_moss

!Vegetation parameters
xlai = GRID_VEG%xlai
xlai_us = GRID_VEG%xlai_us
emiss = GRID_VEG%emiss
zpd = GRID_VEG%zpd
zpd_us = GRID_VEG%zpd_us
z0m = GRID_VEG%z0m
z0h = GRID_VEG%z0h
z0m_us = GRID_VEG%z0m_us
z0h_us = GRID_VEG%z0h_us
rescan = GRID_VEG%rescan
rescan_us = GRID_VEG%rescan_us
emiss_us = GRID_VEG%emiss_us
rsmin = GRID_VEG%rsmin
rsmax = GRID_VEG%rsmax
rsmin_us = GRID_VEG%rsmin_us
rsmax_us = GRID_VEG%rsmax_us
Rpl = GRID_VEG%Rpl
Rpl_us = GRID_VEG%Rpl_us
trefk = GRID_VEG%trefk
trefk_us = GRID_VEG%trefk_us

!COnstants

toleb = GLOBAL%toleb
maxnri = GLOBAL%maxnri
row = GRID_VARS%row
cph2o = GRID_VARS%cph2o
cp = GRID_VARS%cp
roi = GRID_VARS%roi


!Energy balance variables
rib = GRID_VARS%rib

!Water balance variables
rzsm = GRID_VARS%rzsm
tzsm = GRID_VARS%tzsm
rzsm1 = GRID_VARS%rzsm1
tzsm1 = GRID_VARS%tzsm1
r_mossm = GRID_VARS%r_mossm
zrz = GRID_VARS%zrz
smold = GRID_VARS%smold
rzdthetaudtemp = GRID_VARS%rzdthetaudtemp
smpet0 = GLOBAL%smpet0

!DIFF option parameters
iopthermc = GLOBAL%iopthermc
iopgveg = GLOBAL%iopgveg
iopthermc_v = GLOBAL%iopthermc_v
iopstab = GLOBAL%iopstab
iopsmini = GLOBAL%iopsmini

! ====================================================================
! Calculate the incoming solar radiation for the under story and over
! story layers.
! ====================================================================

      call calc_rs(GRID_VEG,i_und,i_moss,Swq_us,&
                   alb_moss,alb_snow,rsd,rs_over,rs_under)

! ====================================================================
! Initialize the liquid rain and snowfall.
! ====================================================================

      snow=0.d0
      rain=0.d0

! ====================================================================
! Initialize soil moisture if first time step.
! ====================================================================

      if ( (i.eq.1).and.(iopsmini.eq.0) ) then

         rzsm1 = smpet0
         rzsm = smpet0
         tzsm = smpet0
         tzsm1 = smpet0

      endif

! ====================================================================
! Initialize soil moisture for the calculation of the thermodynami!
! parameters, as a centered difference.
! ====================================================================

      call sm_cen_dif(iffroz,tkmid,zmid,zrzmax,smtmp,rzsm1,tzsm1,smold,&
                      rzsm,tzsm)

! ====================================================================
! Calculate the soil thermal parameters.
! ====================================================================

      call soiltherm(iopthermc,thermc1,thermc2,rzsm1,smtmp,&
       thetar,thetas,psic,bcbeta,tkmid,quartz,ifcoarse,&
       heatcap1,heatcap2,heatcapold,rocpsoil,row,cph2o,roa,cp,roi,&
       smold,thermc,heatcap,inc_frozen,rzdthetaudtemp)

! ====================================================================
! Initialize potential temperatures.
! ====================================================================

      tkd = tkact
      tkmidd = tkmid
      dshd = dshact
      tkw = tkd
      tkmidw = tkmidd
      dshw = dshd 
      tskinactd_moss = tskinact_moss
      tkactd_moss = tkact_moss
      tkmidactd_moss = tkmid_moss

! --------------------------------------------------------------------
! In case of understory of moss, initialize potential temperatures.
! --------------------------------------------------------------------

      if ( (i_und.gt.0)) then

         tkd_us = tkact_us
         tkmidd_us = tkmid_us
         dshd_us = dshact_us
         tkw_us = tkd_us
         tkmidw_us = tkmidd_us
         dshw_us = dshd_us 

         tk_p_moss = tkact_moss
         tkmid_p_moss = tkmid_moss
         ds_p_moss = dshact_moss
         tskin_p_moss = tskinact_moss

      endif

! ====================================================================
! Solve the energy balance for bare soil.
! ====================================================================
 
      if (ivgtyp.eq.0) then

         call peteb_bs(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       rain,snow,Swq,albd,emiss,ravd,rahd,tkd,tkmidd,rnetd,xled,epetd,hd,&
       gd,dshd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,rld,&
       toleb,maxnri,dt,i,tkw,tkmidw,rnetw,xlew,epetw,hw,gw,dshw,ravw,rahw,&
       PackWater,SurfWater,VaporMassFlux,TPack,TSurf,r_MeltEnergy,&
       Outflow,xleact_snow,hact_snow,rn_snow,dens,za,zpd,albw,&
       z0h,RaSnow,alb_snow,appa,vpsat,uzw,gact,row)

! ====================================================================
! Solve the energy balance for vegetated surfaces.
! ====================================================================

      else

! --------------------------------------------------------------------
! Adapt the soil thermal parameters.
! --------------------------------------------------------------------

         call soiladapt (iopgveg,thermc,iopthermc_v,tcbeta,&
                         xlai,thermc1,heatcap,heatcap1,zero)

         if(i_und.gt.0)then

            call soiladapt (iopgveg,thermc_us,iopthermc_v,tcbeta_us,&
                            xlai_us,thermc1,heatcap_us,heatcap1,zero)

         endif

! --------------------------------------------------------------------
! Calculate the environmental controls on stomatal resistance for
! both the vegetation of the under- and over story, following
! Noilhan and Planton (1989).
! --------------------------------------------------------------------

! ....................................................................
! Solar radiation.
! ....................................................................

         f1par = clcf1par(albd,xlai,rsd,rsmin,rsmax,Rpl)

         if (i_und.gt.0) then

            f1par_us = clcf1par(albd_us,xlai_us,rsd,rsmin_us,rsmax_us,Rpl_us)

         endif

! ....................................................................
! Vapor pressure deficit.
! ....................................................................

         f3vpd = clcf3vpd(tkel,vppa,f3vpdpar)

         if (i_und.gt.0) then

            f3vpd_us = clcf3vpd(tkel_ic,vppa_ic,f3vpdpar_us)

         endif

! ....................................................................
! Air temperature.
! ....................................................................

         f4temp = clcf4temp(tkel,trefk,f4temppar)

         if (i_und.gt.0) then

            f4temp_us = clcf4temp(tkel_ic,trefk_us,f4temppar_us)

         endif

! --------------------------------------------------------------------
! Calculate the potential evaporation for dry vegetation for the
! over story.
! --------------------------------------------------------------------

         call peteb_dv(thermc2,vpsat,heatcap,heatcap2,heatcapold,&
       rs_over,rain,snow,Swq,albd,emiss,thermc,f1par,f3vpd,f4temp,rescan,&
       ravd,rahd,tkd,tkmidd,rnetd,xled,epetd,hd,gd,dshd,&
       tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,&
       rld,toleb,maxnri,dt,i,PackWater,SurfWater,VaporMassFlux,TPack,&
       TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,dens,za,zpd,&
       z0h,RaSnow,appa,uzw,gact,alb_snow,row)

! --------------------------------------------------------------------
! Calculate the potential evaporation for wet vegetation for the
! over story.
! --------------------------------------------------------------------

         call peteb_wv(thermc2,heatcap,heatcap2,heatcapold,&
       rs_over,rain,snow,Swq,albw,emiss,&
       thermc,ravw,rahw,tkw,tkmidw,rnetw,xlew,epetw,hw,&
       gw,dshw,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,&
       zmid,rld,toleb,maxnri,dt,i,iopstab,tkact,zww,&
       za,uzw,zpd,z0m,tkel,press,rib,z0h,PackWater,SurfWater,&
       VaporMassFlux,TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,&
       hact_snow,rn_snow,dens,RaSnow,alb_snow,appa,vpsat,gact,row,ipix)

      endif

!General Vegetation parameters
!canclos = GRID_VEG%canclos
!extinct = GRID_VEG%extinct
GRID_VEG%i_und = i_und
GRID_VEG%i_moss = i_moss
GRID_VEG%ivgtyp = ivgtyp

!Snow Pack variables
GRID_VARS%PackWater = PackWater
GRID_VARS%SurfWater = SurfWater
GRID_VARS%Swq = Swq
GRID_VARS%VaporMassFlux = VaporMassFlux
GRID_VARS%r_MeltEnergy = r_MeltEnergy
GRID_VARS%Outflow = Outflow
GRID_VARS%PackWater_us = PackWater_us
GRID_VARS%SurfWater_us = SurfWater_us
GRID_VARS%Swq_us = Swq_us
GRID_VARS%VaporMassFlux_us = VaporMassFlux_us
GRID_VARS%r_MeltEnergy_us = r_MeltEnergy_us
GRID_VARS%Outflow_us = Outflow_us

!Albedos of the over story, under story, and moss layer
albd_us = GRID_VEG%albd_us
GRID_VEG%alb_moss = alb_moss
GRID_VEG%albd = albd
GRID_VEG%albw = albw
GRID_VEG%albw_us = albw_us
GRID_VARS%alb_snow = alb_snow

!Meteorological data
GRID_MET%rsd = rsd
GRID_MET%rld = rld
GRID_VEG%zww = zww
GRID_VEG%za = za
GRID_MET%uzw = uzw
GRID_MET%press = press

!Temperature variables
GRID_VARS%tkmid = tkmid
GRID_VARS%tkact = tkact

!Energy fluxes and states
GRID_VARS%dshact = dshact
GRID_VARS%gact = gact
GRID_VEG%rnetd = rnetd
GRID_VEG%xled = xled
GRID_VEG%hd = hd
GRID_VEG%gd = gd
GRID_VEG%dshd = dshd
GRID_VEG%tkd = tkd
GRID_VEG%tkmidd = tkmidd
GRID_VEG%rnetw = rnetw
GRID_VEG%xlew = xlew
GRID_VEG%hw = hw
GRID_VEG%gw = gw
GRID_VEG%dshw = dshw
GRID_VEG%tkw = tkw
GRID_VEG%tkmidw = tkmidw
GRID_VARS%epetw = epetw

!Soil Parameters
GRID_SOIL%thetar = thetar
GRID_SOIL%thetas = thetas
GRID_SOIL%psic = psic
GRID_SOIL%bcbeta = bcbeta
GRID_SOIL%quartz = quartz
GRID_SOIL%ifcoarse = ifcoarse
GRID_SOIL%rocpsoil = rocpsoil
GRID_VEG%tcbeta = tcbeta
GRID_VEG%tcbeta_us = tcbeta_us
GRID_SOIL%zdeep = zdeep
GRID_SOIL%zmid = zmid
GLOBAL%zrzmax = zrzmax

!Moss Parameters
GRID_VEG%r_moss_depth = r_moss_depth
GRID_VEG%eps = eps
GRID_VEG%emiss_moss = emiss_moss
GRID_VEG%zpd_moss = zpd_moss
GRID_VEG%z0m_moss = z0m_moss
GRID_VEG%z0h_moss = z0h_moss

!Vegetation parameters
GRID_VEG%xlai = xlai
GRID_VEG%xlai_us = xlai_us
GRID_VEG%emiss = emiss
GRID_VEG%zpd = zpd
GRID_VEG%zpd_us = zpd_us
GRID_VEG%z0m = z0m
GRID_VEG%z0h = z0h
GRID_VEG%z0m_us = z0m_us
GRID_VEG%z0h_us = z0h_us
GRID_VEG%rescan = rescan
GRID_VEG%rescan_us = rescan_us
GRID_VEG%emiss_us = emiss_us
GRID_VEG%rsmin = rsmin
GRID_VEG%rsmax = rsmax
GRID_VEG%rsmin_us = rsmin_us
GRID_VEG%rsmax_us = rsmax_us
GRID_VEG%Rpl = Rpl
GRID_VEG%Rpl_us = Rpl_us

GRID_VEG%trefk = trefk
GRID_VEG%trefk_us = trefk_us

!COnstants

GLOBAL%toleb = toleb
GLOBAL%maxnri = maxnri
GRID_VARS%row = row
GRID_VARS%cph2o = cph2o
GRID_VARS%cp = cp
GRID_VARS%roi = roi


!Energy balance variables
GRID_VARS%rib = rib

!Water balance variables
GRID_VARS%rzsm = rzsm
GRID_VARS%tzsm = tzsm
GRID_VARS%rzsm1 = rzsm1
GRID_VARS%tzsm1 = tzsm1
GRID_VARS%r_mossm = r_mossm
GRID_VARS%zrz = zrz
GRID_VARS%smold = smold
GRID_VARS%rzdthetaudtemp = rzdthetaudtemp
GLOBAL%smpet0 = smpet0

!DIFF option parameters
GLOBAL%iopthermc = iopthermc
GLOBAL%iopgveg = iopgveg
GLOBAL%iopthermc_v = iopthermc_v
GLOBAL%iopstab = iopstab
GLOBAL%iopsmini = iopsmini



      return

      end subroutine peteb

      function esat(Tsk)

! ====================================================================
! Compute saturation vapor pressure in Pa. given temperature
! in Kelvin.
! ====================================================================

      implicit none
      real*8 Tsk,Tsc,esat

      Tsc = Tsk - 273.15

      esat=611.d0*dexp((17.27d0*Tsc)/(237.3d0+Tsc))

      return
      end function esat

! ====================================================================
!
!                         function clcf1par
!
! ====================================================================
!
! Calculate the limiting factor due to Photosynthetically
! Active Radiation (PAR) for canopy resistance following
! Jarvis (1976), Dickenson et al. (1986), Noilhan and Planton (1989)
!
!
!          f1par =          +   f
!                   --------------------
!                       + rsmin/rsmax
!
! where
!                      =   0.55 * (1-alb)*Rg * 2
!                                         ---  ---
!                                         Rgl  LAI
!
!
! The values for rsmax, Rgl are usually from Dickenson et al. (1986).
! Units:
! rsmax [s/m]    = 5000
! Rgl   [w/m^2]  = 30 for forest and 100 for crop/grassland
!
! ====================================================================

      function  clcf1par(alb,LAI,Rg,rsmin,rsmax,Rgl)
      implicit none
      include "help/clcf1par.h"
      !type (GRID_VEG_template) :: GRID_VEG
      data a1/0.19/,a2/1128/,a3/30.8/

      par = 0.55 * (1.d0 - alb)*Rg

      f = 2.d0*par/(Rgl*LAI)

      f1par = (1.d0 + f)/(f + rsmin/rsmax)

      clcf1par = f1par

      return

      end function clcf1par

! ====================================================================
!
!                   function clcf3vpd
!
! ====================================================================
! Calculate the limiting factor due to vapor pressure
! for canopy resistance following
! Jarvis (1976), Dickenson et al. (1986), Noilhan and Planton (1989)
!
!
!         f3vpd =         1
!                 ----------------------
!                 1.d0 - g*(esat(Ts)-ea)
!
!
!
! Revision 06/25/97:  added limits to vpd of 1/g (max) and 0 (min)
!                     by Ted Endreny and Mark Zion.
!
! ====================================================================

      function  clcf3vpd(Ts,ea,g)
      implicit none
      include "help/clcf3vpd.h"

      rmax_vpd = (1.d0/g) - 0.1
      rmin_vpd = 0.0
      vpd = esat(Ts)-ea
      if (vpd.lt.rmin_vpd) vpd=rmin_vpd
      if (vpd.gt.rmax_vpd) vpd=rmax_vpd
      
      if ( (1.d0 - g*(esat(Ts)-ea)).ge.0.25d0) then

         f3vpd =  1.d0/(1.d0 - g*vpd)

      else

         f3vpd = 4.d0

      endif

      clcf3vpd = f3vpd

      return
      end function clcf3vpd

! ====================================================================
!
!                              function clcf4temp
!
! ====================================================================
!
! Calculate the limiting factor due to air temperature
! for canopy resistance following
! Jarvis (1976), Dickenson et al. (1986), Noilhan and Planton (1989)
!
!
!                   f4temp =___________1______________
!                           1.0 - f4par(tref - tair)^2
!
! where
!                    f4par is usually 0.0016
!                    tref and tpar in Kelvin or Celsius
!
! Revision 6/25/97:  added limits on air temperature to:
!                    maximum:   tref + (1/f4par)^.5 - 0.01
!                    minimum:   tref - (1/f4par)^.5 + 0.01
!                    by Ted Endreny and Mark Zion
!
! Revision 7/22/98:  add further limits to air temperature in
!                    order to reduce the maximum f4temp.  0.01
!                    is increased to 5.
!
! ====================================================================

      function clcf4temp(tair,tref,f4par)

      implicit none
      include "help/clcf4temp.h"

! ====================================================================
! Set maximum and minimum values for air temperature so
! factor does not go negative.
! ====================================================================

      rmin_tair = tref - sqrt(1/f4par) + 0.01
      rmax_tair = tref + sqrt(1/f4par) - 0.01

      rmin_tair = tref - sqrt(1/f4par) + 5.
      rmax_tair = tref + sqrt(1/f4par) - 5.

      if (tair.lt.rmin_tair) then

         tair_calc = rmin_tair

      else if (tair.gt.rmax_tair) then

         tair_calc = rmax_tair

      else

         tair_calc = tair

      endif

! ====================================================================
! Calculate temperature factor.
! ====================================================================

      f4temp =  1.d0/(1.d0 - f4par*(tref - tair_calc)**2)

      clcf4temp = f4temp

      return

      end function clcf4temp

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
!       include "SNOW.h"
      include "help/peteb_bs.h"

      if ( (Swq.le.(0.0002d0))) then

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

      end subroutine peteb_bs

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
!       include "SNOW.h"
      include "help/peteb_dv.h"

      if ( (Swq.le.(0.0002d0))) then

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

      end subroutine peteb_dv

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
!       include "SNOW.h"
      include "help/peteb_wv.h"
      integer :: ipix

      if ( (Swq.le.(0.0002d0))) then

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

      end subroutine peteb_wv

! ====================================================================
!
!               subroutine inidum
!
! ====================================================================
!
! Initialize the snow model variables that will be used in the
! snow subroutine but that not will be altered in the main program.
!
! ====================================================================

      subroutine inidum(dum1,PackWater,dum2,SurfWater,dum3,Swq,dum4,&
       VaporMassFlux,dum5,TPack,dum6,TSurf,dum7,r_MeltEnergy,dum8,Outflow,&
       dum9,xleact_snow,dum10,hact_snow,dum11,rn_snow,dum12,dens)

      implicit none
      include "help/inidum.h"

      dum1=PackWater
      dum2=SurfWater
      dum3=Swq
      dum4=VaporMassFlux
      dum5=TPack
      dum6=TSurf
      dum7=r_MeltEnergy
      dum8=Outflow
      dum9=xleact_snow
      dum10=hact_snow
      dum11=rn_snow
      dum12=dens

      return

      end subroutine inidum

END MODULE MODULE_ATMOS

