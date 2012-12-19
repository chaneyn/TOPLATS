MODULE MODULE_ATMOS

USE VARIABLES

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

  subroutine atmos(ipix,i,dt,inc_frozen,i_2l,&

! General vegetation parameters

       canclos,extinct,i_und,i_moss,ivgtyp,&

! Snow pack variables

       PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,PackWater_us,&
       SurfWater_us,Swq_us,VaporMassFlux_us,TPack_us,&
       TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us,&
       hact_snow_us,rn_snow_us,dens,dens_us,&

! Albedos of the over story, under story,&
! and moss layer

       albd_us,alb_moss,alb_snow,albd,albw,albw_us,&

! Meteorological data

       rsd,rld,tcel,vppa,psychr,xlhv,tkel,zww,za,uzw,press,&
       appa,vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,&
       Tslope1,Tint1,Tslope2,Tint2,Tsep,Tincan,tdry,Twslope1,&
       Twint1,Twslope2,Twint2,Twsep,twet_ic,twet,&
       rh,rh_ic,qv,qv_ic,ra,ra_ic,&

! Temperature variables

       tkmid,tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,Tdeepstep,amp,phase,shift,tdeep,&
       tmid0,tmid0_moss,tk0moss,&

! Energy fluxes and states

       dshact,epetd,gact,epetd_us,dshact_moss,xle_act_moss,rnetd,xled,hd,&
       gd,dshd,tkd,tkmidd,rnetw,xlew,hw,gw,dshw,tkw,&
       tkmidw,tskinactd_moss,tkactd_moss,tkmidactd_moss,ds_p_moss,epetw,&
       dshact_us,rnetw_us,xlew_us,hw_us,gw_us,&
       dshw_us,tkw_us,tkmidw_us,epetw_us,&
       rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us,&
       tkmidd_us,rnet_pot_moss,xle_p_moss,&
       h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,&
       tskin_p_moss,eact_moss,ebspot,tsoilold,tkmidpet,tkpet,&
       tkmidpet_us,tkmidpet_moss,dspet,dspet_us,dspet_moss,rnetpn,gbspen,&

! Soil parameters

       thetar,thetas,psic,bcbeta,quartz,ifcoarse,rocpsoil,tcbeta,&
       tcbeta_us,zdeep,zmid,zrzmax,&

! Moss parameters

       r_moss_depth,eps,emiss_moss,zpd_moss,rib_moss,z0m_moss,z0h_moss,&
       epet_moss,&

! Vegetation parameters

       xlai,xlai_us,emiss,zpd,zpd_us,z0m,z0h,z0m_us,z0h_us,&
       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,rescan,&
       rescan_us,f1,f2,f3,emiss_us,rsmin,rsmax,rsmin_us,&
       rsmax_us,Rpl,Rpl_us,f3vpdpar,f3vpdpar_us,trefk,f4temppar,&
       trefk_us,f4temppar_us,&

! Constants

       row,cph2o,roa,cp,roi,toleb,maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,&
       rib,RaSnow,rib_us,ravw,ravw_us,rahw,rahw_us,&

! Water balance variables

       rzsm,tzsm,rzsm1,tzsm1,r_mossm,zrz,&
       smold,rzdthetaudtemp,smpet0,&

! Different option paramters

       iopthermc,iopgveg,iopthermc_v,iopstab,ioppet,iopwv,iopsmini)

    implicit none
    include "help/atmos.h" !take this out when variables are fixed

! ====================================================================
! Define the albedo for the snow layer.
! ====================================================================

    alb_snow = 0.75

   ! call calctempds(amp,phase,shift,Tdeepstep,tdeep) 
! ====================================================================
! Calculate the temperature of the deep soil layer.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : Assuming a constant temperature.
! --------------------------------------------------------------------

      if ( (amp.eq.(0.d0)).and.(phase.eq.(0.d0)).and.(shift.eq.(0.d0)) ) then

         Tdeepstep=tdeep

      else

! --------------------------------------------------------------------
! Option 2 : Assuming a cosine wave form.
! --------------------------------------------------------------------

         rrr=real(i)

         Tdeepstep=tdeep + amp*cos ( rrr*phase - shift )

      endif

! ====================================================================
! Calculate the temperature in the canopy given the temperature
! above the canopy.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : Assuming there is no temperature difference.
! --------------------------------------------------------------------

      if ( (Tslope1.eq.(0.d0)).and.(Tint1.eq.(0.d0)).and.&
           (Tslope2.eq.(0.d0)).and.(Tint2.eq.(0.d0)).and.&
           (Tsep.eq.(0.d0)) ) then

         Tincan=tdry

      else

! --------------------------------------------------------------------
! Option 2 : Assuming the temperature difference depends on the
! total incoming radiation.
! --------------------------------------------------------------------

         rrrr=rld+rsd

         if (rrrr.ge.Tsep) then

            Tincan=tdry+Tint2+Tslope2*rrrr

         endif

         if (rrrr.lt.Tsep) then

            Tincan=tdry+Tint1+Tslope1*rrrr

         endif

      endif

! ====================================================================
! Calculate the dew point temperature in the canopy given the dew point
! temperature above the canopy.
! ====================================================================

! --------------------------------------------------------------------
! Option 1 : Assuming there is no temperature difference.
! --------------------------------------------------------------------

      if ( (Twslope1.eq.(0.d0)).and.(Twint1.eq.(0.d0)).and.&
           (Twslope2.eq.(0.d0)).and.(Twint2.eq.(0.d0)).and.&
           (Twsep.eq.(0.d0)) ) then

         twet_ic=twet

         if (r_moss_depth.lt.0.d0) stop

      else

! --------------------------------------------------------------------
! Option 2 : Assuming the temperature difference depends on the
! total incoming radiation.
! --------------------------------------------------------------------

         rrrr=rld+rsd

         if (rrrr.ge.Twsep) then

            twet_ic=twet+Twint2+Twslope2*rrrr

         endif

         if (rrrr.lt.Twsep) then

            twet_ic=twet+Twint1+Twslope1*rrrr

         endif

      endif

! ====================================================================
! Initialize temperature variables.
! ====================================================================

      tcel=tdry
      tkel=tcel+273.15d0
      tcel_ic=Tincan
      tkel_ic=tcel_ic+273.15d0

! ====================================================================
! If first time step, use air temperature to initialize mid soil temp.
! ====================================================================

      if (i.eq.1) then

         call inittk(tdeep,tmid0,tmid0_moss,tkmid,tkmid_us,tkmid_moss,tkel,&
       tk0moss,tkact,tkact_us,tkact_moss,tskinact_moss,dshact,&
       dshact_us,dshact_moss,tkpet,tkmidpet,tkmidpet_us,tkmidpet_moss,&
       dspet,dspet_us,dspet_moss,TSurf,TPack,TSurf_us,TPack_us)

      endif

      tsoilold=tkmidpet

! ====================================================================
! Vapor pressure variables -- use different method depending
! if input includes wet bulb temperature or relative humidity
!
! If running as a two layer method READ THE HUMIDITY IN AS
! DEW POINT TEMPERATURES !
! ====================================================================

      appa=100.d0*press
      vpsat=611.d0*dexp((17.27d0*tcel)/(237.3d0+tcel))
      vpsat_ic=611.d0*dexp((17.27d0*tcel_ic)/(237.3d0+tcel_ic))

      if (iopwv.eq.0) then

         vppa=611.0d0*dexp((17.27d0*(twet))/(237.3d0+(twet)))
         rh=100.*vppa/vpsat
         vppa_ic=611.0d0*dexp((17.27d0*(twet_ic))/(237.3d0+(twet_ic)))
         rh_ic=100.*vppa_ic/vpsat_ic

      else

         vppa=0.01*rh*vpsat
         vppa_ic=0.01*rh_ic*vpsat_ic

      endif

      qv=0.622d0*(vppa/appa)
      qv_ic=0.622d0*(vppa_ic/appa)

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
      psychr=(cp*appa)/(0.622d0*xlhv)

      ra_ic=287.d0*(one+0.608d0*qv_ic)
      roa_ic=appa/(ra_ic*tkel_ic)
      xlhv_ic=2.501d6-2370.d0*tcel_ic
      psychr_ic=(cp*appa)/(0.622d0*xlhv_ic)



! ====================================================================
! Now, if requested, read in the Richardson Number for stability
! correction for aerodynamic resistance.  If no stability
! correction or first time step then set the Richardson number to zero.
!
! You can only use the stability correction if you are solving
! for a skin temperature i.e. ioppet = 0
! 
! Do this for overstory, understory and moss.
! ====================================================================

      if (iopstab.eq.1.and.i.gt.1.and.ioppet.eq.0) then 
         
         call stabcor(zww,za,uzw,zpd,z0m,tkel,press,tkact,vppa,rib)
         
      else

         rib= zero

      endif


! ====================================================================
! Calculate aerodynami! resistances to heat and mass transfer
! (assumed equal)--including stability correction if rib!=0.
!
! Again, do this for over story, under story and moss.
! Also do this for snow.
! ====================================================================

      rahd=calcra(uzw,zww,za,zpd,z0m,z0h,rib)
      rahw=calcra(uzw,zww,za,zpd,z0m,z0h,rib)

      ravd=rahd
      ravw=rahw

      RaSnow=calcra(uzw,zww,za,zpd,0.005d0,0.0005d0,1.d0)

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

      if(ioppet.eq.0)then

            call peteb(ipix,i,dt,inc_frozen,&

! General vegetation parameters

       canclos,extinct,i_und,i_moss,ivgtyp,&

! Snow pack variables

       PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,PackWater_us,&
       SurfWater_us,Swq_us,VaporMassFlux_us,TPack_us,TSurf_us,r_MeltEnergy_us,&
       Outflow_us,xleact_snow_us,hact_snow_us,rn_snow_us,dens,dens_us,&

! Albedos of the over story, under story,&
! and moss layer

       albd_us,alb_moss,alb_snow,albd,albw,albw_us,&

! Meteorological data

       rsd,rld,tcel,vppa,psychr,xlhv,tkel,zww,za,uzw,press,&
       appa,vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,&

! Temperature variables

       tkmid,tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
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

       thetar,thetas,psic,bcbeta,quartz,ifcoarse,rocpsoil,tcbeta,&
       tcbeta_us,zdeep,zmid,zrzmax,&

! Moss parameters

       r_moss_depth,eps,emiss_moss,zpd_moss,rib_moss,&
       z0m_moss,z0h_moss,epet_moss,&

! Vegetation parameters

       xlai,xlai_us,emiss,zpd,zpd_us,z0m,z0h,z0m_us,z0h_us,&
       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,rescan,&
       rescan_us,f1,f2,f3,emiss_us,rsmin,rsmax,rsmin_us,&
       rsmax_us,Rpl,Rpl_us,f3vpdpar,f3vpdpar_us,trefk,f4temppar,&
       trefk_us,f4temppar_us,&

! Constants

       row,cph2o,roa,cp,roi,toleb,maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,&
       rib,RaSnow,rib_us,ravw,ravw_us,rahw,rahw_us,&

! Water balance variables

       rzsm,tzsm,rzsm1,tzsm1,r_mossm,zrz,smold,rzdthetaudtemp,smpet0,&

! Different option paramters

       iopthermc,iopgveg,iopthermc_v,iopstab,iopsmini)

      else if(ioppet.eq.1)then

        call petpen(tcel,vpsat,vpdef,f1par,albd,&
       xlai,rsd,rsmin,rsmax,Rpl,tkel,vppa,f3vpd,f3vpdpar,f4temp,trefk,&
       f4temppar,rnetpn,gbspen,rnetd,rnetw,gd,gw,rescan,ravd,xlhv,&
       row,epetd,epetw,ravw,psychr,xled,xlew,hd,hw,cp,roa)
 
      endif

      return

  end subroutine atmos

    !subroutine calctempds(amp,phase,shift,Tdeepstep,tdeep) 

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

  subroutine inittk(tdeep,tmid0,tmid0_moss,tkmid,&
       tkmid_us,tkmid_moss,tkel,tk0moss,tkact,tkact_us,&
       tkact_moss,tskinact_moss,dshact,dshact_us,dshact_moss,tkpet,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,dspet,dspet_us,dspet_moss,&
       Tsurf,Tpack,Tsurf_us,Tpack_us)

      implicit none
      include "help/inittk.h"

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

  subroutine petpen(tcel,vpsat,vpdef,f1par,albd,xlai,rsd,rsmin,&
       rsmax,Rpl,tkel,vppa,f3vpd,f3vpdpar,f4temp,trefk,&
       f4temppar,rnetpn,gbspen,rnetd,rnetw,gd,gw,rescan,ravd,xlhv,&
       row,epetd,epetw,ravw,psychr,xled,xlew,hd,hw,cp,roa)

      implicit none
      include "help/petpen.h"

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
      
    

      return

      end subroutine petpen

END MODULE MODULE_ATMOS
