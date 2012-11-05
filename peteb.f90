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

       thetar,thetas,psic,bcbeta,quartz,ifcoarse,&
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

       iopthermc,iopgveg,iopthermc_v,iopstab,iopsmini)

      implicit none
      include "help/peteb.h"

! ====================================================================
! Calculate the incoming solar radiation for the under story and over
! story layers.
! ====================================================================

      call calc_rs(canclos,extinct,i_und,i_moss,Swq_us,&
                   albd_us,alb_moss,alb_snow,rsd,rs_over,rs_under)

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

      if ( (i_und.gt.0).or.&
           (i_moss.gt.0) ) then

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

         if ( (i_und.gt.0).or.&
              (i_moss.gt.0) ) then

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

! --------------------------------------------------------------------
! Calculate the potential energy balance for the moss layer.
! --------------------------------------------------------------------

         if (i_moss.gt.0) then

            call peteb_moss(thermc1,thermc2,vpsat_ic,&
       heatcap_moss,heatcap1,heatcap2,heatcap_us,rs_under,rain,snow,&
       tskin_p_moss,tk_p_moss,tkmid_p_moss,thermc_moss,Tdeepstep,&
       xle_act_moss,h_p_moss,g_p_moss,&
       xle_p_moss,ds_p_moss,tkel_ic,rav_moss,rah_moss,r_moss_depth,&
       zmid,zdeep,eps,dt,toleb,maxnri,rld,alb_moss,&
       rnet_pot_moss,vppa_ic,iopstab,emiss_moss,&
       roa_ic,psychr_ic,xlhv_ic,epet_moss,i,tskinact_moss,zww,&
       za,uzw,zpd_moss,press,z0m_moss,rib_moss,z0h_moss,&
       PackWater_us,SurfWater_us,VaporMassFlux_us,TPack_us,&
       TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us,&
       hact_snow_us,rn_snow_us,dens_us,RaSnow,alb_snow,appa,tcel_ic,gact,&
       epetw_us,row,Swq_us)

         endif

      endif

      return

      end

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
      data a1/0.19/,a2/1128/,a3/30.8/

      par = 0.55 * (1.d0 - alb)*Rg

      f = 2.d0*par/(Rgl*LAI)

      f1par = (1.d0 + f)/(f + rsmin/rsmax)

      clcf1par = f1par

      return

      end

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
      end

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
      end

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

      function  clcf4temp(tair,tref,f4par)

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

      end
