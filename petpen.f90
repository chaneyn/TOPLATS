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
! Calculate potential evapotranpiration rate from dry canopy or
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

      end
