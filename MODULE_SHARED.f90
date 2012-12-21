MODULE MODULE_SHARED

contains


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
!                   function calch!
!
! ====================================================================
!
! Calculate the heat capacity of the soil using the input soil
! moisture requested (theta), and saturated soil moisture or
! porosity (thetasat).  Heat capacity found as density-weighted
! combination of soil, water and air heat capacities.
! From: McCumber and Pielke (JGR vol 86, no c10,&
! pp.9929-9938, 1981
!
! ====================================================================

      function calchc(theta,thetasat,rocp,row,cph2o,roa,cp,&
                      temperature,roi,inc_frozen,rzdthetaudtemp)

      implicit none
      include "help/calchc.h"

! --------------------------------------------------------------------
! Check the soil moisture input.
! --------------------------------------------------------------------

      if ( (theta.le.1.d0).and.(theta.ge.0.d0) ) then

         theta=theta

      else

         write (*,*) 'Soil moisture out of bounds in CALCHC'
         write (*,*) theta,thetasat,rocp,row,cph2o,roa,cp,&
                     temperature,roi,inc_frozen,rzdthetaudtemp
         stop

      endif

! --------------------------------------------------------------------
! Soil contribution.
! --------------------------------------------------------------------

      heatcap = (1. - thetasat) * rocp

! --------------------------------------------------------------------
! Water contribution.
! --------------------------------------------------------------------

      heatcap = heatcap + theta * row * cph2o

! --------------------------------------------------------------------
! Air contribution.
! --------------------------------------------------------------------

      heatcap = heatcap + (thetasat - theta) * roa * cp

! --------------------------------------------------------------------
! Check the size of the resulting heat capacity.
! --------------------------------------------------------------------

      if (heatcap.lt.1000000000.d0) then

         heatcap=heatcap

      else

         write (*,*) 'heatcap too high in CALCHC',heatcap
         write (*,*) theta,thetasat,rocp,row,cph2o,roa,cp,&
                     temperature,roi,inc_frozen,rzdthetaudtemp
         stop

      endif

! --------------------------------------------------------------------
! In case of frozen soil.
! --------------------------------------------------------------------

      if (inc_frozen.eq.1) then

         ttt=temperature-273.15d0

         rlf=333.2+4.995*ttt+0.02987*ttt*ttt
         rlf=rlf*row*1000.d0

         if (ttt.lt.(0.d0)) then

            if (rzdthetaudtemp.gt.0.d0) then

               heatcap=heatcap+rlf*rzdthetaudtemp

            else

               if (heatcap.gt.(-2.d0*rlf*rzdthetaudtemp)) then

                  heatcap=heatcap+rlf*rzdthetaudtemp

               else

                  heatcap=0.5d0*heatcap

               endif

            endif

         endif

      endif

! --------------------------------------------------------------------
! Check the size of the resulting heat capacity.
! --------------------------------------------------------------------

      if (heatcap.le.0.d0) then

         write (*,*) 'HEATCAP lt 0 : ',heatcap,&
                      rzdthetaudtemp,heatcap-rlf*rzdthetaudtemp,&
                      temperature,heatcap,rlf*rzdthetaudtemp
         stop

      endif

      if (heatcap.lt.100000000000000.d0) then

         heatcap=heatcap

      else

!         write (*,*) 'heatcap too high in CALCHC',heatcap
!         write (*,*) theta,thetasat,rocp,row,cph2o,roa,cp,&
!                     temperature,roi,inc_frozen,rzdthetaudtemp
         heatcap=100000000000000.d0

      endif

      calchc = heatcap

      return

      end function calchc


! ====================================================================
!
!		subroutine stabcor
!
! ====================================================================
!
! Subroutine to calculate all variables related to stability
! correction and aerodynami! resistances.
!
! ====================================================================
!
! Parameter description:
!
! zwind		Height of wind speed measurement (m).
! zhum		Height of humidity measurement (m).
! wind_s	Wind speed (m/s).
! zpdis		Zero plane displacement height (m).
! r_m		Roughness length for momentum transfer (m).
! tair		Air temperature (C).
! airp		Air pressure (mbar).
! t_skin	Skin temperature of the soil (C).
! vappres	Vapor pressure (mbar).
! richn		Richardson number (-).
! ====================================================================

      subroutine stabcor(zwind,zhum,wind_s,zpdis,r_m,tair,airp,&
                         t_skin,vappres,richn)

      implicit none
      include "help/stabcor.h"

! ====================================================================
! Calculate the air pressure in Pascals and the absolute
! humidity.
! ====================================================================

      appa=100.d0*airp
      qv=0.622d0*(vappres/appa)
      qskin=qv

! ====================================================================
! Adjust the wind speed (if necessary) for height differences
! between the height of the wind speed measurement and the
! height of the humidity measurement and calculate the Richardson
! number.
! ====================================================================

      if ( ((zwind-zhum).gt.0.01d0).and.(zwind.gt.zhum) ) then

! --------------------------------------------------------------------&
! If the height of the wind speed measurement is higher than the
! height of the humidity measurement adjust the wind speed towards
! the humidity measurement height.
! --------------------------------------------------------------------&

         uza=wind_s * log((zhum-zpdis)/r_m) / (log((zwind-zpdis)/r_m))

         richn=calcrib(tair,qv,airp,uza,zhum,t_skin,qskin)

      else

! --------------------------------------------------------------------&
! If the height of the wind speed measurement is lower than the
! height of the humidity measurement adjust the wind speed towards
! the humidity measurement height.
! --------------------------------------------------------------------&

         if ( ((zhum-zwind).gt.0.01d0).and.(zhum.gt.zwind) ) then

            uza=wind_s*log((zhum-zpdis)/r_m) / (log((zwind-zpdis)/r_m))

            richn=calcrib(tair,qv,airp,uza,zhum,t_skin,qskin)

         else 

! --------------------------------------------------------------------&
! If the height of the wind speed measurement is equal to the
! height of the humidity measurement no adjustment in wind speed has
! to be made.
! --------------------------------------------------------------------&

            richn=calcrib(tair,qv,airp,wind_s,zhum,t_skin,qskin)

         endif

      endif

      return

      end subroutine stabcor


!   ====================================================================
!
!                     function calcrib
!
!   ====================================================================
!
!   Calculate Bulk Richardson number based on input
!   temperature, humidity, pressure, wind profile
!
!   tk2   temperature (K) at 2nd level
!   q2    specifi!   humidity (kg/kg) at 2nd level
!   p2    pressure (hPa) at 2nd level
!   u2    wind speed (m/s) at 2nd level
!   z2    distance (m) between 1st and 2nd level
!   tk1   temperature (K) at 1st level
!   q1    specifi!   humidity (kg/kg) at 1st level
!
!   ====================================================================

      real*8 function calcrib(tk2,q2,p2,u2,z2,tk1,q1)

      implicit none
      include "help/calcrib.h"
      data rdcp,GRAV/-0.286,9.81d0/

!   --------------------------------------------------------------------
!   If wind is below detectable limit, set wind speed to a small number 
!   so ribtmp doesn't divide by zero.
!   --------------------------------------------------------------------

      if (u2.lt.0.1d0) u2=0.1d0 
      
      thta1 = tk1 * (p2/1.d4)**(rdcp)
      thta2 = tk2 * (p2/1.d4)**(rdcp)
      
      thta1v = thta1
      thta2v = thta2
      
      ribtmp = GRAV * z2 *(thta2v-thta1v)/(thta2v*u2*u2) 

!   --------------------------------------------------------------------
!   Check the bounds of the Richardson number.
!   --------------------------------------------------------------------
      
      if ( (ribtmp.ge.-10000.d0).and.(ribtmp.le.10000.d0) ) then

         ribtmp=ribtmp

      else

        write (*,*) 'CALCRIB : Richardson number out of bounds ',ribtmp
        write (*,*) 'A ',tk2,q2
        write (*,*) 'B ',p2,u2,z2
        write (*,*) '!   ',tk1,q1
        write (*,*) 'D ',thta1,thta2
        write (*,*) 'E ',(thta2v-thta1v),(thta2v*u2*u2)
        write (*,*) 'F ',(thta2v-thta1v)/(thta2v*u2*u2),z2,GRAV
        stop

      endif

      calcrib = ribtmp
      
      return
      
      end function calcrib

! ====================================================================
!
!                   function calctc_j
!
! ====================================================================
!
! Calculate the thermal conductivity (thermc) of
! an unsaturated soil using Johansen's method
!
! INPUT:
!
! soil moisture (theta) (% vol),
! residual saturation (thetar) (% vol),
! porosity (thetas) (% vol),
! quartz content (quartz) (%),
! frozen soil flag (iffroz) (1=frozen,0=not).
! coarse soil flag (ifcoarse) (1=coarse,0=fine).
!
! ====================================================================

      function calctc_j(theta,thetar,thetas,quartz,iffrozen,ifcoarse)

      implicit none
      include "help/calctc_j.h"
      data soildens/2700.d0/,ko/2.d0/,kq/7.7d0/

! --------------------------------------------------------------------
! Double check the soil moisture inputs.
! --------------------------------------------------------------------

      if ( (theta.ge.0.d0).and.(theta.le.1.d0) ) then

         theta=theta

      else

         write (*,*) 'Soil moisture out of bounds in CALCTC_J'
         write (*,*) theta,thetar,thetas,quartz,iffrozen,ifcoarse
         stop

      endif

! --------------------------------------------------------------------
! Calculate relative saturation.
! --------------------------------------------------------------------

      satrel=(theta-thetar)/(thetas-thetar)

! --------------------------------------------------------------------
! Calculate bulk density .
! --------------------------------------------------------------------

      bulkdens = soildens*(1.d0 - thetas)

! --------------------------------------------------------------------
! Compute dry thermal conductivity (kdry).
! --------------------------------------------------------------------

      if (ifcoarse.eq.1) then

         kdry = 0.39d0*(thetas)**(-2.2d0)

      else

         kdry = 0.137d0*bulkdens + 64.7d0
         kdry = kdry/(soildens-0.947d0*bulkdens)

      endif
	
! --------------------------------------------------------------------
! Compute Kersten number (Ke).
! --------------------------------------------------------------------

      if (iffrozen.eq.1) then

         Ke = satrel

      else

         if (ifcoarse.eq.1) then

            if (satrel.gt.0.05) then

               Ke= 0.7d0*dlog10(satrel) + 1.d0

            else

               Ke = 0.d0

            endif

         else	

            if (satrel.gt.0.1) then

               Ke= dlog10(satrel) + 1.d0

            else

               Ke = 0.d0

            endif

         endif

      endif

! --------------------------------------------------------------------
! Compute solids thermal conductivity (ks).
! --------------------------------------------------------------------

      if (ifcoarse.eq.1.and.quartz.lt.0.20) then

         ks = kq**quartz * (ko+1.)**(1.d0-quartz)

      else

         ks = kq**quartz * ko**(1.d0-quartz)

      endif

! --------------------------------------------------------------------
! Compute saturated thermal conductivity (ksat).
! --------------------------------------------------------------------

      if (iffrozen.eq.1) then

         ksat = 2.2d0**thetas * ks**(1-thetas) * 0.269d0**theta	

      else

         ksat = 0.57d0**thetas * ks**(1-thetas)

      endif

! --------------------------------------------------------------------
! Compute thermal conductivity (thermc).
! --------------------------------------------------------------------

      thermc = (ksat-kdry)*Ke + kdry

      if ( (thermc.lt.100.d0).or.(thermc.gt.0.d0) ) then

         thermc=thermc

      else

         write (*,*) 'CALCTC_J : therm! unrealisti! ',thermc
         write (*,*) theta,thetar,thetas,quartz,iffrozen,ifcoarse
         stop

      endif

      calctc_j = thermc
	
      return

      end function calctc_j

! ====================================================================
!
!                         function calctc_m
!
!
! ====================================================================
!
! Calculate the the soil water tension (psi) and thermal
! conductivity (thermc) of the soil using the input soil
! moisture requested (soilm).  Thermal conductivity found
! using McCumber and Pielke (JGR vol 86, no c10,
! pp.9929-9938, 1981 corrected to Johansen's method as
! described in Liang et al 1995 and Peters-Lidard et al
! 1995.
!
! ====================================================================

      function calctc_m(soilm,thetar,thetas,psic,bcbeta)
      implicit none
      include "help/calctc_m.h"

      satrel=(soilm-thetar)/(thetas-thetar)
      psi=psic/(satrel**(1.d0/bcbeta))
      psicm=psi*100.d0
      pf=dlog10(dabs(psicm))

      if (pf.le.5.1d0)then

         tmpthermc=dexp(-(pf+2.7d0))

      else

         tmpthermc=0.00041d0

      endif

! ====================================================================
! Convert thermal conductivity of soil from cal/(cm*s*deg)
! to w/m/deg
! ====================================================================

      thermc = tmpthermc*418.6d0

      calctc_m = thermc

      return

      end function calctc_m

END MODULE MODULE_SHARED
