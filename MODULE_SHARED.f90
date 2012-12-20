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

END MODULE MODULE_SHARED
