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

      end
