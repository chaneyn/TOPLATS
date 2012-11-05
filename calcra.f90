! ====================================================================
!
!                         function calcra
!
! ====================================================================
!
! Calculate the aerodynammi! resistance for heat and/or vapor
! transfer following Oregon State Univ. PBL model; Businger et al.
! (1971), Holtslag and Beljaars (1989),&
   ! and Mahrt (1987) for stable conditions.
!
! zw = horizontal wind speed at height zw
! zw = height of wind measurement
! za = height of thermodynami! measurement
! zpd = zero plane displacement
! z0m = roughness length for momentum transfer
! z0h = roughness length for heat/vapor transfer
! rib = bulk Richardson number (used for stability correction)
! vk! = von Karman's constant
!
! ====================================================================

      function calcra(uzw,zw,za,zpd,z0m,z0h,rib)

      implicit none
      include "help/calcra.h"
      data vkc,astab/0.4d0,1.d0/

! --------------------------------------------------------------------
! Calculate the second power of the von Karman constant.
! --------------------------------------------------------------------

      vkctwo=vkc*vkc

! --------------------------------------------------------------------
! Check if uzw is very small or near zero, reset to small number.
! --------------------------------------------------------------------

      if (uzw.lt.0.1d0) uzw=0.1d0

! --------------------------------------------------------------------
! Assuming neutral conditions
! --------------------------------------------------------------------

      resist =(one/((vkc**two)*uzw))*&
              (dlog((za-zpd)/z0h))*&
              (dlog((zw-zpd)/z0m))

! --------------------------------------------------------------------
! Calculate stability correction based on bulk Richardson number
! --------------------------------------------------------------------

      if (rib.lt.zero) then

         stcorr = one-(15.d0*rib /&
                  (one + 7.5d0*(vkctwo)/(dlog((zw-zpd)/z0m)*&
                                         dlog((za-zpd)/z0h))*10&
                               *dsqrt(-1.d0*rib*(zw-zpd)/z0m)))

      else

         stcorr = dexp(-astab*rib)

         if (stcorr.le.0.25d0) stcorr=0.25d0

      endif

! --------------------------------------------------------------------
! Apply stability correction
! --------------------------------------------------------------------

      if (rib.ne.0.d0) then

         resist =resist/stcorr

      endif

! --------------------------------------------------------------------
! Check the bounds for the resistance.
! --------------------------------------------------------------------

      if ( (resist.ge.-100000.d0).and.(resist.le.100000.d0) ) then

         resist=resist

      else

        write (*,*) 'CALCRA : Resistance out of bounds ',resist
        write (*,*) 'A ',uzw,zw,za
        write (*,*) 'B ',zpd,z0m,z0h
        write (*,*) '! ',rib
        write (*,*) 'D ',(one/((vkc**two)*uzw)),&
                            (dlog((za-zpd)/z0h)),&
                            (dlog((zw-zpd)/z0m))
        write (*,*) 'E ',(one/((vkc**two)*uzw))*&
                         (dlog((za-zpd)/z0h))*&
                         (dlog((zw-zpd)/z0m))
        write (*,*) 'F ',stcorr
        stop

      endif

      calcra = resist

      return

      end
