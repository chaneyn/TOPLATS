

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

      real*8 calcrib 
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
