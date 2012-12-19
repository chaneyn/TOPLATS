
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
