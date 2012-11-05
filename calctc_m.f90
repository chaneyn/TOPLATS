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

      end
