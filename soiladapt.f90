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

      end
