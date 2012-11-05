!====================================================================
!
!                   subroutine acttrans
!
! ====================================================================
!
! Set actual transpiration to the minimum of potential
! transpiration or vegetation capacity.
!
! ====================================================================

      subroutine acttrans(Swq,vegcap,epet,evtact,ievcon,zrz)

      implicit none
      include "help/acttrans.h"

      if (Swq.le.(0.d0)) then

! --------------------------------------------------------------------
! If no snow present present on top of the over story the transpiration
! is determined by the plant.
! --------------------------------------------------------------------

         if (vegcap.lt.epet) then

            evtact = vegcap
            ievcon = 1

         else

            evtact = epet

            if (zrz.gt.(0.d0)) then

               ievcon = 2

            else

               ievcon = 3

            endif

         endif

      else

! --------------------------------------------------------------------
! If there is snow the evaporation is determined by the flux out
! of the snow pack.
! --------------------------------------------------------------------

         evtact = epet
         ievcon = 3

      endif

      return

      end
