! ====================================================================
!
!			subroutine rdebal
!
! ====================================================================
!
! Subroutine to read in energy balance calculation specification.
!
! ====================================================================

      subroutine rdebal(ioppet,iopwv,iopstab,iopgveg,iopthermc,&
                        iopthermc_v,maxnri,toleb)

      implicit none
      include "help/rdebal.h"

! ====================================================================
! Read in options for energy balance and calculation specifications. 
! ====================================================================

      ioppet = 0 !Always run in full water and energy balance
      iopwv = 1 !Always read in water vapor using relative humidity
      iopstab = 1 !Always perform stability correction on aero. resis.
      read(1000,*) iopgveg
      read(1000,*) iopthermc
      read(1000,*) iopthermc_v
      read(1000,*) maxnri
      read(1000,*) toleb

      return

      end
