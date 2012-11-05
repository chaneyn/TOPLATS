! ====================================================================
!
!			subroutine rdpet
!
! ====================================================================
!
! Subroutine to read in potential evapotranspiration (m/s) from a file.
! It also sets evaporation from bare soil equal to evaporation
! from wet canopy equal to unstressed transpiration .
!
! ====================================================================

      subroutine rdpet(epetd,epetw,ebspot,xled,xlew,xlhv,row)

      implicit none
      include "help/rdpet.h"

! ====================================================================
! Set potential evapotranspiration for wet and dry canopy.and.&
! bare soil to the input value.
! ====================================================================
      
      epetd = ebspot/(3600*1000)
      epetw = ebspot/(3600*1000)

! ====================================================================
! Calculate mass flux.
! ====================================================================

      epetdmf = epetd*row
      epetwmf = epetw*row

! ====================================================================
! Calculate latent heat flux.
! ====================================================================

      xled = epetdmf*xlhv
      xlew = epetwmf*xlhv

      return

      end





