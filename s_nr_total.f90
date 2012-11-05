! ====================================================================
!
!                       subroutine s_nr_total
!
! ====================================================================
!
! Calculate the catchment average water balance fluxes.
!
! ====================================================================

      subroutine s_nr_total(gwtsum,capsum,etwtsum,rzpsum,tzpsum,&
       s_nr_gwtsum,s_nr_capsum,s_nr_etwtsum,s_nr_rzpsum,s_nr_tzpsum)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/s_nr_total.h"

      gwtsum=gwtsum+s_nr_gwtsum
      capsum=capsum+s_nr_capsum
      etwtsum=etwtsum+s_nr_etwtsum
      rzpsum=rzpsum+s_nr_rzpsum
      tzpsum=tzpsum+s_nr_tzpsum

      return

      end
