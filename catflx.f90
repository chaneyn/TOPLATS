! ====================================================================
!
!			subroutine catflx
!
! ====================================================================
!
!  Calculates the catchment time step totals of evapotranspiration, 
!  runoff, surface energy fluxes and vertical soil-water fluxes.
!
! ====================================================================

      subroutine catflx(i,ic,area,pixsiz,r_lakearea,ettot,&
       etstsum,etwtsum,etlakesum,etbssum,fbs,etdcsum,&
       etwcsum,pptsum,pnetsum,contot,qsurf,sxrtot,xixtot,ranrun,&
       conrun,gwtsum,capsum,tzpsum,rzpsum,fwcat,iprn,&
       s_nr_etwtsum,s_nr_gwtsum,s_nr_capsum,s_nr_tzpsum,s_nr_rzpsum)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/catflx.h"

! ====================================================================
! Calculate the number of pixels in the current catchment.
! ====================================================================

      catpix = area/pixsiz/pixsiz
      catlakpix = r_lakearea/pixsiz/pixsiz
      catvegpix = (area-r_lakearea)/pixsiz/pixsiz

! ====================================================================
! Find catchment average evapotranspiration rates.
! ====================================================================

      ettot = ettot / catpix
      etstsum = etstsum / catvegpix
      etwtsum = etwtsum / catvegpix
      s_nr_etwtsum = s_nr_etwtsum / catvegpix
      etlakesum = etlakesum / catlakpix

      if (fbs.gt.0.) then

         etbssum = etbssum / fbs/catvegpix

      else

         etbssum = 0.

      endif

      if (fbs.lt.1.) then

         etdcsum = etdcsum / (one-fbs)/catvegpix
         etwcsum = etwcsum / (one-fbs)/catvegpix

      else

         etdcsum = 0.
         etwcsum = 0.

      endif

! ====================================================================
! Find catchment precipitation and condensation rate.
! ====================================================================

      pptsum = pptsum / catpix
      pnetsum = pnetsum / catpix
      contot = contot / catpix

! ====================================================================
! Compute total catchment runoff/infiltration rates.
! ====================================================================

      qsurf = qsurf / catvegpix
      sxrtot = sxrtot / catvegpix
      xixtot = xixtot / catvegpix

! ====================================================================
! Compute total runoff due to rainfall and due to condensation.
! ====================================================================

      ranrun = ranrun / catvegpix
      conrun = conrun / catvegpix

! ====================================================================
! Compute water table fluxes.
! ====================================================================

      gwtsum = gwtsum / catvegpix
      s_nr_gwtsum = s_nr_gwtsum / catvegpix
      capsum = capsum / catvegpix
      s_nr_capsum = s_nr_capsum / catvegpix

! ====================================================================
! Compute average available porosity above the water table.
! ====================================================================

      tzpsum = tzpsum / catvegpix 
      s_nr_tzpsum = s_nr_tzpsum / catvegpix 
      rzpsum = rzpsum / catvegpix 
      s_nr_rzpsum = s_nr_rzpsum / catvegpix 

! ====================================================================
! Calculate catchment fractions of wet canopy.
! ====================================================================

      fwcat = fwcat / catvegpix / (one-fbs)

! ====================================================================
! Write evaporation and infiltration/runoff results if required.
! ====================================================================

      if (iprn(80).eq.1)then
        write(80,1000) i,ic,ettot*3600000.,etbssum*3600000.,&
                       etdcsum*3600000.,etwcsum*3600000.,&
                       fwcat,fbs
        endif
      if (iprn(81).eq.1)then
        write(81,1100) i,ic,pptsum*3600000.,pnetsum*3600000.,&
                       contot*3600000.,&
                       (ranrun+contot)*3600000.,&
                       (pnetsum+contot-ranrun-contot)*3600000.,&
                       sxrtot*3600000.,xixtot*3600000.
        endif
! ====================================================================
! Format statements.
! ====================================================================

1000  format(2i5,4f10.5,2f7.3)
1100  format(2i5,7f10.5)

      return

      end
