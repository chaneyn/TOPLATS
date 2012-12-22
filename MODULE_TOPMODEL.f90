MODULE MODULE_TOPMODEL

USE MODULE_VARIABLES
  
!  implicit none

contains

! ====================================================================
!
!		subroutine instep
!
! ====================================================================
!
! Subroutine to initialize regional scale water balance variables.
!
! ====================================================================

  subroutine instep(i,ncatch,djday,dt,CAT_VARS,REG,CAT)

      implicit none
      include "help/instep.h"
      type (CAT_VARS_template) :: CAT_VARS
      type (REGIONAL_template) :: REG
      type (CATCHMENT_template),dimension(:),allocatable :: CAT

! ====================================================================
! Update the decimal Julian day.
! ====================================================================

      djday = djday + 0.0416666667d0*2.777777777d-4*dt

! ====================================================================
! Initialize snow water equivalent sums.
! ====================================================================

      REG%Swqsum = zero
      REG%Swq_ussum = zero
      REG%Sdepthsum = zero
      REG%Sdepth_ussum = zero

! ====================================================================
! Initialize regional state variables.
! ====================================================================

      REG%fwreg = zero
      REG%rzsmav = zero
      REG%tzsmav = zero
      REG%wcsum = REG%wcip1sum
      REG%wcip1sum = zero

! ====================================================================
! Initialize regional evaporation sums.
! ====================================================================

      REG%ettotrg = zero
      REG%etstsumrg = zero
      REG%etwtsumrg = zero
      REG%etbssumrg = zero
      REG%etdcsumrg = zero
      REG%etwcsumrg = zero
      REG%etlakesumrg = zero

! ====================================================================
! Initialize regional condensation and precipitation sums.
! ====================================================================

      REG%pptsumrg = zero
      REG%pnetsumrg = zero
      REG%contotrg = zero

! ====================================================================
! Initialize regional infiltration/runoff and baseflow sums.
! ====================================================================

      REG%sxrtotrg = zero
      REG%xixtotrg = zero
      REG%qsurfrg = zero
      REG%ranrunrg = zero
      REG%conrunrg = zero
      REG%qbreg = zero

! ====================================================================
! Initialize regional vertical moisture flux sums and water table.
! ====================================================================

      REG%capsumrg = zero
      REG%difrzsumrg = zero
      REG%gwtsumrg = zero
      REG%grzsumrg = zero
      REG%gtzsumrg = zero
      REG%zbarrg = REG%zbar1rg
      REG%zbar1rg = zero

! ====================================================================
! Initialize region water balance variables.
! ====================================================================

      REG%dswcsum = zero
      REG%dsrzsum = zero
      REG%dstzsum = zero
      REG%dssum = zero
      REG%wcrhssum = zero
      REG%rzrhssum = zero
      REG%tzrhssum = zero
      REG%svarhssum = zero

! ====================================================================
! Initialize regional actual and potential energy fluxes.
! ====================================================================

      REG%rnsum = zero
      REG%xlesum = zero
      REG%hsum = zero
      REG%gsum = zero
      REG%tksum = zero
      REG%dshsum = zero
      REG%tkmidsum = zero
      REG%tkdeepsum = zero

      REG%rnpetsum = zero
      REG%xlepetsum = zero
      REG%hpetsum = zero
      REG%gpetsum = zero
      REG%tkpetsum = zero
      REG%tkmidpetsum = zero
      REG%dshpetsum = zero

! ====================================================================
! Initialize variables telling percent land coverages of various 
! surface states.
! ====================================================================

      REG%perrg1 = zero
      REG%perrg2 = zero
      REG%pr3sat = zero

      REG%pr2sat = zero
      REG%pr2uns = zero
      REG%pr1sat = zero
      REG%pr1rzs = zero
      REG%pr1tzs = zero
      REG%pr1uns = zero

      REG%persac = zero
      REG%peruac = zero
      REG%perusc = zero

      REG%persxr = zero
      REG%perixr = zero

! ====================================================================
! Initialize variables for catchment average/total values
! ====================================================================

      do kk=1,ncatch

! --------------------------------------------------------------------
! Evaporation and condensation.
! --------------------------------------------------------------------

         CAT_VARS%ettot(kk) = zero
         CAT_VARS%etstsum(kk) = zero
         CAT_VARS%etwtsum(kk) = zero
         CAT_VARS%etbssum(kk) = zero
         CAT_VARS%etdcsum(kk) = zero
         CAT_VARS%etwcsum(kk) = zero
         CAT_VARS%etlakesum(kk) = zero
         CAT_VARS%contot(kk) = zero

! --------------------------------------------------------------------
! Infiltration/runoff/precipitation.
! --------------------------------------------------------------------

         CAT_VARS%pptsum(kk) = zero
         CAT_VARS%pnetsum(kk) = zero
         CAT_VARS%sxrtot(kk) = zero
         CAT_VARS%xixtot(kk) = zero
         CAT_VARS%qsurf(kk) = zero
         CAT_VARS%ranrun(kk) = zero
         CAT_VARS%conrun(kk) = zero

! --------------------------------------------------------------------
! Vertical soil moisture fluxes and water table updating.
! --------------------------------------------------------------------

         CAT_VARS%zbar(kk) = CAT_VARS%zbar1(kk)
         CAT_VARS%capsum(kk) = zero
         CAT_VARS%gwtsum(kk) = zero
         CAT_VARS%rzpsum(kk) = zero
         CAT_VARS%tzpsum(kk) = zero

! --------------------------------------------------------------------
! State variables.
! --------------------------------------------------------------------

         CAT_VARS%fwcat(kk) = zero

! 100   continue

CAT%ettot = CAT_VARS%ettot
CAT%etwtsum = CAT_VARS%etwtsum
CAT%etbssum = CAT_VARS%etbssum
CAT%etdcsum = CAT_VARS%etdcsum
CAT%etwcsum = CAT_VARS%etwcsum
CAT%etlakesum = CAT_VARS%etlakesum
CAT%pptsum = CAT_VARS%pptsum
CAT%pnetsum = CAT_VARS%pnetsum
CAT%sxrtot = CAT_VARS%sxrtot
CAT%xixtot = CAT_VARS%xixtot
CAT%qsurf = CAT_VARS%qsurf
CAT%ranrun = CAT_VARS%ranrun
CAT%conrun = CAT_VARS%conrun
CAT%zbar = CAT_VARS%zbar
CAT%zbar1 = CAT_VARS%zbar1
CAT%capsum = CAT_VARS%capsum
CAT%gwtsum = CAT_VARS%gwtsum
CAT%rzpsum = CAT_VARS%rzpsum
CAT%tzpsum = CAT_VARS%tzpsum
CAT%fwcat = CAT_VARS%fwcat
        enddo
        return
      
      end subroutine instep
      
! ====================================================================
!
!                       subroutine catflx
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

      end subroutine catflx

END MODULE MODULE_TOPMODEL
