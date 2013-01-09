MODULE MODULE_TESTS

USE FRUIT

USE MODULE_VARIABLES_OLD

USE MODULE_VARIABLES

contains

! ====================================================================
!
!			subroutine lswb
!
! ====================================================================
!
!  Check regional scale water balance and sum simulation totals.
!
! ====================================================================

      subroutine lswb(i,REG,GLOBAL,GRID,CAT)

      implicit none
      include "help/lswb.h"
      type (REGIONAL_template),intent(inout) :: REG
      type (GLOBAL_template),intent(in) :: GLOBAL
      type (REGIONAL_template) :: REG_OLD
      type (GRID_template),dimension(:),intent(in) :: GRID
      type (CATCHMENT_template),dimension(:),intent(in) :: CAT
      real*8 rest
      real*8 :: tmp
      integer :: isoil,icatch
      ncatch = GLOBAL%ncatch
      pixsiz = GLOBAL%pixsiz
      npix = GLOBAL%npix
      ettotrg = REG%ettotrg
      etlakesumrg = REG%etlakesumrg
      etstsumrg = REG%etstsumrg
      etwtsumrg = REG%etwtsumrg
      fbsrg = REG%fbsrg
      etbssumrg = REG%etbssumrg
      etdcsumrg = REG%etdcsumrg
      etwcsumrg = REG%etwcsumrg
      pptsumrg = REG%pptsumrg
      pnetsumrg = REG%pnetsumrg
      qsurfrg = REG%qsurfrg
      sxrtotrg = REG%sxrtotrg
      xixtotrg = REG%xixtotrg 
      contotrg = REG%contotrg
      ranrunrg = REG%ranrunrg
      conrunrg = REG%conrunrg
      qbreg = REG%qbreg
      gwtsumrg = REG%gwtsumrg
      grzsumrg = REG%grzsumrg
      gtzsumrg = REG%gtzsumrg
      capsumrg = REG%capsumrg
      difrzsumrg = REG%difrzsumrg
      dswcsum = REG%dswcsum
      wcrhssum = REG%wcrhssum
      dsrzsum = REG%dsrzsum
      rzrhssum = REG%rzrhssum
      dstzsum = REG%dstzsum
      tzrhssum = REG%tzrhssum
      dssum = REG%dssum
      svarhssum = REG%svarhssum
      rzsmav = REG%rzsmav
      tzsmav = REG%tzsmav
      rnpetsum = REG%rnpetsum
      xlepetsum = REG%xlepetsum
      hpetsum = REG%hpetsum
      gpetsum = REG%gpetsum
      dshpetsum = REG%dshpetsum
      tkpetsum = REG%tkpetsum
      tkmidpetsum = REG%tkmidpetsum
      rnsum = REG%rnsum
      xlesum = REG%xlesum
      hsum = REG%hsum
      gsum = REG%gsum
      dshsum = REG%dshsum 
      tksum = REG%tksum
      tkmidsum = REG%tkmidsum
      tkdeepsum = REG%tkdeepsum
      !REG%fwreg = REG%fwreg
      wcip1sum = REG%wcip1sum
      zbar1rg = REG%zbar1rg
      pr3sat = REG%pr3sat 
      perrg2 = REG%perrg2
      pr2sat = REG%pr2sat
      pr2uns = REG%pr2uns
      perrg1 = REG%perrg1
      pr1sat = REG%pr1sat
      pr1rzs = REG%pr1rzs
      pr1tzs = REG%pr1tzs
      pr1uns = REG%pr1uns
      persxr = REG%persxr
      perixr = REG%perixr
      persac = REG%persac
      peruac = REG%peruac
      perusc = REG%perusc
      wcsum = REG%wcsum
      zbarrg = REG%zbarrg
      ivgtyp = GRID%VEG%ivgtyp
      Swqsum = REG%Swqsum
      Swq_ussum = REG%Swq_ussum
      Sdepthsum = REG%Sdepthsum
      Sdepth_ussum = REG%Sdepth_ussum

! ====================================================================
! Validate Results
! ====================================================================

! --------------------------------------------------------------------
! Actual energy fluxes.
! --------------------------------------------------------------------

  read(2091,*) i,REG_OLD%rnsum,REG_OLD%xlesum,REG_OLD%hsum,REG_OLD%gsum,tmp,tmp,&
    REG_OLD%tksum,REG_OLD%tkmidsum,REG_OLD%tkdeepsum

  call set_unit_name ('lswb.f90: rnsum')
  call assert_equals (rnsum,REG_OLD%rnsum)
  call set_unit_name ('lswb.f90: xlesum')
  call assert_equals (xlesum,REG_OLD%xlesum)
  call set_unit_name ('lswb.f90: hsum')
  call assert_equals (hsum,REG_OLD%hsum)
  call set_unit_name ('lswb.f90: gsum')
  call assert_equals (gsum,REG_OLD%gsum)
  call set_unit_name ('lswb.f90: tksum')
  call assert_equals (tksum,REG_OLD%tksum)
  call set_unit_name ('lswb.f90: tkmidsum')
  call assert_equals (tkmidsum,REG_OLD%tkmidsum)
  call set_unit_name ('lswb.f90: tkdeepsum')
  call assert_equals (tkdeepsum,REG_OLD%tkdeepsum)

! --------------------------------------------------------------------
! Canopy water balance.
! --------------------------------------------------------------------

  read(2092,*) i,REG_OLD%wcip1sum,REG_OLD%wcsum,REG_OLD%pptsumrg,REG_OLD%pnetsumrg,&
                  REG_OLD%etwcsumrg,REG_OLD%fwreg,REG_OLD%dswcsum,REG_OLD%wcrhssum,tmp

  call set_unit_name ('lswb.f90: wcip1sum')
  call assert_equals (wcip1sum,REG_OLD%wcip1sum)
  call set_unit_name ('lswb.f90: wcsum')
  call assert_equals (nint(wcsum),nint(REG_OLD%wcsum))
  call set_unit_name ('lswb.f90: pptsumrg')
  call assert_equals (pptsumrg*3600000,REG_OLD%pptsumrg)
  call set_unit_name ('lswb.f90: pnetsumrg')
  call assert_equals (pnetsumrg*3600000,REG_OLD%pnetsumrg)
  call set_unit_name ('lswb.f90: etwcsumrg')
  call assert_equals (etwcsumrg*3600000,REG_OLD%etwcsumrg)
  call set_unit_name ('lswb.f90: REG%fwreg')
  call assert_equals (REG%fwreg,REG_OLD%fwreg)
  call set_unit_name ('lswb.f90: dswcsum')
  call assert_equals (dswcsum,REG_OLD%dswcsum)
  call set_unit_name ('lswb.f90: wcrhssum')
  call assert_equals (wcrhssum,REG_OLD%wcrhssum)

! --------------------------------------------------------------------
! Precipitation/Runoff/Infiltration.
! --------------------------------------------------------------------

  read(2093,*) i,REG_OLD%pptsumrg,REG_OLD%pnetsumrg,REG_OLD%contotrg,REG_OLD%qsurfrg,&
                  tmp,REG_OLD%sxrtotrg,REG_OLD%xixtotrg

  call set_unit_name ('lswb.f90: pptsumrg')
  call assert_equals (pptsumrg*3600000,REG_OLD%pptsumrg)
  call set_unit_name ('lswb.f90: pnetsumrg')
  call assert_equals (pnetsumrg*3600000,REG_OLD%pnetsumrg)
  call set_unit_name ('lswb.f90: contotrg')
  call assert_equals (contotrg*3600000,REG_OLD%contotrg)
  call set_unit_name ('lswb.f90: qsurfrg')
  call assert_equals (qsurfrg*3600000,REG_OLD%qsurfrg)
  call set_unit_name ('lswb.f90: sxrtotrg')
  call assert_equals (sxrtotrg*3600000,REG_OLD%sxrtotrg)
  call set_unit_name ('lswb.f90: xixtotrg')
  call assert_equals (xixtotrg*3600000,REG_OLD%xixtotrg)

! --------------------------------------------------------------------
! Evaporation.
! --------------------------------------------------------------------

  read(2094,*) i,REG_OLD%ettotrg,REG_OLD%etbssumrg,REG_OLD%etdcsumrg,REG_OLD%etwcsumrg,&
               REG_OLD%fwreg,REG_OLD%fbsrg

  rest=ettotrg-fbsrg*etbssumrg-etdcsumrg-etwcsumrg
  etwcsumrg=etwcsumrg+rest
  call set_unit_name ('lswb.f90: ettotrg')
  call assert_equals (ettotrg*3600000,REG_OLD%ettotrg)
  call set_unit_name ('lswb.f90: etbssumrg')
  call assert_equals (etbssumrg*3600000,REG_OLD%etbssumrg)
  call set_unit_name ('lswb.f90: etdcsumrg')
  call assert_equals (etdcsumrg*3600000,REG_OLD%etdcsumrg)
  call set_unit_name ('lswb.f90: etwcsumrg')
  call assert_equals (etwcsumrg*3600000,REG_OLD%etwcsumrg)
  call set_unit_name ('lswb.f90: REG%fwreg')
  call assert_equals (REG%fwreg,REG_OLD%fwreg)
  call set_unit_name ('lswb.f90: fbsrg')
  call assert_equals (fbsrg,REG_OLD%fbsrg)

! --------------------------------------------------------------------
! Root and Transmission Zone Balance Checks.
! --------------------------------------------------------------------

  read(2095,*) i,REG_OLD%rzsmav,REG_OLD%dsrzsum,REG_OLD%rzrhssum,REG_OLD%tzsmav,REG_OLD%dstzsum,&
               REG_OLD%tzrhssum

  call set_unit_name ('lswb.f90: rzsmav')
  call assert_equals (rzsmav,REG_OLD%rzsmav)
  call set_unit_name ('lswb.f90: dsrzsum')
  call assert_equals (dsrzsum*1000,REG_OLD%dsrzsum)
  call set_unit_name ('lswb.f90: rzrhssum')
  call assert_equals (rzrhssum*1000,REG_OLD%rzrhssum)
  call set_unit_name ('lswb.f90: tzsmav')
  call assert_equals (tzsmav,REG_OLD%tzsmav)
  call set_unit_name ('lswb.f90: dstzsum')
  call assert_equals (dstzsum*1000,REG_OLD%dstzsum)
  call set_unit_name ('lswb.f90: tzrhssum')
  call assert_equals (tzrhssum*1000,REG_OLD%tzrhssum)


! --------------------------------------------------------------------
! Water table balance.
! --------------------------------------------------------------------

  read(2096,*) i,REG_OLD%zbar1rg,REG_OLD%zbarrg,REG_OLD%gwtsumrg,&
                 REG_OLD%etwtsumrg,REG_OLD%qbreg,&
                 REG_OLD%grzsumrg,REG_OLD%gtzsumrg,REG_OLD%difrzsumrg,&
                 REG_OLD%capsumrg,REG_OLD%pptsumrg,REG_OLD%pnetsumrg,&
                 REG_OLD%ettotrg

  call set_unit_name ('lswb.f90: zbar1rg')
  call assert_equals (zbar1rg*1000,REG_OLD%zbar1rg)
  call set_unit_name ('lswb.f90: zbarrg')
  call assert_equals (nint(zbarrg*1000),nint(REG_OLD%zbarrg))
  call set_unit_name ('lswb.f90: gwtsumrg')
  call assert_equals (gwtsumrg*3600000,REG_OLD%gwtsumrg)
  call set_unit_name ('lswb.f90: etwtsumrg')
  call assert_equals (etwtsumrg*3600000,REG_OLD%etwtsumrg)
  call set_unit_name ('lswb.f90: qbreg')
  call assert_equals (qbreg/npix/pixsiz/pixsiz*3600000,REG_OLD%qbreg)
  call set_unit_name ('lswb.f90: grzsumrg')
  call assert_equals (grzsumrg*3600000,REG_OLD%grzsumrg)
  call set_unit_name ('lswb.f90: gtzsumrg')
  call assert_equals (gtzsumrg*3600000,REG_OLD%gtzsumrg)
  call set_unit_name ('lswb.f90: difrzsumrg')
  call assert_equals (difrzsumrg*3600000,REG_OLD%difrzsumrg)
  call set_unit_name ('lswb.f90: capsumrg')
  call assert_equals (capsumrg*3600000,REG_OLD%capsumrg)
  call set_unit_name ('lswb.f90: pptsumrg')
  call assert_equals (pptsumrg*3600000,REG_OLD%pptsumrg)
  call set_unit_name ('lswb.f90: pnetsumrg')
  call assert_equals (pnetsumrg*3600000,REG_OLD%pnetsumrg)
  call set_unit_name ('lswb.f90: ettotrg')
  call assert_equals (ettotrg*3600000,REG_OLD%ettotrg)

! --------------------------------------------------------------------
! Write the fractional areas in different regions.
! --------------------------------------------------------------------

  read(2097,*) i,REG_OLD%pr3sat,REG_OLD%perrg2,REG_OLD%pr2sat,&
                 REG_OLD%pr2uns,REG_OLD%perrg1,REG_OLD%pr1sat,&
                 REG_OLD%pr1tzs,REG_OLD%pr1rzs,REG_OLD%pr1uns

  call set_unit_name ('lswb.f90: pr3sat')
  call assert_equals (pr3sat,REG_OLD%pr3sat)
  call set_unit_name ('lswb.f90: perrg2')
  call assert_equals (perrg2,REG_OLD%perrg2)
  call set_unit_name ('lswb.f90: pr2sat')
  call assert_equals (pr2sat,REG_OLD%pr2sat)
  call set_unit_name ('lswb.f90: pr2uns')
  call assert_equals (pr2uns,REG_OLD%pr2uns)
  call set_unit_name ('lswb.f90: perrg1')
  call assert_equals (perrg1,REG_OLD%perrg1)
  call set_unit_name ('lswb.f90: pr1sat')
  call assert_equals (pr1sat,REG_OLD%pr1sat)
  call set_unit_name ('lswb.f90: pr1tzs')
  call assert_equals (pr1tzs,REG_OLD%pr1tzs)
  call set_unit_name ('lswb.f90: pr1rzs')
  call assert_equals (pr1rzs,REG_OLD%pr1rzs)
  call set_unit_name ('lswb.f90: pr1uns')
  call assert_equals (pr1uns,REG_OLD%pr1uns)

! --------------------------------------------------------------------
! Write fractional runoff mechanisms and evapotranspiration control.
! --------------------------------------------------------------------

  read(2098,*) i,REG_OLD%ettotrg,REG_OLD%persac,REG_OLD%peruac,&
                 REG_OLD%perusc,REG_OLD%pnetsumrg,REG_OLD%persxr,&
                 REG_OLD%perixr

  call set_unit_name ('lswb.f90: ettotrg')
  call assert_equals (ettotrg*3600000,REG_OLD%ettotrg)
  call set_unit_name ('lswb.f90: persac')
  call assert_equals (persac,REG_OLD%persac)
  call set_unit_name ('lswb.f90: peruac')
  call assert_equals (peruac,REG_OLD%peruac)
  call set_unit_name ('lswb.f90: perusc')
  call assert_equals (perusc,REG_OLD%perusc)
  call set_unit_name ('lswb.f90: pnetsumrg')
  call assert_equals (pnetsumrg*3600000,REG_OLD%pnetsumrg)
  call set_unit_name ('lswb.f90: persxr')
  call assert_equals (persxr,REG_OLD%persxr)
  call set_unit_name ('lswb.f90: perixr')
  call assert_equals (perixr,REG_OLD%perixr)

! --------------------------------------------------------------------
! Write snow cover resilts.
! --------------------------------------------------------------------

  read(2099,*) i,REG_OLD%Swq_ussum,REG_OLD%Swqsum,tmp,&
                 REG_OLD%Sdepth_ussum,REG_OLD%Sdepthsum,tmp

  call set_unit_name ('lswb.f90: Swq_ussum')
  call assert_equals (Swq_ussum*1000,REG_OLD%Swq_ussum)
  call set_unit_name ('lswb.f90: Swqsum')
  call assert_equals (Swqsum*1000,REG_OLD%Swqsum)
  call set_unit_name ('lswb.f90: Sdepth_ussum')
  call assert_equals (Sdepth_ussum*1000,REG_OLD%Sdepth_ussum)
  call set_unit_name ('lswb.f90: Sdepthsum')
  call assert_equals (Sdepthsum*1000,REG_OLD%Sdepthsum)

      return

      end subroutine lswb


END MODULE MODULE_TESTS
