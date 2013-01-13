MODULE MODULE_CATCHMENT_TEST

USE FRUIT

USE MODULE_VARIABLES

contains

! ====================================================================
!
!                       subroutine Run_Tests
!
! ====================================================================
!
!  Check regional scale water balance and sum simulation totals.
!
! ====================================================================

      subroutine Run_Tests()

      implicit none
      type (REGIONAL_template) :: REG_NEW
      type (REGIONAL_template) :: REG_OLD
      integer :: i

! ====================================================================
! Read in the old results
! ====================================================================

  read(2001,*) i,REG_OLD

! ====================================================================
! Read in the new results
! ====================================================================

  read(2000,*) i,REG_NEW

! --------------------------------------------------------------------
! Actual energy fluxes.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: rnsum')
  call assert_equals (REG_OLD%rnsum,REG_NEW%rnsum)
  call set_unit_name ('lswb.f90: xlesum')
  call assert_equals (REG_OLD%xlesum,REG_NEW%xlesum)
  call set_unit_name ('lswb.f90: hsum')
  call assert_equals (REG_OLD%hsum,REG_NEW%hsum)
  call set_unit_name ('lswb.f90: gsum')
  call assert_equals (REG_OLD%gsum,REG_NEW%gsum)
  call set_unit_name ('lswb.f90: tksum')
  call assert_equals (REG_OLD%tksum,REG_NEW%tksum)
  call set_unit_name ('lswb.f90: tkmidsum')
  call assert_equals (REG_OLD%tkmidsum,REG_NEW%tkmidsum)
  call set_unit_name ('lswb.f90: tkdeepsum')
  call assert_equals (REG_OLD%tkdeepsum,REG_NEW%tkdeepsum)

! --------------------------------------------------------------------
! Canopy water balance.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: wcip1sum')
  call assert_equals (REG_OLD%wcip1sum,REG_NEW%wcip1sum)
  call set_unit_name ('lswb.f90: wcsum')
  call assert_equals (nint(REG_OLD%wcsum),nint(REG_NEW%wcsum))
  call set_unit_name ('lswb.f90: dswcsum')
  call assert_equals (REG_OLD%dswcsum,REG_NEW%dswcsum)
  call set_unit_name ('lswb.f90: wcrhssum')
  call assert_equals (REG_OLD%wcrhssum,REG_NEW%wcrhssum)

! --------------------------------------------------------------------
! Precipitation/Runoff/Infiltration.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: pptsumrg')
  call assert_equals (REG_OLD%pptsumrg,REG_NEW%pptsumrg)
  call set_unit_name ('lswb.f90: pnetsumrg')
  call assert_equals (REG_OLD%pnetsumrg,REG_NEW%pnetsumrg)
  call set_unit_name ('lswb.f90: contotrg')
  call assert_equals (REG_OLD%contotrg,REG_NEW%contotrg)
  call set_unit_name ('lswb.f90: qsurfrg')
  call assert_equals (REG_OLD%qsurfrg,REG_NEW%qsurfrg)
  call set_unit_name ('lswb.f90: sxrtotrg')
  call assert_equals (REG_OLD%sxrtotrg,REG_NEW%sxrtotrg)
  call set_unit_name ('lswb.f90: xixtotrg')
  call assert_equals (REG_OLD%xixtotrg,REG_NEW%xixtotrg)

! --------------------------------------------------------------------
! Evaporation.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: ettotrg')
  call assert_equals (REG_OLD%ettotrg,REG_NEW%ettotrg)
  call set_unit_name ('lswb.f90: etbssumrg')
  call assert_equals (REG_OLD%etbssumrg,REG_NEW%etbssumrg)
  call set_unit_name ('lswb.f90: etdcsumrg')
  call assert_equals (REG_OLD%etdcsumrg,REG_NEW%etdcsumrg)
  call set_unit_name ('lswb.f90: etwcsumrg')
  call assert_equals (REG_OLD%etwcsumrg,REG_NEW%etwcsumrg)
  call set_unit_name ('lswb.f90: REG_OLD%fwreg')
  call assert_equals (REG_OLD%fwreg,REG_NEW%fwreg)
  call set_unit_name ('lswb.f90: fbsrg')
  call assert_equals (REG_OLD%fbsrg,REG_NEW%fbsrg)

! --------------------------------------------------------------------
! Root and Transmission Zone Balance Checks.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: rzsmav')
  call assert_equals (REG_OLD%rzsmav,REG_NEW%rzsmav)
  call set_unit_name ('lswb.f90: dsrzsum')
  call assert_equals (REG_OLD%dsrzsum,REG_NEW%dsrzsum)
  call set_unit_name ('lswb.f90: rzrhssum')
  call assert_equals (REG_OLD%rzrhssum,REG_NEW%rzrhssum)
  call set_unit_name ('lswb.f90: tzsmav')
  call assert_equals (REG_OLD%tzsmav,REG_NEW%tzsmav)
  call set_unit_name ('lswb.f90: dstzsum')
  call assert_equals (REG_OLD%dstzsum,REG_NEW%dstzsum)
  call set_unit_name ('lswb.f90: tzrhssum')
  call assert_equals (REG_OLD%tzrhssum,REG_NEW%tzrhssum)

! --------------------------------------------------------------------
! Water table balance.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: zbar1rg')
  call assert_equals (REG_OLD%zbar1rg,REG_NEW%zbar1rg)
  call set_unit_name ('lswb.f90: zbarrg')
  call assert_equals (nint(REG_OLD%zbarrg),nint(REG_NEW%zbarrg))
  call set_unit_name ('lswb.f90: gwtsumrg')
  call assert_equals (REG_OLD%gwtsumrg,REG_NEW%gwtsumrg)
  call set_unit_name ('lswb.f90: etwtsumrg')
  call assert_equals (REG_OLD%etwtsumrg,REG_NEW%etwtsumrg)
  call set_unit_name ('lswb.f90: qbreg')
  call assert_equals (REG_OLD%qbreg,REG_NEW%qbreg)
  call set_unit_name ('lswb.f90: grzsumrg')
  call assert_equals (REG_OLD%grzsumrg,REG_NEW%grzsumrg)
  call set_unit_name ('lswb.f90: gtzsumrg')
  call assert_equals (REG_OLD%gtzsumrg,REG_NEW%gtzsumrg)
  call set_unit_name ('lswb.f90: difrzsumrg')
  call assert_equals (REG_OLD%difrzsumrg,REG_NEW%difrzsumrg)
  call set_unit_name ('lswb.f90: capsumrg')
  call assert_equals (REG_OLD%capsumrg,REG_NEW%capsumrg)

! --------------------------------------------------------------------
! Write the fractional areas in different regions.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: pr3sat')
  call assert_equals (REG_OLD%pr3sat,REG_NEW%pr3sat)
  call set_unit_name ('lswb.f90: perrg2')
  call assert_equals (REG_OLD%perrg2,REG_NEW%perrg2)
  call set_unit_name ('lswb.f90: pr2sat')
  call assert_equals (REG_OLD%pr2sat,REG_NEW%pr2sat)
  call set_unit_name ('lswb.f90: pr2uns')
  call assert_equals (REG_OLD%pr2uns,REG_NEW%pr2uns)
  call set_unit_name ('lswb.f90: perrg1')
  call assert_equals (REG_OLD%perrg1,REG_NEW%perrg1)
  call set_unit_name ('lswb.f90: pr1sat')
  call assert_equals (REG_OLD%pr1sat,REG_NEW%pr1sat)
  call set_unit_name ('lswb.f90: pr1tzs')
  call assert_equals (REG_OLD%pr1tzs,REG_NEW%pr1tzs)
  call set_unit_name ('lswb.f90: pr1rzs')
  call assert_equals (REG_OLD%pr1rzs,REG_NEW%pr1rzs)
  call set_unit_name ('lswb.f90: pr1uns')
  call assert_equals (REG_OLD%pr1uns,REG_NEW%pr1uns)

! --------------------------------------------------------------------
! Write fractional runoff mechanisms and evapotranspiration control.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: ettotrg')
  call assert_equals (REG_OLD%ettotrg,REG_NEW%ettotrg)
  call set_unit_name ('lswb.f90: persac')
  call assert_equals (REG_OLD%persac,REG_NEW%persac)
  call set_unit_name ('lswb.f90: peruac')
  call assert_equals (REG_OLD%peruac,REG_NEW%peruac)
  call set_unit_name ('lswb.f90: perusc')
  call assert_equals (REG_OLD%perusc,REG_NEW%perusc)
  call set_unit_name ('lswb.f90: pnetsumrg')
  call assert_equals (REG_OLD%pnetsumrg,REG_NEW%pnetsumrg)
  call set_unit_name ('lswb.f90: persxr')
  call assert_equals (REG_OLD%persxr,REG_NEW%persxr)
  call set_unit_name ('lswb.f90: perixr')
  call assert_equals (REG_OLD%perixr,REG_NEW%perixr)

! --------------------------------------------------------------------
! Write snow cover resilts.
! --------------------------------------------------------------------

  call set_unit_name ('lswb.f90: Swq_ussum')
  call assert_equals (REG_OLD%Swq_ussum,REG_NEW%Swq_ussum)
  call set_unit_name ('lswb.f90: Swqsum')
  call assert_equals (REG_OLD%Swqsum,REG_NEW%Swqsum)
  call set_unit_name ('lswb.f90: Sdepth_ussum')
  call assert_equals (REG_OLD%Sdepth_ussum,REG_NEW%Sdepth_ussum)
  call set_unit_name ('lswb.f90: Sdepthsum')
  call assert_equals (REG_OLD%Sdepthsum,REG_NEW%Sdepthsum)

      return

  end subroutine Run_Tests

END MODULE MODULE_CATCHMENT_TEST

PROGRAM CATCHMENT_TESTS_DRIVER

USE MODULE_CATCHMENT_TEST

implicit none
integer :: t,nt
character(len=200) :: filename(2)
nt = 250

!Open files

!Regional Variables Output
call get_command_argument(1,filename(1))
open(2000,file=trim(filename(1)))

!Old Regional Variables Output
call get_command_argument(2,filename(2))
open(2001,file=trim(filename(2)))

!####################################################################
! Initialize unit testing
!####################################################################

call init_fruit

!Compare the regional output of the new model to the old model
do t = 1,nt
  call Run_Tests()
enddo

!####################################################################
! Finalize unit testing and print summary
!####################################################################

call fruit_summary
call fruit_finalize

!Close files
close(2000)
close(2001)

END PROGRAM CATCHMENT_TESTS_DRIVER

