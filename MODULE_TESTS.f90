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
!  Write areal average flux rates, check regional scale water balance
!  and sum simulation totals.
!
! ====================================================================

      subroutine lswb(i,REG,GLOBAL,GRID)

      implicit none
      include "help/lswb.h"
      type (REGIONAL_template),intent(inout) :: REG
      type (GLOBAL_template),intent(in) :: GLOBAL
      type (REGIONAL_template) :: REG_OLD
      type (GRID_template),dimension(:),intent(in) :: GRID
      real*8 rest
      real*8 :: tmp
      integer :: isoil
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
! Calculate the number of pixels covered by vegetation
! ====================================================================

      nvegpix=npix

      vegd=1.d0

! ====================================================================
! Calculate regional fraction of wet canopy.
! ====================================================================

      REG%fwreg = sum(GRID%VARS%fw)

! ====================================================================
! Calculate the percent of land surface in various surface states.
! ====================================================================

      do ipix = 1,npix
        isoil = GRID(ipix)%SOIL%isoil
        call sursat(GRID(ipix)%VARS%zw,GRID(isoil)%SOIL%psic,GRID(ipix)%VARS%fw,&
          GLOBAL%mul_fac,REG%pr3sat,GLOBAL%zrzmax,REG%perrg2,GRID(ipix)%VARS%rzsm1,&
          GRID(isoil)%SOIL%thetas,REG%pr2sat,REG%pr2uns,REG%perrg1,GRID(ipix)%VARS%tzsm1,&
          REG%pr1sat,REG%pr1rzs,REG%pr1tzs,REG%pr1uns,GRID(ipix)%VARS%satxr,&
          REG%persxr,GRID(ipix)%VARS%xinfxr,REG%perixr,GRID(ipix)%VARS%ievcon,&
          REG%persac,REG%peruac,REG%perusc,i,ipix)
      enddo
      pr1sat = REG%pr1sat
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

! ====================================================================
! Compute regional water balance fluxes.
! ====================================================================

! --------------------------------------------------------------------
! Evapotranspiration.
! --------------------------------------------------------------------

      ettotrg = ettotrg / real(npix)

      if (nlakpix.gt.0) then

         etlakesumrg=etlakesumrg/nlakpix

      else

         etlakesumrg=0.d0

      endif

      if (nvegpix.gt.0) then

         etstsumrg = etstsumrg / real(nvegpix)
         etwtsumrg = etwtsumrg / real(nvegpix)

         if (fbsrg.gt.0.) then

            etbssumrg = etbssumrg / fbsrg/real(nvegpix)

         else

            etbssumrg = 0.

         endif

         if (npix.gt.1) then

            etdcsumrg = etdcsumrg / (1-fbsrg)/real(nvegpix)
            etwcsumrg = etwcsumrg / (1-fbsrg)/real(nvegpix)

         else

            if (fbsrg.lt.(1.d0)) then

               etdcsumrg = etdcsumrg / (1-fbsrg)/real(nvegpix)
               etwcsumrg = etwcsumrg / (1-fbsrg)/real(nvegpix)

            else

               etdcsumrg=0.d0
               etwcsumrg=0.d0

            endif

         endif

      else

         etstsumrg=0.d0
         etwtsumrg=0.d0
         etbssumrg=0.d0
         etdcsumrg=0.d0
         etwcsumrg=0.d0

      endif

! --------------------------------------------------------------------
! Precipitation, Runoff, Infiltration and condensation.
! --------------------------------------------------------------------

      pptsumrg = pptsumrg / real(npix)
      pnetsumrg = pnetsumrg / real(npix)

      if (npix.gt.0) then

         qsurfrg = qsurfrg / real(npix)
         sxrtotrg = sxrtotrg / real(npix) 
         xixtotrg = xixtotrg / real(npix)
         contotrg = contotrg / real(npix)
         ranrunrg = ranrunrg / real(npix)
         conrunrg = conrunrg / real(npix)

      else

         qsurfrg=0.d0
         sxrtotrg=0.d0
         xixtotrg=0.d0
         contotrg=0.d0
         ranrunrg=0.d0
         conrunrg=0.d0
         qbreg=0.d0

      endif

! --------------------------------------------------------------------
! Vertical soil moisture flux.
! --------------------------------------------------------------------

      if (nvegpix.gt.0) then

         gwtsumrg = gwtsumrg / real(nvegpix)
         grzsumrg = grzsumrg / real(nvegpix)
         gtzsumrg = gtzsumrg / real(nvegpix)
         capsumrg = capsumrg / real(nvegpix)
         difrzsumrg = difrzsumrg / real(nvegpix)

      else

         gwtsumrg=0.d0
         grzsumrg=0.d0
         gtzsumrg=0.d0
         capsumrg=0.d0
         difrzsumrg=0.d0

      endif

! --------------------------------------------------------------------
! Water balance checks.
! --------------------------------------------------------------------

      if (nvegpix.gt.0) then

         dswcsum = dswcsum / real(nvegpix)
         wcrhssum = wcrhssum / real(nvegpix)
         dsrzsum = dsrzsum / real(nvegpix)
         rzrhssum = rzrhssum / real(nvegpix)
         dstzsum = dstzsum / real(nvegpix)
         tzrhssum = tzrhssum / real(nvegpix)    
         dssum = dssum / real(nvegpix)
         svarhssum = svarhssum / real(nvegpix)

      else

         dswcsum=0.d0
         wcrhssum=0.d0
         dsrzsum=0.d0
         rzrhssum=0.d0
         dstzsum=0.d0
         tzrhssum=0.d0
         dssum=0.d0
         svarhssum=0.d0

      endif

! --------------------------------------------------------------------
! Average soil moistures and snow cover.
! --------------------------------------------------------------------

      if (nvegpix.gt.0) then

         Swqsum = vegd * Swqsum / real(nvegpix)
         Swq_ussum = vegd * Swq_ussum / real(nvegpix)
         Sdepthsum = vegd * Sdepthsum / real(nvegpix)
         Sdepth_ussum = vegd * Sdepth_ussum / real(nvegpix)
         rzsmav = vegd * rzsmav / real(nvegpix)
         tzsmav = vegd * tzsmav / real(nvegpix)

      else

         Swqsum=0.d0
         Swq_ussum=0.d0
         Sdepthsum=0.d0
         Sdepth_ussum=0.d0
         rzsmav=0.d0
         tzsmav=0.d0

      endif

! --------------------------------------------------------------------
! Regional average energy fluxes at PET.
! --------------------------------------------------------------------

      rnpetsum = rnpetsum / real(npix)
      xlepetsum = xlepetsum / real(npix)
      hpetsum = hpetsum / real(npix)
      gpetsum = gpetsum / real(npix)
      dshpetsum = dshpetsum / real(npix)
      tkpetsum = tkpetsum / real(npix)
      if (nvegpix.gt.0) then

         tkmidpetsum = tkmidpetsum / real(nvegpix)

      else

         tkmidpetsum = tkmidpetsum

      endif

! --------------------------------------------------------------------
! Regional average actual surface energy fluxes.
! --------------------------------------------------------------------

      rnsum = rnsum / real(npix)
      xlesum = xlesum / real(npix)
      hsum = hsum / real(npix)
      gsum = gsum / real(npix)
      dshsum = dshsum / real(npix)
      tksum = tksum / real(npix)

      if (nvegpix.gt.0) then

         tkmidsum = tkmidsum / real(nvegpix)

      else

         tkmidsum = tksum

      endif

      tkdeepsum = tkdeepsum / real(npix)

! --------------------------------------------------------------------
! Region percentages of moisture states.
! --------------------------------------------------------------------

      if (nvegpix.gt.0) then

         if (fbsrg.lt.(1.d0)) then

            REG%fwreg = REG%fwreg / real(nvegpix)/(one-fbsrg)

         else

            REG%fwreg=0.d0

         endif

         wcip1sum = wcip1sum / real(nvegpix)

      else

         REG%fwreg=0.d0
         wcip1sum=0.d0

      endif

         zbar1rg = zbar1rg / npix

! --------------------------------------------------------------------
! Find percentage of land cover saturation states.
! --------------------------------------------------------------------

      if (nvegpix.gt.0) then

         pr3sat = pr3sat / real(npix)
   
         perrg2 = perrg2 / real(npix)
         pr2sat = pr2sat / real(npix)
         pr2uns = pr2uns / real(npix)

         perrg1 = perrg1 / real(npix)
         pr1sat = pr1sat / real(npix)
         pr1rzs = pr1rzs / real(npix)
         pr1tzs = pr1tzs / real(npix)
         pr1uns = pr1uns / real(npix)

      else

         pr3sat=0.d0
         perrg2=0.d0
         pr2sat=0.d0
         pr2uns=0.d0
         perrg1=0.d0
         pr1sat=0.d0
         pr1rzs=0.d0
         pr1tzs=0.d0
         pr1uns=0.d0

      endif

! --------------------------------------------------------------------
! Fractions of land surface saturation/infiltration excess runoff
! and atmosphere/land surface controled evapotranspiration.
! --------------------------------------------------------------------

      if (npix.gt.0) then

         persxr = persxr / real(npix)
         perixr = perixr / real(npix)

         persac = persac / real(npix)
         peruac = peruac / real(npix)
         perusc = perusc / real(npix)

      else
      
         persxr=0.d0
         perixr=0.d0

         persac=0.d0
         peruac=0.d0
         perusc=0.d0

      endif

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

      REG%wcip1sum = wcip1sum
      REG%zbar1rg = zbar1rg

      return

      end subroutine lswb

! ====================================================================
!
!                       subroutine sursat
!
! ====================================================================
!
! Define land surface saturation states for the region.
!
! ====================================================================

    subroutine sursat(zw,psic,fw,mul_fac,pr3sat,zrzmax,perrg2,&
       rzsm1,thetas,pr2sat,pr2uns,perrg1,tzsm1,pr1sat,pr1rzs,pr1tzs,pr1uns,&
       satxr,persxr,xinfxr,perixr,ievcon,persac,peruac,perusc,i,ipix)

      implicit none
      include "help/sursat.h"
      integer :: i,ipix

      data tolsat / 0.0001d0 /

! ====================================================================
! Define land surface saturation states for the region:
!
! Region 3:  Water Table at surface.
! Region 2:  Water Table in root zone.
! Region 1:  Water Table below root zone.
! ====================================================================

      if ((zw-psic).le.zero) then

! --------------------------------------------------------------------&
! First find areas in region 3.
! --------------------------------------------------------------------&

         pr3sat = pr3sat + one*mul_fac

      else if (((zw-psic).lt.zrzmax).and.((zw-psic).gt.zero)) then

! --------------------------------------------------------------------&
! For all pixels not in area 3 : first see if the water table is
! in the root zone and if the root zone is not saturated (region 2).
! --------------------------------------------------------------------&

         perrg2 = perrg2 + one*mul_fac


         if (rzsm1.ge.thetas-tolsat) then

            pr2sat = pr2sat + one*mul_fac

         else

            pr2uns = pr2uns + one*mul_fac

         endif

      else

! --------------------------------------------------------------------&
! If a pixel is not in in region 3 or 2 it has to be in region 1.
! Ssplit into four possibilities for root and transmission zone
! saturation:
! --------------------------------------------------------------------&

         perrg1 = perrg1 + one

         if ((rzsm1.ge.thetas-tolsat).and.(tzsm1.ge.thetas-tolsat)) then

! ....................................................................
! 1) Boot root and transmsission zone are saturated.
! ....................................................................

            pr1sat = pr1sat + one*mul_fac

         else if ((rzsm1.ge.thetas-tolsat).and.(tzsm1.lt.thetas-tolsat)) then

! ....................................................................
! 2) Root zone is saturated and transmsission zone is not
!    saturated.
! ....................................................................

            pr1rzs = pr1rzs + one*mul_fac

         else if ((rzsm1.lt.thetas-tolsat).and.(tzsm1.ge.thetas-tolsat)) then

! ....................................................................
! 3) Root zone is not saturated and transmsission zone is
!    saturated.
! ....................................................................

            pr1tzs = pr1tzs + one*mul_fac

         else

! ....................................................................
! 4) Both root and transmsission zone are not saturated.
! ....................................................................

            pr1uns = pr1uns + one*mul_fac

         endif

      endif

! ====================================================================
! Determine fractions of land surface contribtuting saturation
! or infiltration excess runoff.
! ====================================================================

      if (satxr.gt.zero) then

         persxr = persxr + one*mul_fac

      else if (xinfxr.gt.zero) then

         perixr = perixr + one*mul_fac

      endif

! ====================================================================
! Determine areal fractions of bare soil evaporation
! controls - check for atmospheri! contolled (saturated),&
! atmospheri! contolled (unsaturated) and soil controlled.
! ====================================================================

      if (ievcon.eq.3) then

         persac = persac + one*mul_fac

      else if (ievcon.eq.2) then

         peruac = peruac + one*mul_fac

      else

         perusc = perusc + one*mul_fac

      endif

      return

    end subroutine sursat
                         

END MODULE MODULE_TESTS
