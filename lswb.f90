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

      subroutine lswb(i,ncatch,r_lakearea,pixsiz,npix,ettotrg,etlakesumrg,&
       etstsumrg,etwtsumrg,fbsrg,etbssumrg,etdcsumrg,etwcsumrg,pptsumrg,&
       pnetsumrg,qsurfrg,sxrtotrg,xixtotrg,contotrg,ranrunrg,conrunrg,qbreg,&
       gwtsumrg,grzsumrg,gtzsumrg,capsumrg,difrzsumrg,dswcsum,wcrhssum,&
       dsrzsum,rzrhssum,dstzsum,tzrhssum,dssum,svarhssum,rzsmav,tzsmav,&
       rnpetsum,xlepetsum,hpetsum,gpetsum,dshpetsum,tkpetsum,tkmidpetsum,&
       rnsum,xlesum,hsum,gsum,dshsum,tksum,tkmidsum,tkdeepsum,fwreg,wcip1sum,&
       zbar1rg,pr3sat,perrg2,pr2sat,pr2uns,perrg1,pr1sat,pr1rzs,pr1tzs,pr1uns,&
       persxr,perixr,persac,peruac,perusc,iprn,wcsum,zbarrg,MODE,f_lake,&
       veg_pdf,nlcs,ivgtyp,veg,Swqsum,Swq_ussum,Sdepthsum,Sdepth_ussum,qb24sum)

      use fruit
      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/lswb.h"
      real*8 rest
      real*8 :: tmp
      real*8 :: rnsum_old,xlesum_old,hsum_old,gsum_old,tksum_old,tkmidsum_old,tkdeepsum_old

! ====================================================================
! Calculate the number of pixels covered by lakes.
! ====================================================================

      nlakpix=0

      tot=0.d0

      if (MODE.eq.1) then

         do ic=1,ncatch

            tot=tot+r_lakearea(ic)/pixsiz/pixsiz

         enddo

         nlakpix=tot

      endif

      if (MODE.eq.2) then

         do ipix=1,npix

            tot=tot+f_lake(ipix)

         enddo

         nlakpix=tot

      endif

      nvegpix=npix-nlakpix

      vegd=0.d0

      if (MODE.eq.1) vegd=1.d0

      if (MODE.eq.2) then

         do ipix=1,npix

            do landc=1,nlcs

               if (ivgtyp(veg(ipix,landc)).ge.0) then

                  vegd=vegd+veg_pdf(ipix,landc)/real(npix)

               endif

            enddo

         enddo

      endif

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

!cw!      if (MODE.eq.1) then

!cw!         qbreg = qbreg / ncatch 

!cw!      endif

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

            fwreg = fwreg / real(nvegpix)/(one-fbsrg)

         else

            fwreg=0.d0

         endif

         wcip1sum = wcip1sum / real(nvegpix)

      else

         fwreg=0.d0
         wcip1sum=0.d0

      endif

      if (MODE.eq.1) then

         zbar1rg = zbar1rg / npix

      endif

      if (MODE.eq.2) then

         zbar1rg = ( zbar1rg / npix )

      endif

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
! Write results.
! ====================================================================

! --------------------------------------------------------------------
! Energy fluxes at PET.
! --------------------------------------------------------------------

      if (iprn(90).eq.1) then

         write(90,1000) i,rnpetsum,xlepetsum,hpetsum,gpetsum,&
                        rnpetsum-xlepetsum-hpetsum-gpetsum,gpetsum+dshpetsum,&
                        tkpetsum,tkmidpetsum,tkdeepsum

      endif

! --------------------------------------------------------------------
! Actual energy fluxes.
! --------------------------------------------------------------------

      if (iprn(91).eq.1) then
         write(91,1000) i,rnsum,xlesum,hsum,gsum,&
                        rnsum-xlesum-hsum-gsum,gsum+(rnsum-xlesum-hsum-gsum),&
                        tksum,tkmidsum,tkdeepsum
	 !Read the old values
	 read(2091,1000) i,rnsum_old,xlesum_old,hsum_old,gsum_old,tmp,tmp,&
			tksum_old,tkmidsum_old,tkdeepsum_old
	 !Use unit testing to compare
	 call assert_equals (nint(rnsum),nint(rnsum_old))
	 call assert_equals (nint(xlesum),nint(xlesum_old))
	 call assert_equals (nint(hsum),nint(hsum_old))
	 call assert_equals (nint(gsum),nint(gsum_old))
	 call assert_equals (nint(tksum),nint(tksum_old))
	 call assert_equals (nint(tkmidsum),nint(tkmidsum_old))
	 call assert_equals (nint(tkdeepsum),nint(tkdeepsum_old))
	
      endif

! --------------------------------------------------------------------
! Canopy water balance.
! --------------------------------------------------------------------

      if (iprn(92).eq.1)&

         write(92,1100) i,wcip1sum,wcsum,pptsumrg*3600000.,&
                        pnetsumrg*3600000.,etwcsumrg*3600000.,fwreg,&
                        dswcsum,wcrhssum,dswcsum-wcrhssum

! --------------------------------------------------------------------
! Precipitation/Runoff/Infiltration.
! --------------------------------------------------------------------

      if (iprn(93).eq.1)&

         write(93,1200) i,pptsumrg*3600000.,pnetsumrg*3600000.,&
                        contotrg*3600000.,qsurfrg*3600000.,&
                        (pnetsumrg-qsurfrg)*3600000.,&
                        sxrtotrg*3600000.,xixtotrg*3600000.

! --------------------------------------------------------------------
! Evaporation.
! --------------------------------------------------------------------

      if (iprn(94).eq.1) then

         rest=ettotrg-fbsrg*etbssumrg-etdcsumrg-etwcsumrg
         etwcsumrg=etwcsumrg+rest

         write(94,1300) i,ettotrg*3600000.,etbssumrg*3600000.,&
                        etdcsumrg*3600000.,etwcsumrg*3600000.,fwreg,fbsrg

      endif

! --------------------------------------------------------------------
! Root and Transmission Zone Balance Checks.
! --------------------------------------------------------------------

      if (iprn(95).eq.1)&

         write(95,1400) i,rzsmav,dsrzsum*1000,rzrhssum*1000,tzsmav,dstzsum*1000,&
                        tzrhssum*1000

! --------------------------------------------------------------------
! Water table balance.
! --------------------------------------------------------------------


      if (iprn(96).eq.1) then
! Output only after every 24 hour time steps  
!         if (i .eq. 1) then
!	         qb24sum = 0
!			 write(*,*) qb24sum
!			 write(*,*) 'hello'
! 	     endif
!		 if (mod(i,24) == 0) then
!		     write(96,1500) qbreg*3600.	 	     
! 			 qb24sum = qbreg
!		     write(*,*) qb24sum,qbreg
!		 else 
!	     write(*,*) qb24sum
!	     qb24sum = qb24sum + qbreg
!	     write(*,*) qb24sum,qbreg
!	     endif
!		 write(*,*) mod(i,24)
		 if (i .eq. 1) then
			 qb24sum = 0
		 endif
		 if (mod(i,8) == 0) then
!			qb24sum = qb24sum + (qbreg + qsurfrg*npix*pixsiz*pixsiz)*3600
			qb24sum = qb24sum + qbreg*3600*3
        		write(96,1500) qb24sum,zbar1rg 
			qb24sum = 0
			else
			qb24sum = qb24sum + qbreg*3600*3
!			qb24sum = qb24sum + (qbreg + qsurfrg*npix*pixsiz*pixsiz)*3600
		 endif   
!		 write(*,*) zbar1rg		
!        write(96,1500) i,qb24sum,qbreg*3600,qsurfrg*npix*pixsiz*pixsiz*3600
! 		 write(*,*) qb24sum,qbreg,qsurfrg*npix*pixsiz*pixsiz        

	  
! NWC 06/13/11
!         write(96,1500) i,zbar1rg*1000.,zbarrg*1000.,gwtsumrg*3600000.,        
!                        etwtsumrg*3600000.,qbreg/npix/pixsiz/pixsiz*3600000.,&
!                        grzsumrg*3600000.,gtzsumrg*3600000.,difrzsumrg*3600000.,&
!                        capsumrg*3600000.,pptsumrg*3600000.,pnetsumrg*3600000.,&
!                        ettotrg*3600000.
!         write(91,1000) i,rnsum,xlesum,hsum,gsum,&
!                        rnsum-xlesum-hsum-gsum,gsum+(rnsum-xlesum-hsum-gsum),&
!                        tksum,tkmidsum,tkdeepsum
	endif	
! --------------------------------------------------------------------
! Write the fractional areas in different regions.
! --------------------------------------------------------------------

      if (iprn(97).eq.1)&
         write(97,1600) i,pr3sat,perrg2,pr2sat,pr2uns,&
                        perrg1,pr1sat,pr1tzs,pr1rzs,pr1uns

! --------------------------------------------------------------------
! Write fractional runoff mechanisms and evapotranspiration control.
! --------------------------------------------------------------------

      if (iprn(98).eq.1)&
         
         write(98,1700) i,ettotrg*3600000.,persac,peruac,perusc,&
                        pnetsumrg*3600000.,persxr,perixr

! --------------------------------------------------------------------
! Write snow cover resilts.
! --------------------------------------------------------------------

      if (iprn(99).eq.1)&

         write(99,1800) i,1000.*Swq_ussum,1000.*Swqsum,1000.*(Swq_ussum+Swqsum),&
                          1000.*Sdepth_ussum,1000.*Sdepthsum,&
                          1000.*(Sdepthsum+Sdepth_ussum)

! ====================================================================
! Format statements.
! ====================================================================

1000  format(i5,9f10.3)
1100  format(i5,5f10.5,f7.3,3f10.5)
1200  format(i5,7f10.5)
1300  format(i5,4f10.5,2f7.3)
1400  format(i5,f7.3,2f10.5,f12.3,4f10.5)
1500  format(2(g11.4,1x))
1600  format(i5,2f10.3,2f7.3,f10.3,4f7.3)
1700  format(i5,f10.3,3f7.3,f13.3,2f7.3)
1800  format(i5,5(f10.4,' '),f10.4)

      return

      end
