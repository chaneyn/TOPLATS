! ====================================================================
!
! 			subroutine sumflx
!
! ====================================================================
!
! Calculates the time step totals of evapotranspiration, runoff, 
! surface energy fluxes and vertical soil-water fluxes.
!
! ====================================================================

      subroutine sumflx(&

! Factor to rescale all local fluxes

       rescale,&

! General vegetation parameters

       ivgtyp,i_und,i_moss,f_und,f_moss,i,iprn,canclos,ilandc,dt,&

! Canopy variables

       dc,fw,dc_us,fw_us,&

! Regional evapotranspiration variables

       etpix,evtact,ettot,ettotrg,epwms,evtact_us,epwms_us,evtact_moss,&
       etstore,etwt,etstsum,etstsumrg,etwtsumrg,etwtsum,etbssum,etbssumrg,&
       etdcsum,etdcsumrg,etwcsum,etwcsumrg,etlakesumrg,etlakesum,&

! Condensation variables

       bsdew,contot,contotrg,tair,&

! Precipitation and runoff variables

       pptsum,pptms,pptsumrg,pnetsum,pnet,pnetsumrg,qsurf,runtot,qsurfrg,&
       sxrtot,satxr,sxrtotrg,xixtot,xinfxr,xixtotrg,ranrun,ranrunrg,&
       conrun,conrunrg,wcip1sum,wcip1,&

! Checks on water balances

       dswcsum,dswc,wcrhssum,wcrhs,dsrzsum,dsrz,&
       rzrhssum,rzrhs,dstzsum,dstz,tzrhssum,tzrhs,&

! Drainage variables

       zrz,gwt,grz,gtz,ztz,gwtsum,gwtsumrg,grzsumrg,gtzsumrg,&

! Capillary rise/diffusion variables

       capsum,capsumrg,difrzsumrg,diftz,dstore,dssum,&
       svarhs,difrz,svarhssum,&

! Soil moisture variables

       inc_frozen,rzsm1_u,tzsm1_u,rzsm1,tzsm1,rzsmav,tzsmav,tzpsum,thetas,&
       rzpsum,r_mossm,rzsm,tzsm,Swqsum,Swq_ussum,Swq,Swq_us,Sdepthsum,Sdepth_ussum,Sdepth,Sdepth_us,&

! Actual energy balance components

       rnact_moss,xleact_moss,hact_moss,gact_moss,dshact_moss,tskinact_moss,&
       tkact_moss,tkmid_moss,rnact_us,xleact_us,hact_us,gact_us,&
       dshact_us,tkact_us,tkmid_us,rnact,xleact,hact,gact,&
       dshact,tkact,tkmid,rnsum,xlesum,hsum,gsum,dshsum,tksum,&
       tkmidsum,&

! Potential energy balance components

       rnpet,rnpet_us,rnpet_moss,xlepet,xlepet_us,xlepet_moss,&
       hpet,hpet_us,hpet_moss,gpet,gpet_us,gpet_moss,&
       dspet,dspet_us,dspet_moss,tkpet,tkpet_us,tkpet_moss,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rnpetsum,xlepetsum,&
       hpetsum,gpetsum,dshpetsum,tkpetsum,tkmidpetsum,tkdeepsum,&
       Tdeepstep)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      
    
      include "help/sumflx.dif.h"

! ====================================================================
! Compute regional average evapotranspiration rate for the time step.    
! ====================================================================
 
  
      if (ivgtyp.eq.(-1)) then

! --------------------------------------------------------------------&
! In case of an open water pixel.
! --------------------------------------------------------------------&

         etpixloc=evtact
         ettot = ettot + etpixloc*rescale
         ettotrg = ettotrg + etpixloc*rescale

      else

! --------------------------------------------------------------------&
! In case of vegetated or bare soil pixels.
! --------------------------------------------------------------------&

         if ( (i_und.eq.0).and.(i_moss.eq.0) ) then

! ....................................................................
! If understory/moss is not represented.
! ....................................................................

            etpixloc = evtact*dc*(one-fw) + epwms*dc
            ettot = ettot + etpixloc*rescale
            ettotrg = ettotrg + etpixloc*rescale

         else

! ....................................................................
! Add the evapotranspiration rate from the understory/moss layer.
! ....................................................................

!            etpixlo! = (1.d0-f_moss-f_und)*&
!                          (evtact*dc*(one-fw) + epwms*dc) +&
!                          f_und*&
!                          (evtact_us*dc_us*(one-fw_us)+epwms_us*dc_us)
!                            f_moss*evtact_moss
            xlhv=1000000000.*(2.501-0.002361*tair)
            dummy = canclos*xleact+ f_und*xleact_us+ f_moss*xleact_moss
            etpixloc = dummy/xlhv
            ettot = ettot + etpixloc*rescale
            ettotrg = ettotrg + etpixloc*rescale

         endif

      endif

! ====================================================================
! Split evapotranspiration by whether the water is supplied from
! below or above the water table.  Also sum evapotranspiration
! from each surface type.  Separate the evaporation from the
! water table by catchment for water table updating.
! ====================================================================

      if (ivgtyp.ge.0) then

         if (((ivgtyp.ne.2).and.(zrz.le.zero)).or.&
             ((ivgtyp.eq.2).and.(ztz.le.zero))) then

! --------------------------------------------------------------------&
! Add the evapotranspiration values up from pixels where the water
! is supplied from above the water table.
! --------------------------------------------------------------------&

            if (i_moss+i_und.eq.0) then

! ....................................................................
! Understory/moss is not represented.
! ....................................................................

               etstore = epwms*dc

            else

! ....................................................................
! Add understory/moss in the totals.
!....................................................................

               etstore = (1.-f_moss-f_und)*epwms*dc+f_und*epwms_us*dc_us

            endif

         else 

! --------------------------------------------------------------------&
! Add the evapotranspiration values up from pixels where the water
! is supplied from below the water table.
! --------------------------------------------------------------------&

            etstore = etpixloc

         endif

! ====================================================================
! Add up the regional values for evapotranspiration, either from
! the water table or from above the water table.
! ====================================================================

         etwt = etpixloc - etstore
         etpix=etpix+etpixloc*rescale
         etstsum = etstsum + etstore*rescale
         etstsumrg = etstsumrg + etstore*rescale
         etwtsum = etwtsum + etwt*rescale
         etwtsumrg = etwtsumrg + etwt*rescale

         if (ivgtyp.eq.0) then

! ====================================================================
! Add up evaporation from bare soil values.
! ====================================================================

            etbssum = etbssum + evtact*rescale
            etbssumrg = etbssumrg + evtact*rescale

         else

! ====================================================================
! Add up evaporation from dry and wet canopy for vegetated pixels.
! ====================================================================

            if ( (i_und.eq.0).and.(i_moss.eq.0) ) then

! --------------------------------------------------------------------&
! Understory/moss is not represented.
! --------------------------------------------------------------------&

               etdcsum = etdcsum + evtact*dc*(one-fw)*rescale
               etdcsumrg = etdcsumrg + evtact*dc*(one-fw)*rescale
               etwcsum = etwcsum + epwms*dc*rescale
               etwcsumrg = etwcsumrg + epwms*dc*rescale

            else

! --------------------------------------------------------------------&
! Add the values from the understory/moss in the totals.
! --------------------------------------------------------------------&

               etdcsum = etdcsum + ((1.-f_moss- f_und)*&
                                       evtact*dc*(one-fw) +&
                                       f_und*evtact_us*dc_us*(one-fw_us) +&
                                       f_moss*evtact_moss)*rescale
               etdcsumrg = etdcsumrg + ((1.-f_moss- f_und)*&
                                       evtact*dc*(one-fw) +&
                                       f_und*evtact_us*dc_us*(one-fw_us) +&
                                       f_moss*evtact_moss)*rescale
               etwcsum = etwcsum + (epwms*dc* (1.-f_moss- f_und)+&
                                       epwms_us*dc_us*f_und)*rescale
               etwcsumrg = etwcsumrg + (epwms*dc*(1.-f_moss-f_und)+&
                                        epwms_us*dc_us*f_und)*rescale

            endif

         endif

      else

         etlakesum=etlakesum+evtact*rescale
         etlakesumrg=etlakesumrg+evtact*rescale

      endif

! ====================================================================
! Compute total, catchment and regional condensation.
! ====================================================================

      if (ilandc.ge.0) then

         if (ivgtyp.eq.0) then

! --------------------------------------------------------------------&
! In case of bare soil the condensation is the dew onto the soil.
! --------------------------------------------------------------------&

            conpix = bsdew 

         else

! --------------------------------------------------------------------&
! In case of vegetated pixels the condensation is the negative
! evaporation onto the wet canopy.
! --------------------------------------------------------------------&

            conpix = epwms*(one-dc)

         endif

         if (ivgtyp.eq.(-1)) then

! --------------------------------------------------------------------&
! In case of a lake the condensation is zero.
! --------------------------------------------------------------------&

!CVAL            if (evtact.lt.(0.d0)) then

!CVAL               conpix=0.d0-evtact

!CVAL            else

!CVAL               conpix=0.d0

!CVAL            endif

            conpix=0.d0

         endif

      endif

      contot = contot + conpix*rescale
      contotrg = contotrg + conpix*rescale

! ====================================================================
! Compute total precipitation and surface runoff for the time step.
! ====================================================================

      pptsum = pptsum + pptms*rescale
      pptsumrg = pptsumrg + pptms*rescale
      pnetsum = pnetsum + pnet*rescale
      pnetsumrg = pnetsumrg + pnet*rescale
      qsurf = qsurf + runtot*rescale
      qsurfrg = qsurfrg + runtot*rescale

! ====================================================================
! Compute total saturation and infiltration excess runoff for the
! time step.
! ====================================================================

      sxrtot = sxrtot + satxr *rescale
      sxrtotrg = sxrtotrg + satxr *rescale
      xixtot = xixtot + xinfxr*rescale
      xixtotrg = xixtotrg + xinfxr*rescale

! ====================================================================
! Compute total runoff due to rainfall and due to condensation.
! ====================================================================

      if (pptms.gt.zero) then

! --------------------------------------------------------------------&
! When under rainfall runoff has to be due to rain.
! --------------------------------------------------------------------&

         ranrun = ranrun + runtot*rescale
         ranrunrg = ranrunrg + runtot*rescale

      else

! --------------------------------------------------------------------&
! When no precipitation runoff has to be due to condensation.
! --------------------------------------------------------------------&

         conrun = conrun + runtot *rescale
         conrunrg = conrunrg + runtot *rescale

      endif  

! ====================================================================
! Compute checks on canopy water balance, root zone water balance,&
! and transmission zone balance.
! ====================================================================

      if (ivgtyp.ge.(0)) then

         dswcsum = dswcsum + dswc*rescale
         wcrhssum = wcrhssum + wcrhs*rescale

         dsrzsum = dsrzsum + dsrz*rescale
         rzrhssum = rzrhssum + rzrhs*rescale

         dstzsum = dstzsum + dstz*rescale
         tzrhssum = tzrhssum + tzrhs*rescale

! ====================================================================
! Compute drainage to the water table for time step.
! ====================================================================

         if (zrz.eq.zero) then
     
! --------------------------------------------------------------------&
! If the root zone is saturated there is no drainage towards the water 
! table.
! --------------------------------------------------------------------&

            gwt = zero
            difwt = difrz

        else if (ztz.eq.zero) then

! --------------------------------------------------------------------&
! If the transmission zone is saturated and the root zone is not
! the drainage and diffusion towards the water table comes from the
! root zone.
! --------------------------------------------------------------------&

            gwt = grz
            difwt = difrz

        else

! --------------------------------------------------------------------&
! If the transmission zone is not saturated the drainage and diffusion
! towards the water table comes from the  transmission zone.
! --------------------------------------------------------------------&
   
            gwt = gtz
            difwt = diftz
         endif

! ====================================================================
! Compute the regional totals of drainage to the water table and
! root and transmission zone drainage terms.
! ====================================================================


         gwtsum = gwtsum + gwt*rescale
         gwtsumrg = gwtsumrg + gwt*rescale
         grzsumrg = grzsumrg + grz*rescale
         gtzsumrg = gtzsumrg + gtz*rescale

! ====================================================================
! Compute the diffusion totals for the time step.
! ====================================================================

         capsum = capsum - difwt*rescale
         capsumrg = capsumrg - difwt*rescale
         difrzsumrg = difrzsumrg - difrz*rescale

! ====================================================================
! Compute change in storage above the water table for the time step 
! and perform water balance.
! ====================================================================

         dstore = dswc + dsrz + dstz
         dssum = dssum + dstore

! ====================================================================
! Add up all the input terms towards the water table and make the
! regional total.
! ======================================================================


         svarhs = dt*(pptms - etstore - runtot - gwt - difwt)

         svarhssum = svarhssum + svarhs*rescale

! ====================================================================
! Compute average soil moistures for the end of the time step.
! ====================================================================

         if (inc_frozen.eq.0) then

! --------------------------------------------------------------------&
! If frozen and liquid soil water are treated as liquid water than the
! liquid water content is the total water content.
! --------------------------------------------------------------------&

            rzsm1_u=rzsm1
            tzsm1_u=tzsm1

         endif

         rzsmav = rzsmav + rzsm1*rescale
         tzsmav = tzsmav + tzsm1*rescale
         Swqsum = Swqsum + Swq*rescale
         Swq_ussum = Swq_ussum + Swq_us*rescale
         Sdepthsum = Sdepthsum + Sdepth*rescale
         Sdepth_ussum = Sdepth_ussum + Sdepth_us*rescale


! ====================================================================
! Make a regional total of canopy interception storage.
! ====================================================================

         wcip1sum = wcip1sum + wcip1*rescale

! ====================================================================
! Compute available porosity in root and transmission zones .or.&
! updating of average water table depth.  Only interested in
! the porosity for the zone above the water table.
! ====================================================================

         if (ztz.gt.zero) then

! --------------------------------------------------------------------&
! If the root and transmission zone are both unsaturated than
! the region of interest is the transmission zone.
! --------------------------------------------------------------------&

            tzpsum = tzpsum + (thetas-tzsm)*rescale

         else if (zrz.gt.zero) then

! --------------------------------------------------------------------&
! If the root zone is unsaturated and the transmission zone
! unsaturated than the region of interest is the root zone.
! --------------------------------------------------------------------&

            rzpsum = rzpsum + (thetas-rzsm)*rescale

         endif

      endif

! ====================================================================
! Write the components of the energy balance out if requested.
! ====================================================================

      if (ivgtyp.eq.(-1)) then

!CVAL         gact=0.d0
!CVAL         gpet=0.d0
!CVAL         dshact=0.d0
!CVAL         dspet=0.d0
!CVAL         tkact=0.d0
!CVAL         tkpet=0.d0
         tkmid=0.d0
         tkmidpet=0.d0
         tskinact_moss=0.d0

      endif

      if (iprn(110).eq.1) then

         write (110,125) i,rnact_moss,xleact_moss,hact_moss,&
                         gact_moss,dshact_moss,tskinact_moss,&
                         tkact_moss,tkmid_moss,r_mossm

      endif

      if (iprn(111).eq.1) then

         write (111,126) i,rnact_us,xleact_us,hact_us,&
                         gact_us,dshact_us,&
                         tkact_us,tkmid_us

      endif

      if (iprn(112).eq.1) then

         write (112,126) i,rnact,xleact,hact,gact,&
                         dshact,&
                         tkact,tkmid

      endif

! ====================================================================
! Compute pixel total energy fluxes at PET.
! ====================================================================

      if (i_moss+i_und.gt.0) then

! --------------------------------------------------------------------&
! If understory/moss is represented than add their fluxes temperatures
! in the total for the pixel.
! --------------------------------------------------------------------&

!         rnpet = (1.-f_moss-f_und)*&
!                       rnpet+&
!                       f_und*rnpet_us+&
!                       f_moss*rnpet_moss
         rnpet = canclos*rnpet+ f_und*rnpet_us+ f_moss*rnpet_moss
         xlepet = canclos*xlepet+ f_und*xlepet_us+ f_moss*xlepet_moss
         hpet = canclos*hpet+ f_und*hpet_us+ f_moss*hpet_moss
!         gpet = (1.-f_moss-f_und)*&
!                      gpet+&
!                      f_und*gpet_us+&
!                      f_moss*gpet_moss
         gpet = canclos*gpet+ f_und*gpet_us+ f_moss*gpet_moss
!         dspet = (1.-f_moss-f_und)*&
!                       dspet+&
!                       f_und*dspet_us+&
!                       f_moss*dspet_moss
         dspet = canclos*dspet+ f_und*dspet_us+ f_moss*dspet_moss
         tkpet = (1.-f_moss-f_und)* tkpet+ f_und*tkpet_us+ f_moss*tkpet_moss
         tkmidpet = (1.-f_moss-f_und)* tkmidpet+ f_und*tkmidpet_us+&
                    f_moss*tkmidpet_moss

      endif

! ====================================================================
! Compute regional average energy fluxes at PET.
! ====================================================================

      rnpetsum = rnpetsum + rnpet*rescale
      xlepetsum = xlepetsum + xlepet*rescale
      hpetsum = hpetsum + hpet*rescale
      gpetsum = gpetsum + gpet*rescale
      dshpetsum = dshpetsum + dspet*rescale
      tkpetsum = tkpetsum + tkpet*rescale
      tkmidpetsum = tkmidpetsum + tkmidpet*rescale
      tkdeepsum = tkdeepsum + Tdeepstep*rescale

! ====================================================================
! Compute pixel total actual surface energy fluxes for the time step.
! ====================================================================

      if (i_moss+i_und.gt.0) then

! --------------------------------------------------------------------&
! If understory/moss is represented than add their fluxes temperatures
! in the total for the pixel.
! --------------------------------------------------------------------&


!         rnact = (1.-f_moss-f_und)*&
!                       rnact+&
!                       f_und*rnact_us+&
!                       f_moss*rnact_moss
         rnact = canclos*rnact+ f_und*rnact_us+ f_moss*rnact_moss
         xleact = canclos*xleact+ f_und*xleact_us+ f_moss*xleact_moss
         hact = canclos*hact+ f_und*hact_us+ f_moss*hact_moss
         gact = canclos*gact+ f_und*gact_us+ f_moss*gact_moss
!         gact = (1.-f_moss-f_und)*&
!                      gact+&
!                      f_und*gact_us+&
!                      f_moss*gact_moss
!         dshact = (1.-f_moss-f_und)*&
!                        dshact+&
!                        f_und*dshact_us+&
!                        f_moss*dshact_moss
         dshact = canclos*dshact+ f_und*dshact_us+ f_moss*dshact_moss
         tkact = (1.-f_moss-f_und)* tkact+ f_und*tkact_us+ f_moss*tkact_moss
         tkmid = (1.-f_moss-f_und)* tkmid+ f_und*tkmid_us+ f_moss*tkmid_moss

      endif

! ====================================================================
! Compute areal average actual surface energy fluxes for the time step.
! ====================================================================

      rnsum = rnsum + rnact*rescale
      xlesum = xlesum + xleact*rescale
      hsum = hsum + hact*rescale
      gsum = gsum + gact*rescale
      dshsum = dshsum + dshact*rescale
      tksum = tksum + tkact*rescale
      tkmidsum = tkmidsum + tkmid*rescale

! ====================================================================
! Format statements.
! ====================================================================

125   format (1i5,9(f11.5," "))
126   format (1i5,7(f11.5," "))

      return

      end
