MODULE MODULE_REGIONAL

USE MODULE_VARIABLES_OLD

USE MODULE_VARIABLES
  
implicit none

contains

! ====================================================================
!
!               subroutine Update_Regional
!
! ====================================================================
!
! Subroutine to update regional variables
!
! ====================================================================

  subroutine Update_Regional(REG,GRID,GLOBAL)

    implicit none
    type (GLOBAL_template),intent(in) :: GLOBAL
    type (GRID_template),dimension(:),intent(in) :: GRID
    type (REGIONAL_template),intent(inout) :: REG
    integer :: icatch,isoil,ilandc

!#####################################################################
! Initialize regional variables
!#####################################################################

    call instep_regional(REG)

!#####################################################################
! Sum the local water and energy balance fluxes.
!#####################################################################

    do ipix = 1,GLOBAL%npix

      isoil = GRID(ipix)%SOIL%isoil
      ilandc = GRID(ipix)%VEG%ilandc
      icatch = GRID(ipix)%VARS%icatch

      call sumflx_regional(REG,GRID(ipix)%VARS,GLOBAL,&
         GRID(ilandc)%VEG,GRID(isoil)%SOIL,GRID(ipix)%MET,&
         ilandc)

    enddo

  end subroutine Update_Regional

! ====================================================================
!
!		subroutine instep_regional
!
! ====================================================================
!
! Subroutine to initialize regional scale water balance variables.
!
! ====================================================================

  subroutine instep_regional(REG)

      implicit none
      include "help/instep.h"
      type (REGIONAL_template),intent(inout) :: REG

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

        return
      
      end subroutine instep_regional
      
! ====================================================================
!
! 			subroutine sumflx_regional
!
! ====================================================================
!
! Calculates the time step totals of evapotranspiration, runoff, 
! surface energy fluxes and vertical soil-water fluxes.
!
! ====================================================================

      subroutine sumflx_regional(REG,GRID_VARS,GLOBAL,&
       GRID_VEG,GRID_SOIL,GRID_MET,ilandc)

      implicit none
    
      include "help/sumflx.dif.h"
      type (REGIONAL_template),intent(inout) :: REG
      type (GRID_VARS_template),intent(in) :: GRID_VARS
      type (GRID_VEG_template),intent(in) :: GRID_VEG
      type (GRID_SOIL_template),intent(in) :: GRID_SOIL
      type (GRID_MET_template),intent(in) :: GRID_MET
      type (GLOBAL_template),intent(in) :: GLOBAL

! TEMPORARY
!GRID
rnact = GRID_VARS%rnact
xleact = GRID_VARS%xleact
hact = GRID_VARS%hact
gact = GRID_VARS%gact
dshact = GRID_VARS%dshact
tkact = GRID_VARS%tkact
tkmid = GRID_VARS%tkmid
rnpet = GRID_VARS%rnpet
xlepet = GRID_VARS%xlepet
hpet = GRID_VARS%hpet
gpet = GRID_VARS%gpet
dspet = GRID_VARS%dspet
tkpet = GRID_VARS%tkpet
tkmidpet = GRID_VARS%tkmidpet
rzsm1 = GRID_VARS%rzsm1
tzsm1  = GRID_VARS%tzsm1
rzsm = GRID_VARS%rzsm
tzsm = GRID_VARS%tzsm
runtot = GRID_VARS%runtot
pnet = GRID_VARS%pnet
evtact = GRID_VARS%evtact
etpix = GRID_VARS%etpix
ivgtyp = GRID_VEG%ivgtyp
canclos = GRID_VEG%canclos
tair = GRID_MET%tdry
pptms = GRID_MET%pptms
Swq = GRID_VARS%Swq
Swq_us = GRID_VARS%Swq_us
Sdepth = GRID_VARS%Sdepth
Sdepth_us = GRID_VARS%Sdepth_us
Tdeepstep = GRID_SOIL%Tdeepstep

!Soil Data
thetas = GRID_SOIL%thetas

!Point Data
zrz = GRID_VARS%zrz
ztz = GRID_VARS%ztz
smold = GRID_VARS%smold
rzsmold = GRID_VARS%rzsmold
tzsmold = GRID_VARS%tzsmold
capflx = GRID_VARS%capflx
difrz = GRID_VARS%difrz
diftz = GRID_VARS%diftz
grz = GRID_VARS%grz
gtz = GRID_VARS%gtz
satxr = GRID_VARS%satxr
xinfxr = GRID_VARS%xinfxr
dc = GRID_VARS%dc
fw = GRID_VARS%fw
dsrz = GRID_VARS%dsrz
rzrhs = GRID_VARS%rzrhs
dstz = GRID_VARS%dstz
tzrhs = GRID_VARS%tzrhs
dswc = GRID_VARS%dswc
wcrhs = GRID_VARS%wcrhs
!Energy fluxes
epwms = GRID_VARS%epwms
wcip1 = GRID_VARS%wcip1

!GLOBAL
rescale = GLOBAL%mul_fac
dt = GLOBAL%dt
inc_frozen = GLOBAL%inc_frozen

! ====================================================================
! Compute regional average evapotranspiration rate for the time step.    
! ====================================================================
 
  
      if (ivgtyp.eq.(-1)) then

! --------------------------------------------------------------------&
! In case of an open water pixel.
! --------------------------------------------------------------------&

         etpixloc=evtact
         ettot = ettot + etpixloc*rescale
         ettotrg = etpixloc*rescale

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
            ettotrg = etpixloc*rescale

         else

! ....................................................................
! Add the evapotranspiration rate from the understory/moss layer.
! ....................................................................

            xlhv=1000000000.*(2.501-0.002361*tair)
            dummy = canclos*xleact+ f_und*xleact_us+ f_moss*xleact_moss
            etpixloc = dummy/xlhv
            ettot = ettot + etpixloc*rescale
            ettotrg = etpixloc*rescale

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
         etstsumrg = etstore*rescale
         etwtsumrg = etwt*rescale

         if (ivgtyp.eq.0) then

! ====================================================================
! Add up evaporation from bare soil values.
! ====================================================================

            etbssumrg = evtact*rescale

         else

! ====================================================================
! Add up evaporation from dry and wet canopy for vegetated pixels.
! ====================================================================

            if ( (i_und.eq.0).and.(i_moss.eq.0) ) then

! --------------------------------------------------------------------&
! Understory/moss is not represented.
! --------------------------------------------------------------------&

               etdcsumrg = evtact*dc*(one-fw)*rescale
               etwcsumrg = epwms*dc*rescale

            else

! --------------------------------------------------------------------&
! Add the values from the understory/moss in the totals.
! --------------------------------------------------------------------&

               etdcsumrg = ((1.-f_moss- f_und)*&
                                       evtact*dc*(one-fw) +&
                                       f_und*evtact_us*dc_us*(one-fw_us) +&
                                       f_moss*evtact_moss)*rescale

               etwcsumrg = (epwms*dc*(1.-f_moss-f_und)+&
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

            conpix=0.d0

         endif

      endif

      contotrg = conpix*rescale

! ====================================================================
! Compute total precipitation and surface runoff for the time step.
! ====================================================================

      pptsumrg = pptms*rescale
      pnetsumrg = pnet*rescale
      qsurfrg = runtot*rescale

! ====================================================================
! Compute total saturation and infiltration excess runoff for the
! time step.
! ====================================================================

      sxrtotrg = satxr *rescale
      xixtotrg = xinfxr*rescale

! ====================================================================
! Compute total runoff due to rainfall and due to condensation.
! ====================================================================

      if (pptms.gt.zero) then

! --------------------------------------------------------------------&
! When under rainfall runoff has to be due to rain.
! --------------------------------------------------------------------&

         ranrunrg = runtot*rescale

      else

! --------------------------------------------------------------------&
! When no precipitation runoff has to be due to condensation.
! --------------------------------------------------------------------&

         conrunrg = runtot *rescale

      endif  

! ====================================================================
! Compute checks on canopy water balance, root zone water balance,&
! and transmission zone balance.
! ====================================================================

      if (ivgtyp.ge.(0)) then

         dswcsum = dswc*rescale
         wcrhssum = wcrhs*rescale

         dsrzsum = dsrz*rescale
         rzrhssum = rzrhs*rescale

         dstzsum = dstz*rescale
         tzrhssum = tzrhs*rescale

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


         gwtsumrg = gwt*rescale
         grzsumrg = grz*rescale
         gtzsumrg = gtz*rescale

! ====================================================================
! Compute the diffusion totals for the time step.
! ====================================================================

         capsumrg =  - difwt*rescale
         difrzsumrg =  - difrz*rescale

! ====================================================================
! Compute change in storage above the water table for the time step 
! and perform water balance.
! ====================================================================

         dstore = dswc + dsrz + dstz
         dssum = dstore

! ====================================================================
! Add up all the input terms towards the water table and make the
! regional total.
! ======================================================================


         svarhs = dt*(pptms - etstore - runtot - gwt - difwt)

         svarhssum = svarhs*rescale

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

         rzsmav = rzsm1*rescale
         tzsmav = tzsm1*rescale
         Swqsum = Swq*rescale
         Swq_ussum = Swq_us*rescale
         Sdepthsum = Sdepth*rescale
         Sdepth_ussum = Sdepth_us*rescale


! ====================================================================
! Make a regional total of canopy interception storage.
! ====================================================================

         wcip1sum = wcip1*rescale

      endif

! ====================================================================
! Write the components of the energy balance out if requested.
! ====================================================================

      if (ivgtyp.eq.(-1)) then

         tkmid=0.d0
         tkmidpet=0.d0
         tskinact_moss=0.d0

      endif

! ====================================================================
! Compute pixel total energy fluxes at PET.
! ====================================================================

      if (i_moss+i_und.gt.0) then

! --------------------------------------------------------------------&
! If understory/moss is represented than add their fluxes temperatures
! in the total for the pixel.
! --------------------------------------------------------------------&

         rnpet = canclos*rnpet+ f_und*rnpet_us+ f_moss*rnpet_moss
         xlepet = canclos*xlepet+ f_und*xlepet_us+ f_moss*xlepet_moss
         hpet = canclos*hpet+ f_und*hpet_us+ f_moss*hpet_moss
         gpet = canclos*gpet+ f_und*gpet_us+ f_moss*gpet_moss
         dspet = canclos*dspet+ f_und*dspet_us+ f_moss*dspet_moss
         tkpet = (1.-f_moss-f_und)* tkpet+ f_und*tkpet_us+ f_moss*tkpet_moss
         tkmidpet = (1.-f_moss-f_und)* tkmidpet+ f_und*tkmidpet_us+&
                    f_moss*tkmidpet_moss

      endif

! ====================================================================
! Compute regional average energy fluxes at PET.
! ====================================================================

      rnpetsum = rnpet*rescale
      xlepetsum = xlepet*rescale
      hpetsum = hpet*rescale
      gpetsum = gpet*rescale
      dshpetsum = dspet*rescale
      tkpetsum = tkpet*rescale
      tkmidpetsum = tkmidpet*rescale
      tkdeepsum = Tdeepstep*rescale

! ====================================================================
! Compute pixel total actual surface energy fluxes for the time step.
! ====================================================================

      if (i_moss+i_und.gt.0) then

! --------------------------------------------------------------------&
! If understory/moss is represented than add their fluxes temperatures
! in the total for the pixel.
! --------------------------------------------------------------------&


         rnact = canclos*rnact+ f_und*rnact_us+ f_moss*rnact_moss
         xleact = canclos*xleact+ f_und*xleact_us+ f_moss*xleact_moss
         hact = canclos*hact+ f_und*hact_us+ f_moss*hact_moss
         gact = canclos*gact+ f_und*gact_us+ f_moss*gact_moss
         dshact = canclos*dshact+ f_und*dshact_us+ f_moss*dshact_moss
         tkact = (1.-f_moss-f_und)* tkact+ f_und*tkact_us+ f_moss*tkact_moss
         tkmid = (1.-f_moss-f_und)* tkmid+ f_und*tkmid_us+ f_moss*tkmid_moss

      endif

! ====================================================================
! Compute areal average actual surface energy fluxes for the time step.
! ====================================================================

      rnsum = rnact*rescale
      xlesum = xleact*rescale
      hsum = hact*rescale
      gsum = gact*rescale
      dshsum = dshact*rescale
      tksum = tkact*rescale
      tkmidsum = tkmid*rescale

! ====================================================================
! Format statements.
! ====================================================================

125   format (1i5,9(f11.5," "))
126   format (1i5,7(f11.5," "))

!$OMP CRITICAL
!Regional 
REG%ettotrg = REG%ettotrg + ettotrg
REG%etstsumrg = REG%etstsumrg + etstsumrg
REG%etwtsumrg = REG%etwtsumrg + etwtsumrg
REG%etbssumrg = REG%etbssumrg + etbssumrg
REG%etdcsumrg = REG%etdcsumrg + etdcsumrg
REG%etwcsumrg = REG%etwcsumrg + etwcsumrg
REG%pptsumrg = REG%pptsumrg + pptsumrg
REG%pnetsumrg = REG%pnetsumrg + pnetsumrg
REG%contotrg = REG%contotrg + contotrg
REG%sxrtotrg = REG%sxrtotrg + sxrtotrg
REG%xixtotrg = REG%xixtotrg + xixtotrg
REG%qsurfrg = REG%qsurfrg + qsurfrg
REG%ranrunrg = REG%ranrunrg + ranrunrg
REG%conrunrg = REG%conrunrg + conrunrg
REG%gwtsumrg = REG%gwtsumrg + gwtsumrg
REG%grzsumrg = REG%grzsumrg + grzsumrg
REG%gtzsumrg = REG%gtzsumrg + gtzsumrg
REG%capsumrg = REG%capsumrg + capsumrg
REG%difrzsumrg = REG%difrzsumrg + difrzsumrg
REG%rnpetsum = REG%rnpetsum + rnpetsum
REG%xlepetsum = REG%xlepetsum + xlepetsum
REG%hpetsum = REG%hpetsum + hpetsum
REG%gpetsum = REG%gpetsum + gpetsum
REG%dshpetsum = REG%dshpetsum + dshpetsum
REG%tkpetsum = REG%tkpetsum + tkpetsum
REG%tkmidpetsum = REG%tkmidpetsum + tkmidpetsum
REG%tkdeepsum = REG%tkdeepsum + tkdeepsum
REG%wcip1sum = REG%wcip1sum + wcip1sum
REG%dswcsum = REG%dswcsum + dswcsum
REG%wcrhssum = REG%wcrhssum + wcrhssum
REG%dsrzsum = REG%dsrzsum + dsrzsum
REG%rzrhssum = REG%rzrhssum + rzrhssum
REG%dstzsum = REG%dstzsum + dstzsum
REG%tzrhssum = REG%tzrhssum + tzrhssum
REG%dssum = REG%dssum + dssum
REG%svarhssum = REG%svarhssum + svarhssum
REG%rzsmav = REG%rzsmav + rzsmav
REG%tzsmav = REG%tzsmav + tzsmav
REG%Swqsum = REG%Swqsum + Swqsum
REG%Swq_ussum = REG%Swq_ussum + Swq_ussum
REG%Sdepthsum = REG%Sdepthsum + Sdepthsum
REG%Sdepth_ussum = REG%Sdepth_ussum + Sdepth_ussum
REG%rnsum = REG%rnsum + rnsum
REG%xlesum = REG%xlesum + xlesum
REG%hsum = REG%hsum + hsum
REG%gsum = REG%gsum + gsum
REG%dshsum = REG%dshsum + dshsum
REG%tksum = REG%tksum + tksum
REG%tkmidsum = REG%tkmidsum + tkmidsum

!$OMP END CRITICAL 

      return

      end subroutine sumflx_regional

END MODULE MODULE_REGIONAL
