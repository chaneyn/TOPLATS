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

  subroutine Update_Regional(REG,GRID,GLOBAL,CAT)

    implicit none
    type (GLOBAL_template),intent(in) :: GLOBAL
    type (GRID_template),dimension(:),intent(in) :: GRID
    type (REGIONAL_template),intent(inout) :: REG
    type (CATCHMENT_template),dimension(:),intent(in) :: CAT
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

    call Compute_Regional(REG,GLOBAL,GRID,CAT)

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
    
      type (REGIONAL_template),intent(inout) :: REG
      type (GRID_VARS_template),intent(in) :: GRID_VARS
      type (GRID_VEG_template),intent(in) :: GRID_VEG
      type (GRID_SOIL_template),intent(in) :: GRID_SOIL
      type (GRID_MET_template),intent(in) :: GRID_MET
      type (GLOBAL_template),intent(in) :: GLOBAL
      integer ivgtyp,i_und,i_moss,i,ilandc,inc_frozen
      real*8 f_und,f_moss,canclos,dc,fw,dc_us,fw_us,etpix
      real*8 evtact,ettot,ettotrg,epwms,evtact_us,epwms_us
      real*8 evtact_moss,etstore,etwt,etstsum,etstsumrg
      real*8 etwtsumrg,etwtsum,etbssum,etbssumrg,etdcsum
      real*8 etdcsumrg,etwcsum,etwcsumrg,etlakesumrg,etlakesum
      real*8 bsdew,contot,contotrg,pptsum,pptms,pptsumrg
      real*8 pnetsum,pnet,pnetsumrg,qsurf,runtot,qsurfrg
      real*8 sxrtot,satxr,sxrtotrg,xixtot,xinfxr,xixtotrg
      real*8 ranrun,ranrunrg,conrun,conrunrg,wcip1sum
      real*8 wcip1,dswcsum,dswc,wcrhssum,wcrhs,dsrzsum,dsrz
      real*8 rzrhssum,rzrhs,dstzsum,dstz,tzrhssum,tzrhs,zrz
      real*8 gwt,grz,gtz,ztz,gwtsum,gwtsumrg,grzsumrg,gtzsumrg
      real*8 capsum,capsumrg,difrzsumrg,diftz,difrz
      real*8 dstore,dssum,svarhs,rzsm1_u,tzsm1_u,rzsm1,tzsm1
      real*8 rzsmav,tzsmav,tzpsum,thetas,rzpsum,r_mossm,rzsm
      real*8 tzsm,rnact_moss,xleact_moss,hact_moss,gact_moss
      real*8 dshact_moss,tskinact_moss,tkact_moss,tkmid_moss
      real*8 rnact_us,xleact_us,hact_us,gact_us,dshact_us
      real*8 tkact_us,tkmid_us,rnact,xleact,hact,gact,dshact
      real*8 tkact,tkmid,rnsum,xlesum,hsum,gsum,dshsum,tksum
      real*8 tkmidsum,rnpet,rnpet_us,rnpet_moss,xlepet,xlepet_us
      real*8 xlepet_moss,hpet,hpet_us,hpet_moss,gpet,gpet_us
      real*8 gpet_moss,dspet,dspet_us,dspet_moss,tkpet,tkpet_us
      real*8 tkpet_moss,tkmidpet,tkmidpet_us,tkmidpet_moss
      real*8 rnpetsum,xlepetsum,hpetsum,gpetsum,dshpetsum
      real*8 tkpetsum,tkmidpetsum,tkdeepsum,Tdeepstep,dt
      real*8 svarhssum,Swqsum,Swq_ussum,Swq,Swq_us
      real*8 rescale,etpixloc
      real*8 conpix
      real*8 difwt
      real*8 tair,xlhv,dummy
      real*8 Sdepthsum,Sdepth_ussum,Sdepth,Sdepth_us


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

! ====================================================================
!
!                       subroutine lswb
!
! ====================================================================
!
!  Check regional scale water balance and sum simulation totals.
!
! ====================================================================

      subroutine Compute_Regional(REG,GLOBAL,GRID,CAT)

      implicit none
      type (REGIONAL_template),intent(inout) :: REG
      type (GLOBAL_template),intent(in) :: GLOBAL
      type (GRID_template),dimension(:),intent(in) :: GRID
      type (CATCHMENT_template),dimension(:),intent(in) :: CAT
      real*8 rest
      real*8 :: tmp
      integer :: isoil,icatch
      integer ncatch,npix,ic
      integer ipix,landc,nlcs,ivgtyp(GLOBAL%nrow*GLOBAL%ncol)
      integer count24
      real*8 nlakpix,nvegpix
      real*8 qb24sum
      real*8 pixsiz,ettotrg,etlakesumrg
      real*8 etstsumrg,etwtsumrg,fbsrg,etbssumrg,etdcsumrg
      real*8 etwcsumrg,pptsumrg,pnetsumrg,qsurfrg,sxrtotrg
      real*8 xixtotrg,contotrg,ranrunrg,conrunrg,qbreg,gwtsumrg
      real*8 grzsumrg,gtzsumrg,capsumrg,difrzsumrg,dswcsum
      real*8 wcrhssum,dsrzsum,rzrhssum,dstzsum,tzrhssum,dssum
      real*8 svarhssum,rzsmav,tzsmav,rnpetsum,xlepetsum,hpetsum
      real*8 gpetsum,dshpetsum,tkpetsum,tkmidpetsum,rnsum,xlesum
      real*8 hsum,gsum,dshsum,tksum,tkmidsum,tkdeepsum,fwreg
      real*8 wcip1sum,zbar1rg,pr3sat,perrg2,pr2sat,pr2uns,perrg1
      real*8 pr1sat,pr1rzs,pr1tzs,pr1uns,persxr,perixr,persac
      real*8 peruac,perusc,wcsum,zbarrg,tot
      real*8 vegd,Swqsum,Swq_ussum
      real*8 Sdepthsum,Sdepth_ussum

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
! Find new average water table depth and baseflow for entire region.
! ====================================================================

      qbreg = zero
      zbar1rg = zero
      do icatch = 1,GLOBAL%ncatch
        qbreg = qbreg + CAT(icatch)%qb
        zbar1rg = zbar1rg + CAT(icatch)%zbar1*CAT(icatch)%area/(GLOBAL%pixsiz*GLOBAL%pixsiz)
      enddo

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

      REG%wcip1sum = wcip1sum
      REG%zbar1rg = zbar1rg
      REG%ettotrg = ettotrg
      REG%hsum = hsum
      REG%etwtsumrg = etwtsumrg
      REG%dswcsum  = dswcsum
      REG%wcrhssum = wcrhssum
      REG%etdcsumrg = etdcsumrg
      REG%etwcsumrg = etwcsumrg
      REG%tkmidsum = tkmidsum
      REG%tksum = tksum
      REG%gsum = gsum
      REG%xlesum = xlesum
      REG%perrg2 = perrg2
      REG%pr3sat = pr3sat
      REG%capsumrg = capsumrg
      REG%difrzsumrg = difrzsumrg
      REG%pr2uns = pr2uns
      REG%gtzsumrg = gtzsumrg
      REG%grzsumrg = grzsumrg
      REG%qbreg = qbreg
      REG%pr2sat = pr2sat
      REG%perrg1 = perrg1
      REG%pr1uns = pr1uns
      REG%persac = persac
      REG%peruac = peruac
      REG%perusc = perusc
      REG%rnsum = rnsum
      REG%tkdeepsum = tkdeepsum
      REG%rzsmav = rzsmav 
      REG%dsrzsum = dsrzsum
      REG%rzrhssum = rzrhssum
      REG%tzsmav = tzsmav
      REG%contotrg = contotrg
      REG%dstzsum = dstzsum
      REG%tzrhssum = tzrhssum
      REG%gwtsumrg = gwtsumrg
      REG%sxrtotrg = sxrtotrg
      REG%pptsumrg = pptsumrg
      REG%pnetsumrg = pnetsumrg 
      REG%persxr = persxr
      REG%qsurfrg = qsurfrg
      REG%pr1tzs = pr1tzs
      REG%Swqsum = Swqsum
      REG%Sdepthsum = Sdepthsum

      return

      end subroutine Compute_Regional


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
      integer :: i,ipix
      integer ievcon
      real*8 fwcat,fwreg,zw,psic,pr3sat,zrzmax,perrg2,rzsm1,thetas
      real*8 pr2sat,pr2uns,perrg1,tzsm1,pr1sat,persac,peruac,perusc
      real*8 pr1rzs,pr1tzs,pr1uns,satxr,persxr,xinfxr,perixr,tolsat
      real*8 fw,mul_fac
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


END MODULE MODULE_REGIONAL
