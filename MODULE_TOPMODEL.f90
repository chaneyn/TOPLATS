MODULE MODULE_TOPMODEL

USE MODULE_VARIABLES_OLD

USE MODULE_VARIABLES
  
implicit none

contains

! ====================================================================
!
!               subroutine Update_Catchment
!
! ====================================================================
!
! Subroutine to run TOPMODEL section
!
! ====================================================================

  subroutine Update_Catchments(GLOBAL,CAT,GRID)

    implicit none
    type (GLOBAL_template),intent(in) :: GLOBAL
    type (GRID_template),dimension(:),intent(inout) :: GRID
    type (CATCHMENT_template),dimension(:),intent(inout) :: CAT
    integer :: icatch,isoil,ilandc

    do ipix = 1,GLOBAL%npix

      isoil = GRID(ipix)%SOIL%isoil
      ilandc = GRID(ipix)%VEG%ilandc
      icatch = GRID(ipix)%VARS%icatch

      call sumflx_catchment(CAT(icatch),GRID(ipix)%VARS,GLOBAL,&
         GRID(ilandc)%VEG,GRID(isoil)%SOIL,GRID(ipix)%MET,&
         ilandc)

    enddo

    do icatch = 1,GLOBAL%ncatch

      call catflx(GLOBAL%pixsiz,CAT(icatch))

      call upzbar(icatch,CAT(icatch),GLOBAL,GRID)
  
    enddo 

  end subroutine Update_Catchments

! ====================================================================
!
!		subroutine instep_catchment
!
! ====================================================================
!
! Subroutine to initialize regional scale water balance variables.
!
! ====================================================================

  subroutine instep_catchment(ncatch,CAT)

      implicit none
      type (CATCHMENT_template),dimension(:),intent(inout) :: CAT
      integer,intent(in) :: ncatch
      integer :: kk

! ====================================================================
! Initialize variables for catchment average/total values
! ====================================================================

      do kk=1,ncatch

! --------------------------------------------------------------------
! Evaporation and condensation.
! --------------------------------------------------------------------

         CAT(kk)%ettot = zero
         CAT(kk)%etstsum = zero
         CAT(kk)%etwtsum = zero
         CAT(kk)%etbssum = zero
         CAT(kk)%etdcsum = zero
         CAT(kk)%etwcsum = zero
         CAT(kk)%etlakesum = zero
         CAT(kk)%contot = zero

! --------------------------------------------------------------------
! Infiltration/runoff/precipitation.
! --------------------------------------------------------------------

         CAT(kk)%pptsum = zero
         CAT(kk)%pnetsum = zero
         CAT(kk)%sxrtot = zero
         CAT(kk)%xixtot = zero
         CAT(kk)%qsurf = zero
         CAT(kk)%ranrun = zero
         CAT(kk)%conrun = zero

! --------------------------------------------------------------------
! Vertical soil moisture fluxes and water table updating.
! --------------------------------------------------------------------

         CAT(kk)%zbar = CAT(kk)%zbar1
         CAT(kk)%capsum = zero
         CAT(kk)%gwtsum = zero
         CAT(kk)%rzpsum = zero
         CAT(kk)%tzpsum = zero

! --------------------------------------------------------------------
! State variables.
! --------------------------------------------------------------------

         CAT(kk)%fwcat = zero

        enddo

        return
      
      end subroutine instep_catchment
      
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

      subroutine catflx(pixsiz,CAT)

      implicit none
      type (CATCHMENT_template),intent(inout) :: CAT
      real*8 area,pixsiz,r_lakearea,ettot,etstsum,etwtsum,etlakesum
      real*8 etbssum,fbs,etdcsum,etwcsum,pptsum,pnetsum,contot
      real*8 qsurf,sxrtot,xixtot,ranrun,conrun,gwtsum,capsum,tzpsum
      real*8 rzpsum,fwcat
      real*8 catpix,catlakpix,catvegpix
       area = CAT%area
       ettot = CAT%ettot
       etstsum = CAT%etstsum
       etwtsum = CAT%etwtsum
       etlakesum = CAT%etlakesum
       etbssum = CAT%etbssum
       fbs = CAT%fbs
       etdcsum = CAT%etdcsum
       etwcsum = CAT%etwcsum
       pptsum = CAT%pptsum 
       pnetsum = CAT%pnetsum
       contot = CAT%contot
       qsurf = CAT%qsurf
       sxrtot = CAT%sxrtot
       xixtot = CAT%xixtot
       ranrun = CAT%ranrun
       conrun = CAT%conrun
       gwtsum = CAT%gwtsum
       capsum = CAT%capsum
       tzpsum = CAT%tzpsum
       rzpsum = CAT%rzpsum
       fwcat = CAT%fwcat

! ====================================================================
! Calculate the number of pixels in the current catchment.
! ====================================================================

      catpix = area/pixsiz/pixsiz
      !catlakpix = r_lakearea/pixsiz/pixsiz
      catvegpix = (area-r_lakearea)/pixsiz/pixsiz

! ====================================================================
! Find catchment average evapotranspiration rates.
! ====================================================================

      ettot = ettot / catpix
      etstsum = etstsum / catvegpix
      etwtsum = etwtsum / catvegpix
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
      capsum = capsum / catvegpix

! ====================================================================
! Compute average available porosity above the water table.
! ====================================================================

      tzpsum = tzpsum / catvegpix 
      rzpsum = rzpsum / catvegpix 

! ====================================================================
! Calculate catchment fractions of wet canopy.
! ====================================================================

      fwcat = fwcat / catvegpix / (one-fbs)

! ====================================================================
! Format statements.
! ====================================================================

1000  format(2i5,4f10.5,2f7.3)
1100  format(2i5,7f10.5)

       CAT%area = area
       CAT%ettot = ettot
       CAT%etstsum = etstsum
       CAT%etwtsum = etwtsum
       CAT%etlakesum = etlakesum
       CAT%etbssum = etbssum
       CAT%fbs = fbs
       CAT%etdcsum = etdcsum
       CAT%etwcsum = etwcsum
       CAT%pptsum = pptsum
       CAT%pnetsum = pnetsum
       CAT%contot = contot
       CAT%qsurf = qsurf
       CAT%sxrtot = sxrtot
       CAT%xixtot = xixtot
       CAT%ranrun = ranrun
       CAT%conrun = conrun
       CAT%gwtsum = gwtsum
       CAT%capsum = capsum
       CAT%tzpsum = tzpsum
       CAT%rzpsum = rzpsum 
       CAT%fwcat = fwcat

      return

      end subroutine catflx

! ====================================================================
!
!                       subroutine upzbar
!
! ====================================================================
!
!  Updates the average water table depth.
!
! ====================================================================

      subroutine upzbar(ic,CAT,GLOBAL,GRID)

      implicit none
      type (GLOBAL_template),intent(in) :: GLOBAL
      type (GRID_template),dimension(:),intent(inout) :: GRID
      type (CATCHMENT_template),intent(inout) :: CAT
      integer :: ic,iopbf,npix,mm
      integer :: ilandc(GLOBAL%nrow*GLOBAL%ncol),ivgtyp(GLOBAL%nrow*GLOBAL%ncol)
      integer :: icatch(GLOBAL%nrow*GLOBAL%ncol),isoil(GLOBAL%nrow*GLOBAL%ncol)
      real*8 q0,ff,zbar,dtil
      real*8 basink,dd,xlength
      real*8 gwtsum,capsum,area
      real*8 r_lakearea,dt,etwtsum
      real*8 rzpsum,tzpsum,psicav
      real*8 zrzmax,zbar1,qbreg,zbar1rg
      real*8 psic(GLOBAL%nrow*GLOBAL%ncol),tzsm1(GLOBAL%nrow*GLOBAL%ncol),rzsm1(GLOBAL%nrow*GLOBAL%ncol)
      real*8 zw(GLOBAL%nrow*GLOBAL%ncol),thetas(GLOBAL%nrow*GLOBAL%ncol)
      real*8 zero,one,two,three,four,five,six
      real*8 pixsiz
      real*8 qb,hbar,zbrflx,zbrpor,qzbar,dzbar

      iopbf = GLOBAL%iopbf
      q0 = CAT%q0
      ff = CAT%ff
      zbar = CAT%zbar
      dtil = CAT%dtil
      basink = CAT%basink
      dd = CAT%dd
      xlength = CAT%xlength
      gwtsum = CAT%gwtsum
      capsum = CAT%capsum
      area = CAT%area 
      dt = GLOBAL%dt
      etwtsum = CAT%etwtsum
      rzpsum = CAT%rzpsum 
      tzpsum = CAT%tzpsum
      psicav = CAT%psicav 
      ivgtyp = GRID%VEG%ivgtyp
      ilandc = GRID%VEG%ilandc
      npix = GLOBAL%npix
      icatch = GRID%VARS%icatch
      zw = GRID%VARS%zw
      psic = GRID%SOIL%psic
      isoil = GRID%SOIL%isoil
      zrzmax = GLOBAL%zrzmax
      tzsm1 = GRID%VARS%tzsm1
      thetas = GRID%SOIL%thetas 
      rzsm1 = GRID%VARS%rzsm1
      zbar1 = CAT%zbar1
      pixsiz = GLOBAL%pixsiz

! ====================================================================
! Chose option for calculating baseflow.
! ====================================================================

      if (iopbf.eq.0) then

! --------------------------------------------------------------------&
! Calculate baseflow according to Sivapalan et al.(1987).
! --------------------------------------------------------------------&

         qb = q0 * dexp(-ff*zbar)

      else
! --------------------------------------------------------------------&
! Calculate baseflow according to Troch et al.(1992).
! --------------------------------------------------------------------&

         hbar = dtil - zbar
         qb = 5.772*basink*hbar*hbar*dd*xlength

      endif

! ====================================================================
! Determine net recharge to water table and available soil 
! water storage.
! ====================================================================

      zbrflx = (gwtsum - capsum - etwtsum - (qb/ area)) * dt
      zbrpor = (rzpsum+tzpsum) * (zbar-psicav)

      if (zbrflx.gt.zbrpor) then

! ====================================================================
! If net recharge exceeds storage assign overflow to streamflow.
! ====================================================================

         qzbar = (zbrflx - zbrpor)/dt
         zbrflx = zbrpor
         qb = qb + qzbar*area

! --------------------------------------------------------------------&
! If water table is falling but storage is full zbar1 will
! blow up. use rzsm1 and tzsm1 (vs rzsm and tzsm) to
! calculate storage so that it is > zero.
! --------------------------------------------------------------------&

      else if ((zbrflx.le.zero).and.(zbrpor.le.zero)) then

! --------------------------------------------------------------------&
! Recalculate rzpsum and tzpsum with new soil moistures.
! --------------------------------------------------------------------&

         rzpsum = zero
         tzpsum = zero

         do 100 mm=1,npix

            if (ivgtyp(ilandc(mm)).ge.0) then

               if ((icatch(mm)).eq.ic) then

                  if ((zw(mm)-psic(isoil(mm))).gt.zrzmax) then

                     tzpsum = tzpsum+(thetas(isoil(mm))-tzsm1(mm))

                  else if ((zw(mm)-psic(isoil(mm)).gt.zero)) then

                     rzpsum = rzpsum+(thetas(isoil(mm))-rzsm1(mm))

                  endif

               endif

            endif

100      continue

      endif

! ====================================================================
! Update zbar by taking the total flux and dividing by
! the average porosity just above the water table.
! ====================================================================

! --------------------------------------------------------------------&
! If the available porosity is nonzero divide the flux by its value.
! --------------------------------------------------------------------&

      if ( (rzpsum+tzpsum).gt.(0.001d0)) then

         zbar1 = zbar - zbrflx/(rzpsum+tzpsum)

      endif

      if ( (rzpsum+tzpsum).le.(0.001d0)) then

         zbar1=zbar

      endif

! ====================================================================
! Find change in water table.
! ====================================================================

      dzbar = zbar1-zbar

! ====================================================================
! Find new average water table depth and baseflow for entire region.
! ====================================================================

      qbreg = qbreg + qb
      zbar1rg = zbar1rg + zbar1*area/(pixsiz*pixsiz)

! ====================================================================
! Format statements.
! ====================================================================

      CAT%zbar1 = zbar1
      CAT%qb = qb

      return

      end subroutine upzbar

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

      subroutine sumflx_catchment(CAT,GRID_VARS,GLOBAL,&
       GRID_VEG,GRID_SOIL,GRID_MET,ilandc)

      implicit none
    
      type (GRID_VARS_template),intent(in) :: GRID_VARS
      type (GRID_VEG_template),intent(in) :: GRID_VEG
      type (GRID_SOIL_template),intent(in) :: GRID_SOIL
      type (GRID_MET_template),intent(in) :: GRID_MET
      type (CATCHMENT_template),intent(inout) :: CAT
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

         else

! ....................................................................
! Add the evapotranspiration rate from the understory/moss layer.
! ....................................................................

            xlhv=1000000000.*(2.501-0.002361*tair)
            dummy = canclos*xleact+ f_und*xleact_us+ f_moss*xleact_moss
            etpixloc = dummy/xlhv
            ettot = ettot + etpixloc*rescale

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
         etstsum = etstore*rescale
         etwtsum = etwt*rescale

         if (ivgtyp.eq.0) then

! ====================================================================
! Add up evaporation from bare soil values.
! ====================================================================

            etbssum = evtact*rescale

         else

! ====================================================================
! Add up evaporation from dry and wet canopy for vegetated pixels.
! ====================================================================

            if ( (i_und.eq.0).and.(i_moss.eq.0) ) then

! --------------------------------------------------------------------&
! Understory/moss is not represented.
! --------------------------------------------------------------------&

               etdcsum = evtact*dc*(one-fw)*rescale
               etwcsum = epwms*dc*rescale

            else

! --------------------------------------------------------------------&
! Add the values from the understory/moss in the totals.
! --------------------------------------------------------------------&

               etdcsum = ((1.-f_moss- f_und)*&
                                       evtact*dc*(one-fw) +&
                                       f_und*evtact_us*dc_us*(one-fw_us) +&
                                       f_moss*evtact_moss)*rescale

               etwcsum = (epwms*dc* (1.-f_moss- f_und)+&
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

      contot = conpix*rescale

! ====================================================================
! Compute total precipitation and surface runoff for the time step.
! ====================================================================

      pptsum = pptms*rescale
      pnetsum = pnet*rescale
      qsurf = runtot*rescale

! ====================================================================
! Compute total saturation and infiltration excess runoff for the
! time step.
! ====================================================================

      sxrtot = satxr *rescale
      xixtot = xinfxr*rescale

! ====================================================================
! Compute total runoff due to rainfall and due to condensation.
! ====================================================================

      if (pptms.gt.zero) then

! --------------------------------------------------------------------&
! When under rainfall runoff has to be due to rain.
! --------------------------------------------------------------------&

         ranrun = runtot*rescale

      else

! --------------------------------------------------------------------&
! When no precipitation runoff has to be due to condensation.
! --------------------------------------------------------------------&

         conrun = runtot *rescale

      endif  

! ====================================================================
! Compute checks on canopy water balance, root zone water balance,&
! and transmission zone balance.
! ====================================================================

      if (ivgtyp.ge.(0)) then

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


         gwtsum = gwt*rescale

! ====================================================================
! Compute the diffusion totals for the time step.
! ====================================================================

         capsum = - difwt*rescale

! ====================================================================
! Compute change in storage above the water table for the time step 
! and perform water balance.
! ====================================================================

         dstore = dswc + dsrz + dstz

! ====================================================================
! Add up all the input terms towards the water table and make the
! regional total.
! ======================================================================


         svarhs = dt*(pptms - etstore - runtot - gwt - difwt)

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

            tzpsum = (thetas-tzsm)*rescale

         else if (zrz.gt.zero) then

! --------------------------------------------------------------------&
! If the root zone is unsaturated and the transmission zone
! unsaturated than the region of interest is the root zone.
! --------------------------------------------------------------------&

            rzpsum = (thetas-rzsm)*rescale

         endif

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


! ====================================================================
! Format statements.
! ====================================================================

125   format (1i5,9(f11.5," "))
126   format (1i5,7(f11.5," "))

!$OMP CRITICAL

!Catchment Variables
CAT%etstsum = CAT%etstsum + etstsum
CAT%etwtsum = CAT%etwtsum + etwtsum
CAT%etbssum = CAT%etbssum + etbssum
CAT%etdcsum = CAT%etdcsum + etdcsum
CAT%etwcsum = CAT%etwcsum + etwcsum
CAT%contot = CAT%contot + contot
CAT%pptsum = CAT%pptsum + pptsum
CAT%pnetsum = CAT%pnetsum + pnetsum
CAT%qsurf = CAT%qsurf + qsurf
CAT%sxrtot = CAT%sxrtot + sxrtot
CAT%xixtot = CAT%xixtot + xixtot
CAT%ranrun = CAT%ranrun + ranrun
CAT%conrun = CAT%conrun + conrun
CAT%gwtsum = CAT%gwtsum + gwtsum
CAT%capsum = CAT%capsum + capsum
CAT%tzpsum = CAT%tzpsum + tzpsum
CAT%rzpsum = CAT%rzpsum + rzpsum
!$OMP END CRITICAL 

      return

      end subroutine sumflx_catchment

END MODULE MODULE_TOPMODEL
