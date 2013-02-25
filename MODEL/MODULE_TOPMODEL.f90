MODULE MODULE_TOPMODEL

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

!> run TOPMODEL in specified catchments
subroutine Update_Catchments(GLOBAL,CAT,GRID)

  implicit none
  type (GLOBAL_template),intent(in) :: GLOBAL
  type (GRID_template),dimension(:),intent(inout) :: GRID
  type (CATCHMENT_template),dimension(:),intent(inout) :: CAT
  integer :: icatch,isoil,ilandc

!#####################################################################
! Initialize water balance variables for the time step.
!#####################################################################

  call instep_catchment(GLOBAL%ncatch,CAT)

!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(isoil,ilandc,icatch) 

  do ipix = 1,GLOBAL%npix

    isoil = GRID(ipix)%SOIL%isoil
    ilandc = GRID(ipix)%VEG%ilandc
    icatch = GRID(ipix)%VARS%icatch

!!$OMP CRITICAL 
    call sumflx_catchment(CAT(icatch),GRID(ipix)%VARS,GLOBAL,&
       GRID(ilandc)%VEG,GRID(isoil)%SOIL,GRID(ipix)%MET,&
       ilandc)
!!$OMP END CRITICAL 

  enddo

!!$OMP END PARALLEL DO
 
  do icatch = 1,GLOBAL%ncatch

    call catflx(GLOBAL%pixsiz,CAT(icatch))

    call upzbar(icatch,CAT(icatch),GLOBAL,GRID)

    CAT%tzpsum = CAT(icatch)%smpsum(2)
    CAT%rzpsum = CAT(icatch)%smpsum(1)
 
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

!> Initialize regional scale water balance variables
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

    CAT(kk)%capsum = zero
    CAT(kk)%gwtsum = zero
    CAT(kk)%smpsum = zero

! --------------------------------------------------------------------
! State variables.
! --------------------------------------------------------------------

    CAT(kk)%fwcat = zero

! --------------------------------------------------------------------
! Others.
! --------------------------------------------------------------------

    CAT(kk)%psicav = zero

  enddo

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

!> Calculates the catchment time step totals of evapotranspiration, 
!!  runoff, surface energy fluxes and vertical soil-water fluxes.

subroutine catflx(pixsiz,CAT)

  implicit none
  type (CATCHMENT_template),intent(inout) :: CAT
  real*8 :: pixsiz,catpix,catlakpix,catvegpix

! ====================================================================
! Calculate the number of pixels in the current catchment.
! ====================================================================

  catpix = CAT%area/pixsiz/pixsiz
  catvegpix = CAT%area/pixsiz/pixsiz

! ====================================================================
! Find catchment average evapotranspiration rates.
! ====================================================================

  CAT%ettot = CAT%ettot / catpix
  CAT%etstsum = CAT%etstsum / catvegpix
  CAT%etwtsum = CAT%etwtsum / catvegpix
  CAT%etlakesum = CAT%etlakesum / catlakpix

  if (CAT%fbs.gt.0.) then

    CAT%etbssum = CAT%etbssum / CAT%fbs/catvegpix

  else

    CAT%etbssum = 0.

  endif

  if (CAT%fbs.lt.1.) then

    CAT%etdcsum = CAT%etdcsum / (one-CAT%fbs)/catvegpix
    CAT%etwcsum = CAT%etwcsum / (one-CAT%fbs)/catvegpix

  else

    CAT%etdcsum = 0.
    CAT%etwcsum = 0.

  endif

! ====================================================================
! Find catchment precipitation and condensation rate.
! ====================================================================

  CAT%pptsum = CAT%pptsum / catpix
  CAT%pnetsum = CAT%pnetsum / catpix
  CAT%contot = CAT%contot / catpix

! ====================================================================
! Compute total catchment runoff/infiltration rates.
! ====================================================================

  CAT%qsurf = CAT%qsurf / catvegpix
  CAT%sxrtot = CAT%sxrtot / catvegpix
  CAT%xixtot = CAT%xixtot / catvegpix

! ====================================================================
! Compute total runoff due to rainfall and due to condensation.
! ====================================================================

  CAT%ranrun = CAT%ranrun / catvegpix
  CAT%conrun = CAT%conrun / catvegpix

! ====================================================================
! Compute water table fluxes.
! ====================================================================

  CAT%gwtsum = CAT%gwtsum / catvegpix
  CAT%capsum = CAT%capsum / catvegpix

! ====================================================================
! Compute average available porosity above the water table.
! ====================================================================

  CAT%smpsum = CAT%smpsum / catvegpix

! ====================================================================
! Calculate catchment fractions of wet canopy.
! ====================================================================

  CAT%fwcat = CAT%fwcat / catvegpix / (one-CAT%fbs)

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

!> Updates the average water table depth
subroutine upzbar(ic,CAT,GLOBAL,GRID)

  implicit none
  type (GLOBAL_template),intent(in) :: GLOBAL
  type (GRID_template),dimension(:),intent(inout) :: GRID
  type (CATCHMENT_template),intent(inout) :: CAT
  integer :: ic,mm,ilandc,isoil,icatch,kk
  real*8 :: hbar,qzbar,zbrflx,zbrpor

! ====================================================================
! Chose option for calculating baseflow.
! ====================================================================

  if (GLOBAL%KS_TYPE.eq.0)then

! --------------------------------------------------------------------&
! Calculate baseflow according to Sivapalan et al.(1987).
! --------------------------------------------------------------------&

    CAT%qb = CAT%q0 * dexp(-CAT%ff*CAT%zbar)

  else if (GLOBAL%KS_TYPE.eq.1)then

! --------------------------------------------------------------------&
! Calculate baseflow according to Chaney et al.(2013).
! --------------------------------------------------------------------&

    call Calculate_Qb(CAT%n,CAT%ff,CAT%zbar,CAT%q0,CAT%qb)

  endif

! ====================================================================
! Determine net recharge to water table and available soil 
! water storage.
! ====================================================================

  zbrflx = (CAT%gwtsum - CAT%capsum - CAT%etwtsum - (CAT%qb/ CAT%area)) * GLOBAL%dt
  zbrpor = sum(CAT%smpsum) * (CAT%zbar-CAT%psicav)
 
  if (zbrflx.gt.zbrpor) then

! ====================================================================
! If net recharge exceeds storage assign overflow to streamflow.
! ====================================================================

    qzbar = (zbrflx - zbrpor)/GLOBAL%dt
    zbrflx = zbrpor
    CAT%qb = CAT%qb + qzbar*CAT%area

! --------------------------------------------------------------------&
! If water table is falling but storage is full zbar1 will
! blow up. use rzsm1 and tzsm1 (vs rzsm and tzsm) to
! calculate storage so that it is > zero.
! --------------------------------------------------------------------&

  else if ((zbrflx.le.zero).and.(zbrpor.le.zero)) then

! --------------------------------------------------------------------&
! Recalculate rzpsum and tzpsum with new soil moistures.
! --------------------------------------------------------------------&

    CAT%smpsum = zero

    do mm=1,GLOBAL%npix

      isoil = GRID(mm)%SOIL%isoil
      ilandc = GRID(mm)%VEG%ilandc
      icatch = GRID(mm)%VARS%icatch

      if (GRID(ilandc)%VEG%ivgtyp.ge.0) then

        if (icatch.eq.ic) then

          if ((GRID(mm)%VARS%zw-GRID(isoil)%SOIL%psic).gt.GLOBAL%zrzmax) then

            CAT%smpsum(GLOBAL%nlayer) = CAT%smpsum(GLOBAL%nlayer)+(GRID(isoil)%SOIL%thetas-GRID(mm)%VARS%sm1(GLOBAL%nlayer))

          else 
  
            if ((GRID(mm)%VARS%zw-GRID(isoil)%SOIL%psic.gt.zero))then

              do kk=GLOBAL%nlayer-1,1,-1
          
                 if (GRID(mm)%VARS%zw-GRID(isoil)%SOIL%psic.le.GLOBAL%zmax_layer(kk)) then
  
                   CAT%smpsum(kk) = CAT%smpsum(kk)+(GRID(isoil)%SOIL%thetas-GRID(mm)%VARS%sm1(kk))

                 endif
               
              enddo
           
            endif

          endif

        endif

      endif

    enddo

  endif

! ====================================================================
! Update zbar by taking the total flux and dividing by
! the average porosity just above the water table.
! ====================================================================

! --------------------------------------------------------------------&
! If the available porosity is nonzero divide the flux by its value.
! --------------------------------------------------------------------&

  if ( sum(CAT%smpsum).gt.(0.001d0)) then

    CAT%zbar1 = CAT%zbar - zbrflx/sum(CAT%smpsum)

  endif

  if ( sum(CAT%smpsum).le.(0.001d0)) then

    CAT%zbar1=CAT%zbar

  endif


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

!> Calculates the time step totals of evapotranspiration, runoff, 
!! surface energy fluxes and vertical soil-water fluxes.

subroutine sumflx_catchment(CAT,GRID_VARS,GLOBAL,&
               GRID_VEG,GRID_SOIL,GRID_MET,ilandc)

  implicit none
  type (GRID_VARS_template),intent(inout) :: GRID_VARS
  type (GRID_VEG_template),intent(in) :: GRID_VEG
  type (GRID_SOIL_template),intent(in) :: GRID_SOIL
  type (GRID_MET_template),intent(in) :: GRID_MET
  type (CATCHMENT_template),intent(inout) :: CAT
  type (GLOBAL_template),intent(in) :: GLOBAL
  integer,intent(in) :: ilandc
  real*8 :: ettot,etstore,etstsum,etwt,etbssum,etdcsum,etwtsum
  real*8 :: etwcsum,bsdew,contot,pptsum,pnetsum,qsurf
  real*8 :: sxrtot,xixtot,conrun,ranrun,gwt,gwtsum
  real*8 :: capsum,tzpsum,etpixloc,conpix,difwt
  real*8 :: xlhv,dummy,svarhs
  integer :: i,kk
  real*8 :: smpsum(GLOBAL%nlayer)


! ====================================================================
! Compute regional average evapotranspiration rate for the time step.    
! ====================================================================
 
  
  if (GRID_VEG%ivgtyp.eq.(-1)) then

! --------------------------------------------------------------------&
! In case of an open water pixel.
! --------------------------------------------------------------------&

    etpixloc=GRID_VARS%evtact
    ettot = ettot + etpixloc

  else

! --------------------------------------------------------------------&
! In case of vegetated or bare soil pixels.
! --------------------------------------------------------------------&

  if ( (GRID_VEG%i_und.eq.0)) then

! ....................................................................
! If understory/moss is not represented.
! ....................................................................

    etpixloc = GRID_VARS%evtact*GRID_VARS%dc*(one-GRID_VARS%fw) + GRID_VARS%epwms*GRID_VARS%dc
    ettot = ettot + etpixloc

  else

! ....................................................................
! Add the evapotranspiration rate from the understory/moss layer.
! ....................................................................

    xlhv=1000000000.*(2.501-0.002361*GRID_MET%tdry)
    dummy = GRID_VEG%canclos*GRID_VARS%xleact
    etpixloc = dummy/xlhv
    ettot = ettot + etpixloc

  endif

endif

! ====================================================================
! Split evapotranspiration by whether the water is supplied from
! below or above the water table.  Also sum evapotranspiration
! from each surface type.  Separate the evaporation from the
! water table by catchment for water table updating.
! ====================================================================

  if (GRID_VEG%ivgtyp.ge.0) then

    if (((GRID_VEG%ivgtyp.ne.2).and.(GRID_VARS%zrz.le.zero)).or.&
        ((GRID_VEG%ivgtyp.eq.2).and.(GRID_VARS%ztz.le.zero))) then

! --------------------------------------------------------------------&
! Add the evapotranspiration values up from pixels where the water
! is supplied from above the water table.
! --------------------------------------------------------------------&

      if (GRID_VEG%i_und.eq.0) then

! ....................................................................
! Understory is not represented.
! ....................................................................

        etstore = GRID_VARS%epwms*GRID_VARS%dc

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
    GRID_VARS%etpix=GRID_VARS%etpix+etpixloc
    etstsum = etstore
    etwtsum = etwt

    if (GRID_VEG%ivgtyp.eq.0) then

! ====================================================================
! Add up evaporation from bare soil values.
! ====================================================================

      etbssum = GRID_VARS%evtact

    else

! ====================================================================
! Add up evaporation from dry and wet canopy for vegetated pixels.
! ====================================================================

      if ( (GRID_VEG%i_und.eq.0)) then

! --------------------------------------------------------------------&
! Understory/moss is not represented.
! --------------------------------------------------------------------&

        etdcsum = GRID_VARS%evtact*GRID_VARS%dc*(one-GRID_VARS%fw)
        etwcsum = GRID_VARS%epwms*GRID_VARS%dc

      else

! --------------------------------------------------------------------&
! Add the values from the understory/moss in the totals.
! --------------------------------------------------------------------&

        etdcsum = GRID_VARS%evtact*GRID_VARS%dc*(one-GRID_VARS%fw)
        etwcsum = GRID_VARS%epwms*GRID_VARS%dc

      endif

    endif

  else

  endif

! ====================================================================
! Compute total, catchment and regional condensation.
! ====================================================================

  if (ilandc.ge.0) then

    if (GRID_VEG%ivgtyp.eq.0) then

! --------------------------------------------------------------------&
! In case of bare soil the condensation is the dew onto the soil.
! --------------------------------------------------------------------&

      conpix = bsdew 

    else

! --------------------------------------------------------------------&
! In case of vegetated pixels the condensation is the negative
! evaporation onto the wet canopy.
! --------------------------------------------------------------------&

      conpix = GRID_VARS%epwms*(one-GRID_VARS%dc)

    endif

    if (GRID_VEG%ivgtyp.eq.(-1)) then

      conpix=0.d0

    endif

  endif

  contot = conpix

! ====================================================================
! Compute total precipitation and surface runoff for the time step.
! ====================================================================

  pptsum = GRID_MET%pptms
  pnetsum = GRID_VARS%pnet
  qsurf = GRID_VARS%runtot

! ====================================================================
! Compute total saturation and infiltration excess runoff for the
! time step.
! ====================================================================

  sxrtot = GRID_VARS%satxr 
  xixtot = GRID_VARS%xinfxr

! ====================================================================
! Compute total runoff due to rainfall and due to condensation.
! ====================================================================

  if (GRID_MET%pptms.gt.zero) then

! --------------------------------------------------------------------&
! When under rainfall runoff has to be due to rain.
! --------------------------------------------------------------------&

    ranrun = GRID_VARS%runtot

  else

! --------------------------------------------------------------------&
! When no precipitation runoff has to be due to condensation.
! --------------------------------------------------------------------&

    conrun = GRID_VARS%runtot

  endif  

! ====================================================================
! Compute checks on canopy water balance, root zone water balance,&
! and transmission zone balance.
! ====================================================================

  if (GRID_VEG%ivgtyp.ge.(0)) then

! ====================================================================
! Compute drainage to the water table for time step.
! ====================================================================

    if (GRID_VARS%zrz.eq.zero) then
     
! --------------------------------------------------------------------&
! If the root zone is saturated there is no drainage towards the water 
! table.
! --------------------------------------------------------------------&

      gwt = zero
      difwt = GRID_VARS%difrz

    else if (GRID_VARS%ztz.eq.zero) then

! --------------------------------------------------------------------&
! If the transmission zone is saturated and the root zone is not
! the drainage and diffusion towards the water table comes from the
! root zone.
! --------------------------------------------------------------------&

      gwt = GRID_VARS%grz
      difwt = GRID_VARS%difrz

    else

! --------------------------------------------------------------------&
! If the transmission zone is not saturated the drainage and diffusion
! towards the water table comes from the  transmission zone.
! --------------------------------------------------------------------&
   
      gwt = GRID_VARS%gtz
      difwt = GRID_VARS%diftz

    endif

! ====================================================================
! Compute the regional totals of drainage to the water table and
! root and transmission zone drainage terms.
! ====================================================================

    gwtsum = gwt

! ====================================================================
! Compute the diffusion totals for the time step.
! ====================================================================

    capsum = - difwt

! ====================================================================
! Add up all the input terms towards the water table and make the
! regional total.
! ======================================================================


    svarhs = GLOBAL%dt*(GRID_MET%pptms - etstore - GRID_VARS%runtot - gwt - difwt)

! ====================================================================
! Compute average soil moistures for the end of the time step.
! ====================================================================

    if (GLOBAL%inc_frozen.eq.0) then

! --------------------------------------------------------------------&
! If frozen and liquid soil water are treated as liquid water than the
! liquid water content is the total water content.
! --------------------------------------------------------------------&

      do i=1,GLOBAL%nlayer
        GRID_VARS%sm1_u(i) = GRID_VARS%sm1(i)
      enddo

    endif

! ====================================================================
! Compute available porosity in root and transmission zones .or.&
! updating of average water table depth.  Only interested in
! the porosity for the zone above the water table.
! ====================================================================

    smpsum = zero

    do i = GLOBAL%nlayer,1,-1

      if (GRID_VARS%z_layer(i).gt.zero) then

        !This is the region of interest
        smpsum(i) = (GRID_SOIL%thetas-GRID_VARS%sm(i))
        exit

      endif

    enddo

  endif

!Update Catchment Sums
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
  CAT%gwtsum = CAT%gwtsum + gwtsum
  CAT%capsum = CAT%capsum + capsum
  CAT%smpsum = CAT%smpsum + smpsum

end subroutine sumflx_catchment

!>Subroutine to redistribute the catchment mean water table to a given grid cell as defined in Chaney et al.,2013
subroutine Redistribute_Zbar(n,ff,lambda,GSTI,zbar,zw)

  implicit none
  real*8,intent(in) :: n,ff,lambda,GSTI,zbar
  real*8,intent(inout) :: zw

  zw = 1/ff - 1/ff*GSTI/lambda + zbar*GSTI/lambda

end subroutine Redistribute_Zbar

!>Subroutine to create the generalized soils topographix index as defined in Chaney et al.,2013
subroutine Calculate_GSTI(GLOBAL,CAT,GRID)

  implicit none
  type (GLOBAL_template),intent(in) :: GLOBAL
  type (GRID_template),dimension(:),intent(inout) :: GRID
  type (CATCHMENT_template),dimension(:),intent(in) :: CAT
  integer :: ipix
  real*8 :: n,T0star,TI

  do ipix=1,GLOBAL%npix
    n = CAT(GRID(ipix)%VARS%icatch)%n
    T0star = GRID(ipix)%VARS%T0/n
    TI = GRID(ipix)%VARS%TI
    GRID(ipix)%VARS%GSTI = (dexp(TI)/T0star)**(one/n)
  enddo

end subroutine Calculate_GSTI

!>Subroutine to calculate the baseflow using the generalized topmodel defined in
!Chaney et al.,2013

subroutine Calculate_Qb(n,ff,zbar,q0,qb)

  implicit none
  real*8,intent(in) :: q0,ff,zbar,n
  real*8,intent(inout) :: qb

  qb = q0 * (1 - ff*zbar)**n

end subroutine Calculate_Qb

END MODULE MODULE_TOPMODEL
