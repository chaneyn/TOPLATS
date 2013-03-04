MODULE MODULE_CELL

USE MODULE_VARIABLES

USE MODULE_LAND

USE MODULE_ATMOS

USE MODULE_CANOPY

USE MODULE_SNOW

contains

!#####################################################################
!
!                        subroutine Update_Cells
!
!#####################################################################
!
! Solve the water and energy budget for all land surface area
!
!#####################################################################

  subroutine Update_Cells(GRID,CAT,GLOBAL,i)
  
    implicit none
    type (GRID_template),dimension(:),intent(inout) :: GRID
    type (CATCHMENT_template),dimension(:),intent(inout) :: CAT
    type (GLOBAL_template),intent(inout) :: GLOBAL
    integer,intent(in) :: i
    integer :: ipix,isoil,icatch,ilandc
    real*8 :: omp_get_wtime,start_time,end_time
    type (GRID_VEG_template) :: GRID_VEG
    type (GRID_SOIL_template) :: GRID_SOIL
    type (GRID_MET_template) :: GRID_MET
    type (GRID_VARS_template) :: GRID_VARS
    type (CATCHMENT_template) :: CAT_INFO
    type (GLOBAL_template) :: GLOBAL_INFO

!#####################################################################
! Update each grid cell
!#####################################################################

    call OMP_SET_NUM_THREADS(GLOBAL%nthreads)
    GLOBAL%mul_fac = 1.0d0
    start_time = omp_get_wtime()

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(isoil,ilandc,icatch,GRID_VEG,&
!$OMP GRID_SOIL,GRID_MET,GRID_VARS,CAT_INFO,GLOBAL_INFO) 

    do ipix=1,GLOBAL%npix

      !Extract data
      isoil = GRID(ipix)%SOIL%isoil
      ilandc = GRID(ipix)%VEG%ilandc
      icatch = GRID(ipix)%VARS%icatch  
      GRID_MET = GRID(ipix)%MET
      GRID_SOIL = GRID(isoil)%SOIL
      GRID_VEG = GRID(ilandc)%VEG
      GRID_VARS = GRID(ipix)%VARS
      CAT_INFO = CAT(icatch)
      GLOBAL_INFO = GLOBAL
      !Update the cell
      call Update_Cell(ipix,i,GRID_MET,GRID_SOIL,&
         GRID_VEG,GRID_VARS,&
         CAT_INFO,GLOBAL_INFO)
      !Output the soil moisture in zrzmax
      GRID_VARS%rzsm_zrzmax = (GRID_VARS%rzsm*GRID_VARS%zrz + &
        GRID_SOIL%thetas*(GLOBAL_INFO%zrzmax-GRID_VARS%zrz))/GLOBAL_INFO%zrzmax

      !Write back
      GRID(isoil)%SOIL = GRID_SOIL
      GRID(ipix)%VARS = GRID_VARS

    enddo

!$OMP END PARALLEL DO
    end_time = omp_get_wtime()
!    print*,end_time - start_time

  end subroutine Update_Cells


! ====================================================================
!
!                        subroutine Update_Cell
!
! ====================================================================
!
! Solve the water and energy budget for a land surface area
!
! ====================================================================


      subroutine Update_Cell(ipix,i,GRID_MET,GRID_SOIL,GRID_VEG,&
               GRID_VARS,CAT,GLOBAL)

      implicit none
      type (GRID_VEG_template) :: GRID_VEG
      type (GRID_SOIL_template) :: GRID_SOIL
      type (GRID_MET_template) :: GRID_MET
      type (GRID_VARS_template) :: GRID_VARS
      type (CATCHMENT_template) :: CAT
      type (GLOBAL_template) :: GLOBAL
      type (CELL_VARS_template) :: CELL_VARS
      integer ipix,i
      real*8 snow,rain,rrr      
      
!Point Data Initializations
!Water Balance
GRID_VARS%zrz = 0.d0
GRID_VARS%ztz = 0.d0
GRID_VARS%smold = 0.d0
GRID_VARS%rzsmold = 0.d0
GRID_VARS%tzsmold = 0.d0
GRID_VARS%capflx = 0.d0
GRID_VARS%difrz = 0.d0
GRID_VARS%diftz = 0.d0
GRID_VARS%grz = 0.d0
GRID_VARS%gtz = 0.d0
GRID_VARS%satxr = 0.d0
GRID_VARS%xinfxr = 0.d0
GRID_VARS%dc = 0.d0!d
GRID_VARS%fw = 0.d0!fw
GRID_VARS%dsrz = 0.d0!dsrz
GRID_VARS%rzrhs = 0.d0!rzrhs
GRID_VARS%dstz = 0.d0!dstz
GRID_VARS%tzrhs = 0.d0!tzrhs
GRID_VARS%dswc = 0.d0!dswc
GRID_VARS%wcrhs = 0.d0!wcrhs
!Energy Fluxes
GRID_VARS%epwms = 0.d0!epwms


! ====================================================================
! If the vegetation type is greater than or equal to zero then
! solve the water and energy balance for a land area.
! ====================================================================

      if (GRID_VEG%ivgtyp.ge.0) then

! ....................................................................
! Calculate the local energy fluxes and set
! up the storm/interstorm event times and flags.
! ..................................................................

       call atmos(ipix,i,GRID_VEG,CELL_VARS,GRID_MET,GRID_VARS,GRID_SOIL,GLOBAL)

! ....................................................................
! Calculate local wet canopy water balance.
! ....................................................................

 
        call canopy(ipix,GRID_VARS,GRID_VEG,GRID_MET,GLOBAL)

! ....................................................................
! Calculate the local land surface water/energy balance.
! ....................................................................

! ....................................................................
! Option 2 : the incoming long wave radiation for both under and over
! story is equal and is the atmospheri! incoming long wave radiation.
! The uncouples the radiation balances for both layers from each
! other.  This option is also used when under story is not represented.
! ....................................................................

       call land(ipix,i,CELL_VARS,GRID_MET,GRID_VEG,GRID_VARS,GRID_SOIL,CAT,GLOBAL)

! ====================================================================
! Calculate the density and depth of the snow layers.
! ====================================================================

         if ( (GRID_VARS%Swq.gt.0.d0) ) then

            call calcrain (CELL_VARS%tcel,snow,rain,GRID_VARS%precip_o,GLOBAL%dt)
            call snow_density(GRID_VARS%dsty,snow,CELL_VARS%tcel,GRID_VARS%Swq,GRID_VARS%Sdepth,GRID_VARS%TSurf,GLOBAL%dt)

         else

           GRID_VARS%Sdepth=0.d0
           GRID_VARS%dsty=0.d0

         endif

      endif

! ====================================================================
! In the vegetation type is lower than zero then solve the open 
! water energy and water balance.
! ====================================================================

      if (GRID_VEG%ivgtyp.eq.(-1)) then

! ....................................................................
! Calculate the deep soil temperature.
! ....................................................................

         if ( (GRID_SOIL%amp.eq.(0.d0)).and.&
              (GRID_SOIL%phase.eq.(0.d0)).and.&
              (GRID_SOIL%shift.eq.(0.d0)) ) then

           GRID_SOIL%Tdeepstep=GRID_SOIL%tdeep

         else

            rrr=real(i)

            GRID_SOIL%Tdeepstep=GRID_SOIL%tdeep + GRID_SOIL%amp*cos ( rrr*GRID_SOIL%phase - GRID_SOIL%shift )

         endif

      endif

      return

      end subroutine Update_Cell

END MODULE MODULE_CELL
