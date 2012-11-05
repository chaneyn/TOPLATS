! ====================================================================
!
!			subroutine rdsoil
!
! ====================================================================
!
! Subroutine to read and initialize time in-variant soil parameters
!
! ====================================================================

      subroutine rdsoil(nsoil,irestype,ikopt,zrzmax,iopsmini,&
       smpet0,isoil,nrow,ncol,ipixnum,bcbeta,psic,thetas,thetar,xk0,zdeep,&
       tdeep,zmid,tmid0,rocpsoil,quartz,ifcoarse,srespar1,srespar2,srespar3,&
       a_ice,b_ice,bulk_dens,amp,phase,shift,inc_frozen,bcgamm,par,corr,&
       idifind,ncatch,icatch,pixsiz,area,npix,psicav,iprn,tc,tw)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/rdsoil.h"
      type SoilDataTemplate
        real*8,dimension(:,:),allocatable :: bcbeta
        real*8,dimension(:,:),allocatable :: psic
        real*8,dimension(:,:),allocatable :: thetas
        real*8,dimension(:,:),allocatable :: thetar
        real*8,dimension(:,:),allocatable :: xk0
        real*8,dimension(:,:),allocatable :: zdeep
        real*8,dimension(:,:),allocatable :: tdeep
        real*8,dimension(:,:),allocatable :: zmid
        real*8,dimension(:,:),allocatable :: tmid0
        real*8,dimension(:,:),allocatable :: rocpsoil
        real*8,dimension(:,:),allocatable :: quartz
        integer,dimension(:,:),allocatable :: ifcoarse
        real*8,dimension(:,:),allocatable :: srespar1
        real*8,dimension(:,:),allocatable :: srespar2
        real*8,dimension(:,:),allocatable :: srespar3
        real*8,dimension(:,:),allocatable :: a_ice
        real*8,dimension(:,:),allocatable :: b_ice
        real*8,dimension(:,:),allocatable :: bulk_dens
        real*8,dimension(:,:),allocatable :: amp
        real*8,dimension(:,:),allocatable :: phase
        real*8,dimension(:,:),allocatable :: shift
        real*8,dimension(:,:),allocatable :: thetaw
        real*8,dimension(:,:),allocatable :: thetac
      end type SoilDataTemplate
      type (SoilDataTemplate) SoilData
      integer :: soilnvars,ipos,jpos
      real,dimension(:,:,:),allocatable :: TempArray
      soilnvars = 23
      allocate(TempArray(ncol,nrow,soilnvars))

! ====================================================================
! Read spatially constant bare soil parameters and options.
! Then read root and transmission zone data.
! ====================================================================

      !read(1000,*)nsoil
      nsoil = nrow*ncol

      read(1000,*)irestype
      read(1000,*)ikopt
      read(1000,*)zrzmax
      read(1000,*)iopsmini

      if (iopsmini.eq.0) read(1000,*)smpet0

      print*,"rdsoil:  Read spatially constant soil pars"

      if (iopsmini.eq.1)&
         print*,"rdsoil:  Will read initial soil moisture images"

! ====================================================================
! Read the binary soil file
! ====================================================================

      allocate(SoilData%bcbeta(ncol,nrow))
      allocate(SoilData%psic(ncol,nrow))
      allocate(SoilData%thetas(ncol,nrow))
      allocate(SoilData%thetar(ncol,nrow))
      allocate(SoilData%xk0(ncol,nrow))
      allocate(SoilData%zdeep(ncol,nrow))
      allocate(SoilData%tdeep(ncol,nrow))
      allocate(SoilData%zmid(ncol,nrow))
      allocate(SoilData%tmid0(ncol,nrow))
      allocate(SoilData%rocpsoil(ncol,nrow))
      allocate(SoilData%quartz(ncol,nrow))
      allocate(SoilData%ifcoarse(ncol,nrow))
      allocate(SoilData%srespar1(ncol,nrow))
      allocate(SoilData%srespar2(ncol,nrow))
      allocate(SoilData%srespar3(ncol,nrow))
      allocate(SoilData%a_ice(ncol,nrow))
      allocate(SoilData%b_ice(ncol,nrow))
      allocate(SoilData%bulk_dens(ncol,nrow))
      allocate(SoilData%amp(ncol,nrow))
      allocate(SoilData%phase(ncol,nrow))
      allocate(SoilData%shift(ncol,nrow))
      allocate(SoilData%thetaw(ncol,nrow))
      allocate(SoilData%thetac(ncol,nrow))

      print*,"rdsoil:  Reading in all soil properties at once"
      read(1001,rec=1)TempArray(:,:,:)
      SoilData%bcbeta(:,:) = dble(TempArray(:,:,1))
      SoilData%psic(:,:) = dble(TempArray(:,:,2))
      SoilData%thetas(:,:) = dble(TempArray(:,:,3))
      SoilData%thetar(:,:) = dble(TempArray(:,:,4))
      SoilData%xk0(:,:) = dble(TempArray(:,:,5))
      SoilData%zdeep(:,:) = dble(TempArray(:,:,6))
      SoilData%tdeep(:,:) = dble(TempArray(:,:,7))
      SoilData%zmid(:,:) = dble(TempArray(:,:,8))
      SoilData%tmid0(:,:) = dble(TempArray(:,:,9))
      SoilData%rocpsoil(:,:) = dble(TempArray(:,:,10))
      SoilData%quartz(:,:) = dble(TempArray(:,:,11))
      SoilData%ifcoarse(:,:) = int(TempArray(:,:,12))
      SoilData%srespar1(:,:) = dble(TempArray(:,:,13))
      SoilData%srespar2(:,:) = dble(TempArray(:,:,14))
      SoilData%srespar3(:,:) = dble(TempArray(:,:,15))
      SoilData%a_ice(:,:) = dble(TempArray(:,:,16))
      SoilData%b_ice(:,:) = dble(TempArray(:,:,17))
      SoilData%bulk_dens(:,:) = dble(TempArray(:,:,18))
      SoilData%amp(:,:) = dble(TempArray(:,:,19))
      SoilData%phase(:,:) = dble(TempArray(:,:,20))
      SoilData%shift(:,:) = dble(TempArray(:,:,21))
      SoilData%thetaw(:,:) = dble(TempArray(:,:,22))
      SoilData%thetac(:,:) = dble(TempArray(:,:,23))

! ====================================================================
! Read the soil classification image.
! ====================================================================

      call rdimgi(isoil,12,nrow,ncol,ipixnum)

      print*,"rdsoil:  Read soil texture image"

! ====================================================================
! Pass the soil properties from the original i,j pos. to the kk pos.
!  ====================================================================
      do kk=1,nsoil

         !Map the kk position to the i,j position
         if(mod(kk,nrow) .ne. 0)then
                ipos = kk/nrow+1
                jpos = mod(kk,nrow)
         else
                ipos = kk/nrow
                jpos = nrow
         endif

         bcbeta(kk) = SoilData%bcbeta(ipos,jpos)
         psic(kk) = SoilData%psic(ipos,jpos)
         thetas(kk) = SoilData%thetas(ipos,jpos)
         thetar(kk) = SoilData%thetar(ipos,jpos)
         xk0(kk) = SoilData%xk0(ipos,jpos)
         zdeep(kk) = SoilData%zdeep(ipos,jpos)
         tdeep(kk) = SoilData%tdeep(ipos,jpos)
         zmid(kk) = SoilData%zmid(ipos,jpos)
         tmid0(kk) = SoilData%tmid0(ipos,jpos)
         rocpsoil(kk) = SoilData%rocpsoil(ipos,jpos)
         quartz(kk) = SoilData%quartz(ipos,jpos)
         ifcoarse(kk) = SoilData%ifcoarse(ipos,jpos)
         srespar1(kk) = SoilData%srespar1(ipos,jpos)
         srespar2(kk) = SoilData%srespar2(ipos,jpos)
         srespar3(kk) = SoilData%srespar3(ipos,jpos)
         a_ice(kk) = SoilData%a_ice(ipos,jpos)
         b_ice(kk) = SoilData%b_ice(ipos,jpos)
         bulk_dens(kk) = SoilData%bulk_dens(ipos,jpos)
         amp(kk) = SoilData%amp(ipos,jpos)
         phase(kk) = SoilData%phase(ipos,jpos)
         shift(kk) = SoilData%shift(ipos,jpos)
         tc(kk) = SoilData%thetac(ipos,jpos)
         tw(kk) = SoilData%thetaw(ipos,jpos)

      enddo
 
      inc_frozen = 1 !THIS MEANS THAT THE FROZEN ALGORITHM IS ALWAYS RUN

! ====================================================================
! Calculate time in-variant soil parameters for each soil class.
! ====================================================================

      do 400 kk=1,nsoil

! --------------------------------------------------------------------&
! Calculate soil parameters based on Brooks-Corey parameters.
! --------------------------------------------------------------------&

         bcgamm(kk) = two + three * bcbeta(kk)

! --------------------------------------------------------------------&
! Calculate constants for bare soil evaporation desorptivity 
! equation used in Famiglietti PhD Thesis, Princetion Univ, 1992.
! --------------------------------------------------------------------&

         par(kk) = one + ((bcgamm(kk)-one)/bcbeta(kk))
         corr(kk)=((one/(one+three*bcbeta(kk)))-&
                   (0.85d0/(one+four*bcbeta(kk)))-&
                   (0.85d0*0.15d0*0.5d0/(one+five*bcbeta(kk)))+&
                   (0.85d0*0.15d0*1.15d0/&
                   (six*(one+six*bcbeta(kk)))))

! --------------------------------------------------------------------&
! Calculate diffusivity index and dimensionless exfiltration
! diffusivity from Eagleson, WRR, 1972.
! --------------------------------------------------------------------&

         idifind(kk) = ((1.0+2.0*bcbeta(kk))/bcbeta(kk))+0.5
         tempsum=0 

         do 300 nn=1,idifind(kk)

            dtaken = exp(factln(idifind(kk))-factln(nn)-&
                       factln(idifind(kk)-nn))
            tempsum = tempsum+(((-1)**nn)*dtaken/(1.85+nn))

300      continue

400   continue

      print*,"rdsoil:  Calculated time-invariant soil pars"

! ====================================================================
! Calculate the fraction of different soil types in each catchment
! and in total area.  First array index for icount and frcov is
! the soil type, second array index is the catchment number.
! Catchment ncatch+1 is total area.
! ====================================================================

      do 450 kk=1,nsoil

         do 440 jj=1,ncatch+1

            icount(kk,jj) = 0

440      continue

450   continue

      do 500 kk=1,npix

        icount(isoil(kk),icatch(kk))=icount(isoil(kk),icatch(kk))+1
        icount(isoil(kk),ncatch+1) = icount(isoil(kk),ncatch+1) + 1

500   continue

      do 550 kk=1,nsoil

         do 540 jj=1,ncatch

            frsoil(kk,jj) = icount(kk,jj)*pixsiz*pixsiz/area(jj)

540      continue

         frsoil(kk,ncatch+1) = icount(kk,ncatch+1)/real(npix)

550   continue

      print*,"rdsoil:  Calculated fractional coverage for soil types"

! ====================================================================
! Calculate average bubbling pressure in each catchment (used
! in updating average water table depths.
! ====================================================================

      do 570 jj=1,ncatch

         psicav(jj) = zero

         do 560 kk=1,nsoil

            psicav(jj) = psicav(jj) + frsoil(kk,jj)*psic(kk)

560      continue

570   continue

      print*,"rdsoil:  Calculated average psi! for each catchment"

! ====================================================================
! Print summary table.
! ====================================================================

      if (iprn(78).eq.1) then

         do 580 kk=1,nsoil

            write(78,1000) kk,bcbeta(kk),psic(kk)*100,thetas(kk),&
                           thetar(kk),xk0(kk)*3600000.

580      continue

      endif

! ====================================================================
! Print fractional coverage of soils types in each catchment.
! ====================================================================

      if (iprn(79).eq.1) then

         do 700 jj=1,ncatch

            write (79,1100) jj,area(jj)/1000000.

            do 650 kk=1,nsoil

               write(79,1110) kk,frsoil(kk,jj)*100

650         continue

            write(79,*)

700      continue

         write(79,1120)

         do 800 kk=1,nsoil

            write(79,1110) kk,frsoil(kk,ncatch+1)*100

800      continue

      endif        

      print*,"rdsoil:  Printed summary tables"

! ====================================================================
! Format statement.
! ====================================================================

1000  format (i5,5f10.3)
1100  format ('Catchment Number',i4,' (Area = ',f8.2,' km^2):')
1110  format ('   Soil Type',i4,':',f7.2,'%')
1120  format ('Total Area:')

      return

      end
