! ====================================================================
!
!			subroutine rdveg
!
! ====================================================================
!
! Subroutine to read and initiailize simulation constant vegetation
!   and land cover parameters.
!
! ====================================================================

      subroutine rdveg(npix,nrow,ncol,ilandc,ipixnum,nlandc,iopveg,ivgtyp,iprn,&
       xlai,xlai_wsc,albd,albw,emiss,za,zww,z0m,z0h,zpd,rsmin,rsmax,Rpl,&
       f3vpdpar,f4temppar,trefk,tcbeta,extinct,canclos,Tslope1,Tint1,&
       Tslope2,Tint2,Twslope1,Twint1,Twslope2,Twint2,Tsep,Twsep,&
       eps,&
       rtact,rtdens,rtres,psicri,&
       rescan,respla,wsc,wcip1,&
       pixsiz,area,fbs,fbsrg,ncatch)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/rdveg.h"
      type VegDataTemplate
        real*8,dimension(:,:),allocatable :: ivgtyp
        real*8,dimension(:,:),allocatable :: xlai
        real*8,dimension(:,:),allocatable :: xlai_wsc
        real*8,dimension(:,:),allocatable :: albd
        real*8,dimension(:,:),allocatable :: albw
        real*8,dimension(:,:),allocatable :: emiss
        real*8,dimension(:,:),allocatable :: za
        real*8,dimension(:,:),allocatable :: zww
        real*8,dimension(:,:),allocatable :: z0m
        real*8,dimension(:,:),allocatable :: z0h
        real*8,dimension(:,:),allocatable :: zpd
        real*8,dimension(:,:),allocatable :: rsmin
        real*8,dimension(:,:),allocatable :: rsmax
        real*8,dimension(:,:),allocatable :: Rpl
        real*8,dimension(:,:),allocatable :: f3vpdpar
        real*8,dimension(:,:),allocatable :: f4temppar
        real*8,dimension(:,:),allocatable :: trefk
        real*8,dimension(:,:),allocatable :: tcbeta
        real*8,dimension(:,:),allocatable :: extinct
        real*8,dimension(:,:),allocatable :: canclos
      end type VegDataTemplate
      type (VegDataTemplate) VegData
      character(len=200) :: filename
      integer :: vegnvars,dvegnvars,ipos,jpos
      real,dimension(:,:,:),allocatable :: TempArray
      vegnvars = 20
      dvegnvars = 2
      allocate(TempArray(ncol,nrow,vegnvars))

! ====================================================================
! Read the image with the land cover clasifications.
! ====================================================================

      call rdimgi(ilandc,11,nrow,ncol,ipixnum)

! ====================================================================
! Read spatially constant vegetation parameters.
! ====================================================================

      nlandc = nrow*ncol
      iopveg = 0
      read(1000,*) iopwc0

      print*,"rdveg:  Read spatially constant veg pars"

! ====================================================================
! Read lookup table with parameters for vegetation/land cover.
! Read either critical and wilting soil moistures or
! plant/root resistance parameter depending on actual
! transpiration option.
! ====================================================================

      print*,"rdveg:  Read lookup table"

! ====================================================================
! Read the binary vegetation binary file
! ====================================================================

      allocate(VegData%ivgtyp(ncol,nrow))
      allocate(VegData%xlai(ncol,nrow))
      allocate(VegData%xlai_wsc(ncol,nrow))
      allocate(VegData%albd(ncol,nrow))
      allocate(VegData%albw(ncol,nrow))
      allocate(VegData%emiss(ncol,nrow))
      allocate(VegData%za(ncol,nrow))
      allocate(VegData%zww(ncol,nrow))
      allocate(VegData%z0m(ncol,nrow))
      allocate(VegData%z0h(ncol,nrow))
      allocate(VegData%zpd(ncol,nrow))
      allocate(VegData%rsmin(ncol,nrow))
      allocate(VegData%rsmax(ncol,nrow))
      allocate(VegData%Rpl(ncol,nrow))
      allocate(VegData%f3vpdpar(ncol,nrow))
      allocate(VegData%f4temppar(ncol,nrow))
      allocate(VegData%trefk(ncol,nrow))
      allocate(VegData%tcbeta(ncol,nrow))
      allocate(VegData%extinct(ncol,nrow))
      allocate(VegData%canclos(ncol,nrow))
        
      print*,"rdveg:  Reading in all the vegetation properties at once"

      read(1002,rec=1)TempArray(:,:,:)

      VegData%ivgtyp(:,:) = dble(TempArray(:,:,1))
      VegData%xlai(:,:) = dble(TempArray(:,:,2))
      VegData%xlai_wsc(:,:) = dble(TempArray(:,:,3))
      VegData%albd(:,:) = dble(TempArray(:,:,4))
      VegData%albw(:,:) = dble(TempArray(:,:,5))
      VegData%emiss(:,:) = dble(TempArray(:,:,6))
      VegData%za(:,:) = dble(TempArray(:,:,7))
      VegData%zww(:,:) = dble(TempArray(:,:,8))
      VegData%z0m(:,:) = dble(TempArray(:,:,9))
      VegData%z0h(:,:) = dble(TempArray(:,:,10))
      VegData%zpd(:,:) = dble(TempArray(:,:,11))
      VegData%rsmin(:,:) = dble(TempArray(:,:,12))
      VegData%rsmax(:,:) = dble(TempArray(:,:,13))
      VegData%Rpl(:,:) = dble(TempArray(:,:,14))
      VegData%f3vpdpar(:,:) = dble(TempArray(:,:,15))
      VegData%f4temppar(:,:) = dble(TempArray(:,:,16))
      VegData%trefk(:,:) = dble(TempArray(:,:,17))
      VegData%tcbeta(:,:) = dble(TempArray(:,:,18))
      VegData%extinct(:,:) = dble(TempArray(:,:,19))
      VegData%canclos(:,:) = dble(TempArray(:,:,20))

! ####################################################################
! Read the vegetation dynamic parameter file
! ####################################################################

      deallocate(TempArray)
      allocate(TempArray(ncol,nrow,dvegnvars))

      print*,"rdveg:  Reading in the dynamic vegetation properties"

      read(1003,rec=1)TempArray(:,:,:)

      VegData%xlai(:,:) = dble(TempArray(:,:,1))
      VegData%albd(:,:) = dble(TempArray(:,:,2))


! ####################################################################
! Convert the 2-d arrays to the model's 1-d arrays
! ####################################################################

        do kk=1,nlandc

                !Map the kk position to the i,j position
                if(mod(kk,nrow) .ne. 0)then
                        ipos = kk/nrow+1
                        jpos = mod(kk,nrow)
                else
                        ipos = kk/nrow
                        jpos = nrow
                endif

                ivgtyp(kk) = int(VegData%ivgtyp(ipos,jpos))
                emiss(kk) = VegData%emiss(ipos,jpos)
                za(kk) = VegData%za(ipos,jpos)
                zww(kk) = VegData%zww(ipos,jpos)
                z0m(kk) = VegData%z0m(ipos,jpos)
                z0h(kk) = VegData%z0h(ipos,jpos)
                zpd(kk) = VegData%zpd(ipos,jpos)
                rsmin(kk) = VegData%rsmin(ipos,jpos)
                rsmax(kk) = VegData%rsmax(ipos,jpos)
                Rpl(kk) = VegData%Rpl(ipos,jpos)
                f3vpdpar(kk) = VegData%f3vpdpar(ipos,jpos)
                f4temppar(kk) = VegData%f4temppar(ipos,jpos)
                trefk(kk) = VegData%trefk(ipos,jpos)
                xlai(kk) = VegData%xlai(ipos,jpos) !dveg
                albd(kk) = VegData%albd(ipos,jpos) !dveg
                tcbeta(kk) = exp(-0.5*xlai(kk)) !dveg
                xlai_wsc(kk) = xlai(kk)
                albw(kk) = albd(kk) !Move to its own file
                extinct(kk) = 0.00!VegData%extinct(ipos,jpos)
                canclos(kk) = 1.00!VegData%canclos(ipos,jpos)
                Tslope1(kk) = 0.00!VegData%Tslope1(ipos,jpos)
                Tint1(kk) = 0.00!VegData%Tint1(ipos,jpos)
                Tslope2(kk) = 0.00!VegData%Tslope2(ipos,jpos)
                Tint2(kk) = 0.00!VegData%Tint2(ipos,jpos)
                Twslope1(kk) = 0.00!VegData%Twslope1(ipos,jpos)
                Twint1(kk) = 0.00!VegData%Twint1(ipos,jpos)
                Twslope2(kk) = 0.00!VegData%Twslope2(ipos,jpos)
                Twint2(kk) = 0.00!VegData%Twint2(ipos,jpos)
                Tsep(kk) = 0.00!VegData%Tsep(ipos,jpos)
                Twsep(kk) = 0.00!VegData%Twsep(ipos,jpos)

        enddo

! ====================================================================
! Calculate parameters for each land cover type.
! ====================================================================

      do kk=1,nlandc

! --------------------------------------------------------------------&
! If not bare soil then calculate the canopy resistance.
! --------------------------------------------------------------------&

         if (ivgtyp(kk).ne.0) then

            rescan(kk) = rsmin(kk)/xlai(kk)

! --------------------------------------------------------------------&
! If bare soil then set canopy resistance to zero.
! --------------------------------------------------------------------&

         else

            rescan(kk) = 0.

         endif

! --------------------------------------------------------------------&
! Calculate canopy storage capacity and initial canopy storage.
! --------------------------------------------------------------------&

         wsc(kk) = 0.0002*xlai_wsc(kk)

      enddo

      print*,"rdveg:  Set minimum st. resist. and and wet can.stor.cap."

! ====================================================================
! Read the initial canopy storage from an image or as a constant
! as requested.
! ====================================================================

      if (iopwc0.eq.0) then

         read(1000,*) wc0

         do kk=1,npix

            wcip1(kk) = wc0

         enddo  

         print*,"rdveg:  Read initial wet canopy storages"

      else 
		 
         call rdimgr(wcip1,13,nrow,ncol,ipixnum)

      endif

! ====================================================================
! Calculate the fraction of different cover types in each catchment
! and in total area.  First array index for icount and frcov is
! the land cover type, second array index is the catchment number.
! Catchment ncatch+1 is total area.
! ====================================================================

      do 550 kk=1,nlandc

         do 540 jj=1,ncatch

            frcov(kk,jj) = pixsiz*pixsiz/area(jj)

540      continue

         frcov(kk,ncatch+1) = 1/real(npix)

550   continue

      print*,"rdveg:  Calculated fractional covers for vegetation"

! ====================================================================
! Find fraction of bare soil in each catchment.
! ====================================================================

      fbsrg = zero

      do 570 jj=1,ncatch+1

         fbs(jj)  = zero

         do 560 kk=1,nlandc

            if (ivgtyp(kk).eq.0) then
               
               if (jj .eq. ncatch+1) then
               
                  fbsrg = fbsrg + frcov(kk,jj)
        
               else

                  fbs(jj) = fbs(jj) + frcov(kk,jj)

               endif

            endif

560      continue

570   continue

      print*,"rdveg:  Calculated fractional covers for bare soil"

! ====================================================================
! Print land cover table summary.
! ====================================================================

      if (iprn(76).eq.1) then

         do 600 kk = 1,nlandc

            write(76,1000) kk,xlai(kk),rescan(kk),wsc(kk)

600      continue

      endif

! ====================================================================
! Print land cover fractions in each catchment.
! ====================================================================

      if (iprn(77).eq.1) then

         do 700 jj=1,ncatch

            write (77,1100) jj,area(jj)/1000000.

            do 650 kk=1,nlandc

               write(77,1110) kk,frcov(kk,jj)*100

650         continue

            write (77,*)

700      continue

         write(77,1120)

         do 800 kk=1,nlandc

            write(77,1110) kk,frcov(kk,ncatch+1)*100

800      continue

      endif        

! ====================================================================
! Format statements.
! ====================================================================

1000  format (i5,2f10.3,f11.7)
1100  format ('Catchment Number',i4,' (Area = ',f8.2,' km^2):')
1110  format ('   Land Cover Type',i4,':',f7.2,'%')
1120  format ('Total Area:')

      return

      end
