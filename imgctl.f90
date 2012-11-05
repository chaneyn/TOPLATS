! ====================================================================
! 
!			subroutine imgctl
!
! ====================================================================
!
! Subroutine to control output of time varying images by deciding
!   what images to print each time step, openning the output
!   file and writing the image
!
! ====================================================================

      subroutine imgctl(tkact,gact,hact,xleact,rnact,etpix,&
       runtot,xinact,irntyp,zw,tzsm1,rzsm1,tkpet,wcip1,gpet,hpet,xlepet,rnpet,&
       r_mossm,pptms,ievcon,ipixnum,icurser,nseries,ioutst,iouten,iprn,ioutsp,&
       pptsumrg,i,ncol,nrow,fnimg,img_opt,iyear,iday,ihour,Swq,Swq_us)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/imgctl.h"

! ====================================================================
! Loop through each of the possible output image unit numbers.
! ====================================================================

      do 100 iu=35,55

! ====================================================================
! Set the default value of the flag to signify image output
! for this time step to not print (zero).
! ====================================================================

         imgprn(iu) = 0

! ====================================================================
! Check if any images of this type are to be printed.
! ====================================================================

         if (iprn(iu).eq.1) then

! ====================================================================
! Check if this time step has image output.  If so, then 
! open the image file and reset the flag to print the 
! image.
! ====================================================================

            if ( (iouten(iu,icurser(iu)).lt.i).and. &
                 (nseries(iu).gt.icurser(iu)) ) icurser(iu) = icurser(iu) + 1 

            if ( (mod(i-ioutst(iu,icurser(iu)),ioutsp(iu,icurser(iu))).eq.0)&
               .and.&
                 (i.ge.ioutst(iu,icurser(iu)))&
               .and.&
                 (i.le.iouten(iu,icurser(iu))) ) then

                   call imgopn(fnimg(iu),i,iu,img_opt,iyear,iday,ihour)
                   imgprn(iu) = 1

            endif

         endif

100   continue

      do 200 iu=151,199

! ====================================================================
! Set the default value of the flag to signify image output
! for this time step to not print (zero).
! ====================================================================

         imgprn(iu) = 0

! ====================================================================
! Check if any images of this type are to be printed.
! ====================================================================

         if (iprn(iu).eq.1) then

! ====================================================================
! Check if this time step has image output.  If so, then 
! open the image file and reset the flag to print the 
! image.
! ====================================================================

            if ( (iouten(iu,icurser(iu)).lt.i).and. &
                 (nseries(iu).gt.icurser(iu)) ) icurser(iu) = icurser(iu) + 1 

            if ( (mod(i-ioutst(iu,icurser(iu)),ioutsp(iu,icurser(iu))).eq.0)&
               .and.&
                 (i.ge.ioutst(iu,icurser(iu)))&
               .and.&
                 (i.le.iouten(iu,icurser(iu))) ) then

                   call imgopn(fnimg(iu),i,iu,img_opt,iyear,iday,ihour)
                   imgprn(iu) = 1

            endif

         endif

200   continue

! ====================================================================
! Print the requested images.
! ====================================================================

      if (imgprn(35).eq.1) then

         if (pptsumrg.gt.(0.d0))&
             call wrimgr(pptms,35,3600000.,nrow,ncol,ipixnum)

      endif

      if (imgprn(36).eq.1) call wrimgr(r_mossm,36,100.,nrow,ncol,ipixnum)
      if (imgprn(37).eq.1) call wrimgr(rnpet,37,1.,nrow,ncol,ipixnum)
      if (imgprn(38).eq.1) call wrimgr(xlepet,38,1.,nrow,ncol,ipixnum)
      if (imgprn(39).eq.1) call wrimgr(hpet,39,1.,nrow,ncol,ipixnum)
      if (imgprn(40).eq.1) call wrimgr(gpet,40,1.,nrow,ncol,ipixnum)
      if (imgprn(41).eq.1) call wrimgr(tkpet,41,1.,nrow,ncol,ipixnum)
      if (imgprn(42).eq.1) call wrimgr(wcip1,42,1000000.,nrow,ncol,ipixnum)
      if (imgprn(43).eq.1) call wrimgr(rzsm1,43,100.,nrow,ncol,ipixnum)
      if (imgprn(43).eq.1) call WRITE_BINARY(rzsm1,100.,nrow,ncol,ipixnum,i)
      if (imgprn(44).eq.1) call wrimgr(tzsm1,44,100.,nrow,ncol,ipixnum)
      if (imgprn(45).eq.1) call wrimgr(zw,45,100.,nrow,ncol,ipixnum)

! ====================================================================
! Close any open image files.
! ====================================================================

      do 203 iu=35,45

         if (imgprn(iu).eq.1) close(iu)

203   continue

      if (imgprn(46).eq.1) call wrimgr(xinact,46,3600000.,nrow,ncol,ipixnum)
      if (imgprn(47).eq.1) call wrimgr(runtot,47,3600000.,nrow,ncol,ipixnum)
      if (imgprn(48).eq.1) call wrimgi(irntyp,48,1,nrow,ncol,ipixnum)
      if (imgprn(49).eq.1) call wrimgr(etpix,49,3600000.,nrow,ncol,ipixnum)
      if (imgprn(50).eq.1) call wrimgi(ievcon,50,1,nrow,ncol,ipixnum)
      if (imgprn(51).eq.1) call wrimgr(rnact,51,1.,nrow,ncol,ipixnum)
      if (imgprn(52).eq.1) call wrimgr(xleact,52,1.,nrow,ncol,ipixnum)
      if (imgprn(53).eq.1) call wrimgr(hact,53,1.,nrow,ncol,ipixnum)
      if (imgprn(54).eq.1) call wrimgr(gact,54,1.,nrow,ncol,ipixnum)
      if (imgprn(55).eq.1) call wrimgr(tkact,55,1.,nrow,ncol,ipixnum)
      if (imgprn(151).eq.1) call wrimgr(Swq,151,1.,nrow,ncol,ipixnum)
      if (imgprn(152).eq.1) call wrimgr(Swq_us,152,1.,nrow,ncol,ipixnum)

! ====================================================================
! Close any open image files.
! ====================================================================

      do 201 iu=46,55

         if (imgprn(iu).eq.1) close(iu)

201   continue

      do 202 iu=151,199

         if (imgprn(iu).eq.1) close(iu)

202   continue

      return

      end
