! ====================================================================
!
!			subroutine wrimgr
!
! ====================================================================
!
! Subroutine to write an image of 4 byte reals from an array
!   of these values indexed by the soils-topographi! index
!   pixel numbers.
!
! ====================================================================
!
!  Parameter definitions:
!
!    a:         values to write to the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    dummy:     dummy value to be writen to row-column location
!                 with no pixel number
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    iu:        unit number to read data
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
!    rmult:     all image output to be multiplied by this amount
!    tmpval:    temporary value used to write the real value
!                 in four byte form when the input array
!                 is eight byte
! ====================================================================

      subroutine wrimgr(a,iu,rmult,nrow,ncol,ipixnum)

      implicit real*8 (a-h,o-z)
      include "SNOW.h"
      include "wgtpar.h"
      include "sun_sgi.h"
      include "help/wrimgr.h"

! ====================================================================
! Loop through the image and write each value in proper location.
! ====================================================================

      dummy = 0.0

      do 200 irow = 1,nrow

         do 100 icol = 1,ncol

! --------------------------------------------------------------------&
! If the location is within the area of interest then
! write the correct value to the image, otherwise 
! write the dummy value.
! --------------------------------------------------------------------&

            if (SUN_SGI.eq.1) then

               if (ipixnum(irow,icol).gt.0) then

                  tmpval = rmult*a(ipixnum(irow,icol))
                  write(iu,rec=((irow-1)*ncol) + icol) tmpval

               else

                  write(iu,rec=((irow-1)*ncol) + icol) dummy

               endif          

            endif

            if (SUN_SGI.eq.2) then

               if (ipixnum(irow,icol).gt.0) then

                  tmpval = rmult*a(ipixnum(irow,icol))
                  write(iu) tmpval

               else

                  write(iu) dummy

               endif

            endif

            if (SUN_SGI.eq.3) then

               if (ipixnum(irow,icol).gt.0) then

                  tmpval = rmult*a(ipixnum(irow,icol))
!	          call swap_r(tmpval,1)
                  write(iu,rec=((irow-1)*ncol) + icol) tmpval

               else
!	          call swap_r(dummy,1)	
                  write(iu,rec=((irow-1)*ncol) + icol) dummy
         
               endif

            endif

            if ( (SUN_SGI.lt.1).or.(SUN_SGI.gt.3) ) then

               write (*,*) 'wrimgr.f : Check the parameter SUN_SGI '
               write (*,*) 'in sun_sgi.h'
               stop

            endif

100      continue

200   continue

      return

      end
