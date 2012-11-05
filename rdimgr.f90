! ====================================================================
!
!			subroutine rdimgr
!
! ====================================================================
!
! Subroutine to read in an image of 4 byte reals and return an array
!   of these values indexed by the soils-topographic index
!   pixel numbers.
!
! ====================================================================
!
!
!  Parameter definitions:
!
!    a:         values read from the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    iu:        unit number to read data
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
!    tmpval:    temporary 4 byte real value read from the image
! ====================================================================

      subroutine rdimgr(a,iu,nrow,ncol,ipixnum)

      implicit real*8 (a-h,o-z)
      include "SNOW.h"
      include "wgtpar.h"
      include "help/rdimgr.h"

! ====================================================================
! Loop through the image and read each value.
! ====================================================================
      do 200 irow = 1,nrow
         do 100 icol = 1,ncol
               read(iu,rec=((irow-1)*ncol) + icol) tmpval

! --------------------------------------------------------------------
! If the location is within the area of interest then
! assign to array 'a', otherwise read next value.
! --------------------------------------------------------------------

            if (ipixnum(irow,icol).gt.0) then

               a(ipixnum(irow,icol)) = tmpval
            
            endif

100      continue

200   continue

      return
      end
