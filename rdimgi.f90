! ====================================================================
!
!			subroutine rdimgi
!
! ====================================================================
!
! Subroutine to read in an image of integers and return an array
!   of these values indexed by the soils-topographi! index
!   pixel numbers.
!
! ====================================================================
!
!
!  Parameter definitions:
!
!    ia:        values read from the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    itmpval:   temporary 4 byte integer value read from the image
!    iu:        unit number to read data
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
! ====================================================================

      subroutine rdimgi(ia,iu,nrow,ncol,ipixnum)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "sun_sgi.h"
      include "help/rdimgi.h"

! ====================================================================
! Loop through the image and read each value.
! ====================================================================

      do irow = 1,nrow

         do icol = 1,ncol
                
             if (iu .eq. 12 .or. iu .eq. 11)then
                itmpval = (irow-1)*ncol + icol
             else
                read(iu,rec=((irow-1)*ncol) + icol) itmpval
             endif
! --------------------------------------------------------------------
! If the location is within the area of interest then
! assign to array 'ia', otherwise read next value.
! --------------------------------------------------------------------

            if (ipixnum(irow,icol).gt.0) then

               ia(ipixnum(irow,icol)) = itmpval

            endif

          enddo
        
      enddo

      return

      end
