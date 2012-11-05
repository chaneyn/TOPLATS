! ====================================================================
!
!			subroutine rdatb
!
! ====================================================================
!
! Subroutine to read in the topographi! index image and set up
! at translation between the pixel number and the image row.and.&
! column.
!
! ====================================================================

      subroutine rdatb(atb,nrow,ncol,ipixnum,ixpix,iypix,npix)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "sun_sgi.h"
      include "help/rdatb.h"

! ====================================================================
! Initialize pixel number counter.
! ====================================================================

      ip = 1

! ====================================================================
! Loop through image row by row.
! ====================================================================

      do 500 irow=1,nrow

         do 400 icol=1,ncol

! --------------------------------------------------------------------
! Read image value at irow,icol.
! --------------------------------------------------------------------


               read(9,rec=((irow-1)*ncol) + icol) tempatb

! --------------------------------------------------------------------
! Check for atb value at current location.  If there is
! a value then record the value and location.  If not
! then mark location as outside the area of interest.
! --------------------------------------------------------------------

            if (tempatb.gt.0.0) then

               atb(ip) = tempatb
               ipixnum(irow,icol) = ip
               ixpix(ip) = icol
               iypix(ip) = irow
               ip = ip + 1

            else

               ipixnum(irow,icol) = 0

            endif

400      continue

500   continue

! ====================================================================
! Set the total number of pixels.
! ====================================================================

      npix = ip - 1

      if (npix.gt.MAX_PIX) then

         write (*,*) 'rdatb : npix greater then MAX_PIX ',npix,MAX_PIX
         stop

      endif

      return

      end
