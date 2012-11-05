! ====================================================================
!
!			subroutine wrimgi
!
! ====================================================================
!
! Subroutine to write an image of 4 byte integers from an array
!   of these values indexed by the soils-topographi! index
!   pixel numbers.
!
! ====================================================================
!
!  Parameter definitions:
!
!    ia:        values to write to the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    idummy:    dummy value to be writen to row-column location
!                 with no pixel number
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    itmpvl:    temporary value used to hold the output value
!    iu:        unit number to read data
!    mult:      all image output to be multiplied by this amount
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
! ====================================================================

      subroutine wrimgi(ia,iu,mult,nrow,ncol,ipixnum)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "sun_sgi.h"
      include "help/wrimgi.h"

! ====================================================================
! Loop through the image and write each value in proper location.
! ====================================================================

      do 200 irow = 1,nrow

         do 100 icol = 1,ncol

! --------------------------------------------------------------------&
! If the location is within the area of interest then
! write the correct value to the image, otherwise 
! write the dummy value.
! --------------------------------------------------------------------&

            if (SUN_SGI.eq.1) then

               if (ipixnum(irow,icol).gt.0) then

                  itmpvl = mult*ia(ipixnum(irow,icol))
                  write(iu,rec=((irow-1)*ncol) + icol) itmpvl

               else

                  write(iu,rec=((irow-1)*ncol) + icol) idummy

               endif          

            endif

            if (SUN_SGI.eq.2) then

               if (ipixnum(irow,icol).gt.0) then

                  itmpvl = mult*ia(ipixnum(irow,icol))
                  write(iu) itmpvl

               else

                  write(iu) idummy

               endif

            endif

            if (SUN_SGI.eq.3) then

               if (ipixnum(irow,icol).gt.0) then

                  itmpvl = mult*ia(ipixnum(irow,icol))
!	          call swap_i(itmpvl,1)
                  write(iu,rec=((irow-1)*ncol) + icol) itmpvl

               else
!	          call swap_i(idummy,1)	
                  write(iu,rec=((irow-1)*ncol) + icol) idummy
         
               endif

            endif

            if ( (SUN_SGI.lt.1).or.(SUN_SGI.gt.3) ) then

               write (*,*) 'wrimgi.f : Check the parameter SUN_SGI '
               write (*,*) 'in sun_sgi.h'
               stop

            endif

100      continue

200   continue

      return

      end
