! ====================================================================
!
!  Parameter definitions:
!
!    icol:      loop index for image column
!    dummy:     dummy value to be writen to row-column location
!                 with no pixel number
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
!    rmult:     all image output to be multiplied by this amount
! ====================================================================

      subroutine WRITE_BINARY(datain,rmult,nrow,ncol,ipixnum,i)

      implicit none
      real :: dummy 
      real :: rmult
      real*8 :: datain(nrow*ncol)
      real :: dataout(ncol,nrow)
      integer :: ipixnum(nrow,ncol)
      integer :: irow,icol,nrow,ncol,i,x,y

! ====================================================================
! Loop through the image and write each value in proper location.
! ====================================================================

      dummy = -999.0
       x = 1
       y = 0
       do irow = 1,nrow

         do icol = 1,ncol

! --------------------------------------------------------------------&
! If the location is within the area of interest then
! write the correct value to the image, otherwise 
! write the dummy value.
! --------------------------------------------------------------------&

               if (y .eq. nrow)then
                        y = 0
                        x = x + 1
               endif

               y = y + 1

               if (ipixnum(irow,icol).gt.0) then

                  dataout(x,y) = rmult*datain(ipixnum(irow,icol))

               else
                  
                  dataout(x,y) = -999

               endif

        enddo

      enddo

      write (1005,rec=i) dataout

      return

      end
