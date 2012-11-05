! ====================================================================
!
!			subroutine imgopn
!
! ====================================================================
!
! Subroutine to set the name of and open a file time varying
!   image input or output.
!
! ====================================================================

      subroutine imgopn(prefix,i,iu,img_opt,iyear,iday,ihour)

      implicit none
      include "help/imgopn.h"

! ====================================================================
! Initialize the file name to null and the file name suffix to 
! the time step, or if requested to the year_day_hour string.
! ====================================================================

      fname=''
      write (suffix,1001) iyear,iday,ihour

! ====================================================================
! Find the length of the prefix.
! ====================================================================

      lgth = index(prefix,' ')

! ====================================================================
! Build the file name one character at a time until a blank space
! in the prefix is encountered then add a period and the suffix.
! ====================================================================

      do 100 k=1,lgth-1

        fname(k:k)=prefix(k:k)

100   continue

      fname(lgth:lgth)='_'
      fname(lgth+1:lgth+13)=suffix

! ====================================================================
! Open the file for binary I/O.
! ====================================================================

         open(unit=iu,file=fname,form='unformatted',&
              access='direct',recl=4)

! ====================================================================
! Format statment.
! ====================================================================

1000  format (I5.5)
1001  format (I2.2,'_',I3.3,'_',I2.2,'.bin')

      return

      end
