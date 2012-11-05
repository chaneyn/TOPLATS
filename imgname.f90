! ====================================================================
!
!			subroutine imgopn
!
! ====================================================================
!
! Subroutine to set the name of a file time varying
!   image input or output.
!
! ====================================================================

      subroutine imgname(prefix,iyear,iday,ihour,fname)

      implicit none
      character*100 prefix
      character*105 fname
      character*13 suffix
      integer i,iu,lgth,k,iyear,iday,ihour

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
! Format statment.
! ====================================================================

1000  format (I5.5)
1001  format (I2.2,'_',I3.3,'_',I2.2,'.bin')

      return

      end subroutine
