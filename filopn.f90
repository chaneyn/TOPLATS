! ====================================================================
!
!			subroutine filopn
!
! ====================================================================
!
! Subroutine to open input/output files and set up variable to 
! control the printing of output files.
!
! ====================================================================

      subroutine filopn(iprn,ioutst,ioutsp,iouten,nseries,icurser,fnimg)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "sun_sgi.h"
      include "help/filopn.h"

! ====================================================================
! Initialize variable which controls output file printing.
! ====================================================================
     

      do 40 jj = 1,MAX_FIL

        iprn(jj) = 0

40    continue 
     
! ====================================================================
! Open the file which lists the input/output file names.  This
! file should be given as the second word in the command line.
! ====================================================================

       call getarg(1,topfil)
     
       open(unit=4,file=topfil)

! ====================================================================
! Read the filenames from file containing i/o filenames.
! ====================================================================

50    read(4,1000) iu,fname
    
      if (iu.gt.MAX_FIL) then

         write (*,*) 'filopn.f : File unit number is greater than '
         write (*,*) 'MAX_FIL ',iu,MAX_FIL
         stop

      endif

! ====================================================================
! If unit number is between 35 and 55, then read the 
! image output start time and spacing and the image
! name prefix for each series of image output.
! ====================================================================

      if ( ((iu.ge.35).and.(iu.le.55)).or.&
           ((iu.ge.151).and.(iu.le.199)) ) then

         iseries = 0
60       iseries = iseries + 1 

         if (iseries.gt.MAX_SER) then

            write (*,*) 'filopn.f : Number of input/output series'
            write (*,*) 'is greater MAX_SER ',iseries,MAX_SER
            stop

         endif
 
! --------------------------------------------------------------------
! As long as there are parts of the timeseries for which this variable
! is to be printed out keep on reading these parts in.
! --------------------------------------------------------------------

         read(4,*) ioutst(iu,iseries),ioutsp(iu,iseries),&
                   iouten(iu,iseries),icont

         if (icont.eq.1) goto 60

         nseries(iu) = iseries
         icurser(iu) = 1
         fnimg(iu) = fname

! ====================================================================
! If unit number is between 20 and 34, then set the
! input image prefix.
! ====================================================================

      else if ((iu.ge.20).and.(iu.le.34)) then

         fnimg(iu) = fname      

! ====================================================================
! If the unit number is between 7 and 19 then open file for
! binary image.
! ====================================================================

      else if ((iu.ge.7).and.(iu.le.19)) then

         open(unit=iu,file=fname,form='unformatted',&
                 access='direct',recl=4)

! ====================================================================
! If the unit number is 55 or greater then
! open that file for ASCII sequential input/output.
! ====================================================================

      else if ((iu.gt.55)) then

         open(unit=iu,file=fname)

! ====================================================================
! If the unit number is negative then done reading and return
! to calling routine.
! ====================================================================

      else

         return

      endif

! ====================================================================
! Set 'iprn' so that the output file will print.
! ====================================================================

      iprn(iu) = 1

! ====================================================================
! Read the next filename.
! ====================================================================
   
      goto 50

! ====================================================================
! Format statement
! ====================================================================

1000  format(i5,a100)        

      end
