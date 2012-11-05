!   ====================================================================
!
!                     function calcrib
!
!   ====================================================================
!
!   Calculate Bulk Richardson number based on input
!   temperature, humidity, pressure, wind profile
!
!   tk2   temperature (K) at 2nd level
!   q2    specifi!   humidity (kg/kg) at 2nd level
!   p2    pressure (hPa) at 2nd level
!   u2    wind speed (m/s) at 2nd level
!   z2    distance (m) between 1st and 2nd level
!   tk1   temperature (K) at 1st level
!   q1    specifi!   humidity (kg/kg) at 1st level
!
!   ====================================================================

      function calcrib(tk2,q2,p2,u2,z2,tk1,q1)

      implicit none
      include "help/calcrib.h"
      data rdcp,GRAV/-0.286,9.81d0/

!   --------------------------------------------------------------------
!   If wind is below detectable limit, set wind speed to a small number 
!   so ribtmp doesn't divide by zero.
!   --------------------------------------------------------------------

      if (u2.lt.0.1d0) u2=0.1d0 
      
      thta1 = tk1 * (p2/1.d4)**(rdcp)
      thta2 = tk2 * (p2/1.d4)**(rdcp)
      
      thta1v = thta1
      thta2v = thta2
      
      ribtmp = GRAV * z2 *(thta2v-thta1v)/(thta2v*u2*u2) 

!   --------------------------------------------------------------------
!   Check the bounds of the Richardson number.
!   --------------------------------------------------------------------
      
      if ( (ribtmp.ge.-10000.d0).and.(ribtmp.le.10000.d0) ) then

         ribtmp=ribtmp

      else

        write (*,*) 'CALCRIB : Richardson number out of bounds ',ribtmp
        write (*,*) 'A ',tk2,q2
        write (*,*) 'B ',p2,u2,z2
        write (*,*) '!   ',tk1,q1
        write (*,*) 'D ',thta1,thta2
        write (*,*) 'E ',(thta2v-thta1v),(thta2v*u2*u2)
        write (*,*) 'F ',(thta2v-thta1v)/(thta2v*u2*u2),z2,GRAV
        stop

      endif

      calcrib = ribtmp
      
      return
      
      end
