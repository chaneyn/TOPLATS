      integer SUN_SGI

      parameter (SUN_SGI=3)

! ====================================================================
! SUN_SGI	1 if you are running the program on a SUN.
! 		2 if you are running the program on an SGI.
! 		3 if you are running the program on an linux machine.
!
!  Note: Only use "3" if you want bit-order of i/o binary images to be
!  correct for sun or sgi (solaris or irix) machines.  If you want
!  bit order to follow DEC convention use "1".  
! 
!
! ====================================================================
