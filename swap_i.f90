!==================================================
!
!       subroutine swap_i
!       
!===================================================
! Byte swapping for new machine - intergers
!==================================================

       subroutine swap_i(A,N)

       implicit none
       include "help/swap_i.h"

       EQUIVALENCE (JTEMP(1),ITEMP)

       SAVE

       DO 10 I = 1,N
         ITEMP    = A(I)
         KTEMP    = JTEMP(4)
         JTEMP(4) = JTEMP(1)
         JTEMP(1) = KTEMP
         KTEMP    = JTEMP(3)
         JTEMP(3) = JTEMP(2)
         JTEMP(2) = KTEMP
         A(I)   =    ITEMP
 10    CONTINUE
       RETURN
       END
