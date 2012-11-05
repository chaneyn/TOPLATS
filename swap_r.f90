!==================================================
!
!       subroutine swap_r
!       
!===================================================
! Byte swapping for new machine - reals
!==================================================

       subroutine swap_r(A,N)

       implicit none
       include "help/swap_r.h"

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
         A(I)     = ITEMP
 10    CONTINUE
       RETURN
       END
