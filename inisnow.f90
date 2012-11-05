! ====================================================================
!
!                     subroutine inisnow
!
! ====================================================================
!
! Initialize the snow pack variables to zero.
!
! ====================================================================

      subroutine inisnow(pw,sw,swq,vm,me,o,np,nl,na)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/inisnow.h"

      integer np,nl,na,t,u,v

      do t=1,na

         do u=1,nl

            do v=1,np

               pw(t,u,v)=0.d0
               sw(t,u,v)=0.d0
               swq(t,u,v)=0.d0
               vm(t,u,v)=0.d0
               me(t,u,v)=0.d0
               o(t,u,v)=0.d0

            enddo

         enddo

      enddo

      return

      end
