! ====================================================================
!
!                         function clcddif
!
! ====================================================================
!
! Calculate the derivative of the diffusive flux out of the surface zone
! and transmission zone with respect to soil moisture.
!
! ====================================================================

      function clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,&
         ikopt,xksrz,xkstz,ff,&
         zrz,ztz,bcbeta,thetas,thetar,psic,tzsm)
        
      implicit none
      include "help/clcddif.h"

! ====================================================================
! Calculate the derivatives of the diffusive flux with respect to
! soil moisture
! ====================================================================
    
      if(zrz.eq.0.and.ztz.eq.0)then
       ddifrzdth1=0.0
       ddifrzdth2=0.0
      else     
        ddifrzdth1 = (.75*bcbeta*(2+bcbeta)*psic*(.5*rzsm - .5*tzsm)&
               *((.75*rzsm + .25*tzsm - thetar)/(-thetar + thetas))**(1+bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas)**2) 

        ddifrzdth1 = ddifrzdth1 + (.5*bcbeta*psic*((.75*rzsm - thetar + &
                .25*tzsm)/(- thetar + thetas))**(2+ bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas))

        ddifrzdth2 = (.25*bcbeta*(2+bcbeta)*psic*(.5*rzsm - .5*tzsm)&
               *((.75*rzsm + .25*tzsm - thetar)/(-thetar + thetas))**(1+bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas)**2)

        ddifrzdth2 = ddifrzdth2 - (.5*bcbeta*psic*&
               ((.75*rzsm - thetar + .25*tzsm)/(-thetar + thetas))**(2+ bcbeta))&
               /((.5/xkstz + .5/xksrz)*(.5*zrz + .5*ztz)*(-thetar + thetas))

      endif

      if(ztz.eq.0.)then
        ddiftzdth2 = 0.0
      else
        ddiftzdth2 = (.75*bcbeta*(2+bcbeta)*xkstz*psic*(.5*tzsm - .5*thetas)&
               *((.75*tzsm + .25*thetas - thetar)/(-thetar + thetas))**(1+bcbeta))&
               /(ztz*(-thetar + thetas)**2)  
      
        ddiftzdth2 = ddiftzdth2 + (.5*bcbeta*xkstz*psic*&
               ((.75*tzsm - thetar + .25*thetas)/(-thetar + thetas))**(2 + bcbeta))&
               /(ztz*(-thetar + thetas))
      endif

      return

      end
