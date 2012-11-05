! ====================================================================
!
!            subroutine zero_snowvar
!
! ====================================================================
!
! If there is no snow all the snow pack variablesa are set to zero.
!
! ====================================================================

      subroutine zero_snowvar(PackWater,SurfWater,Swq,VaporMassFlux,&
       TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow,dens)

      implicit none
      include "help/zero_snowvar.h"

      PackWater=0.d0
      SurfWater=0.d0
      Swq=0.d0
      VaporMassFlux=0.d0
      TPack=0.d0
      TSurf=0.d0
      r_MeltEnergy=0.d0
      Outflow=0.d0
      xleact_snow=0.d0
      hact_snow=0.d0
      dens=0.d0

      return

      end
