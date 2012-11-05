! ====================================================================
!
!               subroutine inidum
!
! ====================================================================
!
! Initialize the snow model variables that will be used in the
! snow subroutine but that not will be altered in the main program.
!
! ====================================================================

      subroutine inidum(dum1,PackWater,dum2,SurfWater,dum3,Swq,dum4,&
       VaporMassFlux,dum5,TPack,dum6,TSurf,dum7,r_MeltEnergy,dum8,Outflow,&
       dum9,xleact_snow,dum10,hact_snow,dum11,rn_snow,dum12,dens)

      implicit none
      include "help/inidum.h"

      dum1=PackWater
      dum2=SurfWater
      dum3=Swq
      dum4=VaporMassFlux
      dum5=TPack
      dum6=TSurf
      dum7=r_MeltEnergy
      dum8=Outflow
      dum9=xleact_snow
      dum10=hact_snow
      dum11=rn_snow
      dum12=dens

      return

      end
