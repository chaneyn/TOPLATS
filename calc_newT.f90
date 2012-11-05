!****************************************************************
!****************************************************************
! Subroutine that solves the energy balance 2x2 matrix system.	*
! Is only called from nreb_moss().				*
!****************************************************************
!****************************************************************

      subroutine calc_newT(c1,c2,cn1,cn2,dz0,dz1,dz2,dz01,dz12,&
       eps,tc0,tc1,tc2,tcn0,tcn1,tcn2,delt,Tn0,T0,T1,T2,T3,Tn1,Tn2,dGdT,ddsdT,&
       heatcapn1,dsold,dsnew,Told)

      implicit none
      include "help/calc_newT.h"

!----------------------------------------------------------------
! Set up the different matrix coefficients for the matrix of	-
! the new timestep.						-
!----------------------------------------------------------------

      fn00=-eps*tcn0/(dz0*dz01)
      f00=-(1.d0-eps)*tc0/(dz0*dz01)
      fn01=cn1/delt+eps*tcn0/(dz0*dz01)+eps*tcn1/(dz1*dz01)
      f01=-c1/delt+(1.d0-eps)*tc0/(dz0*dz01)+(1.d0-eps)*tc1/(dz1*dz01)
      fn02=-eps*tcn1/(dz1*dz01)
      f02=-(1.d0-eps)*tc1/(dz1*dz01)

!----------------------------------------------------------------
! Set up the different matrix coefficients for the matrix of	-
! the old timestep.						-
!----------------------------------------------------------------

      fn10=-eps*tcn1/(dz1*dz12)
      f10=-(1.d0-eps)*tc1/(dz1*dz12)
      fn11=cn2/delt+eps*tcn1/(dz1*dz12)+eps*tcn2/(dz2*dz12)
      f11=-c2/delt+(1.d0-eps)*tc1/(dz1*dz12)+(1.d0-eps)*tc2/(dz2*dz12)
      fn12=-eps*tcn2/(dz2*dz12)
      f12=-(1.d0-eps)*tc2/(dz2*dz12)

!----------------------------------------------------------------
! Calculate the right hand side elements of the system.		-
!----------------------------------------------------------------

      rhs1=-fn00*Tn0-f00*T0
      rhs2=-fn12*T3-f12*T3

      rhs1=rhs1-f01*T1-f02*T2
      rhs2=rhs2-f10*T1-f11*T2

!----------------------------------------------------------------
! Solve for the soil temperatures.				-
!----------------------------------------------------------------

      Tn2=rhs1-fn01*rhs2/fn10
      Tn2=Tn2/(fn02-fn11*fn01/fn10)
      Tn1=(rhs1-fn02*Tn2)/fn01

!----------------------------------------------------------------
! Calculate the derivative of ground heat flux with respect to	-
! skin temperature.						-
!----------------------------------------------------------------

      Gold=tc0*(T0-T1)/dz0
      Gnew=tcn0*(Tn0-Tn1)/dz0

      if (Tn0.eq.T0) then

         dGdt=0.d0

      else

         dGdT=(Gnew-Gold)/(Tn0-T0)

      endif

!----------------------------------------------------------------
! Calculate the derivative of the heat storage with respect to	-
! skin temperature.						-
!----------------------------------------------------------------

      dsnew=heatcapn1*(Tn1-T1)*dz0/delt

      if (Tn0.eq.Told) then

         ddsdT=0.d0

      else

         ddsdT=(dsnew-dsold)/(Tn0-Told)

      endif

      return

      end
