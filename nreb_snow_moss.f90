! ====================================================================
! Calculate the soil temperatures given that a snow pack is on top
! of the soil with a moss layer.
! ====================================================================

      subroutine nreb_snow_moss(tc0,tc1,tc2,tcn0,tcn1,tcn2,heatcap1,heatcap2,&
       heatcap3,heatcapn1,heatcapn2,heatcapn3,T0,T1,T2,T3,Tn0,Tn1,Tn2,Ta,&
       zmoss,z1,z2,eps,delt,gnew)

      implicit none
      include "help/nreb_snow_moss.h"

! --------------------------------------------------------------------
! Calculate the distance between the nodes and the distance between
! the centers of the layers.
! --------------------------------------------------------------------

      d1=zmoss
      d2=zmoss+z1
      d3=zmoss+z2

      dz0=d1
      dz1=d2-d1
      dz2=d3-d2

      dz01=0.5d0*(d2-d1)+0.5d0*d1
      dz12=0.5d0*(d3-d2)+0.5d0*(d2-d1)

      c1=0.5d0*(heatcap1+heatcap2)
      c2=0.5d0*(heatcap2+heatcap3)

      cn1=0.5d0*(heatcapn1+heatcapn2)
      cn2=0.5d0*(heatcapn2+heatcapn3)

! --------------------------------------------------------------------
! Initialize the ground heat storage terms.
! --------------------------------------------------------------------

      dsold=0.d0
      dsnew=0.d0

! --------------------------------------------------------------------
! The moss skin temperature is assumed to be equal to the temperature
! at the bottom of the snow pack.
! --------------------------------------------------------------------

      Tn0=Ta

      dsold=0.d0
      Told=Tn0

! --------------------------------------------------------------------
! Calculate the skin temperature of the soil given the skin temperature
! of the moss and the deep soil temperature.
! --------------------------------------------------------------------

      call calc_newT(c1,c2,cn1,cn2,dz0,dz1,dz2,dz01,dz12,&
       eps,tc0,tc1,tc2,tcn0,tcn1,tcn2,delt,Tn0,T0,T1,T2,T3,Tn1,Tn2,dGdT,ddsdT,&
       heatcapn1,dsold,dsnew,Told)

      Gnew=tcn0*(Tn0-Tn1)/dz0

      if ( (Tn0.le.(-50.)).or.(Tn0.ge.(543.)).or.&
           (Tn1.le.(10.)).or.(Tn1.ge.(503.)).or.&
           (Tn2.le.(10.)).or.(Tn2.ge.(503.)).or.&
           (dsnew.le.(-10000.)).or.(dsnew.ge.(10000.)) ) then

        write (*,*) 'Error in the energy balance under snow/moss '
        write (*,*) tc0,tc1,tc2,tcn0,tcn1,tcn2,&
                                heatcap1,heatcap2,heatcap3,&
                                heatcapn1,heatcapn2,heatcapn3,&
                                T0,T1,T2,T3,Tn0,Tn1,Tn2,Ta,&
                                zmoss,z1,z2,eps,delt,gnew
        stop

      endif

      return

      end
