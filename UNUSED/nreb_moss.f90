!*****************************************************************&
!*****************************************************************&
! This subroutine solves the energy balance equations for the	*&
! skin temperature of the soil, the skin temperature of the	*&
! moss layer, the mid level soil temperature and the latent,	*&
! sensible and ground heat fluxes.				*&
!								*&
! The skin temperature of the soil is solved for using a	*&
! Newton-Rapsody iteration scheme.				*&
!								*&
! Programmer : Valentijn Pauwels, July 1997.			*&
!*****************************************************************&
!*****************************************************************&

      subroutine nreb_moss(T0,T1,T2,T3,Tn0,Tn1,Tn2,tc0,tc1,tc2,tcn0,tcn1,tcn2,&
       heatcap1,heatcap2,heatcap3,heatcapn1,heatcapn2,heatcapn3,ipetopt,&
       r_LE_act,H,G,r_LE,ds,Ta,rav,rah,zmoss,z1,z2,eps,delt,toleb,maxiter,rsd,&
       rld,alb,Rn,vappres,emiss,rhoa,psychr,r_lhv,epot,itime)

!*****************************************************************&
! Inputs (stay unchanged) :					*&
!								*&
! z1		Depth of the top soil layer (m).		*&
! z2		Depth of the deep soil layer (m).		*&
! Ta		Air temperature (K).				*&
! T0		Skin temperature of the moss at end of previous	*&
!		timestep (K).					*&
! T1		Skin temperature of the soil at end of previous	*&
!		timestep (K).					*&
! T2		Mid level soil temperature of the soil at the	*&
!		end of previous timestep (K).			*&
! T3		Deep soil temperature (K).			*&
! rah		Aerodynamis resistance to heat transport (s/m).	*&
! rav		Aerodynamis resistance to mass transport (s/m).	*&
! tc0		Thermal conductivity of the moss layer of the	*&
!		previous timestep (W/mK).			*&
! tc1		Thermal conductivity of the top soil layer of	*&
!		previous timestep (W/mK).			*&
! tc2		Thermal conductivity of the lower soil layer of	*&
!		previous timestep (W/mK).			*&
! eps		0 for fully explicit scheme, 1 for fully	*&
!		implicit.					*&
! rhoa		Air density (kg/m3).				*&
! delt		Timestep (s).					*&
! itime		Number of the current timestep (-).		*&
! zmoss		Depth of the moss layer (m).			*&
! emiss		Emissivity of the moss (-).			*&
! toleb		Maximum deviation of the radiation balance	*&
!		(W/m2).						*&
! psychr	Psychometri! constant (Pa/C).			*&
! maxiter	Maximum number of iterations (-).		*&
! vappres	Air vapor pressure (Pa).			*&
! ipetopt	0 if calculations are made at potential rate,	*&
!		1 if at actual rate.				*&
! r_LE_act	Actual latent heat flux if ipetopt = 1.		*&
! heatcap1	Heat capacity of the moss layer of the previous	*&
!		timestep (J/m3K).				*&
! heatcap2	Heat capacity of the top soil layer of the	*&
!		previous timestep (J/m3K).			*&
! heatcap3	Heat capacity of the lower soil layer of the	*&
!		previous timestep (J/m3K).			*&
!								*&
! Used variables (can be changed in the code) :			*&
!								*&
! tc0		Thermal conductivity of the moss layer of the	*&
!		present timestep (W/mK).			*&
! tc1		Thermal conductivity of the top soil layer of	*&
!		present timestep (W/mK).			*&
! tc2		Thermal conductivity of the lower soil layer of	*&
!		present timestep (W/mK).			*&
! heatcapn1	Heat capacity of the moss layer of the present	*&
!		timestep (J/m3K).				*&
! heatcapn2	Heat capacity of the top soil layer of the	*&
!		present timestep (J/m3K).			*&
! heatcapn3	Heat capacity of the lower soil layer of the	*&
!		present timestep (J/m3K).			*&
!								*&
! Calculated variables :					*&
!								*&
! G		Ground heat flux (W/m2).			*&
! H		Sensible heat flux (W/m2).			*&
! ds		Energy balance error (W/m2).			*&
! Rn		Net radiation (W/m2).				*&
! Tn0		Skin temperature of the moss at end of present	*&
!		timestep (K).					*&
! Tn1		Skin temperature of the soil at end of present	*&
!		timestep (K).					*&
! Tn2		Mid level soil temperature of the soil at the	*&
!		end of present timestep (K).			*&
! r_LE		Latent heat flux (W/m2).			*&
! epot		Potential evaportranspiration (m/s).		*&
!								*&
! Intermediate variables :					*&
!								*&
! vpsat		Saturated vapor pressure (Pa).			*&
! vpdef		Vapor pressure deficit (Pa).			*&
! r_lhv		Latent heat of vaporization (J/kg).		*&
!****************************************************************&

      implicit none
      include "help/nreb_moss.h"

      parameter (SIGMA=5.6695d-8)
      parameter (row=997.d0)
      parameter (Cp=1013.)

!----------------------------------------------------------------
! Calculate the geometry parameters.				-
!----------------------------------------------------------------

      d1=zmoss
      d2=zmoss+z1
      d3=zmoss+z2

      dz0=d1
      dz1=d2-d1
      dz2=d3-d2

      dz01=0.5d0*(d2-d1)+0.5d0*d1
      dz12=0.5d0*(d3-d2)+0.5d0*(d2-d1)

!----------------------------------------------------------------
! Calculate the heat capacity at the different nodes.		-
!----------------------------------------------------------------

      c1=0.5d0*(heatcap1+heatcap2)
      c2=0.5d0*(heatcap2+heatcap3)

      cn1=0.5d0*(heatcapn1+heatcapn2)
      cn2=0.5d0*(heatcapn2+heatcapn3)

!----------------------------------------------------------------
! Start the iteration for the skin temperature.			-
!----------------------------------------------------------------

      Tn0=Ta

      iter=0

      dsold=0.d0
      dold=10.d0**(10.d0)
      Told=Tn0

200   call calc_newT(c1,c2,cn1,cn2,dz0,dz1,dz2,dz01,dz12,&
       eps,tc0,tc1,tc2,tcn0,tcn1,tcn2,delt,Tn0,T0,T1,T2,T3,Tn1,Tn2,dGdT,ddsdT,&
       heatcapn1,dsold,dsnew,Told)

!................................................................
! Update the heat storage parameters and the skin temperature	.
!................................................................

      dsold=dsnew
      dstmp=dsnew
      Told=Tn0

!................................................................
! Keep track of the number of iterations.			.
!................................................................

      iter=iter+1

!................................................................
! Calculate the saturated vapor pressure and the vapor		.
! pressure deficit.						.
!................................................................

      vpsat=611.0d0*dexp(17.27d0*(Tn0-273.15d0)/&
                         (237.3d0+Tn0-273.15d0))
      vpdef=vpsat-vappres

!................................................................
! Make the atmospheri! radiation balance to find the net	.
! radiation.							.
!................................................................

      Rntmp=rsd*(1.d0-alb)+rld-emiss*SIGMA*(Tn0**4.d0)

!................................................................
! Calculate the land components of the radiation balance.	.
!................................................................

      if (ipetopt.eq.0) then

         r_LEtmp=((rhoa*Cp)/(psychr*rav))*vpdef

      else

         r_LEtmp=r_LE_act

      endif

      Gtmp=tcn0*(Tn0-Tn1)/dz0
      Htmp=((rhoa*Cp)/rah)*(Tn0-Ta)

!................................................................
! Calculate the derivatives of the radiation balance with	.
! respect to skin temperature.					.
!................................................................

      dRndT=-4.d0*emiss*SIGMA*(Tn0**3.d0)
      dftfac=(237.d0*17.27d0)/((237.3d0+Tn0-273.15d0)**2.d0)

      if (ipetopt.eq.0) then

         dLEdT=((rhoa*Cp)/(psychr*rav))*vpsat*dftfac

      else

         dLEdT=0.d0

      endif

      dHdT=(rhoa*Cp)/rah

!................................................................
! Calculate the change in skin temperature.			.
!................................................................

      fT=Rntmp-Gtmp-r_LEtmp-Htmp-dstmp
      dfTdT=dRndT-dLEdT-dHdT-dGdT-ddsdT
      deltnr=-fT/dfTdT

      if (iter.ge.(maxiter-2)) then

         dtest=dold/deltnr

         if (dtest.le.(0.d0)) then

            dtest=0.d0-dtest

         endif

         dtest=dtest-1.d0

         if (dtest.le.(10.d0**(-2.d0))) then

!CVAL            write (*,*) 'NREB_MOSS extra solution method ',T0,Tn0
            deltnr=toleb/1.1d0

         endif

      endif

      dold=deltnr

!................................................................
! Update the skin temperature and iterate again if necessary.	.
!................................................................

      Tn0=Tn0+deltnr

      if (deltnr.le.(0.d0)) deltnr=0.d0-deltnr

      if (iter.gt.maxiter) then

         write (*,*) 'No convergence ',H,G,r_le,ds
         write (*,*) 'Temps : ',T0,Tn0,T1,Tn1,T2,Tn2
         write (*,*) 'Vars : ',heatcap1,heatcap2,heatcap3
         write (*,*) heatcapn1,heatcapn2,heatcapn3
         write (*,*) tc0,tc1,tc2,tcn0,tcn1,tcn2
!CVAL         stop
            deltnew=delt/2.d0
            if (deltnew.le.0.5) stop
            write (*,*) 'Convergence time backstepping algorithm ',deltnew
            T0new=T0
            T1new=T1
            T2new=T2
!cwc            call nreb_moss(T0new,T1new,T2new,T3,Tn0,Tn1,Tn2,&
!cwc                           tc0,tc1,tc2,tcn0,tcn1,tcn2,&
!cwc                           heatcap1,heatcap2,heatcap3,&
!cwc                           heatcapn1,heatcapn2,heatcapn3,&
!cwc                           ipetopt,r_LE_act,&
!cwc                           H,G,r_LE,ds,Ta,rav,&
!cwc                           rah,zmoss,z1,z2,eps,deltnew,toleb,maxiter,&
!cwc                           rsd,rld,alb,Rn,vappres,emiss,rhoa,&
!cwc                           psychr,r_lhv,epot,itime)
            T0new=Tn0
            T1new=Tn1
            T2new=Tn2
!cwc            call nreb_moss(T0new,T1new,T2new,T3,Tn0,Tn1,Tn2,&
!cwc                           tc0,tc1,tc2,tcn0,tcn1,tcn2,&
!cwc                           heatcap1,heatcap2,heatcap3,&
!cwc                           heatcapn1,heatcapn2,heatcapn3,&
!cwc                           ipetopt,r_LE_act,&
!cwc                           H,G,r_LE,ds,Ta,rav,&
!cwc                           rah,zmoss,z1,z2,eps,deltnew,toleb,maxiter,&
!cwc                           rsd,rld,alb,Rn,vappres,emiss,rhoa,&
!cwc                           psychr,r_lhv,epot,itime)

      endif

      if (deltnr.le.toleb) then

         goto 100

      else

         goto 200

      endif

!................................................................
! Recalculate the saturated vapor pressure and the vapor	.
! pressure deficit.						.
!................................................................

100   vpsat=611.0d0*dexp(17.27d0*(Tn0-273.15d0)/&
                         (237.3d0+Tn0-273.15d0))
      vpdef=vpsat-vappres

!................................................................
! Recalculate all components of the radiation balance.		.
!................................................................

      Rn=rsd*(1.d0-alb)+rld-emiss*SIGMA*(Tn0**4.d0)

      if (ipetopt.eq.0) then

         r_LE=((rhoa*Cp)/(psychr*rav))*vpdef

      else

         r_LE=r_LE_act

      endif

      G=tcn0*(Tn0-Tn1)/dz0
      H=((rhoa*Cp)/rah)*(Tn0-Ta)
      ds=Rn-r_LE-G-H

      if ( (G.le.(-50000.)).or.(G.ge.(50000.)).or.&
           (H.le.(-50000.)).or.(H.ge.(50000.)).or.&
           (r_LE.le.(-50000.)).or.(r_LE.ge.(50000.)).or.&
           (Tn0.le.(0.)).or.(Tn0.ge.(593.)).or.&
           (Tn1.le.(20.)).or.(Tn1.ge.(523.)).or.&
           (Tn2.le.(20.)).or.(Tn2.ge.(523.)).or.&
           (ds.le.(-50000.)).or.(ds.ge.(50000.)) ) then

            write (*,*) 'No solution ',H,G,r_le,ds
            write (*,*) 'Temps : ',T0,Tn0,T1,Tn1,T2,Tn2
            write (*,*) 'Vars : ',heatcap1,heatcap2,heatcap3
            write (*,*) heatcapn1,heatcapn2,heatcapn3
            write (*,*) tc0,tc1,tc2,tcn0,tcn1,tcn2
            write (*,*) ipetopt,r_LE_act,Ta,rav,&
                           rah,zmoss,z1,z2,eps,delt,toleb,maxiter,&
                           rsd,rld,alb,Rn,vappres,emiss,rhoa,&
                           psychr,r_lhv,epot,itime
            deltnew=delt/2.d0
            if (deltnew.le.0.5) stop
            write (*,*) 'Solution time backstepping algorithm ',deltnew
            T0new=T0
            T1new=T1
            T2new=T2
!cwc            call nreb_moss(T0new,T1new,T2new,T3,Tn0,Tn1,Tn2,&
!cwc                           tc0,tc1,tc2,tcn0,tcn1,tcn2,&
!cwc                           heatcap1,heatcap2,heatcap3,&
!cwc                           heatcapn1,heatcapn2,heatcapn3,&
!cwc                           ipetopt,r_LE_act,&
!cwc                           H,G,r_LE,ds,Ta,rav,&
!cwc                           rah,zmoss,z1,z2,eps,deltnew,toleb,maxiter,&
!cwc                           rsd,rld,alb,Rn,vappres,emiss,rhoa,&
!cwc                           psychr,r_lhv,epot,itime)
            T0new=Tn0
            T1new=Tn1
            T2new=Tn2
!cwc            call nreb_moss(T0new,T1new,T2new,T3,Tn0,Tn1,Tn2,&
!cwc                           tc0,tc1,tc2,tcn0,tcn1,tcn2,&
!cwc                           heatcap1,heatcap2,heatcap3,&
!cwc                           heatcapn1,heatcapn2,heatcapn3,&
!cwc                           ipetopt,r_LE_act,&
!cwc                           H,G,r_LE,ds,Ta,rav,&
!cwc                           rah,zmoss,z1,z2,eps,deltnew,toleb,maxiter,&
!cwc                           rsd,rld,alb,Rn,vappres,emiss,rhoa,&
!cwc                           psychr,r_lhv,epot,itime)

      endif

      epot=r_LE/(r_lhv*row)

      return

      end
