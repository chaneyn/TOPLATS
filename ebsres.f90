! ====================================================================
!
!			subroutine ebsres
!
! ====================================================================
!
! Subroutine ebares calculates the actual rate of evaporation 
! from bare soils using a soil resistance parameterization.
!
! ====================================================================

      subroutine ebsres(inc_frozen,irestype,rsoil,rzsm,srespar1,tkact,&
       srespar2,rzsm_u,srespar3,ravd,iffroz,thetas,tkmid,&
       zmid,zrzmax,smtmp,tzsm,smold,rzsmold,tzsmold,&
       iopthermc,thermc1,thermc2,thetar,heatcapold,psic,bcbeta,&
       quartz,heatcap1,ifcoarse,heatcap2,rocpsoil,row,cph2o,roa,cp,&
       roi,thermc,heatcap,rzdthetaudtemp,dshact,albd,emiss,rahd,ebscap,&
       tcel,vppa,psychr,xlhv,zdeep,Tdeepstep,rsd,rld,toleb,maxnri,dt,i,tkel,&
       zww,za,uzw,zpd,z0m,press,rib,rnetpn,gbspen,epetd,evtact,ievcon,&
       bsdew,z0h,ioppet)

      implicit none
      include "help/ebsres.h"

      data TOLSTAB/0.1/
      data MAXITER/10/

! ====================================================================
! Calculate the bare soil resistance to evaporation.
! ====================================================================

      if (inc_frozen.eq.0) then

! --------------------------------------------------------------------
! No distinction between ice and unfrozen water is made.
! --------------------------------------------------------------------

         call calcrsoil(irestype,rsoil,srespar1,&
                        srespar2,srespar3,thetas,rzsm,tkact)

      else

! --------------------------------------------------------------------
! Treat the ice as mineral soil in the calculation of soil
! resistance.
! --------------------------------------------------------------------

         call calcrsoil(irestype,rsoil,srespar1,&
                        srespar2,srespar3,thetas,rzsm_u,tkact)

      endif

! --------------------------------------------------------------------
! Total resistance = soil resistance plus aerodynami! resistance.
! --------------------------------------------------------------------

      raveff = rsoil + ravd

! ====================================================================
! Initialize soil moisture for the calculation of the thermodynamic
! parameters, as a centered difference.
! ====================================================================

      call sm_cen_dif(iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm,smold,&
                      rzsmold,tzsmold)

! ====================================================================
! Calculate the soil thermal parameters.
! ====================================================================

      call soiltherm(iopthermc,thermc1,thermc2,rzsm,smtmp,&
       thetar,thetas,psic,bcbeta,tkmid,quartz,ifcoarse,&
       heatcap1,heatcap2,heatcapold,rocpsoil,row,cph2o,roa,cp,roi,&
       smold,thermc,heatcap,inc_frozen,rzdthetaudtemp)

! ====================================================================
! Initialize temperatures for the solution of the energy balance
! equations.
! ====================================================================

      dumtk = tkact
      dumtkmid = tkmid
      dumds = dshact

! ====================================================================
! If energy balance method is used then call subroutine to 
! do newton-raphson iterations to find the energy balance.
! ====================================================================

      if (ioppet.eq.0) then

         ttemp = dumtk
         tmidtemp = dumtkmid

! --------------------------------------------------------------------
! Solve the energy balance at potential evaporation rate.
! --------------------------------------------------------------------

         call nreb(1,albd,emiss,thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
       raveff,rahd,dumtk,dumtkmid,dumrn,dumxle,ebscap,dumh,dumg,dumds,tcel,&
       vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,rld,toleb,maxnri,dt,i)

! --------------------------------------------------------------------
! Check for large change in skin temperature which affects
! magnitude of stability correction.
! --------------------------------------------------------------------

         iter = 0

200      if ((abs(ttemp-dumtk)).gt.TOLSTAB.and.iter.lt.MAXITER) then

            iter = iter + 1

! --------------------------------------------------------------------
! As long as there is no convergence in skin temperature of the
! bare soil, keep on iterating for the skin temperature.
! --------------------------------------------------------------------

            ttemp = dumtk
            tacttemp = tkact
            tkact = dumtk

! --------------------------------------------------------------------
! Recalculate the stability correction.
! --------------------------------------------------------------------

            call stabcor(zww,za,uzw,zpd,z0m,tkel,press,tkact,vppa,rib)
            tkact = tacttemp
            dumtkmid = tmidtemp

! --------------------------------------------------------------------
! Recalculate the aerodynami! resistances.
! --------------------------------------------------------------------

            rahd = calcra(uzw,zww,za,zpd,z0m,z0h,rib)

            ravd = rahd
            raveff = rsoil + ravd

! --------------------------------------------------------------------
! Resolve the energy balance at potential evaporation rate.
! --------------------------------------------------------------------

            call nreb(1,albd,emiss,thermc1,thermc2,heatcap1,heatcap2,&
       heatcapold,raveff,rahd,dumtk,dumtkmid,dumrn,dumxle,ebscap,dumh,dumg,&
       dumds,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,rld,toleb,&
       maxnri,dt,i)

! --------------------------------------------------------------------
! Check again whether convergence has been achieved.
! --------------------------------------------------------------------

            goto 200

         endif

! ====================================================================
! If unreasonable skin temperatures print the variables out and
! abort the program.
! ====================================================================

         if (dumtk.lt.100.d0.or.dumtk.gt.400.d0) then

            print*,"ebsres:ERROR:skin temp outside reas. bounds Tskin="&
                ,dumtk
            print*,'thermc1 =             ',thermc1
            print*,'thermc2 =             ',thermc2
            print*,'heatcap1 =             ',heatcap1
            print*,'heatcap2 =             ',heatcap2
            print*,'dt =                  ',dt
            print*,'tmidk =             ',dumtkmid
            print*,'rav =             ',raveff
            print*,'roa =             ',roa
            stop

         endif

      else

! ====================================================================
! If penman-montieth is used then just recalculate evaporation
! using rsoil like canopy resistance.
! ====================================================================

         vpsat = 611.d0*dexp((17.27d0*tcel)/(237.3d0+tcel))
         vpdef = vpsat-vppa
         dvpsdt = 4098.d0*vpsat/((237.3d0+tcel)**two)
         pstar = psychr*(one+rsoil/ravd)
         ebsmf = ((dvpsdt*(rnetpn-gbspen))+&
                  ((cp*roa*vpdef)/ravd)) /&
                 ((dvpsdt+pstar)*xlhv)
         ebscap = ebsmf/row

      endif

! ====================================================================
! Make sure evaporative capacity is not negative.
! ====================================================================

      if (ebscap.lt.zero) ebscap=zero

! ====================================================================
! Set actual evaporation to minimum of potential evaportation from
! the bare soil and the soil controlled evaporation rate.
! Also add the land area to the soil controlled or atmospheric
! controlled percentage.
! ====================================================================

      if (ebscap.le.epetd) then

         evtact = ebscap
         ievcon = 1

      else

         evtact = epetd
         ievcon = 2

      endif

! ====================================================================
! If actual evaporation is negative then condensation occuring.
! ====================================================================

      if (evtact.lt.zero) then

         bsdew = -evtact
         evtact = zero

      else

         bsdew = zero

      endif

      return

      end
