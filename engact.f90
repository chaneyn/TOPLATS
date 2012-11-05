! ====================================================================
!
!		subroutine engact
!
! ====================================================================
!
! Subroutine to calculate the actual surface energy fluxes.
!
! ====================================================================

      subroutine engact(canclos,ievcon,xlhv,row,ivgtyp,xleactd,evtact,&
       bsdew,ioppet,iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,&
       tzsm,smold,rzsmold,tzsmold,iopthermc,thermc1,thermc2,thetar,&
       thetas,psic,bcbeta,quartz,ifcoarse,heatcap1,heatcap2,&
       heatcapold,rocpsoil,cph2o,roa,cp,roi,thermc,heatcap,rzdthetaudtemp,&
       iopgveg,iopthermc_v,tcbeta,xlai,tkact,tkactd,&
       tkmidactd,i_2l,f1,f2,f3,emiss,rescan,ravd,rahd,rnactd,&
       hactd,gactd,dshactd,tcel,vppa,psychr,zdeep,Tdeepstep,&
       rsd,r_lup,rld,toleb,maxnri,dt,i,albd,r_sdn,rnetpn,&
       gbspen,rnetd,xled,hd,gd,dshd,tkd,tkmidd,rnact,xleact,hact,&
       gact,dshact,rnetw,xlew,hw,gw,&
       dshw,tkw,tkmidw,dc,fw,tdiff,inc_frozen,ipix,initer)

      implicit none
      include "help/engact.h"

      ccc=canclos

! ====================================================================
! If evaporation is land surface controled then recalculate
! actual fluxes.
! ====================================================================

      if (ievcon.eq.1) then

! --------------------------------------------------------------------
! Calculate actual latent heat flux from dry canopy or bare soil.
! --------------------------------------------------------------------

         if (ivgtyp.eq.0) then

            xleactd = xlhv * row * (evtact-bsdew)

         else

            xleactd = xlhv * row * evtact

         endif

! --------------------------------------------------------------------
! Solve energy balance again if this option is specified.
! --------------------------------------------------------------------

         if (ioppet.eq.0) then

            call sm_cen_dif(iffroz,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm,smold,&
                            rzsmold,tzsmold)

! ====================================================================
! Calculate the soil thermal parameters.
! ====================================================================

            call soiltherm(iopthermc,thermc1,thermc2,rzsm,smtmp,&
       thetar,thetas,psic,bcbeta,tkmid,quartz,ifcoarse,&
       heatcap1,heatcap2,heatcapold,rocpsoil,row,cph2o,roa,cp,roi,&
       smold,thermc,heatcap,inc_frozen,rzdthetaudtemp)

! --------------------------------------------------------------------
! Correct the thermal parameters for soils under vegetation.
! --------------------------------------------------------------------

            if (ivgtyp.ne.0) then

               call soiladapt(iopgveg,thermc,iopthermc_v,tcbeta,&
                              xlai,thermc1,heatcap,heatcap1,zero)

            endif

! ====================================================================
! Initialize actual temperatures.
! ====================================================================

            tkactd = tkact
            tkmidactd = tkmid

! ====================================================================
! Solve the energy balance for actual temperatures and fluxes.
! ====================================================================

            tdum1=tkactd
            tdum2=tkmidactd

            if (i_2l.eq.1) then

! --------------------------------------------------------------------
! Solve assuming that the radiation budgets for over story and under
! story are related with each other through the long wave radiation
! term.
! --------------------------------------------------------------------

               call nreb(0,f1+f2-f3,2.*emiss,thermc,thermc2,heatcap,heatcap2,&
       heatcapold,rescan+ravd,rahd,tdum1,tdum2,rnactd,xleactd,evtact,hactd,&
       gactd,dshactd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,rsd,&
       r_lup+rld,toleb,maxnri,dt,i)

            else

! --------------------------------------------------------------------
! Solve  assuming that the radiation budget for over and under story
! are independent from each other, which means that the incoming long
! wave radiation for both layers are equal.
! --------------------------------------------------------------------

               call nreb(0,albd,emiss,thermc,thermc2,heatcap,heatcap2,&
       heatcapold,rescan+ravd,rahd,tdum1,tdum2,rnactd,xleactd,evtact,hactd,&
       gactd,dshactd,tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid,r_sdn,&
       rld,toleb,maxnri,dt,i)

            endif

! ====================================================================
! For penmen-montieth case just solve for sensible heat flux
! since Rnet, G and latent heat are known.
! ====================================================================

         else if(ioppet.eq.1) then

! --------------------------------------------------------------------
! Solve for sensible heat.
! --------------------------------------------------------------------

            hactd = rnetpn - xleactd - gbspen

! --------------------------------------------------------------------
! Set values for Rnet and Ground Heat Flux.
! --------------------------------------------------------------------

            rnactd = rnetpn
            gactd = gbspen

         endif

! ====================================================================
! If the evapotranspiration is atmosphere controled then set
! the fluxes equal to potential fluxes.
! ====================================================================

      else

         rnactd = rnetd
         xleactd = xled
         hactd = hd
         gactd = gd
         dshactd = dshd
         tdum1 = tkd
         tdum2 = tkmidd

      endif

! ====================================================================
! Compute average actual surface energy fluxes including any
! wet canopy fluxes.
! ====================================================================

      if (initer.eq.2) then

         tkactd=tdum1
         tkmidactd=tdum2

         rnact = rnactd*dc*(1-fw) + rnetw*(one-dc*(1-fw))
         xleact = xleactd*dc*(1-fw) + xlew*(one-dc*(1-fw))
         hact = hactd*dc*(1-fw) + hw*(one-dc*(1-fw))
         gact = gactd*dc*(1-fw) + gw*(one-dc*(1-fw))
         dshact = dshactd*dc*(1-fw) + dshw*(one-dc*(1-fw))
         tkact = tkactd*dc*(1-fw) + tkw*(one-dc*(1-fw))
         tkmid = tkmidactd*dc*(1-fw) + tkmidw*(one-dc*(1-fw))

      endif

      tdiff=tdum1*dc*(1-fw) + tkw*(one-dc*(1-fw))

! ====================================================================
! In case of unreasonable soil temperatures, which will lead to
! unreasonable ground heat flux values, abort the program and
! print the most important variables out.
! ====================================================================

      if(gactd.gt.15000) then

         write (*,*) 'ENGACT error '
         write (*,*) 'tkactd =',tkactd
         write (*,*) 'tkw =',tkw
         write (*,*) 'tkmidactd =',tkmidactd
         write (*,*) 'tkmidw =',tkmidw
         write (*,*) 'gactd =',gactd
         write (*,*) 'gw =',gw
         write (*,*) 'gact    =',gact
         write (*,*) 'dc    =',dc
         write (*,*) 'fw    =',fw 
         write (*,*) 'gactd    =',gactd
         write (*,*) 'gw    =',gw
         write (*,*) 'dc    =',dc
         write (*,*) 'fw    =',fw
         write (*,*) 'gact = ',gact
         write (*,*) 'therm! =',thermc
         write (*,*) 'thermc1 =',thermc1
         write (*,*) 'thermc2 =',thermc2
         write (*,*) 'heatcap1 =',heatcap1
         write (*,*) 'heatcap2 =',heatcap2
         write (*,*) 'heatcapold =',heatcapold
         write (*,*) 'smtmp = ',smtmp
         write (*,*) 'smold = ',smold
         write (*,*) 'ipix = ',ipix
         write (*,*) 'rzsm = ',rzsm
         write (*,*) 'tzsm = ',tzsm
         write (*,*) 'thetar = ',thetar
         write (*,*) 'thetas = ',thetas
         write (*,*) 'quartz = ',quartz
         write (*,*) 'iffroz = ',iffroz
         write (*,*) 'ifcoarse = ',ifcoarse
         write (*,*) 'rocpsoil = ',rocpsoil
         write (*,*) 'row = ',row
         write (*,*) 'cph2o = ',cph2o
         write (*,*) 'roa = ',roa
         write (*,*) 'cp = ',cp
         write (*,*) 'zmid = ',zmid
         write (*,*) 'zrzmax = ',zrzmax
         write (*,*) 'rzsmold = ',rzsmold
         write (*,*) 'tzsmold = ',tzsmold

         stop

      endif

      return

      end
