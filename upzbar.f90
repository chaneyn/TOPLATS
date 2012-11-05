! ====================================================================
!
!			subroutine upzbar
!
! ====================================================================
!
!  Updates the average water table depth.
!
! ====================================================================

      subroutine upzbar(i,ic,iopbf,q0,ff,zbar,dtil,basink,dd,xlength,&
       gwtsum,capsum,area,r_lakearea,dt,etwtsum,rzpsum,tzpsum,psicav,ivgtyp,&
       ilandc,npix,icatch,zw,psic,isoil,zrzmax,&
       tzsm1,thetas,rzsm1,zbar1,qbreg,zbar1rg,iprn,pixsiz)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/upzbar.h"

! ====================================================================
! Chose option for calculating baseflow.
! ====================================================================

      if (iopbf.eq.0) then

! --------------------------------------------------------------------&
! Calculate baseflow according to Sivapalan et al.(1987).
! --------------------------------------------------------------------&

         qb = q0 * dexp(-ff*zbar)

      else
! --------------------------------------------------------------------&
! Calculate baseflow according to Troch et al.(1992).
! --------------------------------------------------------------------&

         hbar = dtil - zbar
         qb = 5.772*basink*hbar*hbar*dd*xlength

      endif

! ====================================================================
! Determine net recharge to water table and available soil 
! water storage.
! ====================================================================

      zbrflx = (gwtsum - capsum - etwtsum - (qb/ (area-r_lakearea) )) * dt
      zbrpor = (rzpsum+tzpsum) * (zbar-psicav)

      if (zbrflx.gt.zbrpor) then

! ====================================================================
! If net recharge exceeds storage assign overflow to streamflow.
! ====================================================================

         qzbar = (zbrflx - zbrpor)/dt
         zbrflx = zbrpor
         qb = qb + qzbar*(area-r_lakearea)

! --------------------------------------------------------------------&
! If water table is falling but storage is full zbar1 will
! blow up. use rzsm1 and tzsm1 (vs rzsm and tzsm) to
! calculate storage so that it is > zero.
! --------------------------------------------------------------------&

      else if ((zbrflx.le.zero).and.(zbrpor.le.zero)) then

! --------------------------------------------------------------------&
! Recalculate rzpsum and tzpsum with new soil moistures.
! --------------------------------------------------------------------&

         rzpsum = zero
         tzpsum = zero

         do 100 mm=1,npix

            if (ivgtyp(ilandc(mm)).ge.0) then

               if ((icatch(mm)).eq.ic) then

                  if ((zw(mm)-psic(isoil(mm))).gt.zrzmax) then

                     tzpsum = tzpsum+(thetas(isoil(mm))-tzsm1(mm))

                  else if ((zw(mm)-psic(isoil(mm)).gt.zero)) then

                     rzpsum = rzpsum+(thetas(isoil(mm))-rzsm1(mm))

                  endif

               endif

            endif

100      continue

      endif

! ====================================================================
! Update zbar by taking the total flux and dividing by
! the average porosity just above the water table.
! ====================================================================

! --------------------------------------------------------------------&
! If the available porosity is nonzero divide the flux by its value.
! --------------------------------------------------------------------&

      if ( (rzpsum+tzpsum).gt.(0.001d0)) then

         zbar1 = zbar - zbrflx/(rzpsum+tzpsum)

      endif

      if ( (rzpsum+tzpsum).le.(0.001d0)) then

         zbar1=zbar

      endif

! ====================================================================
! Find change in water table.
! ====================================================================

      dzbar = zbar1-zbar

! ====================================================================
! Find new average water table depth and baseflow for entire region.
! ====================================================================

      qbreg = qbreg + qb
      zbar1rg = zbar1rg + zbar1*area/(pixsiz*pixsiz)

! ====================================================================
! Write results.
! ====================================================================

      if(iprn(82).eq.1)&
        write(82,1000) i,ic,zbar1,zbar,capsum*3600000.,&
                       gwtsum*3600000.,etwtsum*3600000.,&
                       qb/area*3600000.,rzpsum,tzpsum

! ====================================================================
! Format statements.
! ====================================================================

1000  format(2i5,2f7.3,4f10.5,2f7.3)

      return

      end
