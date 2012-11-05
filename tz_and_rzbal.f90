! ====================================================================
!
!			subroutine tz_and_rzbal
!
! ====================================================================
!
! Calculate the root zone and transmission zone water balances.
!
! ====================================================================
!
! Created by Wade Crow 3/1/1999
! 
! ====================================================================
 
      subroutine tz_and_rzbal(i,newstorm,inc_frozen,ikopt,ivgtyp,dt,&

! Surface zone and transmission zone soil moistures 

       rzsm,tzsm,rzsm_u,tzsm_u,rzsm1,tzsm1,&

! Layer geometry

       zrz,ztz,zw0,zrzmax,&

! Water fluxes

       evtact,evtact_us,bsdew,dewrun,grz,gtz,diftz,difrz,&
       satxr,runtot,xinact,cuminf,&

! Soil parameters

       ff,thetar,thetas,bcbeta,xk0,psic,&

! Snow
       Swq,Swq_us,&

! Understory/moss

       dc,i_und,i_moss,fw,dc_us,fw_us,evrz_moss,f_und,&

! Storage changes

       dstz,dsrz,&

! Summations of fluxes

       tzrhs,rzrhs)

      implicit none
      include "help/tz_and_rzbal.h"
      integer :: test_flag
      
 
      rzsm0 = rzsm
      tzsm0 = tzsm
      rzsm1old = rzsm
      rzsm1new = rzsm
      tzsm1new = tzsm
      tzsm1old = tzsm
      evtran_rz = zero
      evtran_tz = zero
      cor_flx_tz = zero
      cor_flx_rz = zero
      ddifrzdth1 = zero
      ddifrzdth2 = zero
      dgrzdth1 = zero
      dgrzdth2 = zero
      dgtzdth1 = zero
      dgtzdth2 = zero
      ddiftzdth1 = zero
      ddiftzdth2 = zero
      dewrun = zero
      iter = 1
      max_iter = 10
      num_exp_iter = 60

! ====================================================================
! Calculate change of hydrauli! conductivity with depth
! ====================================================================

      if(ikopt.eq.1)then

        xksrz=xk0
        xkstz=xk0

      else

        xksrz=xk0*dexp(-ff*(zrz/two))
        xkstz=xk0*dexp(-ff*(zrz+(ztz/two)))

      end if

! ====================================================================
! Calculate the second derivative of soil moisture with 
! respect to time.
! ====================================================================

      call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u) 

      call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm, &
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

      call clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,ikopt,&
       xksrz,xkstz,ff,zrz,ztz,bcbeta,thetas,thetar,psic,tzsm)

      dgrzdth1= clcdg(rzsm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)
      dgtzdth2= clcdg(tzsm,ikopt,xkstz,ff,ztz,bcbeta,thetas,thetar)

      if (zrz.gt.zero) then
	
        d2rzsmdt2 = (-dgrzdth1 - ddifrzdth1)*(-grz - difrz + xinact - evtran_rz)*(dt*dt)/(zrz*zrz)

      else

        d2rzsmdt2 = zero

      endif

      if (ztz.gt.zero) then

        d2tzsmdt2 = (ddifrzdth2 -dgtzdth2- ddiftzdth2)*(grz + difrz - gtz - diftz - evtran_tz)*(dt*dt)/(ztz*ztz)

      else

        d2tzsmdt2 = zero

      endif

! ====================================================================
! Modify the numerical approach based on linearity of problem.
! ====================================================================
  
      if (abs(d2rzsmdt2).lt..15d0.and.abs(d2tzsmdt2).lt..15d0) then
 
        explicit_weight = .5d0
        implicit_weight = .5d0
        non_linear_flag = 0
        tol = .001d0

      else if (abs(d2rzsmdt2).gt.5.or.abs(d2tzsmdt2).gt.5) then
        explicit_weight = zero
        implicit_weight = one
        tol = .0001d0
        non_linear_flag = 1

      else

        explicit_weight = zero
        implicit_weight = one
        tol = .0001d0
        non_linear_flag = 0
      endif

! ====================================================================
! Turn on and off numerical testing procedure
        test_flag = 0
! ====================================================================

! ====================================================================
! Identify created and destroyed surface and transmission zones
! - define correction fluxes to account for created water
! don't change the order of these.
! ====================================================================

! ====================================================================
! New surface zone created.
! ====================================================================


      if (zw0.le.zero.and.zrz.gt.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_rz = -zrz*(thetas-thetar) 

      endif

! ====================================================================
! New transmission zone created
! ====================================================================

      if (zw0.le.zrzmax.and.ztz.gt.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_tz = -ztz*(thetas-thetar)

      endif

! ====================================================================
! Transmission zone destroyed
! ====================================================================

      if (zw0.gt.zrzmax.and.ztz.le.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_rz = (zw0-zrzmax)*(thetas-thetar)

      endif

! ====================================================================
! Surface zone or surface and transmission zone destroyed.
! ====================================================================

      if (zw0.gt.zero.and.zrz.le.zero.and.i.gt.1.and.newstorm.eq.0) then

        cor_flx_rz = zw0*(thetas-thetar)

      endif


!================================================
! shut off correction fluxes
!=================================================

      cor_flx_rz = 0
      cor_flx_tz = 0

      
! ====================================================================
! Calculate evaporation and transpiration from surface and transmission
! zone.
! ====================================================================

      call clc_evrz(evtran_rz,Swq,Swq_us,ivgtyp,evtact,dc,i_und,&
       i_moss,fw,evtact_us,dc_us,fw_us,evrz_moss,dummy,f_und)

      if (ivgtyp.eq.2) then

        evtran_tz=evtact*dc*(1-fw)

      else

        evtran_tz  = 0

      endif

! ====================================================================
! Decide if dew goes into surface surface.
! ====================================================================

      if (ivgtyp.eq.0) then 

        dewrz = bsdew

      else

        dewrz = zero

      endif

! ====================================================================
! Case I - water table is at surface
! ====================================================================
 
      if (zrz.eq.zero) then

        case_flag = 1

        rzsm1 = thetas
        tzsm1 = thetas
        difrz = zero
        diftz = zero
        grz = zero
        gtz = zero
        evtran_rz = zero
        evtran_tz = zero
        dewrun = dewrz
        dewrz = zero 
        satxr = satxr + dewrun
        runtot = runtot + dewrun 
        
      endif


! ====================================================================
! Case II - water table is in transmission zone
! ====================================================================

      if (ztz.gt.zero.and.non_linear_flag.eq.0) then

        case_flag = 2.1

500     dgrzdth1= clcdg(rzsm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)
        dgtzdth2= clcdg(tzsm,ikopt,xkstz,ff,ztz,bcbeta,thetas,thetar)
        dgrzdth2 = zero
        dgtzdth1 = zero 

        call clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,ikopt,xksrz,xkstz,ff,&
       zrz,ztz,bcbeta,thetas,thetar,psic,tzsm)

        ddiftzdth1 = zero 

        call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)
          
        call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

        srzflx = dt * (dewrz + xinact - difrz - grz - evtran_rz)
        stzflx = dt * (difrz - diftz  + grz - gtz - evtran_tz)
     
        F11 = rzsm0-rzsm1new+srzflx/zrz
        F22 = tzsm0-tzsm1new+stzflx/ztz

        dF1dtheta1 = dt/zrz*(-ddifrzdth1 - dgrzdth1) - one
        dF1dtheta2 = dt/zrz*(-ddifrzdth2 - dgrzdth2) 
        dF2dtheta1 = dt/ztz*(ddifrzdth1 + dgrzdth1 - dgtzdth1 - ddiftzdth1)
        dF2dtheta2 = dt/ztz*(ddifrzdth2 + dgrzdth2 - dgtzdth2 - ddiftzdth2) - one
        
! --------------------------------------------------------------------&
! Solve 2 x 2 matrix */&
!
! dF1dtheta1(rzsm)   dF1dtheta2(tzsm)           del_rzsm          F1(rzsm)
! dF2dtheta1(rzsm)   dF2dtheta2(tzsm)           del_tzsm          F2(tzsm)
! --------------------------------------------------------------------&

        del_rzsm = one/(dF2dtheta1-(dF1dtheta1*dF2dtheta2/dF1dtheta2))*(F22 - (dF2dtheta2*F11/dF1dtheta2))
        del_tzsm = (F11 - dF1dtheta1*del_rzsm)/dF1dtheta2
        rzsm1new = -(del_rzsm) + rzsm1old
        tzsm1new = -(del_tzsm) + tzsm1old

        if (rzsm1new.gt.thetas) then

          rzsm1new=thetas

        endif

        if (rzsm1new.lt.thetar) then

          rzsm1new=thetar
!         write(*,*) rzsm1new

        endif

        if (tzsm1new.gt.thetas) then

          tzsm1new=thetas

        endif

        if (tzsm1new.lt.thetar) then

          tzsm1new=thetar

        endif

        if ((((rzsm1new.lt.rzsm1old-tol).or.(rzsm1new.gt.rzsm1old+tol)).and.(iter.le.max_iter)).or.&
       (((tzsm1new.lt.tzsm1old-tol).or.(tzsm1new.gt.tzsm1old+tol)).and.(iter.le.max_iter)))   then

          rzsm=implicit_weight*rzsm1new + explicit_weight*rzsm0
          tzsm=implicit_weight*tzsm1new + explicit_weight*tzsm0
          iter=iter+1
          rzsm1old = rzsm1new
          tzsm1old = tzsm1new
          goto 500

        endif

        rzsm1 = rzsm0 + dt/zrz*(dewrz -difrz - grz + xinact - evtran_rz)
        tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

      endif

! ====================================================================
! Case III - water table is in surface zone.
! ====================================================================

      if (zrz.gt.zero.and.zrz.lt.zrzmax.and.non_linear_flag.eq.0) then

        case_flag = 3.1

        diftz = zero
        gtz = zero
        evtran_tz  = zero
        tzsm = thetas

600     dgrzdth1= clcdg(rzsm,ikopt,xksrz,ff,zrz,bcbeta,thetas,thetar)

        call clcddif(ddifrzdth1,ddifrzdth2,ddiftzdth2,rzsm,ikopt,xksrz,xkstz,ff,zrz,ztz,bcbeta,thetas,&
       thetar,psic,tzsm)

        call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

        call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

        srzflx = dt * (dewrz + xinact - difrz - grz - evtran_rz)
        F11 = rzsm0-rzsm1new+srzflx/zrz    
        dF1dtheta1 = dt/zrz*(-ddifrzdth1 - dgrzdth1) - one
        rzsm1new = rzsm1old - F11/dF1dtheta1       

        if (rzsm1new.gt.thetas) then

          rzsm1new=thetas

        endif

        if (rzsm1new.lt.thetar) then

          rzsm1new=thetar

        endif

        if ((((rzsm1new.lt.rzsm1old-tol).or.(rzsm1new.gt.rzsm1old+tol)).and.(iter.le.max_iter))) then

          rzsm=implicit_weight*rzsm1new + explicit_weight*rzsm0
          iter=iter+1
          rzsm1old = rzsm1new
          goto 600

        endif

        rzsm1 = rzsm0 + dt/zrz*(dewrz -difrz - grz + xinact - evtran_rz)
        tzsm1 = thetas

      endif

! ====================================================================
! Back-up Numerics.
! ====================================================================
    
      if (tzsm1.gt.thetas.or.tzsm1.lt.thetar.or.rzsm1.lt.thetar.or.rzsm1.gt.thetas.or.&
       iter.gt.max_iter.or.non_linear_flag.eq.1.or.test_flag.eq.1) then
        
        grz_sum = zero
        difrz_sum = zero
        srzflx = zero
        gtz_sum = zero
        diftz_sum = zero
        stzflx = zero
        tzsm = tzsm0
        rzsm = rzsm0

! --------------------------------------------------------------------&
! Case II - water table is in tranmission zone
! --------------------------------------------------------------------&

        if (ztz.gt.zero) then

          case_flag = 2.2

          do 700 num=1,num_exp_iter

            call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

            call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

            dtsrzflx = dt * (dewrz + xinact - difrz  - grz - evtran_rz)
            srzflx = srzflx + dtsrzflx*(num_exp_iter**(-one))
            difrz_sum = difrz_sum + difrz*(num_exp_iter**(-one))
            grz_sum = grz_sum + grz*(num_exp_iter**(-one))
            rzsm = rzsm0 + srzflx/zrz
            dtstzflx = dt * (difrz  + grz - diftz - gtz - evtran_tz)
            stzflx = stzflx + dtstzflx*(num_exp_iter**(-one))
            diftz_sum = diftz_sum + diftz*(num_exp_iter**(-one))
            gtz_sum = gtz_sum + gtz*(num_exp_iter**(-one))            
            tzsm = tzsm0 + stzflx/ztz

            if (rzsm.gt.thetas) then

              rzsm=thetas

            endif

            if (rzsm.lt.thetar) then

              rzsm=thetar

            endif

            if (tzsm.gt.thetas) then

              tzsm=thetas

            endif
            if (tzsm.lt.thetar) then

              tzsm=thetar

            endif

 700      continue
          difrz =  difrz_sum
          grz =  grz_sum
          diftz = diftz_sum
          gtz = gtz_sum
          rzsm1 = rzsm0 + dt/zrz*(dewrz - difrz - grz + xinact - evtran_rz)
          tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

        endif
     
! --------------------------------------------------------------------&
! Case III - water table is in surface zone.
! --------------------------------------------------------------------&

        if (zrz.gt.zero.and.zrz.lt.zrzmax) then

          case_flag = 3.2

          diftz = zero
          gtz = zero
          evtran_tz  = zero 
          tzsm = thetas

         do 800 num=1,num_exp_iter

            call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,rzsm,tzsm,&
       thetas,thetar,bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz)

            call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
       thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,tzsm,tzsm_u)

            dtsrzflx = dt * (dewrz + xinact - difrz  - grz - evtran_rz)
            srzflx = srzflx + dtsrzflx*(num_exp_iter**(-one))
            difrz_sum = difrz_sum + difrz*(num_exp_iter**(-one))
            grz_sum = grz_sum + grz*(num_exp_iter**(-one))
            rzsm = rzsm0 + srzflx/zrz

            if (rzsm.gt.thetas) then

              rzsm=thetas

            endif

            if (rzsm.lt.thetar) then

              rzsm=thetar

            endif

 800      continue

          difrz =  difrz_sum
          grz =  grz_sum
          diftz = zero
          gtz = zero
          rzsm1 = rzsm0 + dt/zrz*(dewrz - difrz - grz + xinact - evtran_rz)
          tzsm1 = thetas

        endif

      endif
 
! ====================================================================
! Balance checks and corrections.
! ====================================================================

      if (rzsm1.gt.thetas) then
      
        grz = grz + (rzsm1-thetas)*zrz/dt

        if (ztz.gt.zero) then

          tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

        endif

        rzsm1 = thetas

      endif

      if (tzsm1.gt.thetas) then

        gtz = gtz + (tzsm1-thetas)*ztz/dt
        tzsm1=thetas

      endif
   
      if (rzsm1.lt.thetar) then

        difrz = difrz + (rzsm1-thetar)*zrz/dt

        if (ztz.gt.zero) then

           tzsm1 = tzsm0 + dt/ztz*(difrz + grz - gtz - diftz - evtran_tz)

        endif

        rzsm1=thetar

      endif

      if (tzsm1.lt.thetar) then

        diftz = diftz + (tzsm1-thetar)*ztz/dt
        tzsm1=thetar

      endif

! ====================================================================
! Calculate storage changes (used to check water balance in unit 95).
! ====================================================================

      dstz = ztz*(tzsm1-tzsm0)
      dsrz = zrz*(rzsm1-rzsm0)

! ====================================================================
! Sum fluxes (used to check water balance in unit 95).
! ====================================================================

      if (ztz.gt.zero) then

        rzrhs = dt*(xinact-difrz + dewrz-evtran_rz-grz)
        tzrhs = dt*(grz-gtz-evtran_tz+ difrz - diftz )

      else
      
        rzrhs = dt*(xinact-difrz + dewrz-evtran_rz-grz)
        tzrhs = dt*(-gtz-evtran_tz-diftz )
 
      endif

! ======================================================================
! Corrections for "created" and "destroyed" layers
! ======================================================================

      difrz = difrz + cor_flx_rz/dt
      diftz = diftz + cor_flx_tz/dt

! ====================================================================
! Update cumulative infiltration.
! ====================================================================

      cuminf = cuminf + xinact * dt
 
      rzsm = rzsm0
      tzsm = tzsm0

      return
      end


