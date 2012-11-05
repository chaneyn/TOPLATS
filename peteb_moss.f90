! ====================================================================
!
!               subroutine peteb_moss
!
! ====================================================================
!
! Solve the energy balance at potential rate for the moss layer.
!
! ====================================================================

      subroutine peteb_moss(thermc1,thermc2,vpsat_ic,&
       heatcap_moss,heatcap1,heatcap2,heatcap_us,rs_under,rain,snow,&
       tskin_p_moss,tk_p_moss,tkmid_p_moss,thermc_moss,Tdeepstep,&
       xle_act_moss,h_p_moss,g_p_moss,&
       xle_p_moss,ds_p_moss,tkel_ic,rav_moss,rah_moss,r_moss_depth,&
       zmid,zdeep,eps,dt,toleb,maxnri,rld,alb_moss,&
       rnet_pot_moss,vppa_ic,iopstab,emiss_moss,&
       roa_ic,psychr_ic,xlhv_ic,epet_moss,i,tskinact_moss,zww,&
       za,uzw,zpd_moss,press,z0m_moss,rib_moss,z0h_moss,&
       PackWater_us,SurfWater_us,VaporMassFlux_us,TPack_us,&
       TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us,&
       hact_snow_us,rn_snow_us,dens_us,RaSnow,alb_snow,&
       appa,tcel_ic,gact,epetw_us,row,Swq_us)

      implicit none
      include "SNOW.h"
      include "help/peteb_moss.h"

      return

      end
