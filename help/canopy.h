      integer i,ipix,i_und,i_moss,intstp_moss,istmst_moss
      integer istorm_moss,intstm_moss
      integer intstp,istmst,istorm,intstm

      real*8 wc,wc_us,wcip1,wcip1_us,Swq,fw,Swq_us,fw_us
      real*8 wsc,wsc_us,dc
      real*8 epetw,dc_us,epetw_us,dc_moss
      real*8 epet_moss,epwms,epwms_us,pnet,pptms,precip_o
      real*8 dswc,wcrhs,canclos,precip_u,dswc_us,wcrhs_us,Outflow_us
      real*8 PackWater_us,SurfWater_us,xintst_moss,endstm,zmoss
      real*8 xintst,Outflow,PackWater,SurfWater,rnpet,xlepet,hpet
      real*8 gpet,rnetd,xled,hd,gd,rnetw,xlew,hw,gw,tkpet,tkmidpet
      real*8 dspet,tkd,tkmidd,dshd,tkw,tkmidw,dshw,rnpet_us,xlepet_us
      real*8 hpet_us,gpet_us,rnetd_us,rnetw_us,xled_us,xlew_us,hd_us
      real*8 hw_us,gd_us,gw_us,tkpet_us,tkmidpet_us,dspet_us,tkd_us,tkw_us
      real*8 tkmidd_us,tkmidw_us,dshd_us,dshw_us,rnpet_moss,gpet_moss
      real*8 hpet_moss,rnet_pot_moss,g_p_moss,h_p_moss,xlepet_moss
      real*8 xle_p_moss,tskinpet_moss,tskin_p_moss,tkpet_moss,tk_p_moss
      real*8 tkmidpet_moss,tkmid_p_moss,dspet_moss,ds_p_moss
      real*8 zero,one,two,three,four,five,six,dt,dummy
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
