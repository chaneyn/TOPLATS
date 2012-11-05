      integer iprn(MAX_FIL),ivgtyp,i_und,i_moss,i,ilandc,inc_frozen

      real*8 f_und,f_moss,canclos,dc,fw,dc_us,fw_us,etpix
      real*8 evtact,ettot,ettotrg,epwms,evtact_us,epwms_us
      real*8 evtact_moss,etstore,etwt,etstsum,etstsumrg
      real*8 etwtsumrg,etwtsum,etbssum,etbssumrg,etdcsum
      real*8 etdcsumrg,etwcsum,etwcsumrg,etlakesumrg,etlakesum
      real*8 bsdew,contot,contotrg,pptsum,pptms,pptsumrg
      real*8 pnetsum,pnet,pnetsumrg,qsurf,runtot,qsurfrg
      real*8 sxrtot,satxr,sxrtotrg,xixtot,xinfxr,xixtotrg
      real*8 ranrun,ranrunrg,conrun,conrunrg,wcip1sum
      real*8 wcip1,dswcsum,dswc,wcrhssum,wcrhs,dsrzsum,dsrz
      real*8 rzrhssum,rzrhs,dstzsum,dstz,tzrhssum,tzrhs,zrz
      real*8 gwt,grz,gtz,ztz,gwtsum,gwtsumrg,grzsumrg,gtzsumrg
      real*8 capsum,capsumrg,difrzsumrg,diftz,difrz
      real*8 dstore,dssum,svarhs,rzsm1_u,tzsm1_u,rzsm1,tzsm1
      real*8 rzsmav,tzsmav,tzpsum,thetas,rzpsum,r_mossm,rzsm
      real*8 tzsm,rnact_moss,xleact_moss,hact_moss,gact_moss
      real*8 dshact_moss,tskinact_moss,tkact_moss,tkmid_moss
      real*8 rnact_us,xleact_us,hact_us,gact_us,dshact_us
      real*8 tkact_us,tkmid_us,rnact,xleact,hact,gact,dshact
      real*8 tkact,tkmid,rnsum,xlesum,hsum,gsum,dshsum,tksum
      real*8 tkmidsum,rnpet,rnpet_us,rnpet_moss,xlepet,xlepet_us
      real*8 xlepet_moss,hpet,hpet_us,hpet_moss,gpet,gpet_us
      real*8 gpet_moss,dspet,dspet_us,dspet_moss,tkpet,tkpet_us
      real*8 tkpet_moss,tkmidpet,tkmidpet_us,tkmidpet_moss
      real*8 rnpetsum,xlepetsum,hpetsum,gpetsum,dshpetsum
      real*8 tkpetsum,tkmidpetsum,tkdeepsum,Tdeepstep,dt
      real*8 svarhssum,Swqsum,Swq_ussum,Swq,Swq_us
      real*8 zero,one,two,three,four,five,six,rescale,etpixloc
      real*8 conpix
      real*8 difwt
      real*8 tair,xlhv,dummy
      real*8 Sdepthsum,Sdepth_ussum,Sdepth,Sdepth_us
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/

