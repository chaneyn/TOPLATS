      integer i,ncatch

      real*8 djday,dt
      real*8 fwreg,rzsmav,tzsmav,wcsum,wcip1sum
      real*8 ettotrg,etstsumrg,etwtsumrg,etbssumrg,etdcsumrg
      real*8 etwcsumrg,etlakesumrg
      real*8 pptsumrg,pnetsumrg,contotrg
      real*8 sxrtotrg,xixtotrg,qsurfrg,ranrunrg,conrunrg,qbreg
      real*8 capsumrg,difrzsumrg,gwtsumrg,grzsumrg,gtzsumrg
      real*8 zbarrg,zbar1rg
      real*8 dswcsum,dsrzsum,dstzsum,dssum,wcrhssum,rzrhssum
      real*8 tzrhssum,svarhssum
      real*8 rnsum,xlesum,hsum,gsum,tksum,dshsum,tkmidsum,tkdeepsum
      real*8 rnpetsum,xlepetsum,hpetsum,gpetsum,tkpetsum
      real*8 tkmidpetsum,dshpetsum
      real*8 perrg1,perrg2,pr3sat,pr2sat,pr2uns,pr1sat,pr1rzs
      real*8 pr1tzs,pr1uns,persac,peruac,perusc,persxr,perixr
      real*8 ettot(MAX_CAT),etstsum(MAX_CAT),etwtsum(MAX_CAT)
      real*8 etbssum(MAX_CAT),etdcsum(MAX_CAT),etwcsum(MAX_CAT)
      real*8 etlakesum(MAX_CAT),sxrtot(MAX_CAT),conrun(MAX_CAT)
      real*8 contot(MAX_CAT),pptsum(MAX_CAT),pnetsum(MAX_CAT)
      real*8 xixtot(MAX_CAT),qsurf(MAX_CAT),ranrun(MAX_CAT)
      real*8 zbar(MAX_CAT),zbar1(MAX_CAT),capsum(MAX_CAT),gwtsum(MAX_CAT)
      real*8 rzpsum(MAX_CAT),tzpsum(MAX_CAT),fwcat(MAX_CAT)
      real*8 Swqsum,Swq_ussum,Sdepthsum,Sdepth_ussum
      real*8 zero,one,two,three,four,five,six
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/


      integer kk
