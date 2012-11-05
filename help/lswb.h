      integer iprn(MAX_FIL),i,ncatch,npix,ic,MODE
      integer ipix,landc,nlcs,ivgtyp(MAX_VEG)
      integer veg(MAX_PP1,MAX_VST)
	  integer count24

      real*8 nlakpix,nvegpix
      real*8 qb24sum
      real*8 r_lakearea(MAX_CAT),pixsiz,ettotrg,etlakesumrg
      real*8 etstsumrg,etwtsumrg,fbsrg,etbssumrg,etdcsumrg
      real*8 etwcsumrg,pptsumrg,pnetsumrg,qsurfrg,sxrtotrg
      real*8 xixtotrg,contotrg,ranrunrg,conrunrg,qbreg,gwtsumrg
      real*8 grzsumrg,gtzsumrg,capsumrg,difrzsumrg,dswcsum
      real*8 wcrhssum,dsrzsum,rzrhssum,dstzsum,tzrhssum,dssum
      real*8 svarhssum,rzsmav,tzsmav,rnpetsum,xlepetsum,hpetsum
      real*8 gpetsum,dshpetsum,tkpetsum,tkmidpetsum,rnsum,xlesum
      real*8 hsum,gsum,dshsum,tksum,tkmidsum,tkdeepsum,fwreg
      real*8 wcip1sum,zbar1rg,pr3sat,perrg2,pr2sat,pr2uns,perrg1
      real*8 pr1sat,pr1rzs,pr1tzs,pr1uns,persxr,perixr,persac
      real*8 peruac,perusc,wcsum,zbarrg,f_lake(MAX_PIX),tot
      real*8 veg_pdf(MAX_PIX,MAX_VST),vegd,Swqsum,Swq_ussum
      real*8 zero,one,two,three,four,five,six
      real*8 Sdepthsum,Sdepth_ussum
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
