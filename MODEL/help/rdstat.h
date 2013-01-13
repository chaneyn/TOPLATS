      integer natb,nlandc,npix,ii,jj,ivgtyp(MAX_VEG)
      integer icatch(MAX_PIX),ncatch,icount(MAX_PIX),nlcs
      integer veg(MAX_PP1,MAX_VST)

      real*8 atb(MAX_PP1,MAX_ATB),atb_pdf(MAX_PP1,MAX_ATB)
      real*8 veg_pdf(MAX_PP1,MAX_VST),f_lake(MAX_PIX)
      real*8 fbs(MAX_CAT+1),fbsrg,xlamda(MAX_CAT)
      real*8 pixlamda(MAX_PIX)

      real*8 fin(MAX_CAT),ntot
