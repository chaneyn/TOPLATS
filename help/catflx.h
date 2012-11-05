      integer iprn(MAX_FIL),i,ic

      real*8 area,pixsiz,r_lakearea,ettot,etstsum,etwtsum,etlakesum
      real*8 etbssum,fbs,etdcsum,etwcsum,pptsum,pnetsum,contot
      real*8 qsurf,sxrtot,xixtot,ranrun,conrun,gwtsum,capsum,tzpsum
      real*8 rzpsum,fwcat
      real*8 s_nr_etwtsum,s_nr_gwtsum,s_nr_capsum
      real*8 s_nr_tzpsum,s_nr_rzpsum
      real*8 zero,one,two,three,four,five,six

      real*8 catpix,catlakpix,catvegpix
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
            3.d0,4.d0,5.d0,6.d0/
