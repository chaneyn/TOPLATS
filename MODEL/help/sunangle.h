      integer idaylight,idaysperyear
      integer minpdeg,minphour
 
      real*8 waspt,wslop
      real*8 albgen,B,cellaspct,cellslope,declination,degprad
      real*8 adjsol,adjsol2,adjsol3,anglehour
      real*8 adjdif,difskyview,eqnoftime,endhour,halfdaylength
      real*8 radpdeg,radphour,rcosncidenangl,rcosihlfdaylng
      real*8 reflect,reflctsky,rlatitude,rlng_merid
      real*8 rlongitude,rlongitudadjst,rnoonhour,rsinsolaralt
      real*8 solaraltitude,solarazimuth,solartimestp,solarzenith
      real*8 starthour,sunearthdist,sunrise,sunset,sunmax,timeadjust

      integer iihour,iiday,iimonth,iiyear,jday,kk
     
      real*8 averag,dt,timestep


!      real*8 PI
!      data PI, minpdeg, minphour/3.1415927,4, 60/,&
!           albgen /0.3/

!
!
! Descriptions:
!      albgen:          general value for the albedo of a pixel
!      B:               coefficient for equation of time
!      cellaspct:       slope aspect of pixel in radians
!      cellslope:       slope angle of pixel in radians
!      declination:     solar declination (radians)
!      degprad:         degrees per radian (180/PI)
!      adjsol:          adjusted direct solar radiation striking pixel
!      adjdif:          adjusted diffuse solar radiation striking pixel
!      difskyview:      sky view factor for diffuse radiation (0-1)
!      eqnoftime:       adjustment for equation of time (min)
!      endhour:         mid-point of current solar hour (hr)
!      halfdaylength:   half day length (rad)
!      anglehour:       angle of current 'halfhr' from solar noon (rad)
!      idaylight:       flag for daylight:
!                        1 = true
!                        0 = false
!      lat_deg:         2 digit central latitude of site in degree
!      lat_min:         2 digit central latitude of site in minutes
!      lng_deg:         2 digit central longitude of site in degree
!      lng_min:         2 digit central longitude of site in minutes
!      lng_mer:         input standard meridian for central time zone of region
!      minpdeg:        minutes per degree of rotation (min)
!      minphour:       minutes per hour (min = 60)
!      PI:              trigonmetri! constant (3.1415927)
!      radpdeg:         radians per degreee (PI/180)
!      radphour:        radians per hour (radpdeg*degphour)
!      rcosncidenangl:    cosine of the incidence angle between solar
!                        rays and the normal to the surface
!      rcosihlfdaylng:    cosine of half day length
!      reflect:         incoming reflected solar radiation (W/m^2)
!      reflctsky:       view factor for incoming reflected solar radiation
!                        (0-1)
!      rlatitude:       central latitude for the basin (radians)
!      rlng_merid:      calculated standard meridian for central time zone
!                        of region
!      rlongitude:      central longitude for the basin (radians)
!      rlongitudadjst:  adjustment for longitude (min)
!      rnoonhour:       true solar noon (hr)
!      rsinsolaralt:    sine of sun's solar altitude
!      solaraltitude:   solar altitude of sun from the horizon (radians)
!      solarazimuth:    solar azimuth angle of sun from east (radians)
!      solartimestp:    fraction of the timestep the sun is above the horizon
!      solarzenith:     solar zenith angle of sun from east (radians)
!      starthour:       correct hour in solar time (hr)
!      sunearthdist:    distance from the sun to the earth
!      sunrise:         time of sunrise (hr)
!      sunset:          time of sunset (hr)
!      sunmax:          calculated solar radiation at the top of atmosphere
!                        (W/m^2)
!      timeadjust:      required time adjustment to convert local time to solar
!                        time
!
!
!
! ***************************************************************
