      integer ibeginhour,ibeginday,ibeginmonth,ibeginyear
      integer iihour,iiday,iimonth,iiyear,ileapyear
      integer idaysperyear,imonthsinyear
      integer inewday,inewmonth,inewyear,jday,loop
      integer iendofmonth(12)
      real*8 timestep

      data imonthsinyear /12/,&
           iendofmonth /31,59,90,120,151,181,212,243,273,304,334,365/
 



 !Descriptions:
 !
 !     i:          	current time step
 !     ibeginday:   	day simulation begins 
 !     ibeginhour:  	hour simulation begins
 !     ibeginmonth: 	month simulation begins
 !     ibeginyear:  	year simulation begins
 !     iiday:        	day of year in simulation
 !     iendofmonth:     the julian day at end of each month
 !     itimestep:       timestep for simulation in hours
 !     iihour:       	hour of day in simulation (0-23)
 !     ileapyear:   	flag for leapyear dates
 !     imonthsinyear:   number of months in year (12)
 !     iimonth:      	month of year (1-12) in simulation
 !     inewday:         flag for newday
 !     inewyear:        flag for newyear
 !     iiyear:       	year of simulation
 !     jday:       	julian day
 !
 !***************************************************************
