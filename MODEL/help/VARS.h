! ==========================================================
! ==========================================================
! VARIABLE	DESCRIPTION				UNIT
! ==========================================================
! ==========================================================
!
! a_ice		Multiplying coefficient in the
!		calculation of the soil frozen
!		water content				-
! a_ice_moss	Multiplying coefficient in the
!		calculation of the moss frozen
!		water content				-
! adjdif        Adjusted diffuse radiation		W/m2
! adjsol	Adjusted direct solar radiation		W/m2
! albd		Dry surface albedo			-
! albgen        General albedo value in solar adjust 	- 
! alb_moss	Moss albedo				-
! albd_us	Dry surface albedo for under story	-
! alb_snow	Snow albedo				-
! albw		Wet surface albedo			-
! albw_us	Wet surface albedo for under story	-
! area		Area of the current catchment		m2
! amp		Amplitude of the deep soil temperaure
!		cosine wave				K
! anglehour	Angle of 'halfhr' from solar noon	radians
! appa		Air pressure				Pa
! atanb		Soil topographi! index for each pixel
!		in the region				-
! atb		Topographi! index for each pixel in
!		the region				-
! avpx		Air vapor pressure			Pa
! averag	Value for met. value (missch.f)		-
! B		Coefficient for time eqn. (sunangle.f)	-
! basink	Basin scale hydrauli! conductivity	m/s
! b_ice		Exponential coefficient in the
!		calculation of the soil frozen
!		water content				-
! b_ice_moss	Exponential coefficient in the
!		calculation of the moss frozen
!		water content				-
! bcbeta	Brooks-Corey '64 pore size
!		distribution index			-
! bcgamm	2.0+3.0*bcbeta				-
! bsdew		Condensation on bare soil		m/s
! bsdew_moss	Condensation onto moss			m/s
! bulk_dens	Bulk density of the soil		kg/m3
! bulk_dens_moss	Bulk density of the moss	kg/m3
! canclos	Canopy closure of the over story	-
! capflx	Rate of capillary rise from the water
!		table to the root zone			m/s
! capmax	Maximum rate of capillary rise from
!		the water table to the root zone	m/s
! capsum	Sum of capillary flux for each
!		catchment				m/s
! capsumrg	Sum of capillary flux for the region	m/s
! catlakpix	Number of lake pixels in the current
!		catchment				-
! catpix	Number of pixels in the current
!		catchment				-
! catvegpix	Number of vegetated pixels in the
!		catchment				-
! c! 		Gravity term in philip's infiltration
!		equation				-
! cc!		Canopy closure, temporally variable
! 		to make the equations shorter (ccc
!		instead of canclos(ilandc(ipix))	-
! cellaspct	Slope aspect of pixel 			radians
! cellslope	Slope angle of pixel 			radians
! conpix	Condensation rate for the pixel		m/s
! conrun	Rate of runoff due to condensation for
!		each catchment				m/s
! conrunrg	Rate of runoff due to condensation for
!		the region				m/s
! contot	Average condensation rate for the
!		current catcment			m/s
! contotrg	Average condensation rate for the
!		region					m/s
! corr		Constant in the bare soil
!		desorptivity equation
! cp		Specifi! heat of air			J/kgK
! cph2o		Specifi! heat of water			J/kgK
! cumdep	Initial cumulative infiltration		m
! cumexf	Cumulative exfiltration			m
! cuminf	Cumulative infiltration			m
! d!		Zero during condensation and 1 during
!		no condensation				-
! dc_moss	Zero during condensation and 1 during
!		no condensation, for the moss layer	-
! dc_us		Zero during condensation and 1 during
!		no condensation, for the under story	-
! dd		Drainage density			1/m
! ddifrzdth	Derivative of the diffusive flux out
!		of the root zone with respect to soil
!		moisture				m/s
! declination   Solar declination			radians
! degprad	Degrees per radian (180/PI)		(Deg/radian)
! deltrz	Difference between soil moisture at
!		beginning of storm event and
!		residual soil moisture			-
! dens		Density of the snow pack		kg/m3
! dewrun	Runoff rate due to condensation on
!		bare soil				m/s
! dewrz		Condensation forming on bare soil	m/s
! dgrzdth	Derivative of the downward flux out of
!		the root zone with respect to soil
!		moisture				m/s
! djday		Decimal Julian day			-
! difskyview	Sky view factor for diffuse radiation	-
! difrz		Diffusive flux from root zone into
!		transmission zone			m/s
! difrzsum	Catcment average diffusive flux		m/s
! difrzsumrg	Regional average diffusive flux		m/s
! diftz		Diffusive flux from transmission zone
!		towards the water table			m/s
! dshact	Ground heat storage flux at actual
!		evaporation rate			W/m2
! dshactd	Ground heat storage flux at actual
!		evaporation rate for dry canopy		W/m2
! dshactd_moss	Ground heat storage flux at actual
!		evaporation rate for the moss		W/m2
! dshact_moss	Ground heat storage flux at actual
!		evaporation rate for the moss		W/m2
! dshact_us	Ground heat storage flux at actual
!		evaporation rate for the under
!		story					W/m2
! dshactd_us	Ground heat storage flux at actual
!		evaporation rate for dry canopy for
!		the under story				W/m2
! dshd		Ground heat flux at potential
!		evapotranspiration for dry vetation	W/m2
! dshd_us	Ground heat storage at potential
!		evapotranspiration for dry vetation for
!		the under story				W/m2
! dshsum	Average ground heat storage for the
!		region					W/m2
! dshpetsum	Regional average ground heat storage
!		at potential rate			W/m2
! dspet		Ground heat storage at potential
!		evapotranspiration			W/m2
! dspet_moss	Ground heat storage at potential
!		evapotranspiration for the moss layer	W/m2
! dspet_us	Ground heat storage at potential
!		evapotranspiration for the under story	W/m2
! ds_p_moss	Ground heat storage at potential
!		evapotranspiration for the moss layer	W/m2
! dshw		Ground heat storage at potential
!		evapotranspiration for wet vegetation	W/m2
! dshw_us	Ground heat storage at potential
!		evapotranspiration for wet vegetation
!		for the under story			W/m2
! dsrz		Change in storage in the root zone	m
! dsrzsum	Change in storage in the root zone for
!		the region				m
! dssum		Total catchment change in storage
!		above the water table			m
! dstz		Change in storage in the transmission
!		zone					m
! dstzsum	Change in storage in the transmission
!		zone for the region			m
! dsw!		Change of water amount in interception
!		storage					m/s
! dswc_us	Change of water amount in interception
!		storage for the under story		m/s
! dswcsum	Change of water in interception
!		store for the region			m/s
! dtil		Depth to impervious layer		m
! dumds		Dummy variable, used in the calculation
!		of the heat storage			K
! dumtk		Dummy variable, used in the calculation
!		of the soil skin temperature		K
! dumtkmid	Dummy variable, used in the calculation
!		of the soil mid layer temperature	K
! dumtskin	Dummy variable, used in the calculation
!		of the moss skin temperature		
! dt		Time step size				s
! dvpsdt	Slope of vpsat versus air temperature	Pa/C
! dzbar		Change in average water table depth	m
! ebscap	Bare soil exfiltration capacity		m/s
! ebsmf		Actual evaporation mass flux for bare
!		soil					kg/sm2
! ebspot	Potential evaporation forcing for
!		bare soil				m/s
! eqnoftime	Adjustment of eqn. of time		-
! emiss		Emissivity				-
! emiss_moss	Emissivity of the moss layer		-
! emiss_us	Emissivity of the under story		-
! endhour	Midpoint of current solar hour		-
! endstm	Time after a storm period which
!		defines the end of a storm period	s
! eps		Degree of explicitliness of the
!		solution sheme solving for the soil
!		and moss tempertures			-
! eta_a		Decline of solar radiation with depth
!		in the lake				1/m
! etbssum	Total bare soil evaporation for the
!		current catchment			m/s
! etbssumrg	Total bare soil evaporation for the
!		region 					m/s
! etdcsum	Total evaporation form the dry
!		canopy for the current catchment	m/s
! etdcsumrg	Total evaporation form the dry
!		canopy for the region			m/s
! etlakesum	Total evaporation from the lakes for
!		the current catchment			m/s
! etlakesumrg	Total evaporation from the lakes for
!		the region				m/s
! etpix		Total evapotranspiration rate for the
!		time step for the pixel			m/s
! etstore	Eapotranspiration rate from the
!		soil-vegetation-atmosphere column
!		above the water table			m/s
! etstsum	Total evapotranspiration rate from the
!		soil-vegetation column for the current
!		catchment				m/s
! etstsumrg	Total evapotranspiration rate from the
!		soil-vegetation column for the region	m/s
! ettot		Evapotranspiration in the current
!		catchment				m/s
! ettotrg	Evapotranspiration in the region	m/s
! etwcsum	Total evaporation form the wet
!		canopy for the current catchment	m/s
! etwcsumrg	Total evaporation form the wet
!		canopy for the region			m/s
! etwt		Evapotranspiration rate from the
!		water table				m/s
! etwtsum	Total evapotranspiration rate from
!		water table for the current catchment	m/s
! etwtsumrg	Total evapotranspiration rate from
!		water table for the region		m/s
! epetd		Potential evapotranspiration for dry	m/s
!		canopy
! epetd_us	Potential evapotranspiration for dry	m/s
!		canopy for the under story
! epet_moss	Potential evaporation for the moss
!		layer					m/s
! epetw		Potential evapotranspiration for wet	m/s
!		canopy
! epetw_us	Potential evapotranspiration for wet	m/s
!		canopy for the under story
! epwms		Actual evaporation rate from wet canopy	m/s
! epwms_us	Actual evaporation rate from wet canopy
! evrz		Evaporative demand out of the root
!		zone layer				m/s
! evrz_moss	Evaporative demand out of the moss
!		layer					m/s
!		for the under story			m/s
! evtact	Actual local evaporation/transpiration
!		from bare soil or dry canopy		m/s
! evtact_us	Actual local evaporation/transpiration
!		from bare soil or dry canopy for the
!		under story				m/s
! evtact_moss	Actual local evaporation from the
!		moss layer				m/s
! evtran	Actual local evaporation/transpiration
!		from bare soil or dry canopy		m/s
! evtz		Evaporative demand out of the
!		transmission zone layer			m/s
! exfdif	Exfiltration diffusiviuty		-
! extinct	Extinction of solar radiation in
!		the canopy				-
! f1		Amount of solar radiation reflected
!		back into the atmosphere		-
! f1par		Limiting factor due to PAR to
!		canopy resistance			-
! f1par_us	Limiting factor due to PAR to
!		canopy resistance for under story	-
! f2		Amount of solar radiation reaching
!		the under story				-
! f3		Amount of solar radiation reflected
!		back from the under story towards
!		the over story				-
! f3vpd		Limiting factor due to vapor pressure
!		deficit to canopy resistance		-
! f3vpdpar	Parameter used in the calculation of
!		f3vpd					1/Pa
! f3vpdpar_us	Parameter used in the calculation of
!		f3vpd_us				1/Pa
! f3vpd_us	Limiting factor due to vapor pressure
!		deficit to canopy resistance for under
!		story					-
! f4temp	Limiting factor due to air temperature
!		to canopy resistance			-
! f4temppar	Parameter used in the calculation of
!		f4temp					1/K2
! f4temppar_us	Parameter used in the calculation of
!		f4temp_us				1/K2
! f4temp_us	Limiting factor due to air temperature
!		to canopy resistance for under story	-
! fbs		Fraction of bare soil in the catchment	-
! fbsrg		Fraction of bare soil in the region	-
! ff		Topmodel parameter describing
!		exponential decay of saturated
!		hydrauli! conductivity with depth	1/m
! fnimg		Prefix name for output image		character
! fraci_a	Fraction of lake covered with ice	-
! frcbeta	Fractional coverage paramter relating	h/mm
!		fractional coverage with rain rate :
!		Fr. Cov. = 1 - exp (beta*Rainrate)
!		Rainrate is in mm/h.
! frsoil	Fractional coverage of each soil type
!		in each catchment and in the total area	-
! fw		Fraction of wet canopy			-
! fwcat		Catchment average wet canopy fraction	-
! fwreg		Regional average wet canopy fraction	-
! fw_us		Fraction of wet canopy for the under
!		story					-
! gact		Ground heat flux at actual evaporation
!		rate					W/m2
! gactd		Ground heat flux at actual evaporation
!		rate for dry canopy			W/m2
! gactd_moss	Ground heat flux at actual evaporation
!		rate for the moss			W/m2
! gact_moss	Ground heat flux at actual evaporation
!		rate for the moss			W/m2
! gact_us	Ground heat flux at actual evaporation
!		rate for the under story		W/m2
! gactd_us	Ground heat flux at actual evaporation
!		rate for dry canopy for the under story	W/m2
! gb_max	Maximum acceptable value for ground
!		heat flux when using the
!		Penman-Moneith method			W/m2
! gb_min	Minimum acceptable value for ground
!		heat flux when using the
!		Penman-Moneith method			W/m2
! gbspen	Ground heat flux input into the Penman-
!		Monteith equation			-
! gd		Ground heat flux at potential
!		evapotranspiration for dry vetation	W/m2
! gd_us		Ground heat flux at potential
!		evapotranspiration for dry vetation for
!		the under story				W/m2
! gpet		Ground heat flux at potential
!		evapotranspiration			W/m2
! gpet_moss	Ground heat flux at potential
!		evapotranspiration for the moss layer	W/m2
! gpetsum	Regional average ground heat flux
!		at potential rate			W/m2
! gpet_us	Ground heat flux at potential
!		evapotranspiration for the under story	W/m2
! g_p_moss	Ground heat flux at potential
!		evapotranspiration for the moss layer	W/m2
! grz		Rate of drainage out of the root zone 	m/s
! grzsum	Average rate of drainage out of the
!		root zone for each catchment		m/s
! grzsumrg	Average rate of drainage out of the
!		root zone for the region		m/s
! gsum		Total ground heat flux for all
!		surfaces for a time step for the region W/m2
! gtz		Rate of drainage out of the
!		transmission zone			m/s
! gtzsum	Catchment average rate of drainage
!		out of the transmission zone		m/s
! gtzsumrg	Regional average rate of drainage
!		out of the transmission zone		m/s
! gw		Ground heat flux at potential
!		evapotranspiration for wet vegetation	W/m2
! gw_us		Ground heat flux at potential
!		evapotranspiration for wet vegetation
!		for the under story			W/m2
! gwt		Rate of drainage to water table for
!		each time step				m/s
! gwtsum	Rate of drainage to water table for
!		each time step for each catchment	m/s
! gwtsumrg	Rate of drainage to water table for
!		each time step for the region		m/s
! hact		Sensible heat flux at actual
!		evaporation rate			W/m2
! hactd		Sensible heat flux at actual
!		evaporation rate for dry canopy		W/m2
! hactd_moss	Sensible heat flux at actual
!		evaporation rate for the moss		W/m2
! hact_moss	Sensible heat flux at actual
!		evaporation rate for the moss		W/m2
! hact_us	Sensible heat flux at actual
!		evaporation rate for the under story	W/m2
! hact_snow	Sensible heat flux out of a snow pack	W/m2
! hact_snow_us	Sensible heat flux out of a snow pack on
! 		top of the under story			W/m2
! hactd_us	Sensible heat flux at actual
!		evaporation rate for dry canopy	for
!		the under story				W/m2
! halfdaylength Half day length				radians
! hbar		Areal average water table height from
!		impervious layer			m
! hbar0		Initial depth of ground water table
!		from impervious layer			m
! hd		Sensible heat flux at potential
!		evapotranspiration for dry vetation	W/m2
! hd_us		Sensible heat flux at potential
!		evapotranspiration for dry vetation for
!		the under story				W/m2
! heatcap	Soil heat capacity			J/kgK
! heatcap1	Soil heat capacity of the upper layer	J/kgK
! heatcap2	Soil heat capacity of the lower layer	J/kgK
! heatcap_moss	Heat capacity of the moss layer		J/kgK
! heatcapold	Soil heat capacity of the previous
!		time step				J/kgK
! heatcap_us	Heat capacity of the under story	J/kgK
! hice_a	Ice depth in the lake			m
! hice_in	Initial ice depth in the lake		m
! hpet		Sensible heat flux at potential
!		evapotranspiration			W/m2
! hpet_moss	Sensible heat flux at potential
!		evapotranspiration for the moss layer	W/m2
! hpetsum	Regional average sensible heat flux
!		at potential rate			W/m2
! hpet_us	Sensible heat flux at potential
!		evapotranspiration for the under story	W/m2
! h_p_moss	Sensible heat flux at potential
! hu_max	Maximum acceptable value for water
!		vapor content				depends on inputs
! hu_min	Minimum acceptable value for water
!		vapor content				depends on inputs
! hsnw_a	Snow depth in the lake			m
! hsnw_in	Initial snow depth in the lake		m
!		evapotranspiration for the moss layer	W/m2
! hsum		Total sensible heat flux for all
!		surfaces for a time step for the region	W/m2
! hw		Sensible heat flux at potential
!		evapotranspiration for wet vegetation	W/m2
! hw_us		Sensible heat flux at potential
!		evapotranspiration for wet vegetation
!		for the under story			W/m2
! i		Time step				-
! i_2l		Option to solve the under story
!		equations:
!		0) Incoming long wave radiation for
!		   both layers is equal
!		1) Both layers interact through the
!		   long wave radiation terms
! ibeginday     Day simulation begins			-
! ibeginhour    Hour simulation begins			-
! ibeginmonth   Month simulation begins			-
! ibeginyear    Year simulation begins			-
! i!		Current catchment number		-
! icount	Counter used for calculating
!		fractional coverages for each soil
!		type in each catchment			-
! iday		Day of the timestep being solved	-
! idaylight	Flag for daylight			-
!		0) False
!		1) True
! idaysperyear	Number of days in year (365/366)	-
! idifind	Dffusivity index =
!		(1 + 2 * bcbeta)/bcbeta			-
! iendofmonth   1D array of end of month values		-
! ievcon 	Limiting factor in the actual
!		evaporation/transpiration:
!		1) Soil controlled
!		2) Atmosphere controlled - unsaturated
!		   surface
!		3) Atmosphere controlled - saturated
!		   surface				-
! ievcon_moss 	Limiting factor in the actual
!		evaporation/transpiration for the moss
!		layer
!		1) Soil controlled
!		2) Atmosphere controlled - unsaturated
!		   surface
!		3) Atmosphere controlled - saturated
!		   surface				-
! ievcon_us 	Limiting factor in the actual
!		evaporation/transpiration for the under
!		story
!		1) Soil controlled
!		2) Atmosphere controlled - unsaturated
!		   surface
!		3) Atmosphere controlled - saturated
!		   surface				-
! ifcoarse	Describe the coarse state of the soil:
!		0) No coarse material
!		1) Coarse material			-
! iffroz	Describes thermal state of soil:
!		1) Frozen
!		0) Not frozen				-
! iffroz_us	Describes thermal state of soil under
!		the under story:
!		1) Frozen
!		0) Not frozen				-
! ikopt		Specify the variation of hydraulic
!		conductivity with depth:
!		0) Declines exponentially with depth
!		1) Stays uniform with depth
! ihour		Hour of the time step being solved	-
! iiday         Day of month in the simulation (1-31)	-
! iihour	Hour of day in the simulation (0-23)	-
! iimonth       Month of year in the simulation (1-12)	-
! iiyear        Year of the simulation (1997...)	-
! iland!	Land cover type for the pixel		-
! ileapyear	Flag of leapyear dates			-
! i_moss	Option for a moss layer under the
!		canopy:
!		1) A moss layer is under the canopy
!		0) No moss layer
! img_opt	Option of the image input format:
!		0) Listed with time step
!		1) Listed with year_day_hour.bin
! imgprn	Flag for printing of an output image
!		1) Print it
!		0) Do not print it
! imonthsinyear Number of months in year (12)		-
! inc_frozen	Specify the method of treating frozen
!		soil water
!		0) Treat is as if it were not frozen
!		1) Treat is as solid soil particles
! inewday	Flag for new day			-
! inewmonth	Flag for new month			-
! inewyear	Flag for new year			-
! intstm	State of infiltration of the soil:
!		0) Storm
!		1) Interstorm				-
! intstm_moss	State of infiltration of the soil under
!		a moss cover
!		0) Storm
!		1) Interstorm				-
! intstp	Number of steps into interstorm period	-
! intstp_moss	Number of steps into interstorm period
!		under a moss layer			-
! iopbf		Method of calculating baseflow:
!		0) Sivaplan et al. [1987]
!		1) Troch et al. [1992]
! iopebs	Specify the method of calculation of
!		actual evaporation
!		0) Desorptivity from Famiglietti
!		1) Desorptivity from Eagleson
!		2) Soil resistance			-
! iopflg	Option for storm/interstorm
!		initialization:
!		0) Constant value everywhere
!		1) Image is read in for storm/
!		   interstorm flags, initial soil
!		   moisture and cumulative infiltration/
!		   exfiltration				-
! iopgveg	Calculate ground heat flux under
!		vegation assuming:
!		0) No ground heat flux under vegetation
!		1) Exponential decline with LAI		-
! ioppet	Method of energy balance calculation:
!		2) PET values are externally input
!		1) Using the Penman-Monteith equations
!		0) Solve the energy balance through
!		   itereation for the skin temperature	-
! iopslp	Option for slope angle & aspect adjust	-
!		1) Use weights for adjustment
!		0) Use no weights for adjustment
! iopsmini	Method of initalizing of soil moisture:
!		1) User inputs initial root zone and
!		   transmission zone soil moisture
!		2) Use Brooks/Corey and the initial
!		   water table depth			-
! iopstab	Stabitily correction for aerodynamic
!		resistance:
!		1) Use the Richardson number
!		2) No stability correction		-
! ioptherm!	Method of calculating thermal
!		conductivity:
!		0) Johanssen
!		1) McCumber and Pielke
! iopthermc_v	Calculate ground heat flux under
!		vegation assuming:
!		0) Exponential decline with LAI
!		1) 7 * exponential decline with LAI	-
! ioptlr	Option for temperature Lapse Rate	-
!		0) Use method of Lapse Rate
!		1) Use no Lapse Rate
! iopwc0	Initial canopy storage input:
!		0) single value used for whole area
!		1) An image of initial canopy storage
!		   is used				-
! iopwt0	Input of initial conditions:
!		0) Input the initial water table depth
!		1) Input initial baseflow, derive
!		   initial water table depth
! iopwv		Input form of humidity:
!		1) Relative humidity
!		2) Wet bulb temperature			-
! iouten	Time step for last output of image	-
! ioutsp	Time step spacing for output of image	-
! ioutst	Time step for first output of image	-
! ipix		Pixel number				-
! ipixnum	Pixel number for each row and column
!		in the image (can be 0 if the pixel is
!		not in the catchment)			-
! iprn		Flag wether or not to print an output
!		file
!		1) Print it
!		0) Do not print it			-
! irestype	Method of calculating the soil
!		resistance to evaporation:
!		2) Sun (1982)
!		3) Kondo (1990)
!		4) Camillo and Gurney (1986)
!		5) Passerat (1986)			-
! irntyp	Type of runoff occuring at the pixel:
!		0) No runoff
!		1) Infiltration excess runoff
!		2) Saturation excess runoff		-
! isoil		Soil type of the current pixel		-
! istep		Number of time steps into the current
!		storm/interstorm event			-
! istmst	Number of steps into storm period	-
! istmst_moss	Number of steps into storm period for
!		the soil under a moss layer		-
! istorm	Progress of storm:
!		0) Interstorm period
!		1) Storm period				-
! istorm_moss	Progress of storm under a moss layer:
!		0) Interstorm period
!		1) Storm period				-
! i_und		Parameter describing the presence of
!		normal under story vegetation (NOT
!		MOSS) under the canopy			-
! ivgtyp	Vegetation ype of the current pixel:
!		0) Bare soil
!		1) Vegetation with upper layer roots
!		2) Vegetation with lower layer roots	-
! iwel		1D array of pixel weights for temp.
!		 lapse rate				degree/m
! iyear		Year of the timestep being solved	-
! ixpix		Column in the input image for each
!		pixel in the catchment			-
! iypix		Row in the input image for each
!		pixel in the catchment			-
! jday          Julian day timestep being solved	-
! lakpix	Pixel number of the lake		-
! lat_deg	2 digit central latitude for site	degrees
! lat_min	2 digit central latitude for site	minutes
! loop		Looping variable			-
! lng_deg	2 digit central longitude for site	degrees
! lng_min	2 digit central longitude for site	minutes
! lng_mer	2 digit standard meridian for region	degrees
! lw_max	Maximum acceptable value for longwave
!		radiation				W/m2
! lw_min	Minimum acceptable value for longwave
!		radiation				W/m2
! maxnri	Maximum number of iterations in the
!		solution of the energy balance
!		equations				-
! minpdeg	Minutes per degree of rotation		minutes
! minphour	Minutes per hour (60)			-
! mixmax	Maximal depth of convective mixing	node number
! MODE		1 if the model is run in GIS mode, 2
!		if in statistical mode			-
! ncatch	Number of catchments			- 
! nlakpix	Number of lake pixels in the catchment	-
! ngb		List of the stations to be used for
!		each pixel in the calculation of
!		ground heat flux (in case of solution
!		using the Penman-Moneith equations)	-
! nhu		List of the stations to be used for
!		each pixel in the calculation of
!		relative humidity			-
! nlw		List of the stations to be used for
!		each pixel in the calculation of
!		incoming long wave radiation
! npa		List of the stations to be used for
!		each pixel in the calculation of
!		air pressure				-
! npet		List of the stations to be used for
!		each pixel in the calculation of
!		potential evapotranspiration (in
!		case the evapotranspiration values
!		are externally input)			-
! nppt		List of the stations to be used for
!		each pixel in the calculation of
!		precipitation				-
! nrn		List of the stations to be used for
!		each pixel in the calculation of
!		net radiation (in case of using the
!		Penman-Monteith equations)		-
! nsta_gb	Number of stations measuring ground
!		heat flux				-
! nsta_hu	Number of stations measuring air
!		humidity				-
! nsta_lw	Number of stations measuring incoming
!		long wave radiation			-
! nsta_pa	Number of stations measuring air
!		pressure				-
! nsta_pet	Number of stations measuring potential
!		evapotranspiration			-
! nsta_ppt	Number of stations measuring
!		precipitation				-
! nsta_rn	Number of stations measuring net
!		radiation				-
! nsta_sw	Number of stations measuring
!		incoming short wave radiation		-
! nsta_ta	Number of stations measuring air
!		temperature				-
! nsta_ws	Number of stations measuring wind speed	-
! nsw		List of the stations to be used for
!		each pixel in the calculation of
!		incoming short wave radiation		-
! nta		List of the stations to be used for
!		each pixel in the calculation of
!		air temperature				-
! numnod	Number of nodes in the lake		-
! nws		List of the stations to be used for
!		each pixel in the calculation of
!		wind speed				-
! nvegpix	Number of vegetated or bare soil pixels
!		in the catchment			-
! Outflow	Liquid water outflow out of the
!		snow pack				m
! Outflow_us	Liquid water outflow out of the
!		snow pack on top of the under story	m
! PI		Circumference divided by diameter	radians
! Packwater	Liquid water content of snow pack	m
! Packwater_us	Liquid water content of snow pack on
!		top of the under story			m
! pa_max	Maximum acceptable value for air
!		pressure				mbar
! pa_min	Minimum acceptable value for air
!		pressure				mbar
! par		Constant in the bare soil
!		desorptivity equation			-
! perixr	Percent of land surface contributing
!		to saturation excess runoff		-
! perrg1	Percent of land in region 1
!		(zw-psi! > zrzmax)			-
! perrg2	Percent of land surface in region 2
!		(zw-psi! < zrzmax, zw-psi! > 0)		-
! persa!	Percent of bare soil land surface that
!		is saturated and under atmosphere
!		controlled evaporation			-
! persxr	Percent of land surface contributing
!		saturation excess runoff		-
! perua!	Percent of bare soil land surface that
!		is unsaturated and under atmosphere
!		controlled evaporation			-
! perus!	Percent of bare soil land surface that
!		is unsaturated and under soil
!		controlled evaporation			-
! pet_max	Maximum acceptable value for potential
!		evapotranspiration if externally input	m/s
! pet_min	Minimum acceptable value for potential
!		evapotranspiration if externally input	m/s
! phase		Period of the deep soil temperature
!		cosine wave				1/timesteps
! pixsiz	Pixel resolution			m
! pnet		Net precipitation reaching the soil	m/s
! pnetsum	Average net precipitation for the
!		catchment				m/s
! pnetsumrg	Avrage net precipitation for the
!		region					m/s
! ppt_max	Maximum acceptable value for
!		precipitation				m/s
! ppt_min	Minimum acceptable value for
!		precipitation				m/s
! pptms		Precipitation coming in from atmosphere	m/s
! pptsum	Average preciptiation for the catchment	m/s
! pptsumrg	Average preciptiation for the region	m/s
! pr1rzs	Percent of total land surface in
!		region 1 with saturated root zone	-
! pr1sat	Percent of total land surface in
!		region 1 with saturated root zone
!		and transmission zone			-
! pr1tzs	Percent of total land surface in
!		region 1 with saturated
!		transmission zone			-
! pr1uns	Percent of total land surface in
!		region 1 which is unsaturated		-
! pr2sat	Percent of total land surface in
!		region 2 with saturated root zone	-
! pr2uns	Percent of total land surface in
!		region 2 which is unsaturated		-
! pr3sat	Percent of land surface in region 3
!		(zw-psi! <= 0)				-
! preca_a	Precipitation rate			m/s
! precac!	Accumulated precipitation in the lake	m/s
! precip_o	Precipitation that reaches the over
!		story and does not fall through to
!		the under story				m/s
! precip_u	Precipitation that reaches the under
!		story and that does not reach
!		the soil beneath the under story	m/s
! press		Air pressure				mbar
! psi!		Bubbling pressure			m
! psicav	Average bubbling pressure for each
!		catchment				m
! psicri	Critical leaf water potential		m
! psicri_us	Critical leaf water potential for
!		under story				m
! pstar		Psychrometri! constant adjusted for
!		soil resistance				Pa/K
! psurfx	Air pressure				Pa
! psychr	Psychrometri! constant			Pa/K
! psychr_i!	Psychrometri! constant in the canopy	Pa/K
! q0		Subsurface flow at complete saturation	m3/s
! qax		Specifi! humidity			-
! qb0		Initial base flow			m3/s
! qb		Baseflow				m3/s
! qbreg		Regional total baseflow			m3/s
! qsurf		Total surface runoff for the catchment	m/s
! qsurfrg	Total surface runoff for the region	m/s
! quartz	Quartz content of the soil		-
! qv		Specifi! humidity			-
! qv_i!		Specifi! humidity in the canopy		-
! qzbar		Excess net recharge to water table
!		over available soil water storage
!		allocated to runoff			m
! ra		Dry air gas constant			J/kgK
! radpdeg	Radians per degree (PI/180)		radians/degree
! radphour	Radians per hour (radpdeg * degphour)	radians/hour
! rahd		Aerodynami! resistance to heat flow for
!		dry surface				s/m
! rahd_us	Aerodynami! resistance to heat flow for
!		dry surface for under story		s/m
! rahd_us	Aerodynami! resistance to heat flow for
!		dry surface for the under story		s/m
! rah_moss	Aerodynami! resistance to heat flow for
!		the moss layer				s/m
! rahw		Aerodynami! resistance to heat flow for
!		wet surface				s/m
! rahw_us	Aerodynami! resistance to heat flow for
!		wet surface for the under story		s/m
! ra_i!		Dry air gas constant in the canopy	J/kgK
! rain		Liquid precipitation			m/s
! ranrun	Rate of runoff due to rainfall for
!		each catchment				m/s
! ranrunrg	Rate of runoff due to rainfall for
!		the region				m/s
! ravd		Aerodynami! resistance to vapor flow
!		for dry surface				s/m
! ravd_us	Aerodynami! resistance to vapor flow
!		for dry surface for under story		s/m
! raveff	Resulting total resistance for bare
!		soil (rsoil+ravd)			s/m
! rav_moss	Aerodynami! resistance to vapor flow
!		for the moss layer			s/m
! ravw		Aerodynami! resistance to vapor flow
!		for wet surface				s/m
! ravw_us	Aerodynami! resistance to vapor flow
!		for wet surface for under story		s/m
! RaSnow	Aerodynami! resistance to both heat
!		and vapor flow for the snow layer	s/m
! rcosncidenangl Cosine of the incidence angle between
!		 solar rays and the normal to the 
! 		 surface				-
! rcosihlfdaylng Cosine of the half day length		-
! r_diff	Difference between original skin
!		temperature and skin temperature
!		calculated by solution of the
!		energy balance				K
! reflect	Incoming reflected soloar radiation	W/m2
! reflectsky	View factor for variable reflect	
!		 Range between 0 and 1			-
! respla	Plant resistance			s/m
! respla_us	Plant resistance for under story	s/m
! rh		Relative humidity			%
! rh_i!		Relative humidity in canopy		%
! rha		Relative humidity for lake model	-
! rhostp	Water density 				kg/m3
! relrze	Relative saturation for desorptivity
!		calculation				-
! relsrz	Relative saturation in the root zone	-
! relstz	Relative saturation in the
!		transmission zone			-
! rescan	Canopyn resistance to transpiration	s/m
! rescan_us	Canopyn resistance to transpiration
!		for under story				s/m
! rib		Richardson number			-
! rib_moss	Richardson number for the moss layer	-
! rib_us	Richardson number for the under story	-
! r_lakearea	Area of the lakes in the current
!		catchment				m2
! rlatitude	Central latitude for basin		radians
! rld		Incoming long wave radiation		W/m2
! r_ldn		Incoming long wave radiation used in
!		the solution of the 2-layer energy
!		balance equations			W/m2
! rlng_merid	Calculated standard meridian for 
!		 central time zone			-
! rlongitude	Central longitude for the basin		radians
! rlongitudadjst Adjustment for longitude 		minutes
! rlwdx		Downwelling long wave radiation		W/m2
! rlwdn		Downwelling long wave radiation from
!		over story to under story		W/m2
! rlwup		Upwelling long wave radiation from
!		under story to over story		W/m2
! r_MeltEnergy	Energy used for melting and heating of
!		snow pack				W/m2
! r_MeltEnergy_us	Energy used for melting and
!		heating of snow pack on top of under
!		story					W/m2
! r_mindiff	Difference between input over story
!		skin temperature and the new
!		calculated skin temperature for the
!		optimal temperature			K
! r_moss_depth	Depth of the moss layer			m
! r_mossm	Moisture content in the moss layer	-
! r_mossm1	New calculated moisture content
!		in the moss layer			-
! r_mossm_u	Unfrozen water content in the moss
!		layer					-
! r_mossm1_u	New calculated unfrozen water
!		content in the moss layer		-
! r_mossm_f	Frozen water content in the moss
!		layer					-
! r_mossm1_f	New calculated rozen water content
!		in the moss layer			-
! rnact		Net radiation at actual
!		evaporation rate			W/m2
! rnactd	Net ratiation at actual
!		evaporation rate for dry canopy		W/m2
! rnactd_moss	Net ratiation at actual
!		evaporation rate for the moss		W/m2
! rnact_moss	Net ratiation at actual
!		evaporation rate for the moss		W/m2
! rnact_us	Net radiation at actual
!		evaporation rate for the under story	W/m2
! rnactd_us	Net ratiation at actual
!		evaporation rate for dry canopy for
!		the under story				W/m2
! rnetd		Net radiation at potential
!		evapotranspiration for dry vegetation	W/m2
! rnetd_us	Net radiation at potential
!		evapotranspiration for dry vegetation
!		for the under story			W/m2
! rnet_pot_moss	Net radiation at potential
!		evapotranspiration for the moss layer	W/m2
! rnetpn	Net radiation input into the Penman-
!		Monteith equation			W/m2
! rnetw		Net radiation at potential
!		evapotranspiration for wet vegetation	W/m2
! rnetw_us	Net radiation at potential
!		evapotranspiration for wet vegetation
!		for the under story			W/m2
! rn_max	Maximum acceptable value for net
!		radiation when the Penman-Monteith
!		method is used				W/m2
! rn_min	Minimum acceptable value for net
!		radiation when the Penman-Monteith
!		method is used				W/m2
! rnoonhour	True solar noon				hour
! rnpet		Net radiation at potential
!		evapotranspiration			W/m2
! rnpet_moss	Net radiation at potential
!		evapotranspiration for the moss layer	W/m2
! rnpetsum	Regional average net radiation at
!		potential rate				W/m2
! rnpet_us	Net radiation at potential
!		evapotranspiration for the under story	W/m2
! rn_snow	Net radiation of a snow pack		W/m2
! rn_snow_us	Net radiation of a snow pack on top
! 		of the under story			W/m2
! rnsum		Total net radiation for all surfaces
!		for a time step for the region		W/m2
! roa		Air density				kg/m3
! rn_snow	Net radiation of a snow pack		W/m2
! roa_i!	Air density in the canopy		kg/m3
! rocpsoil	Heat capacity of the soil (without
!		water and air)				J/kgK
! roi		Ice density				kg/m3
! row		Water density				kg/m3
! Rpl		Parameter used in the calculation of
!		f1par					W/m2
! Rpl_us	Parameter used in the calculation of
!		f1par_us				W/m2
! rsd		Incoming solar radiation		W/m2
! r_sdn		Incoming solar radiation used if the
!		solution of the 2-layer energy balance
!		equations				W/m2
! rsinsolaralt	Sine of sun's solar altitude		radians
! rsmax		Parameter used in the calculation of
!		f1par					s/m
! rsmax_us	Parameter used in the calculation of
!		f1par_us				s/m
! rsmin		Parameter used in the calculation of
!		f1par					s/m
! rsmin_us	Parameter used in the calculation of
!		f1par_us				s/m
! rsoil		Soil resistance to evaporation		s/m
! rs_over	Incoming solar radiation for the over
!		story					W/m2
! rs_under	Incoming solar radiation for the under
!		story					W/m2
! rtact		Root activity factor			-
! rtact_us	Root activity factor for under story	-
! rtdens	Length of roots per unit volume soil	1/m2
! rtdens_us	Length of roots per unit volume soil
!		for under story				1/m2
! rtres		Root resistivity			s/m
! rtres_us	Root resistivity for under story	s/m
! runtot	Runoff rate for each pixel		m/s
! rzdthetaudtemp	Change of unfrozen water
!		content in the soil versus temperature
!		in the root zone			1/K
! rzdthetaidt	Change of ice content over time in
!		the root zone				1/s
! rzrhs		Sum of fluxes into the root zone	m
! rzrhssum	Sum of fluxes into the root zone for
!		the region				m
! rzpore	Available pore space in root zone	m
! rzpsum	Sum of available porosity in root
!		zone just above water table for
!		each catchment				-
! rzsm		Root zone soil moisture			-
! rzsmav	Regional average root zone soil
!		moisture				-
! rzsm1		New calculated root zone soil moisture	-
! rzsmis	Root zone soil moisture at the
!		beginning of an interstorm event	-
! rzsm1_f	New calculatd ice content of the root
!		zone					-
! rzsmold	Root zone soil moisture at
!		previous time step			-
! rzsm_u	Liquid water content of the root zone	-
! rzsm1_u	New calculated liquid water content
!		of the root zone			-
! rzsmst	Value of root zone soil moisture
!		at the beginning of a storm event	-
! rzsm_tol	Tolerance in root zone soil moisture
!		in the iteration scheme			-
! rzwat		Available water in root zone		m
! satxr		Saturation excess runoff		m/s
! sesq		Desorptivity squared			-
! shift		Phase shift of the deep soil
!		temperature cosine wave			timesteps
! smbeg		Soil moisture at the beginning of
!		the current event			-
! smold		Soil moisture of previous timestep	-
! smcond	Soil moisture conductance for
!		vegetation				-
! smcond_us	Soil moisture conductance for
!		vegetation for under story		-
! smpet0	Initial soil moisture used to estimate
!		thermal parameters			-
! smtmp		Value of soil moisture used to
!		calculate thermal parameters in the
!		transmission zone			-
! snow		Snow fall				m/s
! solaraltitude	Solar altitude of sun from horizon	radians
! solarazimuth	Solar azimuth angle of sun from east	radians
! solartimestp	Fraction of timestep sun is above hrzn	hours
! solarzenith	Solar zenith angle of sun from east	radians
! sorp		Sorptivity
! srespar1	First parameter in the calculation of
!		the soil resitance to evaporation	depends on method
! srespar1_moss	First parameter in the calculation of
!		the soil resitance to evaporation,
!		for the moss layer			depends on method
! srespar2	Second parameter in the calculation of
!		the soil resitance to evaporation	depends on method
! srespar2_moss	Second parameter in the calculation of
!		the soil resitance to evaporation,
!		for the moss layer			depends on method
! srespar3	Third parameter in the calculation of
!		the soil resitance to evaporation	depends on method
! srespar3_moss	Third parameter in the calculation of
!		the soil resitance to evaporation,
!		for the moss layer			depends on method
! srzflx	Sum of the fluxes in and out of the
!		root zone				m
! starthour	Correct hour in solar time		hours
! sunearthdist	Distance from sun to earth		m
! sunrise	Time of sunrise				hours
! sunset	Time of sunset				hours
! sunmax	Calculated solar radiation at top
!		 of atmosphere. This is read in
!		 as the station measured value.		W/m2
! surface	Area of each lake node			m2
! SurfWater	Liquid water content of surface layer	m
! SurfWater_us	Liquid water content of surface layer
!		on top of under storu			m
! stzflx	Sum of the fluxes in and out of the
!		transmission zone			m
! svarhs	Sum of fluxes into the soil-vegetation-
!		atmosphere column			m/s
! svarhssm	Sum of fluxes into the soil-vegetation-
!		atmosphere column for the region	m/s
! sw_max	Maximum acceptable value for shortwave
!		radiation				W/m2
! sw_min	Minimum acceptable value for shortwave
!		radiation				W/m2
! Swq		Snow water equivalent			m
! Swq_us	Snow water equivalent on top of under
!		story		 			m
! swx		Incoming short wave radiation		W/m2
! sxrtot	Saturation excess runoff rate for the
!		catchment				m/s
! sxrtotrg	Saturation excess runoff rate for the
!		region					m/s
! ta_max	Maximum acceptable value for air
!		temperature				C
! ta_min	Minimum acceptable value for air
!		temperature				C
! tax		Air temperature				K
! tcbeta	Rate of exponential decline with LAI
! t!		Critical value of soil moisture below
!		which vegetation falls below
!		unstressed transpiration		-
!		of thermal conductivity under
!		vegetation				-
! tc_us		Critical value of soil moisture below
!		which vegetation falls below
!		unstressed transpiration		-
!		of thermal conductivity under
!		vegetation for under story		-
! tcbeta_us	Rate of exponential decline with LAI
!		of thermal conductivity under
!		vegetation for under story		-
! Tcutoff	Temperature at which water freezes	K
! tdiff		Difference between input
!		and newly calculated over story
!		skin temperature, used in the solution
!		of the 2-layer energy balance		K
! tcel		Air temperature				C
! tcel_i!	Air temperature in the canopy		C
! tdeep		Average deep soil temperature for the
!		entire year				K
! Tdeepstep	Deep soil temperature for this time
!		step					K
! tdry		Air temperature				C
! tempi_a	Ice temperature of each lake node	C
! temp_a	Temperature of each lake node		C
! therm!	Soil thermal conductivity		W/mK
! thermc1	Soil thermal conductivity of
!		upper layer				W/mK
! thermc2	Soil thermal conductivity of
!		lower layer				W/mK
! thermc2_us	Soil thermal conductivity of
!		lower layer for the under story		W/mK
! thermc_moss	Thermal conductivity of the moss
!		layer					W/mK
! thermc_us	Thermal conductivity of the under
!		story					W/mK
! thetar	Residual moisture content		-
! thetas	Saturated moisture content		-
! thetas_moss	Saturated moisture content of the moss
!		layer					-
! tfin		Optimal over story skin temperature	K
! timeadjust	Required time adjustment to convert
!		 local time to solar time		-
! timestep	Timestep for simulation in hours
! Tincan	Air temperature in the canopy		C
! Tint1		Intercept of the first part of the
!		regression of temperature difference
!		above and in the canopy versus total
!		incoming radiation			C
! Tint2		Intercept of the second part of the
!		regression of temperature difference
!		above and in the canopy versus total
!		incoming radiation			C
! tkact		Skin temperature			K
! tkactd	Skin temperature at actual
!		evaporation rate for dry canopy		W/m2
! tkactd_us	Skin temperature at actual
!		evaporation rate for dry canopy for
!		the under story				W/m2
! tkactd_moss	Skin temperature at actual
!		rate for the moss layer			K
! tkact_moss	Skin temperature of the moss layer	K
! tkact_us	Skin temperature of the under story	K
! tkd		Skin temperature at potential
!		evapotranspiration for dry vetation	K
! tkdeepsum	Regional average deep soil temperature	K
! tkd_us	Skin temperature at potential
!		evapotranspiration for dry vetation for
!		the under story				K
! tkel		Air temperature				K
! tkel_i!	Air temperature in the canopy		K
! tkmid		Mid soil temperature of the soil	K
! tkmidactd	Mid soil temperature at actual
!		evaporation rate for dry canopy		K
! tkmidactd_us	Mid soil temperature at actual
!		evaporation rate for dry canopy for the
!		under story				K
! tkmidactd_moss	Mid soil temperature at actual
!		evaporation rate for the moss		K
! tkmidd	Mid soil temperature at potential
!		evapotranspiration for dry vetation	K
! tkmidd_us	Mid soil temperature at potential
!		evapotranspiration for dry vetation for
!		the under story				K
! tkmidpet	Skin temperature at potential
!		evapotranspiration			K
! tkmidpet_moss	Skin temperature at potential
!		evapotranspiration for the moss layer	K
! tkmidpetsum	Regional average mid soil temperature
!		at potential rate			K
! tkmidpet_us	Skin temperature at potential
!		evapotranspiration for the under story	K
! tkmid_p_moss	Skin temperature at potential
!		evapotranspiration for the moss layer	K
! tkmid_us	Mid soil temperature of the soil for
!		the under story				K
! tkmidsum	Regional average mid soil temperature	K
! tkmidw	Skin temperature at potential
!		evapotranspiration for wet vegetation	K
! tkmidw_us	Skin temperature at potential
!		evapotranspiration for wet vegetation
!		for the under story			K
! tkpet		Skin temperature at potential
!		evapotranspiration			K
! tkpetsum	Regional average skin temperature
!		at potential rate			K
! toleb		Tolerance for skin temperature in the
!		solution of the energy balance
!		equations				K
! tolinf	Small number to be subtracted from
!		rzsmst so that sorptivity expression
!		does not blow up
! tkpet_moss	Skin temperature at potential
!		evapotranspiration for the moss layer	K
! tkpet_us	Skin temperature at potential
!		evapotranspiration for the under story	K
! tk_p_moss	Skin temperature at potential
!		evapotranspiration for the moss layer	K
! tksum		Average skin temperature for all
!		surfaces for a time step for the
!		region					K
! tkw		Skin temperature at potential
!		evapotranspiration for wet vegetation	K
! tkw_us	Skin temperature at potential
!		evapotranspiration for wet vegetation
!		for the under story			K
! tmid0		Mid layer initial soil temperature	K
! tmid0_moss	Mid layer initial soil temperature
!		under a moss layer			K
! Tpack		Temperature of snow pack		C
! Tpack_us	Temperature of snow pack on top of
!		the under story				C
! tp_in		Initial lake temperature		C
! trefk		Parameter used in the calculation of
!		f4temp					1/K2
! trefk_us	Parameter used in the calculation of
!		f4temp_us				1/K2
! trlup		Temperature of the under story layer
!		used to calculate the incoming long
!		wave radiation for the over story
!		layer					K
! tskinact_moss	Skin temperature of the moss layer	K
! tskinactd_moss	Skin temperature of the moss
!		layer					K
! tskinpet_moss	Skin temperature of the moss layer at
!		potential evapotranspiration		K
! tskin_p_moss	Skin temperature of the moss layer at
!		potential evapotranspiration		K
! Tslope1	Slope of the first part of the
!		regression of temperature difference
!		above and in the canopy versus total
!		incoming radiation			Cm2/W
! Tslope2	Slope of the second part of the
!		regression of temperature difference
!		above and in the canopy versus total
!		incoming radiation			Cm2/W
! tsoilold	Soil temperature at beginning of time
!		step					K
! Tsurf		Temperature of snow pack surface layer	C
! Tsurf_us	Temperature of snow pack surface layer
!		on top of the under story		C
! tw		Value of soil moisture below which
!		wilting begins				-
! tw_us		Value of soil moisture below which
!		wilting begins for under story		-
! twet		Dew point temperature			C
! twet_i!	Dew point temperature in the canopy	C
! Twint1	Intercept of the first part of the
!		regression of dew temperature difference
!		above and in the canopy versus total
!		incoming radiation			C
! Twint2	Intercept of the second part of the
!		regression of dew temperature difference
!		above and in the canopy versus total
!		incoming radiation			C
! Twslope1	Slope of the first part of the
!		regression of dew temperature difference
!		above and in the canopy versus total
!		incoming radiation			Cm2/W
! Twslope2	Slope of the second part of the
!		regression of dew temperature difference
!		above and in the canopy versus total
!		incoming radiation			Cm2/W
! tzdthetaidt	Change of ice content over time in
!		the transmission zone			1/s
! tzpore	Available pore space in transmission
!		zone					m
! tzpsum	Sum of available porosity in
!		transmission zone just above water
!		table for each catchment		-
! tzrhs		Sum of fluxes into the transmission
!		zone					m
! tzrhssum	Sum of fluxes into the transmission
!		zone for the region			m
! tzsm		Transmission zone soil moisture		-
! tzsmav	Regional average transmission zone soil
!		moisture				-
! tzsm1		New calculated transmission zone soil
!		moisture				-
! tzsm1_f	New calculated ice content of
!		the transmission zone			-
! tzsmold	Transmissiob zone soil moisture at
!		previous time step			-
! tzsm1_u	New calculated liquid water content
!		of the transmission zone		-
! tzwat		Available moisture in transmission zone	m
! uzw		Wind speed at height zww		m/s
! VaporMassFlux	Mass flux of water vapor to or from
!		the intercepted snow			m
! VaporMassFlux_us	Mass flux of water vapor to or
!		from the intercepted snow on top of
!		the under story				m
! vegcap	Maximum flux of water vapor
!		sustainable by plants			m/s
! vegcap_us	Maximum flux of water vapor
!		sustainable by under story plants	m/s
! vpdef		Vapor pressure deficit			Pa
! vpdef_i!	Vapor pressure deficit inside the
!		canopy					Pa
! vppa		Vapor pressure at level za		Pa
! vppa_i!	Vapor pressure at level za inside the
!		canopy					Pa
! vpsat		Vapor pressure in the air		Pa
! vpsat_i!	Vapor pressure in the air in the canopy	Pa
! waspt		Weight of slope aspect 			degrees
!		 0/360 = North
!		 90 = East 
! w!		Canopy water storage			m
! wcsum		Areal average canopy water storage	m
! wc_us		Canopy water storage for the under
!		story					m
! wcip1		Canopy water storage from previous
!		time step				m
! wcip1sum	Areal average canopy water storage
!		from previous time step			m
! wcip1_us	Canopy water storage from previous
!		time step for the under story		m
! wcrhs		Sum of fluxes in and out of the
!		canopy					m/s
! wcrhssum	Sum of fluxes in and out of the
!		canopy for the region			m/s
! wcrhs_us	Sum of fluxes in and out of the
!		canopy for the under story		m/s
! wgb		Weight of each station in the
!		calculation of ground heat flux		-
! whu		Weight of each station in the
!		calculation of air humidity		-
! wlw		Weight of each station in the
!		calculation of long wave radiation	-
! wpa		Weight of each station in the
!		calculation of air pressure		-
! wpet		Weight of each station in the
!		calculation of potential
!		evapotranspiration			-
! wppt		Weight of each station in the
!		calculation of precipitation		-
! wrn		Weight of each station in the
!		calculation of net radiation		-
! ws!		Water storage capacity			m
! wsc_us	Water storage capacity for under story	m
! wslop		2D array of weights for slope		degrees
!		 0 = flat
!		 90 = right angle with surface
! ws_max	Maximum acceptable value for wind speed	m/s
! ws_min	Minimum acceptable value for wind speed	m/s
! wsw		Weight of each station in the
!		calculation of wind speed		-
! wta		Weight of each station in the
!		calculation of air temperature		-
! wws		Weight of each station in the
!		calculation of wind speed		-
! xinact	Actual infiltration rate		m/s
! xinfcp	Infiltration capacity			m/s
! xinfxr	Infiltration excess runoff rate		m/s
! xintst	Cumulative time after ppt ends: when
!		greater than threshold endstm reset
!		rzsm for infiltration calculations	s
! xintst_moss	Cumulative time after ppt ends: when
!		greater than threshold endstm reset
!		rzsm for infiltration calculations,
!		for moss layer				s
! xixtot	Infiltration excess runoff rate for the
!		catchment				m/s
! xixtotrg	Infiltration excess runoff rate for the
!		region					m/s
! xk0		Saturated hydrauli! conductivity at
!		the soil surface			m/s
! xkscap	Saturated hydrauli! conductivity at
!		the midpoint of the surface and the
!		water table depth			m/s
! xksrz		Saturated hydrauli! conductivity in
!		the root zone				m/s
! xkstz		Saturated hydrauli! conductivity in
!		the transmission zone			m/s
! xlai		Leaf area index				-
! xlai_us	Leaf area index of the under story	-
! xlai_ws!	Leaf area index for canopy
!		interception storage			-
! xlai_wsc_us	Leaf area index for canopy
!		interception storage for the under
!		story					-
! xlamda	Areal average topography-soils index	-
! xlat_a	Latitude of the lake			degrees
! xleact	Latent heat flux at actual
!		evaporation rate			W/m2
! xleactd	Actual latent heat flux from dry
!		canopy					W/m2
! xleact_us	Latent heat flux at actual
!		evaporation rate for the under story	W/m2
! xleact_snow	Latent heat flux out of a snow pack	W/m2
! xleact_snow_us	Latent heat flux out of a snow
!		pack on top of the under story		W/m2
! xleactd_us	Actual latent heat flux from dry
!		canopy for the under story		W/m2
! xleactd_moss	Actual latent heat flux from the
!		moss layer				W/m2
! xleact_moss	Actual latent heat flux out of the
!		moss layer				W/m2
! xle_act_moss	Actual latent heat flux out of the
!		moss layer				W/m2
! xle_act_moss	Actual latent heat flux out of the
!		moss layer				W/m2
! xled		Latent heat flux at potential
!		evapotranspiration for dry vetation	W/m2
! xled_us	Latent heat flux at potential
!		evapotranspiration for dry vetation for
!		the under story				W/m2
! xlength	Catchment total length of channel
!		network					m
! xlepet	Latent heat flux at potential
!		evapotranspiration			W/m2
! xlepet_moss	Latent heat flux at potential
!		evapotranspiration for the moss layer	W/m2
! xlepetsum	Regional average latent heat flux
!		at potential rate			W/m2
! xlepet_us	Latent heat flux at potential
!		evapotranspiration for the under story	W/m2
! xle_p_moss	Latent heat flux at potential
!		evapotranspiration for the moss layer	W/m2
! xlesum	Total latent heat flux for all surfaces
!		for a time step for the region		W/m2
! xlew		Latent heat flux at potential		
!		evapotranspiration for wet vegetation	W/m2
! xlew_us	Latent heat flux at potential
!		evapotranspiration for wet vegetation
!		for the under story			W/m2
! xlhv		Latent heat of vaporization		J/kg
! xlhv_i!	Latent heat of vaporization in the
!		canopy					J/kg
! xlong_a	Longitude of the lake			degrees
! z0hd		Roughness length for dry surface
!		for heat transfer			m
! z0h_us	Roughness length for heat transfer
!		for the under story			m
! z0hw		Roughness length for wet surface
!		for heat transfer			m
! z0m		Roughness length for momentum transfer	m
! z0m_moss	Roughness length for momentum transfer
!		for the moss layer			m
! z0m_us	Roughness length for momentum transfer
!		for the under story			m
! z0vd		Roughness length for dry surface
!		for momentum transfer			m
! z0vw		Roughness length for wet surface
!		for momentum transfer			m
! za		Height of measurements of
!		meteorological data but not wind speed	m
! zbar		Catchment average water table depth	m
! zbarrg	Regional average water table depth	m
! zbar0		Initial water table depth		m
! zbar1		Catchment average water table depth
!		at end of previous time step		m
! zbar1rg	Regional average water table depth at
!		end of previous time step		m
! zbrflx	Net recharge to the water table		m
! zbrpor	Available soil water storage for water
!		table recharge				m
! zdeep		Depth of the deep soil temperature	m
! zmid		Depth of the mid soil temperature	m
! zmoss		Water depth in the moss layer		m
! zpd		Zero plane dislacement height		m
! zpd_us	Zero plane dislacement height for under
!		story					m
! zrz		Depth of the root zone			m
! zrzmax	Maximum value of root zone depth	m
! ztz		Transmission zone length		m
! zw		Local water table depth (m)		m
! zww		Height of the wind speed measuremets	m
