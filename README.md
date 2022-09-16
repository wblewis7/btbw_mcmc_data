# BTBW_metadata
Data, JAGS code, R code, and tables of parameter and realized estimates for Bayesian hierarchical models estimating population dynamics and demographic rates of black-throated blue warblers (Setophaga caerulescens) breeding across a range of elevations at the trailing edge of the range in North Carolina and the core of the range in New Hampshire.
---
authors: William B. Lewis, Robert J. Cooper, Richard B. Chandler, Ryan W. Chitwood, Mason H. Cline, Michael T. Hallworth, Joanna L. Hatt, Jeff Hepinstall-Cymerman, Sara A. Kaiser, Nicholas L. Rodenhouse, T. Scott Scillett, Kirk W. Stodola, Michael S. Webster, and Richard T. Holmes


---

# Metadata

# BTBW_rawdata_Lewis_etal.gzip


The data for the black-throated blue warbler (Setophaga caerulescens, BTBW) project are stored in the 'BTBW_rawdata' gzip file. Climate and BTBW
mark-recapture data were collected from study sites two range positions: at the trailing edge of the range near
the Coweeta LTER in North Carolina (CWT) and at the range core at the Hubbard Brook Experimental Forest 
in New Hampshire (HB).


There are five data sources to describe, all contained within the gzip file:

- capData contains capture histories and individual data.
- countData contains counts of unmarked individuals on each study plot in each year.
- plotData contains plot-level data.
- tempData contains daily minimum, maximum, and average temperature (C) data from climate stations at both range positions.
- precipData contains daily precipitation (mm) data from weather stations at both range positions.

These first 3 data files only contain information on birds that were within
the plot boundaries. Some birds were captured or first-detected off
the plots. This is indicated by discrepancies between
"cap_data > cap_year" and encounter histories.

The last 2 data files contain climate data from 2002-2019 for weather
stations operated by the USDA forest service. Stations were not generally located on the actual BTBW study plots,
but were located adjacent to plots in the nearby area. Temperature data can be found at https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0042
and https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hbr.59.9 while precipitation data can
be found at https://www.fs.usda.gov/rds/archive/catalog/RDS-2017-0031 and 
https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hbr.13.13.
Only stations at elevations which were within 175m of the elevation range of the study plots at each range position were included.
Both temperature and precipitation are correlated with elevation at both range positions, but this relationship does not differ
over time at either site. 





## plotData

## plot

Codes denoting study plot. Variable names match other datasets. Note that the 
mid-elevation plot at the trailing edge was expanded by 11 ha in 2006. 
The original 18ha section is denoted 'bs' and is used for all cross-plot 
comparisons. The expanded section is denoted 'bs_exp' and is used only 
for comparison with the original section.

### hectares

The area in hectares of each plot. Important for calculating population density.

### site

The range postion of each plot. "CWT" = Coweeta LTER, NC (trailing edge) and 
"HB" = Hubbard Brook Experimental Forest, NH (range core).

### elev

A three-level factor describing each plot's relative elevation within each range postion.

### elev_val

The middle elevation of each study plot (m ASL).

### first_year

First year used in analysis

### last_year

Last year used in analysis. The low-elevation plot at CWT ('rk') was not
intensively sampled after 2008, but follow-up surveys were performed in 2017-2018.





## capData

### Detections

The first 18 columns are detections named yYEAR. A value of `1` indicates that
that individuals was not detected in that year. A value of `2` indicates that
bird was detected as an SY. A value of '3' indicates that the bird was detected
as an ASY. See "age_cap" below for appreviations. A value of `NA` indicates that
the study plot was not sampled in that year. The low-elevation plot at CWT ('rk') 
was not intensively sampled after 2008, but follow-up surveys were performed in 
2017-2018.

### plot

A seven-level factor name for each plot. These correspond
directly to the further plot information in plotData.

### sex

Subset down to females only.

### age_cap

A three-level factor describing the age at first capture, where "SY" =
second-year (first-time breeder), "ASY" = after-second-year (at least 
second-time breeder), and "AHY" = after-hatch-year (age unknown). 

### cap_year

The year of capture. IMPORTANT NOTE: This will not always correspond directly
with the first detection. If an individual was first captured off plot, it is
not included in the detection history.

Using age_cap and cap_year it is possible to construct a matrix of individual
age over time.





## countData

The number of unbanded female BTBW breeding on each study plot in each
year. Calculated from PDF maps of spot-mapping territory data. 
Manually transcribed. Note the plot variable matches directly with both
other data sets.





## tempData

Used to determine average daily temperatures during the early-breeding period. This is the
period of maximal breeding activity, and is calculated as the mean of average daily temperatures
between the average first laying date (day of year 127 (trailing edge) and 143 (range core)) and
the average fledge date of first nesting attempts (day of year 159 (trailing edge) and 175 (range
core)).

### Station

Name of climate station

### YEAR/MONTH/DAY

Date of the temperature readings

### TMIN

Daily minimum temperature in degrees C recorded from the climate station. This column is not used
in analysis.

### TMAX

Daily maximum temperature in degrees C recorded from the climate station. This column is not used
in analysis.

### TAVG

Daily average temperature in degrees C recorded from the climate station.

### Lat/Long/Elevation

The latitude, longitude, and elevation (m ASL) of the rain gaige

### DOY

Calendar day of the year of the reading

### Site

Variable matches other datasets





## precipData

Used to determine total annual precipitation. Note that for the CWT measurements many days
do not have readings. These indicate days where precipitation was too low to detect (assumed 0).

### Gaige

Name of rain gaige

### YEAR/MONTH/DAY

Date of the precipitation reading

### DailyPrecipmm

Daily precipitation in mm from the rain gaige. Values at CWT were reported in inches and 
converted.

### Lat/Long/Elevation

The latitude, longitude, and elevation (m ASL) of the rain gaige

### DOY

Calendar day of the year of the reading

### Site

Variable matches other datasets






# BTBW_JAGS_sample_script.jag


Sample JAGS code for running the Bayesian hierarchical population models is contained in the 'BTBW_JAGS_sample_script.jag' file, which
can be opened in a text editor. This script is called in the file 'SampleCode_BayesianModels_BTBW.R'. The file contains sample code for
modeling trend effects of average daily early-breeding temperatures on per-capitura recruitment and apparent survival. For temporal models,
clim.act.stand[t-1] would be replaced with (t-7)/4 and the for loop under 'Incorporating uncertainty in climate measures' would be ommitted.
For annual precipitation models, the prior distribution for clim.act would be changed to clim.act[n] ~ dunif(1000,3500). The model for the 
trailing-low study plot was modified slightly so that η in 2017 and 2018 was multiplied by a measure of detection probability with an 
informative prior to account for imperfect detection during surveys in those years. This informative prior was based on the mean (0.61) 
and sd (0.27) of daily detection probabilities of females at the higher-elevation study plots at the trailing edge from 2011 - 2019. 
In the model, 'ebtemp', 'ebtempse', 'ebtemp.mean', 'ebtemp.sd', 'u', 'ycap', and 'ydet' are supplied as data.






# SampleCode_BayesianModels_BTBW.R


Sample R code for running the Bayesian hierarchical population models is contained in the 'SampleCode_BayesianModels_BTBW' R file. This
code calls 'BTBW_rawdata_Lewis_etal.gzip' and 'BTBW_JAGS_sample_script.jag'. Sample code is provided for running the JAGS model, which
estimates trend effects of average daily early-breeding temperature on per-capita recruitment and apparent survival. The R code is the 
same when calling models with effects of time or annual precipitation. Code is shown for running population models at the mid-elevation
plot at the trailing edge of the range in North Carolina, though code is similar for most other study plots. The exception was the 
low-elevation plot at the trailing edge, where the plot was sampled from 2002 - 2008 and then was revisited once each year in 2017 and 2018.
The capture data from line 59 at this plot was reformatted with NAs from 2009 - 2016 and with 1s from 2017 - 2018. Similarly, the count 
data from line 98 was reformatted with NAs from 2009 - 2016 and 0s from 2017 - 2018. MCMC chains are saved in the file 'mcmc_out_BTBW_ebtemp_bs.gzip',
which is called by the file 'SampleCode_forecast_BTBW.R'.






# SampleCode_forecast_BTBW.R


Sample R code for performing statistical forecasting to assess population viability through 2040, contained in the 'SampleCode_forecast_BTBW'
R file. This code calls 'BTBW_rawdata_Lewis_etal.gzip' the mcmc output from 'SampleCode_forecast_BTBW.R', saved in the 'mcmc_out_BTBW_ebtemp_bs.gzip'
file. The file shows sample code for forecasting dynamics in response to projected trends in average daily early-breeding temperature at 
the mid-elevation plot at the trailing edge of the range in North Carolina, but code is similar at other study plots or when forecasting
dynamics in response to temporal or precipitation trends. For climate models, we projected climate at the study plots in future years 
(last year of study at the study plot - 2040) based on the observed 2002 - 2019 trend in climate variables at each range position (trailing
edge or range core). We assumed that future changes in climate would be similar across elevations within each range position.






# ParameterEstimatesfromBayesianHierarchicalModel.pdf


Parameter estimates of population dynamics of black-throated blue warblers breeding at the trailing edge of the range in North Carolina and
core of the range in New Hampshire. Mean, SD, median, lower, and upper 95% credible intervals are provided for estimates of the bounding
parameter on recruitment, per-capita recrutiment intercept, temporal/climate trend effect on per-capita recruitment, density-dependent effect
on recruitment, first-time breeder (SY) apparent survival intercept, after-first-time breeder (ASY) apparent survival intercept, temporal/climate
trend effect on apparent survival, age ratio of ASYs to SYs, probability of classifying birds into a specific age class upon capture, capture 
probability intercept, temporal trend in capture probability, and probability of being detected while breeding on the study plot. Models were
run separatley at each study plot, allowing for trend effects of year, average daily early-breeding temperatures, or annual precipitation on
per-capita recruitment and apparent survival.





# RealizedEstimatesfromBayesianHierachicalModel.pdf


Yearly estimates of abundance and realized demographic rates of black-throated blue warblers breeding at the trailing edge of the range in
North Carolina and core of the range in New Hamphire. Median and 95% credible intervals are provided for yearly estimates of abundance,
population growth rate, per-capita recruitment, and apparent survival. Estimates at each plot incorporated temporal trends on per-capita
recruitment and apparent survival.
