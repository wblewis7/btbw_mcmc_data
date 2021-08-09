# btbw_mcmc_data
Data and JAGS code for running mcmc analyses for Lewis et al. Demographic drivers of trailing edge range contractions in a migratory songbird

---
authors: William B. Lewis, Robert J. Cooper, Richard B. Chandler, T. Scott Sillett, Ryan W. Chitwood, Kirk W. Stodola, Mason H. Cline, Joanna L. Hatt, Michael T. Hallworth, Sara A. Kaiser, and Nicholas L. Rodenhouse.

Data for publication: Demographic drivers of trailing edge range contractions in a migratory songbird. Submitted to Global Change Biology
---

# Metadata for btbw-mcmc and climate data

The data for the btbw-mcmc project are stored in the 'BTBW_markre_clim_data' gzip file. Climate and black-throated blue warbler
(Setophaga caerulescens, BTBW) mark-recapture data were collected from the trailing edge of the range near
the Coweeta LTER in North Carolina (CWT) and at the range core at the Hubbard Brook Experimental Forest 
in New Hampshire (HB).

Sample JAGS code for running the Bayesian hierarchical population models is contained in the 'BTBW_JAGS_code' text file. The
file contains code for temporal trends on per-capita and recruitment; for climate models (t-7)/4 would be substituted with
climate variable. All climate variables were standardized prior to running models. The model for the trailing-low study plot
was modified slightly so that η in 2017 and 2018 was multiplied by a measure of detection probability with an informative prior
to account for imperfect detection during surveys in those years. 

There are five data sources to describe, all contained within the gzip file:

- capData contains capture histories and individual data.
- countData contains counts of unmarked individuals on each plot in each year.
- plotData contains plot-level data.
- precipData contains daily precipitation (mm) data from weather stations at both sites.
- thermalsumsData contains daily degree days over 4C for calculating thermal sums at both sites.

These first 3 data files only contain information on birds that were within
the plot boundaries. Some birds were captured or first-detected off
the plots. This is indicated by discrepancies between
"cap_data > cap_year" and encounter histories.

The last 2 data files contain climate data from 2001-2018 for weather
stations operated by the USDA forest service. Stations are not located on the actual BTBW study plots,
but are in the nearby area. Temperature data can be found at https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0042
and https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hbr.59.9 while precipitation data can
be found at https://www.fs.usda.gov/rds/archive/catalog/RDS-2017-0031 and 
https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hbr.13.13.
Only stations which were located within 175m in elevation of the study plots were included.
Both temperature and precipitation are correlated with elevation, but this relationship does not differ
over time at either range postions. 





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

The mid elevation at each plot (m ASL).

### first_year

First year used in analysis

### last_year

Last year used in analysis. The low-elevation plot at CWT ('rk') was not
intensively sampled after 2008, but follow-up surveys were performed in 2017-2018.





## capData

### Detections

The first 18 columns are detections named yYEAR. A value of `1` indicates that
that individuals was detected in that year. A value of `0` indicates the
opposite. A value of `NA` indicates that the study plot was not sampled in that
year. The low-elevation plot at CWT ('rk') was not intensively sampled after 
2008, but follow-up surveys were performed in 2017-2018.

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

### alum_band

This is the unique Federal USGS aluminum band number for each
individual. Some values are left blank, indicating that the bird
was banded with plastic colored leg bands but not an aluminum
band. This column is not used for analysis.





## countData

The number of unbanded female BTBW breeding on each study plot in each
year. Calculated from PDF maps of spot-mapping territory data. 
Manually transcribed. Note the plot variable matches directly with both
other data sets.





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





## thermalsumsData

Daily degree days calculated based on daily minimum and maximum values from temperature data,
as well as site-specific hours of sunrise and sunset (rounded down) as in Cesaraccio 
et al. 2001,Lany et al. 2016. Method provides measure of heat accumulation (and hence plant growth
and phenology) by estimating hourly temps and then calulating growing degree days based on
daytime temps greater than 4C. Degree days summed together over specific time periods to 
calculate thermal sum. Important determinant of seasonal progression and caterpillars 
(Lany et al. 2016, Reynolds et al. 2007). 
Used to determine total thermal sums during first broods at both sites, as
temperature could cause thermal stress to eggs or nestlings. Using DOY 127-
159 at CWT (average day of first egg and fledge date of first broods, 
respectively) and DOY 143-175 at HB. Roughly corresponds to May temps,
which have been increasing at both sites.

### Site

Variable matches other datasets

### Year

Year of measurement

### DOY

Calendar day of the year of the measurement

### DegreeDays
Daily values of growing degree days over 4C calculated via Cesaraccio et al. 2001

### Elevation
Elevation (m ASL) of the weather station
