# btbw_mcmc_data
Data for running mcmc analyses for Lewis et al. Demographic drivers of trailing edge range contractions in a migratory songbird

---
authors: William B. Lewis, Robert J. Cooper, Richard B. Chandler, T. Scott Sillett, Ryan W. Chitwood, Kirk W. Stodola, Mason H. Cline, Joanna L. Hatt, Michael T. Hallworth, Sara A. Kaiser, and Nicholas L. Rodenhouse
Data for publication: Demographic drivers of trailing edge range contractions in a migratory songbird. Submitted to Global Change Biology
---

# Metadata for btbw-mcmc and climate data

This document describes all of the data for the btbw-mcmc project. Climate and black-throated blue warbler
(Setophaga caerulescens, BTBW) mark-recapture data were collected from the trailing edge of the range near
the Coweeta LTER in North Carolina (CWT) and at the range core at the Hubbard Brook Experimental Forest 
in New Hampshire (HB).

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
##############################################################################

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




