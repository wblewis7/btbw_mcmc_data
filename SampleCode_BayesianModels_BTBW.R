###############################################################################
## Sample code for analyzing long-term demographic trend of local populations
##   of black-throated blue warblers (Setophaga caerulescens) data using 
##   Bayesian hierarchical models.
## Data were collected at the Coweeta Hydrologic Laboratory in North Carolina (CWT)
##   and the Hubbard Brook Experimental Forest (HB)
##
## Lewis, W. B., R. J. Cooper, R. B. Chandler, R. W. Chitwood, 
##   M. H. Cline, M. T. Hallworth, J. L. Hatt, J. Hepinstall-Cymerman,
##   S. A. Kaiser, N. L. Rodenhouse, T. S. Sillett, K. W. Stodola, 
##   M. S. Webster, and R. T. Holmes. Climate-mediated population dynamics
##   of a migratory songbird differ between the trailing edge and range core.
##   Submitted to Ecological Monographs.
################################################################################


library(rjags)
library(parallel)

# Loading in data
load("BTBW_rawdata_Lewis_etal.gzip")
ls()

# capData contains capture histories
# countData contains counts of unmarked females
# plotData contains plot-level data
# tempData contains temperature data from climate stations
# precipData contains precipitation data from climate stations


yrs <- 2002:2019
nyears <- length(yrs)
ch.cols <- paste("y", yrs, sep="")

# Models are run separately at each study plot
# This example is for the mid-elevation plot at the trailing edge (bs), though
#   code is similar for other plots.

(firstYear <- plotData[plotData$plot=="bs","first_year"]-2001)
(lastYear <- plotData[plotData$plot=="bs","last_year"]-2001)




###############################################################################
## Formatting Data
################################################################################



## Capture Data
###############################################################################

str(capData)
# 1 represents not captured/observed
# 2 represents captured/observed as SY (first-time breeder)
# 3 represents captured/observed as ASY (at least second-time breeder)

ch.bs <- data.matrix(capData[capData$plot=="bs", ch.cols])
str(ch.bs)

(nind.bs <- nrow(ch.bs))

(yrs.bs <- yrs[!is.na(ch.bs[1,])]) ##ch.bs[1,]
yri.bs <- which(!is.na(ch.bs[1,]))
yri.bs

## Augment the capture histories
M.bs <- 600
ch.bs.aug <- matrix(1L, M.bs, nyears)
ch.bs.aug[1:nrow(ch.bs),] <- ch.bs


## Detections (resight data)
ydet.bs.aug <- ifelse(ch.bs.aug>1, 1L, 0L)
head(ydet.bs.aug, 20)


## Is the individual marked?
marked <- matrix(0L, M.bs, nyears)
## Year of first capture or nyears if not captured
firstcap <- rep(lastYear, M.bs)
for(i in 1:nrow(ch.bs)) {
  firstcap[i] <- min(which(ch.bs[i,]>1))
  marked[i,firstcap[i]:lastYear] <- 1L
}

head(marked, 20)





## Counts of unbanded birds
################################################################################

countData.bs <- countData[countData$plot=="bs",]
countData.bs

uF.bs <- rep(NA, nyears)
names(uF.bs) <- yrs
uF.bs[yrs.bs-2001] <- countData.bs$unbanded_females#[-15]
uF.bs





## Climate Data
################################################################################

# Need to predict climate variables to the elevation of the study plots
# Predicting for t-1 years

tempData <- tempData[tempData$YEAR < yrs[length(yrs)],]
precipData <- precipData[precipData$YEAR < yrs[length(yrs)],]

# Climate variables are normalized in the model

# Early-breeding temperatures. Incorporating uncertainty in mean temperature estimates.
EBcutoff <- matrix(c(127,159,143,175),byrow=T,ncol=2)

EBtemp <- aggregate(TAVG~YEAR*Site*Elevation,data=tempData[tempData$Site=="CWT"&tempData$DOY>=EBcutoff[1,1]&tempData$DOY<=EBcutoff[1,2],], FUN=mean)
EBtemp_sd <- aggregate(TAVG~YEAR*Site*Elevation,data=tempData[tempData$Site=="CWT"&tempData$DOY>=EBcutoff[1,1]&tempData$DOY<=EBcutoff[1,2],], FUN=sd)
EBtemp_lm <- lm(TAVG~as.factor(YEAR)*Elevation,data=EBtemp,weights=1/EBtemp_sd$TAVG)
EBtemp_data <- data.frame(Elevation=plotData$elev_val[plotData$plot=="bs"],YEAR=as.factor(min(EBtemp$YEAR):max(EBtemp$YEAR)))
EBtemp_pred <-predict(EBtemp_lm,newdata=EBtemp_data,se.fit=T)
EBtemp_data$fit <- EBtemp_pred$fit
EBtemp_data$se <- EBtemp_pred$se.fit
# Only need climate data from start of study 
eb.temp <- EBtemp_data$fit[(as.integer(as.character(EBtemp_data$YEAR))-2001) %in% firstYear:(lastYear-1)]
eb.temp.se <- EBtemp_data$se[(as.integer(as.character(EBtemp_data$YEAR))-2001) %in% firstYear:(lastYear-1)]
mean.eb.temp <- mean(eb.temp,na.rm=T)
sd.eb.temp <- sd(eb.temp,na.rm=T)


# Annual precipitation
APrecip <- aggregate(DailyPrecipmm~Gaige*YEAR*Elevation*Site,data=precipData[precipData$Site=="CWT",],FUN=sum)
APrecip_lm <- lm(DailyPrecipmm~as.factor(YEAR)*Elevation,data=APrecip)
APrecip_data <- data.frame(Elevation=plotData$elev_val[plotData$plot=="bs"],YEAR=as.factor(min(APrecip$YEAR):max(APrecip$YEAR)))
APrecip_pred <-predict(APrecip_lm,newdata=APrecip_data,se.fit=T)
APrecip_data$fit <- APrecip_pred$fit
APrecip_data$se <- APrecip_pred$se.fit
# Only need climate data from start of study
an.precip <- APrecip_data$fit[(as.integer(as.character(APrecip_data$YEAR))-2001) %in% firstYear:(lastYear-1)]
an.precip.se <- APrecip_data$se[(as.integer(as.character(APrecip_data$YEAR))-2001) %in% firstYear:(lastYear-1)]
mean.an.precip <- mean(an.precip,na.rm=T)
sd.an.precip <- sd(an.precip,na.rm=T)




jdat.bs <- list(ycap=ch.bs.aug, ydet=ydet.bs.aug,
                M=M.bs,
                firstYear=as.integer(firstYear),
                lastYear=as.integer(lastYear),
                marked=marked, firstcap=firstcap,
                u=uF.bs,
                ebtemp=eb.temp, ebtempse=eb.temp.se,ebtemp.mean=mean.eb.temp, ebtemp.sd=sd.eb.temp,
                precip=an.precip, precipse=an.precip.se, precip.mean=mean.an.precip, precip.sd=sd.an.precip,
                Area=plotData[plotData$plot=="bs","hectares"])

str(jdat.bs)
################################################################################




## Initial Values
################################################################################

# Initializing year bird enters population
# To start, initializing birds observed on the plot with the year of first capture
#   if first captured as an SY or the year before first capture if first captured
#   as an ASY. Sets at last year of study + 1 if bird was never captured.
bi.bs <- rep(lastYear+1, M.bs) 
for(i in 1:nrow(ch.bs)) {
  first.det.i <- min(which(ch.bs[i,]>1)) 
  first.age.i <- ch.bs[i,first.det.i]
  known.yrs.i <- first.det.i:ncol(ch.bs)
  nknown.yrs.i <- length(known.yrs.i)
  if(first.age.i==2) {
    bi.bs[i] <- first.det.i
  } else if(first.age.i==3) {
    bi.bs[i] <- firstYear 
    if(first.det.i>firstYear)
      bi.bs[i] <- first.det.i-1 
  } else if(first.age.i==4) {
    bi.bs[i] <- 1
    if(first.det.i>firstYear)
      bi.bs[i] <- first.det.i-1
  } else stop("problem")
}
# Staggering initial entry year of uncaptured individuals
bi.bs[(nind.bs+1):(nind.bs+15)] <- 1
bi.bs[(nind.bs+16):(nind.bs+30)] <- 2
bi.bs[(nind.bs+31):(nind.bs+40)] <- 3
bi.bs[(nind.bs+41):(nind.bs+50)] <- 18
set.seed(35609)
bi.bs[351:550] <- sample(firstYear:lastYear, 200, replace=TRUE)
str(bi.bs)


# Initializing ages of captured birds in first year of study. Setting age (SY=1, ASY=2) 
#   for birds captured in the first year of study. Setting age at 1 in first year of study
#   for birds captured as ASY in second year of study.
agei <- matrix(NA, M.bs, ncol(ch.bs))
for(i in 1:nrow(ch.bs)) {
  if(ch.bs[i,firstYear]==2) {
    agei[i,firstYear] <- 1
  } else if(ch.bs[i,firstYear]==3) {
    agei[i,firstYear] <- 2
  } else if(ch.bs[i,1]==1 & ch.bs[i,2]==3) {
    agei[i,firstYear] <- 1 # Putting in a bird as SY in 2002 if detected as ASY in 2003 since must have been undetected SY in 2002
  }
}
agei[,1]


# zi.bs is matrix of initial values for when birds are known to be alive (first-last capture) beyond the first year
# zi0.bs is matrix of initial values for years in which birds were known to be alive and available for re-sighting beyond the first year
zi.bs <- matrix(0L, M.bs, ncol(ch.bs))
zi0.bs <- matrix(0L, M.bs, ncol(ch.bs))
for(i in 1:M.bs) {
  if(bi.bs[i]>lastYear) {
    next 
  }
  last.det.i <- bi.bs[i] 
  if(i<=nind.bs) {
    last.det.i <- max(which(ch.bs[i,]>1))
  } 
  zi.bs[i,bi.bs[i]:last.det.i] <- 1L
  if(last.det.i>bi.bs[i]) {
    zi0.bs[i,(bi.bs[i]+1):last.det.i] <- 1L
  }
}
zi.bs[,1] <- zi0.bs[,1] <- NA
zi.bs[,nyears] <- zi0.bs[,nyears] <- NA




# Compiling
ji.bs <- function() {
  out <- list(b=bi.bs,
              age=agei[,1:lastYear],
              z0=zi0.bs[,1:lastYear],
              psi=runif(1, 0.2, 0.3),
              phi0sy=runif(1),
              phi0asy=runif(1),
              phi1=0,
              gamma00=runif(1, 0.9, 1),
              gamma0=log(runif(1, 0, 0.01)),
              gamma1=0, 
              gamma2=0.001,
              tau=0.5, 
              pcap0=runif(1, 0.5, 0.9),
              pcap1=0,
              pdet0=runif(1, 0.5, 0.9),
              pdet1=0,
              kappa=runif(1, 0.5, 1))
  out$".RNG.name" <- "base::Mersenne-Twister" 
  out$".RNG.seed" <- sample.int(1e9, 1)
  return(out)    
}
##################################################################################





## Model
################################################################################

## Parameters to monitor
jp.bs <- c("phi0sy", "phi0asy", "phi1", "phi", "phiSY", "phiASY",
           "gamma0", "gamma1",  "gamma", "gamma2","gamma00",
           "tau", "pcap0", "pcap1", "pdet0",  "kappa",
           "Nsuper", "deviance", "N", "R", "S", "Q",
           "lambda", "Elam", "ER", "clim.act","b","EN",paste("state[", 1:jdat.bs$M, ",", jdat.bs$lastYear,
                                                             "]", sep=""))



cl1 <- makeCluster(3)



clusterExport(cl1, c("ji.bs", "bi.bs", "agei", "zi0.bs", "lastYear",
                     "jp.bs", "jdat.bs"))

clusterSetRNGStream(cl=cl1, iseed=3499)


system.time({
  out.ebtemp.bs <- clusterEvalQ(cl1, {
    library(rjags)
    load.module("dic")
    jm1.bs <- jags.model("BTBW_jags_sample_script.jag",
                         data=jdat.bs, inits=ji.bs,
                         n.adapt=1000)
    jcDD.bs1 <- coda.samples(jm1.bs, c(jp.bs), n.iter=35000)
    return(as.mcmc(jcDD.bs1))
  })
}) 
attr(out.ebtemp.bs, "Area") <- jdat.bs$Area
save(out.ebtemp.bs, file="mcmc_out_BTBW_ebtemp_bs.gzip")