###############################################################################
## Sample code for forecasting population dynamics and assessing sensitivity to
##   vital rates of black-throated blue warblers (Setophaga caerulescens).
## Data were collected at the Coweeta Hydrologic Laboratory in North Carolina (CWT)
##   and the Hubbard Brook Experimental Forest (HB).
## Using mcmc output from Bayesian hierarchical models run in SampleCode_BTBW
##
## Lewis, W. B., R. J. Cooper, R. B. Chandler, R. W. Chitwood, 
##   M. H. Cline, M. T. Hallworth, J. L. Hatt, J. Hepinstall-Cymerman,
##   S. A. Kaiser, N. L. Rodenhouse, T. S. Sillett, K. W. Stodola, 
##   M. S. Webster, and R. T. Holmes. Climate-mediated population
##   dynamics of a migratory songbird differ between the trailing
##   edge and range core. Submitted to Ecological Monographs.
################################################################################

library(coda)

## Years with data
yrs <- 2002:2019 




# Forecasting population dynamics through 2030. Forecasting with a temporal trend
#   in recruitment and apparent survival, or with a climate trend.
# For the climate trend, predicting temperatures in future years at each site 
#    (CWT,HB) based on the observed trend in temperature during the course of 
#    the study. Assuming that future temperature change at each site will be
#    similar across elevations.




####################################################################################
## Climate Trends
###################################################################################

load("BTBW_rawdata_Lewis_etal.gzip")

# Early-breeding temperature
tempData <- tempData[tempData$YEAR %in% yrs,]
EBT_CWT <- aggregate(TAVG~YEAR*Elevation,data=tempData[tempData$DOY>=127&tempData$DOY<=159&tempData$Site=="CWT",], FUN=mean)
EBT_CWT_sd <- aggregate(TAVG~YEAR*Elevation,data=tempData[tempData$DOY>=127&tempData$DOY<=159&tempData$Site=="CWT",], FUN=sd)
EBT_CWT_lm <- lm(TAVG~YEAR + Elevation, data=EBT_CWT,weights=1/EBT_CWT_sd$TAVG)
EBT_HB <- aggregate(TAVG~YEAR*Elevation,data=tempData[tempData$DOY>=143&tempData$DOY<=175&tempData$Site=="HB",], FUN=mean)
EBT_HB_sd <- aggregate(TAVG~YEAR*Elevation,data=tempData[tempData$DOY>=143&tempData$DOY<=175&tempData$Site=="HB",], FUN=sd)
EBT_HB_lm <- lm(TAVG~YEAR + Elevation, data=EBT_HB,weights=1/EBT_HB_sd$TAVG)

# Annual precipitation
precipData <- precipData[precipData$YEAR %in% yrs,]
precip_CWT <- aggregate(DailyPrecipmm~Gaige*YEAR*Elevation,data=precipData[precipData$Site=="CWT",],FUN=sum)
precip_CWT_lm <- lm(DailyPrecipmm~YEAR + Elevation, data=precip_CWT)
precip_HB <- aggregate(DailyPrecipmm~Gaige*YEAR*Elevation,data=precipData[precipData$Site=="HB",],FUN=sum)
precip_HB_lm <- lm(DailyPrecipmm~YEAR + Elevation, data=precip_HB)


# Putting estimates in a dataframe
clim_beta <- data.frame(Site=rep(c("CWT","HB"),times=2),
                        Variable=rep(c("ebtemp","precip"),each=2),
                        b1=c(coef(EBT_CWT_lm)[names(coef(EBT_CWT_lm))=="YEAR"],coef(EBT_HB_lm)[names(coef(EBT_HB_lm))=="YEAR"],
                             coef(precip_CWT_lm)[names(coef(precip_CWT_lm))=="YEAR"],coef(precip_HB_lm)[names(coef(precip_HB_lm))=="YEAR"]),
                        se=c(coefficients(summary(EBT_CWT_lm))[rownames(coefficients(summary(EBT_CWT_lm)))=="YEAR",2],coefficients(summary(EBT_HB_lm))[rownames(coefficients(summary(EBT_HB_lm)))=="YEAR",2],
                             coefficients(summary(precip_CWT_lm))[rownames(coefficients(summary(precip_CWT_lm)))=="YEAR",2],coefficients(summary(precip_HB_lm))[rownames(coefficients(summary(precip_HB_lm)))=="YEAR",2]))


# Also need to calculate mean and sd of data for standardization
EBT_CWT$sd <- EBT_CWT_sd$TAVG
EBTpred <- matrix(NA, ncol=3, nrow=length(yrs))
for (i in 1:length(yrs)){
  EBTlm <- lm(TAVG~Elevation,data=EBT_CWT[EBT_CWT$YEAR==yrs[i],],weights=1/sd)
  EBTpred[i,] <- predict(EBTlm,newdata=data.frame(Elevation=plotData$elev_val[plotData$site=="CWT"]))
}
CWTpredparams <- cbind(plotData[plotData$site=="CWT",colnames(plotData) %in% c("site","elev")],
                       Mean=colMeans(EBTpred),
                       SD=apply(EBTpred,2,sd))
EBT_HB$sd <- EBT_HB_sd$TAVG
EBTpred <- matrix(NA, ncol=3, nrow=length(yrs))
for (i in 1:length(yrs)){
  EBTlm <- lm(TAVG~Elevation,data=EBT_HB[EBT_HB$YEAR==yrs[i],],weights=1/sd)
  EBTpred[i,] <- predict(EBTlm,newdata=data.frame(Elevation=plotData$elev_val[plotData$site=="HB"]))
}
HBpredparams <- cbind(plotData[plotData$site=="HB",colnames(plotData) %in% c("site","elev")],
                      Mean=colMeans(EBTpred),
                      SD=apply(EBTpred,2,sd))
EBTpredparams <- rbind(CWTpredparams,HBpredparams)
EBTpredparams$Variable <- "ebtemp"


precippred <- matrix(NA, ncol=3, nrow=length(yrs))
for (i in 1:length(yrs)){
  preciplm <- lm(DailyPrecipmm~Elevation,data=precip_CWT[precip_CWT$YEAR==yrs[i],])
  precippred[i,] <- predict(preciplm,newdata=data.frame(Elevation=plotData$elev_val[plotData$site=="CWT"]))
}
CWTpredparams <- cbind(plotData[plotData$site=="CWT",colnames(plotData) %in% c("site","elev")],
                       Mean=colMeans(precippred),
                       SD=apply(precippred,2,sd))
precippred <- matrix(NA, ncol=3, nrow=length(yrs))
for (i in 1:length(yrs)){
  preciplm <- lm(DailyPrecipmm~Elevation,data=precip_HB[precip_HB$YEAR==yrs[i],])
  precippred[i,] <- predict(preciplm,newdata=data.frame(Elevation=plotData$elev_val[plotData$site=="HB"]))
}
HBpredparams <- cbind(plotData[plotData$site=="HB",colnames(plotData) %in% c("site","elev")],
                      Mean=colMeans(precippred),
                      SD=apply(precippred,2,sd))
precippredparams <- rbind(CWTpredparams,HBpredparams)
precippredparams$Variable <- "precip"

predparams <- rbind(EBTpredparams,precippredparams)






################################################################################
## Forecasting Function
################################################################################

# Forecasting is performed for each iteration of the mcmc output. 
# For climate predictions, predicting climate variable in future years based on
#    the estimate of the climate variable in the lastyear - 1 

forecast <- function(mclist, nFuture=10, M.new, Area, thin=1, 
                     sensitivity=c(cgamma=1, cphiSY=1, cphiASY=1),
                     method="",
                     mean.inc=0, sd.inc=0,
                     mean.clim=0, sd.clim=0,
                     returnState=TRUE, report=0) {
  if(class(mclist) != "mcmc.list")
    stop("mclist should be a 'mcmc.list'")
  if(any(sensitivity<0))
    stop("sensitivity factors should be >0")
  cgamma <- sensitivity["cgamma"]
  cphiSY <- sensitivity["cphiSY"]
  cphiASY <- sensitivity["cphiASY"]
  vn <- varnames(mclist)
  state.names <- grep("state\\[", vn, value=TRUE)
  b.names <- grep("b\\[", vn, value=TRUE)
  if(length(state.names)<1 || length(b.names)<1)
    stop("'mclist' should have posterior samples for 'state' and 'b'")
  N.names <- grep("^N\\[", vn, value=TRUE)
  yrs <- as.integer(gsub("[^0-9]", "", N.names))
  first.yr <- min(yrs)
  last.yr <- max(yrs)
  yrs.all <- 1:(last.yr+nFuture)
  n.yrs.past <- length(N.names)
  n.yrs <- length(yrs.all) ##n.yrs.past+nFuture
  mcl <- window(mclist, thin=thin)
  n.chains <- nchain(mcl)
  n.iters <- niter(mcl)
  n.samples <- n.chains*n.iters
  M.old <- length(state.names)
  
  state.names <- paste("state[", 1:M.old, ",", last.yr, "]", sep="")
  b.names <- paste("b[", 1:M.old, "]", sep="")
  state <- array(3L, ## Initialize in dead/unborn state
                 c(M.new, n.yrs, n.samples))
  state[1:M.old,last.yr,] <- t(as.matrix(mcl[,state.names]))
  b <- matrix(NA_integer_, M.new, n.samples)
  b[1:M.old,] <- t(as.matrix(mcl[,b.names]))
  EN <- array(NA_real_, c(2, n.yrs))
  ED <- array(NA_real_, n.yrs)
  mcm <- as.matrix(mcl)
  
  gamma00 <- mcm[,"gamma00"]
  gamma0 <- mcm[,"gamma0"]
  gamma1 <- mcm[,"gamma1"]
  gamma2 <- mcm[,"gamma2"] 
  phi0sy <- mcm[,"phi0sy"]
  phi0asy <- mcm[,"phi0asy"]
  phi1 <- mcm[,"phi1"]
  
  # Climate data in last year - 1
  if(method=="clim"){
    clim.last <- mcm[,paste("clim.act[",last.yr-1,"]",sep="")]
  }
  
  clim.data <- matrix(NA,nrow=n.samples,ncol=n.yrs)
  
  # Simulating dynamics
  ER <- phiSY <- phiASY <- gamma <- numeric(n.yrs)
  pi.b <- numeric(n.yrs+1)
  W <- array(NA_real_, c(2, 2, n.yrs))
  A <- M.new-M.old
  bad.samples <- rep(0L, n.samples)
  reportit <- report>0
  for(s in 1:n.samples) {
    
    # Predicting future climate for each iteration
    if(method=="clim"){
      clim.data[s,last.yr-1] <- clim.last[s]
      for (t in (last.yr):n.yrs){
        clim.data[s,t] <- clim.data[s,t-1] + rnorm(1,mean.inc,sd.inc)
      }
      clim <- (clim.data - mean.clim)/sd.clim
    }
    
    
    if(reportit) {
      if(s %% report == 0)
        cat("Doing sample", s, "of", n.samples, "\t", format(Sys.time()), "\n")
    }
    EN[1,last.yr-1] <- mcm[s, paste("EN[1,", last.yr-1, "]", sep="")]
    EN[2,last.yr-1] <- mcm[s, paste("EN[2,", last.yr-1, "]", sep="")]
    ED[last.yr-1] <- sum(EN[1:2,last.yr-1])/Area
    if(method=="time"){
      gamma[last.yr-1] <- gamma00[s] /
        (1+exp(-(gamma0[s] + gamma1[s]*(last.yr-7)/4 -
                   gamma2[s]*(ED[last.yr-1]^2-0.3)/0.2))) 
      phiSY[last.yr-1] <- plogis(phi0sy[s] + phi1[s]*(last.yr-7)/4)
      phiASY[last.yr-1] <- plogis(phi0asy[s] + phi1[s]*(last.yr-7)/4)
    }
    if(method=="clim"){
      gamma[last.yr-1] <- gamma00[s] /
        (1+exp(-(gamma0[s] + gamma1[s]*clim[s,last.yr-1] -
                   gamma2[s]*(ED[last.yr-1]^2-0.3)/0.2)))
      phiSY[last.yr-1] <- plogis(phi0sy[s] + phi1[s]*clim[s,last.yr-1])
      phiASY[last.yr-1] <- plogis(phi0asy[s] + phi1[s]*clim[s,last.yr-1])
    }
    if(gamma[last.yr-1]>1) {
      warning("gamma's high")
    }
    W[1,,last.yr] <- gamma[last.yr-1]
    W[2,,last.yr] <- c(phiSY[last.yr-1],phiASY[last.yr-1])
    EN[1,last.yr] <- mcm[s, paste("EN[1,", last.yr, "]", sep="")]
    EN[2,last.yr] <- mcm[s, paste("EN[2,", last.yr, "]", sep="")]
    ED[last.yr] <- sum(EN[1:2,last.yr])/Area
    
    for(t in (last.yr+1):n.yrs) {
      if(method=="time"){
        gamma[t-1] <- gamma00[s] /
          (1+exp(-(gamma0[s] + gamma1[s]*(t-7)/4 -
                     gamma2[s]*(ED[t-1]^2-0.3)/0.2)))
        phiSY[t-1] <- plogis(phi0sy[s] + phi1[s]*(t-7)/4)
        phiASY[t-1] <- plogis(phi0asy[s] + phi1[s]*(t-7)/4)
      }
      if(method=="clim"){
        gamma[t-1] <- gamma00[s] /
          (1+exp(-(gamma0[s] + gamma1[s]*clim[s,t-1] -
                     gamma2[s]*(ED[t-1]^2-0.3)/0.2)))
        phiSY[t-1] <- plogis(phi0sy[s] + phi1[s]*clim[s,t-1])
        phiASY[t-1] <- plogis(phi0asy[s] + phi1[s]*clim[s,t-1])
      }
      
      if(gamma[t-1]>1) {
        warning("gamma's high")
        ## browser()
        ## gamma[t-1] <- 1
      }
      
      W[1,,t] <- gamma[t-1]
      W[2,,t] <- c(phiSY[t-1],phiASY[t-1])
      EN[1:2,t] <- W[,,t] %*% EN[,t-1]  ## Expected number of SYs and ASYs
      ED[t] <- sum(EN[1:2,t])/Area
      ER[t] <- EN[1,t]
      pi.b[t] <- ER[t]/A
    }
    
    ENsuper <- sum(ER) ## number of recruits during future years
    if(is.infinite(ENsuper)) {
      warning("Infinite superpop")
      bad.samples[s] <- 1
      next
    }
    omega <- ENsuper/A
    if(omega>1) {
      warning("Increase M.new")
      omega <- 1
      bad.samples[s] <- 2
      next
    }
    pi.b[n.yrs+1] <- 1-omega
    z <- ifelse(state[,,s]==3L, 0L, 1L)
    age <- matrix(NA_integer_, M.new, n.yrs)
    for(i in 1:M.new) {
      if(i <= M.old)
        if(b[i,s]==(last.yr+1))
          b[i,s] <- sample(n.yrs+1, 1, prob=pi.b)
      if(i > M.old)
        b[i,s] <- sample(n.yrs+1, 1, prob=pi.b)
      if(b[i,s]>n.yrs)
        next
      age[i,last.yr] <- ifelse(b[i,s]<last.yr, 2L, 1L)
      for(t in (last.yr+1):n.yrs) {
        b.year <- b[i,s]==t
        age[i,t] <- ifelse(b.year, 1L, 2L)
        z[i,t] <- rbinom(1, 1, b.year +
                           (1-b.year)*z[i,t-1]*W[2,age[i,t-1],t-1])
        state[i,t,s] <- z[i,t]*age[i,t] + (1-z[i,t])*3
      }
    }
  }
  N <- apply(state<3, c(2, 3), sum)
  N[first.yr:last.yr,] <- t(as.matrix(mcl[,N.names]))
  out <- list(N=N, data.years=first.yr:last.yr,
              bad.samples=bad.samples)
  if(method=="clim"){
    out <- list(N=N, data.years=first.yr:last.yr,
                bad.samples=bad.samples, clim=clim.data.inc)
  }
  if(returnState)
    out$state <- state
  return(out)
}
################################################################################






################################################################################
## Forecasting dynamics to 2030
################################################################################

# Forecasts are run separately for each trend effect (temporal, early-breeding 
#    average temperature, annual precipitation) at each study plot
# This example is for early-breeding temperatures at the mid-elevation plot at 
#   the trailing edge (bs), though code is similar for other plots/trend effects.


# Loading in mcmc output
load("mcmc_out_BTBW_ebtemp_bs.gzip")

mc.bs.ebtemp <- as.mcmc.list(out.ebtemp.bs)

area.bs <- attr(out.ebtemp.bs, "Area")

# Forecasting to year 2030 so nFuture is from end of study period
#   on each plot to 2030
# M.new is number of new individuals which could potentially enter the population
#   in future years.
# method is either "clim" for early-breeding temp and annual precip or "time"
#   for temporal trend
# mean.inc and sd.inc are the beta1 parameter for the mean and se estimate for
#   the b1 paramater for the change in the climate variable in the study area over the
#   course of the study. This is used to predict future temperatures and should
#   be set to NA if method=time.
# mean.clim and sd.clim are the mean and sd of the climate variables at the study
#   plot over the course of the study.
f1.bs.ebtemp <- forecast(mc.bs.ebtemp, nFuture=11, Area=area.bs, M.new=3000,
                         method="clim",
                         mean.inc=clim_beta$b1[clim_beta$Site=="CWT"&clim_beta$Variable=="ebtemp"],
                         sd.inc=clim_beta$se[clim_beta$Site=="CWT"&clim_beta$Variable=="ebtemp"],
                         mean.clim=predparams$Mean[predparams$site=="CWT"&predparams$elev=="Mid"&predparams$Variable=="ebtemp"],
                         sd.clim=predparams$SD[predparams$site=="CWT"&predparams$elev=="Mid"&predparams$Variable=="ebtemp"],
                         thin=10, report=50) 

# 1 represents iterations with infinite superpopulations, 2 represents iterations
#    where M.new is too low
table(f1.bs.ebtemp$bad.samples)

apply(f1.bs.ebtemp$N, 1, quantile, prob=c(0.025, 0.5, 0.975))

save(f1.bs.ebtemp, file="forecastebtemp-bs.RData")
###################################################################################
