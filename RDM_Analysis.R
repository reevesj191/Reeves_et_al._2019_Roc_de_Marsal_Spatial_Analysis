#
# This code applies the moving window functions and saves the results in a data file
# that is then used by the main paper to make the figures.
#
# AUTHOR: Jonathan S Reeves and Shannon P. McPherron
#
# Last modified November 14, 2018
#
############################

## This code requires the following libraries or packages
##    sp, spdep
## Additionally, the code is optimized for multiple threads using these packages
##    doParallel and foreach
## If threading causes problems, they can easily be written out of the moving_window_par function


## The following three functions are the core of the moving window analysis.
## The first computes the test statistic for the window around each artifact.
## The second determines what spatial scale gives the most statistically signicant clustering results.
## The third one computes the clustering at that scale.

moving_window_par <- function(dat, eq, neighborhood, colname, ncores = 1) {
  ## dat = the spatial dataframe that holds the data to be tested
  ## eq = the function to be applied to each neighborhood
  ## neighborhood = the size of the moving window
  ## colname = the name given to the column that holds the result
  ## ncores = the number of CPU cores to be used (careful - a number larger than the available cores with crash the program)  

  ## The par in the function name refers to the parallelization done here.
  library(doParallel)
  library(foreach)
  
  print(paste("  Starting moving_window_par using",ncores,"CPU cores."))
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  results <- foreach(i = 1:nrow(dat@data), .combine = c, .packages= c('rgeos', "sp")) %dopar% {
    buff <- gBuffer(spgeom = dat[i,], width = neighborhood, capStyle = "ROUND", joinStyle = "ROUND")
    eq(subset(dat, over(dat, buff) == 1))
  }
  dat$eq = results
  names(dat)[names(dat)=="eq"] <- colname
  stopCluster(cl)
  return(dat)
}  

## This function is not used in this version of our code.  It was an earlier version
## that is now more optimized with the mi.dist.par() function that follows.  However,
## we have left it here because it returns a data.frame that is more useful for
## examining the effects of scale on the Moran I's index.

mi.dist.plot <- function(dat, var, lower, icr, n) {
  ## dat = the spatial dataframe that holds the data that you want to test 
  ## var = the variable from the moving window results that you are testing
  ## lower = the lower bound for your neighborhood
  ## icr = the increment that you want your neighborhood to increase by 
  ## n = the number of increments to run
  
  q <- data.frame("band" = 1:n,
                  "neighborhood" = 1:n,
                  "morans.I" = 1:n,
                  "p.value" = 1:n,
                  "class" = 1:n)
  dmax <- lower
  for (i in 1:n){
    print(i) ## a progress bar to make sure its actually doing something and not crashed. 
    dmax <- dmax + icr
    q$neighborhood[i] <- dmax
    nb <- dnearneigh(dat, d1 = lower, d2 = dmax)
    wt <- nb2listw(nb, zero.policy = TRUE)
    mt <- moran.test(var, listw = wt, zero.policy = TRUE, na.action = na.exclude)
    q$morans.I[i] <- mt$estimate[1]
    q$p.value[i] <- mt$p.value
  }
  f <- list(correlo.dat = q, scale.band = q$neighborhood[q$morans.I == max(q$morans.I)])
  return(f)
}

## This function is the same as the one just above except that it has been optimized
## to use multiple cores and to return just the statistic that we need for the 
## analysis that follows after.
mi.dist.par <- function(dat, var, lower, icr, n, ncores = 1) {
  ## dat = the spatial dataframe that holds the data that you want to test 
  ## var = the variable from the moving window results that you are testing
  ## lower = the lower bound for your neighborhood
  ## icr = the increment that you want your neighborhood to increase by 
  ## n = the number of increments to run

  ## The par in the function name refers to the parallelization done here.
  library(doParallel)
  library(foreach)

  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  print(paste("  Starting mi.dist.par using",ncores,"CPU cores."))
  
  results <- foreach(dmax = seq(lower + icr, by = icr, length.out = n),
                     .combine = c, .packages= c('spdep')) %dopar% {
                       nb <- dnearneigh(dat, d1 = lower, d2 = dmax)
                       wt <- nb2listw(nb, zero.policy = TRUE)
                       mt <- moran.test(var, listw = wt, zero.policy = TRUE, na.action = na.exclude)
                       mt$estimate[1]
                     }
  stopCluster(cl)
  return(seq(lower + icr, by = icr, length.out = n)[which.max(results)])
}



local.moransI <- function(dat, lower, upper, var, colname) {
  ## dat <- must be a matrix of coordinates. SpatialPointDataFrames work on their own. If you are using a SpatialPolygon then you need to the coordinates() function retrieve its coordinate matrix. 
  ## lower <- minimum neighborhood distance. 
  ## upper <- maximum neighborhood distance.
  ## var <- a continuous variabile you are interested in.
  ## colname <- the name of the column for the results

  print("  Computing local.moransI.")
  
  nb <- dnearneigh(x = dat, d1 = lower, d2 = upper)
  wt <- nb2listw(nb, zero.policy = TRUE)
  lm <- localmoran(var, wt, zero.policy = TRUE, p.adjust.method = "bonferroni", na.action = na.exclude)
  lmI <- as.data.frame(lm)
  dat$LMI_Z <- lmI$Z.Ii 
  dat$S_COT <- "Nonsig" #### Where Cluster Outlier Data will be stored
  dat$S_COT[dat$LMI_Z >= 2 & var > mean(var)] <- "HH"
  dat$S_COT[dat$LMI_Z >= 2 & var < mean(var)] <- "LL"
  dat$S_COT[dat$LMI_Z <= -2 & var < mean(var)] <- "LH"
  dat$S_COT[dat$LMI_Z <= -2 & var > mean(var) ] <- "HL"
  names(dat)[names(dat)=="LMI_Z"] <- paste0(colname,"_Z")
  names(dat)[names(dat)=="S_COT"] <- paste0(colname,"_COT")
  return(dat)
} 

############################
#
# What follows are a series of functions that are used to characterize each neighborhood.
# These functions are passed to the moving.window.par() function listed above.
#
############################

# Density of Artifacts
GPF.Density <- function(dat) {
  nrow(dat)  
} 

# Burning Ratio (Standardized by flakes)
# Name stands for burned flakes to unburned flakes
bflk2ubflk <- function(dat) {
  (length(subset(dat, DATACLASS == "PROXFLAKE" & BURNED == "YES")) + (length(subset(dat, DATACLASS == "COMPFLAKE" & BURNED == "YES"))))/ 
    (length(subset(dat, DATACLASS == "PROXFLAKE")) + (length(subset(dat, DATACLASS == "COMPFLAKE"))))
}

# Scraper to flake Ratio
scraper2flake2 <- function(dat) {
  length(subset(dat, scraper == "TRUE")) / 
    (length(subset(dat, DATACLASS == "PROXFLAKE")) + (length(subset(dat, DATACLASS == "COMPFLAKE"))))
}

# Flake to core ratio
flk2core <- function(dat) {
  length(subset(dat, core == TRUE)) / (length(subset(dat, DATACLASS == "PROXFLAKE")) + (length(subset(dat, DATACLASS == "COMPFLAKE"))))
}

# Cortex To Mass ratio
cortex2mass <- function(dat){
  dat <- subset(dat, is.na(cortical_SA) == FALSE & is.na(WEIGHT) == FALSE)
  sum(dat$cortical_SA) / sum(dat$WEIGHT)
}

# Breakage ratio (Name stands for unburned flake breakage)
ubflkbreak <- function(dat) {
  length(subset(dat, DATACLASS == "PROXFLAKE" & BURNED == "NO")) / 
    (length(subset(dat, DATACLASS == "PROXFLAKE" & BURNED == "NO")) + (length(subset(dat, DATACLASS == "COMPFLAKE" & BURNED == "NO"))))
}

# Flake Weight
med.weight.comp.ub.flakes <- function(dat){
  dat <- subset(dat, BURNED == "NO" & DATACLASS == "COMPFLAKE")
  median(dat$WEIGHT, na.rm = TRUE)
}

############################
#
# Here begins the main code.
#
############################

#### Load libraries ####
library(sp)
library(spdep)

#### Load data ####
## Minimally the data consist of XY coordinates, layers, and columns with attributes on each 
## artifact that go into characterizing each neighborhood using the functions list above.
xdata <- readRDS("RDM_MovingWindow_Data.RDS")
xdata$LEVEL <- as.factor(xdata$LEVEL)

#### Script Parameters #####
moranIstartband <- .01
neighborhood.size = .3 # Set the size of the neighborhood
ncores = 4  # The moving window function is set up to parallelize work. This should be set to the specifications of your computer. It is safest to level this at 1. Inappropriate allocation of proccessor cores can result in catastrophic crashes.
sample.size <- 1000 # This is for subsampling when using the moran's I for levels that are very big.

# The name of the moving window functions that you wish to have applied to the data set.
fun.names <- c( 
   GPF.Density,
   bflk2ubflk,
   scraper2flake2,
   flk2core,
   cortex2mass,
   ubflkbreak,
   med.weight.comp.ub.flakes
)

# The names of the columns where the results from the above functions will be stored.
result.cols <- c( 
  "density",
  "burning",
  "Scraper2flake",
  "core2flake",
  "cortex2mass",
  "breakage",
  "flake.weight"
)

# The archaeological layers for which the analysis will be applied
levels <- c("09",
            "08",
            "07",
            "05",
            "04")

#### Analysis Script ####

RESULT = NA

for(i in levels) {
  level <- i
  print(paste('Working on Layer',level))
  RDM <- subset(xdata, LEVEL == level) ## Get just the level we are working on
  level.size = nrow(RDM)
  
  coordinates(RDM) <- c("X", "Y") ## Turns it into a spatial data set

  ## Loop through each of the functions
  for(j in 1:length(fun.names)){

    print(paste(' Looking at',result.cols[j]))

    ## Apply the moving window analysis
    RDM <- moving_window_par(dat = RDM,
                             eq = fun.names[[j]],
                             neighborhood = neighborhood.size,
                             colname = result.cols[j],
                             ncores = ncores )
    
    ## Replace infinte results (like dividing by zero) with zero
    RDM@data[!is.finite(RDM@data[,result.cols[j]]), result.cols[j]] <- 0
    
    ## Adjust neighborhood flake weight by the layer median flake weight
    if (result.cols[j]=="flake.weight") {
      RDM@data[,"flake.weight"] = RDM@data[,"flake.weight"] / median(RDM@data$WEIGHT, na.rm = TRUE)
    }
    
    ## Before computing spatial scale for the analysis,
    ## if the sample is large, subsample it.
    if (level.size > sample.size) {
      RDM.sample <- RDM[sample(nrow(RDM), sample.size),] }
    else {
      RDM.sample <- RDM }

    ## Get the scale at which there is the most clustering
    correlation  <- mi.dist.par(dat = RDM.sample,
                                var = RDM.sample@data[,result.cols[j]],
                                lower = moranIstartband,
                                icr = .025,
                                n = 20,
                                ncores = ncores)

    #print(correlation)

    ## Computer MoransI at the scale determined in the previous step
    RDM <- local.moransI(dat = RDM,
                         lower = moranIstartband,
                         upper = correlation,
                         var = RDM@data[,result.cols[j]], 
                         colname = result.cols[j])
  }

  ## After repeating each of the test functions for a particular layer
  ## save the results.  Then repeat for the next layer and append
  ## the new results to the of the results output.
  if (length(RESULT)<2) {
    print("Creating Results table")
    RESULT <- RDM }
  else {
      print("Appending Results Table")
      RESULT <- rbind(RESULT,RDM) }
}

## After doing all of the test functions for all of the layers
## write the results in an ASCII file that will be read and plotted
## by the main paper.
write.csv(RESULT, "RDM_MovingWindow_Results.csv")
