  # source("../../media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/downloadNCDF_byPix.R")
# download some ncdf4 files, in a more efficient way: build a function to read each
remove(list=ls())
library(ncdf4)
library(raster)
library(maptools)
library(rgdal)
library(rgeos)
library(rasterVis) # rasterTheme, levelplot
library(parallel)
library(rEDM)
#library(smooth)
#source("/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_entropy/func_entropy_v1.R")
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/"
#repos <- paste("/home/emiliano/causaLearner/", sep="")
setwd(repos)
source("./pkg_entropy/func_entropy_v1.R")


#block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


#print(block)
#q('no')
block <- 38

pmAll <- proc.time()
print(paste("start block: ", block))

#server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"
reposData <- "/home/emiliano/Documents/ISP/proyectos/causality/causalCarbon/carboncausality/data/ccmData/" 
#reposData <- paste("/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/data/world", sep="")
#reposData <- "/home/emiliano/CCM/data/"
# dir(reposData)

variable <- c("root_moisture","soil_moisture","surface_moisture",
              "air_temperature_2m","gross_primary_productivity", 
              "radiation", "latent_energy", "precipitation")[1]

print(variable)
years <- 2001:2011
fileNames <- paste(paste(years, variable, sep="_"), "nc", sep=".")
files <- paste(reposData, variable, fileNames, sep="/")

pm <- proc.time()
dat <- stack(files)
proc.time() - pm # 30 seconds
#plot(dat[[1:4]])

# get timeseries for pixels in a square block specified by a starting row and column and the edge size
getPixData <- function(stck, row, col, edgeSize){
  #data <- getValuesBlock(stck, row=row, nrows=min(edgeSize, nrow(stck)-row), col=col, ncols=min(edgeSize, ncol(stck)-col), format="")
  data <- crop(stck, extent(stck, r1=row, r2=min(nrow(stck), row+edgeSize-1), c1=col, c2=min(ncol(stck), col+edgeSize-1)))
  return(data)
}

# test getPixData
# pm <- proc.time()
# miniStck <- getPixData(dat, 50, 587, 4)
# proc.time() - pm # 30 seconds

# functions to apply to each pixel-timeseries

# max number of consecutive non-empty observations 
cumOnes <- function(x){ 
  tmp1 <- cumsum(x)
  tmp2 <- (!x)*tmp1
  return(tmp1-cummax(tmp2))
}
consecValid <- function(x, paramsFunc) max(cumOnes(!(is.na(x)|is.nan(x))))

# get % of time windows with no missing data
getUsefulWindows <- function(x, paramsFunc){
  params <- paramsFunc$params
  
  
  noMiss <- apply(params, 1, function(row){
    train <- c(row["trainStart"], row["trainStart"] + row["trainLength"] -1)
    test <- c(row["testStart"], row["testStart"] + row["testLength"])
    resTr <- any(is.na(x[train[1]:train[2]])|is.nan(x[train[1]:train[2]]))
    resTe <- any(is.na(x[test[1]:test[2]])|is.nan(x[test[1]:test[2]]))
    return((!resTr)&(!resTe))
  })
  return(sum(noMiss)/length(noMiss)*100)
}


# obtain the time embedding of a time series, paramsFunc is a list with params: a data-frame with 
# trainStart, testStart, trainLength and testLength and minConsec, the minimum number of consecutive
# non-empty observations necessary to produce a result (otherwise NA is given)
getTimeDim_medDiff <- function(x, paramsFunc){
  params <- paramsFunc$params
  minConsec <- paramsFunc$minConsec
  # ts <- getValuesBlock(dat, row=100, nrows=1, col=100, ncols=1); ts <- ts[1,]; summary(ts)
  
  
  consecVal <- cumOnes(!(is.na(x)|is.nan(x))) 
  
  if(max(consecVal)<minConsec) return(NA)
  
  indxMax <- which.max(consecVal)
  maxConsecVal <- consecVal[indxMax]
  indxFirst <- indxMax-maxConsecVal+1
  tsConsec <- x[indxFirst:indxMax]
  
  rhoDist <- apply(params, 1, function(row){
    train <- c(row["trainStart"], row["trainStart"] + row["trainLength"] -1)
    test <- c(row["testStart"], row["testStart"] + row["testLength"])
    simp <-simplex(tsConsec, train, test)
    return(simp$rho)
  })
  # boxplot(t(rhoDist))
  dif_rhoDist <- diff(rhoDist)
  # boxplot(t(dif_rhoDist))
  med_rho <- apply(dif_rhoDist, 1, median)
  # plot(med_rho); abline(h=0, col="red)
  indx_dec  <- which(med_rho  < 0)[1]
  return(indx_dec)
}


calcSNR_ent <- function(x, paramsFunc){
  
  minConsec <- paramsFunc$minConsec
  # ts <- getValuesBlock(dat, row=100, nrows=1, col=100, ncols=1); ts <- ts[1,]; summary(ts)
  
  
  consecVal <- cumOnes(!(is.na(x)|is.nan(x))) 
  
  if(max(consecVal)<minConsec) return(rep(NA,8))
  
  indxMax <- which.max(consecVal)
  maxConsecVal <- consecVal[indxMax]
  indxFirst <- indxMax-maxConsecVal+1
  tsConsec <- x[indxFirst:indxMax]
  if(sd(tsConsec)>0) tsConsec <- (tsConsec-mean(tsConsec))/sd(tsConsec)
  
  # snr1 - n = diff(x); snr = 20*log10(std(n)/std(x))
  ns <- diff(tsConsec)
  snr <- 20*log(sd(tsConsec)/sd(ns),10)
  
  # snr2 - xp=smooth(x); n = x-xp; snr = 20*log10(std(n)/std(x))
  
  pm <- proc.time()
  x_sm <- movavg(tsConsec, n=3, type="s")
  #x_sm <- sma(tsConsec, h=0, silent=TRUE, order=3)
  proc.time()-pm
  
  
  
  #plot(tsConsec, type="l")
  #lines(x_sm, col="red")
  ns_sm <- tsConsec-x_sm#$fitted[,1]
  snr_sm <- 20*log(sd(tsConsec)/sd(ns_sm),10)
  
  
  # entropy - Shannon_MaxEnt1, Shannon_PSD_SzegoT, Shannon_KDP
  
  h_ts_max <- Shannon_MaxEnt1(tsConsec)
  #h_ts_psd <- Shannon_PSD_SzegoT(tsConsec)
  h_ts_kdp <- Shannon_KDP(tsConsec)
  h_ns_max <- Shannon_MaxEnt1(ns)
  #h_ns_psd <- Shannon_PSD_SzegoT(ns)
  h_ns_kdp <- Shannon_KDP(ns)
  h_ns_sm_max <- Shannon_MaxEnt1(ns_sm)
  #h_ns_sm_psd <- Shannon_PSD_SzegoT(ns_sm)
  h_ns_sm_kdp <- Shannon_KDP(ns_sm)
  
  res <- c(snr, snr_sm, h_ts_max, h_ts_kdp,h_ns_max, h_ns_kdp, h_ns_sm_max, h_ns_sm_kdp)
  
  return(res)
}


# apply a pixel by pixel function ina hierarchical way first by block (one per slurm process)
# and then within each block also parallelize by row and col
applyFuncStack <- function(stck,  func, paramsFunc, mcCoresRow, mcCoresCol){
  
  rows <- 1:nrow(stck)
  cols <- 1:ncol(stck)
  vals <- as.array(stck)
  # mcCoresRow <- 4; mcCoresCol=1
  pm <- proc.time()
  res <- mcmapply(function(row){
    #print(paste("rowNum", which(row==rows)))
    res <- mcmapply(function(col){
      # row <- 2; col <- 19
      #print(paste("row= ", row, " ; col= ", col))
      res <- do.call(func, list(x=vals[row, col,], paramsFunc=paramsFunc))
      return(res)
    }, col=cols, SIMPLIFY="array", mc.cores=mcCoresCol)
    return(res)
  }, row=rows, SIMPLIFY="array", mc.cores=mcCoresRow)
  proc.time() - pm 
  
  res <- lapply(1:dim(res)[1], function(i) raster(t(res[i,,]), template=stck))
  res <- stack(res)
  
  return(res)
}

# lets figure out how many blocks we need to divide the world into to parallelize
# Each node has 28 cpus so each block will have 28 cores, lets use 25, 5 for rows 5 for columns within each block
# lets say we only want to process 0.1gb of info at any one block
# how many pixels per block
numPix <- round(1e8/(506*8)) # 0.1e9 bytes/block /  506 timeobs/pixels * 8 bytes /timeobs = 24704 pixels / block
edgeSize <- round(sqrt(numPix))
rows <- c(seq(1, nrow(dat), by=edgeSize))
cols <- c(seq(1, ncol(dat), by=edgeSize))
numSuperRows <- length(rows)
numSuperCols <- length(cols)
numSuperCells <- numSuperRows*numSuperCols # so we would need 50 tasks
numPix*numSuperCells
nrow(dat)*ncol(dat)



# lets go by column as usual in R
superRow <- c(seq(4),5)[match(block %% numSuperRows, c(seq(4), 0))]
superCol <- ceiling(block/numSuperRows) 

row <- (superRow-1)*edgeSize+1
col <- (superCol-1)*edgeSize+1
 
pm <- proc.time() 
miniStck <- getPixData(dat, row, col, edgeSize)
proc.time() - pm # 88 seconds
# plot(miniStck[[1:4]])

# number of consecutive non-empty observations
#pm <- proc.time()
#fconsec_vals <- applyFuncStack(stck=miniStck, func=consecValid, paramsFunc=list(), mcCoresRow=5, mcCoresCol=5)
#proc.time() - pm # 15 seconds
#as.array(fconsec_vals)
# plot(fconsec_vals)

# time embedding
#trainStart <- seq(1,100, length.out=10) 
#testStart <- round(seq(277, 376, length.out=10))
#trainLength <- 46
#testLength <- 46

#params <- expand.grid(trainStart=trainStart, testStart=testStart, trainLength=trainLength, testLength=testLength)
minConsec <- 376+46

# % of useful windows
#pm <- proc.time()
#fusefulWinds_vals <- applyFuncStack(stck=miniStck, func=getUsefulWindows, paramsFunc=list(params=params), mcCoresRow=5, mcCoresCol=5)
#proc.time() - pm # 76 seconds
# plot(stack(fconsec_vals, fusefulWinds_vals, fconsec_vals>300 ,fusefulWinds_vals>10))

# mcCoresRow <- mcCoresCol <- 1
mcCoresRow <- mcCoresCol <- 5
 
pm <- proc.time()
fsnrEnt_vals <- applyFuncStack(stck=miniStck, func="calcSNR_ent", paramsFunc=list(minConsec=minConsec), mcCoresRow=mcCoresRow, mcCoresCol=mcCoresCol)
proc.time() - pm # 36 mins locally, 8.6 mins in erc
#as.array(fsnrEnt_vals)
# plot(fsnrEnt_vals)

#fconsec_vals_nm <- paste("fconsec_vals", variable, min(years), max(years), block, sep="_")
fsnrEnt_vals_nm <- paste("fsnrEnt_vals", variable, block, sep="_")
#assign(fconsec_vals_nm, fconsec_vals)
assign(fsnrEnt_vals_nm, fsnrEnt_vals)

# save results

reposResults <- "/home/emiliano/CCM/results_snr_ent/"
#reposResults <- paste("/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps_byBlock/",variable, "/", sep="")
#reposResults <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps/", sep="") 
fileName <- paste(paste(variable, min(years), max(years), block, "snrEnt", sep="_"), "RData", sep=".")
save(list=c(fsnrEntDim_vals_nm), file=paste(reposResults, fileName, sep=""))

proc.time()-pmAll #

print(paste("finish block: ", block))

