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
library(reshape)


block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#block <- 7

numBlocksPerCombo <- 50
block2 <- block%%numBlocksPerCombo
block2[which(block2==0)] <- numBlocksPerCombo 

#print(block)
#q('no')
#block <- 1

pmAll <- proc.time()
print(paste("start block: ", block))

# server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"



# reposDataRaw <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/data/world"
# reposDataE <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps"
reposDataRaw <- paste("/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/data/world", sep="")
reposDataE <- paste("/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps", sep="")
# dir(reposDataRaw); dir(reposDataE)


variables <- c("soil_moisture","surface_moisture",
              "air_temperature_2m","gross_primary_productivity", 
              "radiation", "latent_energy", "precipitation")

combos <- as.data.frame(t(combn(variables, 2)))
colnames(combos) <- c("V1","V2")

combos$blockIni <- ((1:nrow(combos))-1)*numBlocksPerCombo + 1
combos$blockFin <- ((1:nrow(combos)))*numBlocksPerCombo 

indx <- findInterval(block, combos$blockIni)
V1 <- combos$V1[indx]
V2 <- combos$V2[indx]

print("variables")
print(V1); print(V2)

years <- 2001:2011
fileNamesRawV1 <- paste(paste(years, V1, sep="_"), "nc", sep=".")
fileNamesRawV2 <- paste(paste(years, V2, sep="_"), "nc", sep=".")
filesRawV1 <- paste(reposDataRaw, V1, fileNamesRawV1, sep="/")
filesRawV2 <- paste(reposDataRaw, V2, fileNamesRawV2, sep="/")

fileNamesE_V1 <- paste(paste(c("fconsec","timeDim"), V1, min(years), max(years), sep="_"), "nc", sep=".")
fileNamesE_V2 <- paste(paste(c("fconsec","timeDim"), V2, min(years), max(years), sep="_"), "nc", sep=".")
filesE_V1 <- paste(reposDataE,  fileNamesE_V1, sep="/")
filesE_V2 <- paste(reposDataE,  fileNamesE_V2, sep="/")

pm <- proc.time()
rawV1 <- stack(filesRawV1)
rawV2 <- stack(filesRawV2)
E_V1 <- stack(filesE_V1)
E_V2 <- stack(filesE_V2)
proc.time() - pm # 62 seconds
# plot(rawV1[[1:4]])
# plot(rawV2[[1:4]])
# plot(E_V1)
# plot(E_V2)

# get timeseries for pixels in a square block specified by a starting row and column and the edge size
getPixData <- function(stck, row, col, edgeSize){
  #data <- getValuesBlock(stck, row=row, nrows=min(edgeSize, nrow(stck)-row), col=col, ncols=min(edgeSize, ncol(stck)-col), format="")
  data <- crop(stck, extent(stck, r1=row, r2=min(nrow(stck), row+edgeSize-1), c1=col, c2=min(ncol(stck), col+edgeSize-1)))
  return(data)
}



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

getCause <- function(x, y, Ex, Ey){
  
  # lib- causa, target- efecto
  # E- time embedding de causa
  # convention for xmap names: cause_xmap_effect
  ts <- cbind(x, y)
  colnames(ts) <- c("x","y")
  
  lib_sizes <- c(25, 50, 100, 250, 500) #seq(100, 500, by=100)
  N <- length(lib_sizes)
  #pm <- proc.time()
  x_xmap_y <- ccm(ts, E=Ex, lib_column="x", target_column="y", lib_sizes=lib_sizes, num_samples = 50, random_libs = TRUE, replace = TRUE, silent = TRUE)
  #proc.time() - pm # 0.517 for 5 lib sizes
  x_xmap_y_means <- ccm_means(x_xmap_y)
  y_xmap_x <- ccm(ts, E=Ey, lib_column="y", target_column="x", lib_sizes=lib_sizes, num_samples = 50, random_libs = TRUE, replace = TRUE, silent = TRUE)
  y_xmap_x_means <- ccm_means(y_xmap_x)
  # plot(x_xmap_y_means$lib_size,  pmax(0, x_xmap_y_means$rho), type = "l", col = "red", xlab = "Library size", ylab = "Cross Map Skill (rho)",ylim=range(0.0,1.0))
  # lines(y_xmap_x_means$lib_size, pmax(0, y_xmap_x_means$rho), col = "cyan")
  
  
  #cause     <- which.max(c(x_xmap_y_means$rho[N], y_xmap_x_means$rho[N]))
  #cause_dif <-  x_xmap_y_means$rho[N] - y_xmap_x_means$rho[N]
  params <- expand.grid(lib_column= c("x","y"), target_column = c("x", "y"), tp = round(seq(-15,15, length.out=31)))
  params <- params[which(params$lib_column != params$target_column), ]
  params$E <- c(Ex, Ey)[match(params$lib_column, c("x","y"))]
  output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
    ccm(ts, E = params$E[i], lib_sizes = NROW(ts),
        random_libs = FALSE, lib_column = params$lib_column[i], target_column = params$target_column[i],
        tp = params$tp[i], silent = TRUE)
  }))
  output_ts <- cast(output, tp~lib_column, value="rho")
  
  indx_x <- which.max(output_ts[,"x"])
  indx_y <- which.max(output_ts[,"y"])
  #dif_x_y <- indx_x - indx_y
  
  maxTPx <- output_ts[indx_x,"tp"]
  maxTPy <- output_ts[indx_y,"tp"]
  if(length(indx_x)==0) maxTPx <- NA
  if(length(indx_y)==0) maxTPy <- NA
  
  # plot(output_ts[,"tp"], output_ts[,"x"], type="l", ylim=range(output_ts[,c("x","y")]))
  # lines(output_ts[,"tp"], output_ts[,"y"], col="red")
  # lines(output_ts[c(indx_x, indx_y),1], c(output_ts[indx_x, "x"], output_ts[indx_y, "y"]), col=c("black","red"), type="p", cex=2)
  
  return(c(x_xmap_y_means$rho, y_xmap_x_means$rho,  maxTPx, maxTPy))
  
}

# apply a pixel by pixel function ina hierarchical way first by block (one per slurm process)
# and then within each block also parallelize by row and col
applyFuncStack <- function(stck,  func, paramsFunc, mcCoresRow, mcCoresCol){
  
  rows <- 1:nrow(stck)
  cols <- 1:ncol(stck)
  vals <- as.array(stck)
  res <- mcmapply(function(row){
    #print(paste("rowNum", which(row==rows)))
    res <- mcmapply(function(col){
      # row <- 1; col <- 1
      res <- do.call(func, list(x=vals[row, col,], paramsFunc=paramsFunc))
      return(res)
    }, col=cols, SIMPLIFY="array", mc.cores=mcCoresCol)
    return(res)
  }, row=rows, SIMPLIFY="array", mc.cores=mcCoresRow)
  
  res <- raster(t(res), template=stck)
  return(res)
}

applyFuncStack2 <- function(rawV1, rawV2, E_V1, E_V2,  func, mcCoresRow, mcCoresCol){
  
  rows <- 1:nrow(rawV1)
  cols <- 1:ncol(rawV1)
  rawV1_vals <- as.array(rawV1)
  rawV2_vals <- as.array(rawV2)
  E_V1_vals <- as.array(E_V1)
  E_V2_vals <- as.array(E_V2)
  pm <- proc.time()
  res <- mcmapply(function(row){
    # row <- 74
    #print(paste("rowNum", which(row==rows)))
    res1 <- mcmapply(function(col){
      # row <- 97; col <- 97
      # print(paste("row-col:", row, "-", col))
      x <- rawV1_vals[row, col,]
      y <- rawV2_vals[row, col,]
      Ex <- max(E_V1_vals[row, col, 2],2)
      Ey <- max(E_V2_vals[row, col, 2],2)
      if(E_V1_vals[row,col,1]==506 & E_V2_vals[row,col,1]==506 & !is.na(Ex) & !is.na(Ey)){
        # plot(rawV1_vals[row, col,])
        # plot(rawV2_vals[row, col,])
        res2 <- do.call(func, list(x, y, Ex, Ey))
      } else{
        res2 <- rep(NA,12)
      }
      #print(paste("res2", paste(res2, sep="-")))
      return(res2)
    }, col=cols, SIMPLIFY="array", mc.cores=mcCoresCol)
    #print(paste("dim(res1): ", paste(dim(res1), collapse="-")))
    return(res1)
  }, row=rows, SIMPLIFY="array", mc.cores=mcCoresRow)
  proc.time() - pm 
  # 38 seconds for edgeSize=4 for mcCoresRow = mcCoresCol = 5
  
  res <- lapply(1:dim(res)[1], function(i) raster(t(res[i,,]), template=E_V1))
  res <- stack(res)
  return(res)
}

# lets figure out how many blocks we need to divide the world into to parallelize
# Each node has 28 cpus so each block will have 28 cores, lets use 25, 5 for rows 5 for columns within each block
# lets say we only want to process 0.1gb of info at any one block
# how many pixels per block
numPix <- round(1e8/(506*8)) # 0.1e9 bytes/block /  506 timeobs/pixels * 8 bytes /timeobs = 24704 pixels / block
edgeSize <- round(sqrt(numPix))
rows <- c(seq(1, nrow(rawV1), by=edgeSize))
cols <- c(seq(1, ncol(rawV1), by=edgeSize))
numSuperRows <- length(rows)
numSuperCols <- length(cols)
numSuperCells <- numSuperRows*numSuperCols # so we would need 50 tasks
numPix*numSuperCells
nrow(rawV1)*ncol(rawV1)



# lets go by column as usual in R
superRow <- c(seq(4),5)[match(block2 %% numSuperRows, c(seq(4), 0))]
superCol <- ceiling(block2/numSuperRows) 

row <- (superRow-1)*edgeSize+1
col <- (superCol-1)*edgeSize+1

pm <- proc.time() 
rawV1 <- getPixData(rawV1, row, col, edgeSize)
rawV2 <- getPixData(rawV2, row, col, edgeSize)
E_V1 <- getPixData(E_V1, row, col, edgeSize)
E_V2 <- getPixData(E_V2, row, col, edgeSize)
proc.time() - pm # 3.3 mins
# plot(rawV1[[1:4]])
# plot(rawV2[[1:4]])
# plot(E_V1)
# plot(E_V2)


# rawV1 <- getPixData(rawV1, 1, 122, 4)
# rawV2 <- getPixData(rawV2, 1, 122, 4)
# E_V1 <- getPixData(E_V1, 1, 122, 4)
# E_V2 <- getPixData(E_V2, 1, 122, 4)
# plot(rawV1[[1:4]])
# plot(rawV2[[1:4]])
# plot(E_V1)
# plot(E_V2)

pm <- proc.time()
cause <- applyFuncStack2(rawV1, rawV2, E_V1, E_V2,  func="getCause", mcCoresRow=5, mcCoresCol=5)
proc.time() - pm # 7.8 mins for block 14, mins for block 1 
# as.array(cause)
# plot(cause)
# plot(stack(is.na(cause[[1]]),cause[[1]]-cause[[2]], cause[[3]]-cause[[4]]))


cause_nm <- paste("ccm_cause_vals", V1, V2, min(years), max(years), block2, sep="_")
assign(cause_nm, cause)


# save results

reposResults <- paste("/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses_byBlock/", paste(V1, "vs", V2, sep="_"), "/",sep="")
# reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses_byBlock/", paste(V1, "vs", V2, sep="_"),"/", sep="") 
dir.exists(reposResults)
dir.create(reposResults)
fileName <- paste(paste(V1, "vs", V2, min(years), max(years), block2, "ccmCauses", sep="_"), "RData", sep=".")
save(list=c(cause_nm), file=paste(reposResults, fileName, sep=""))

proc.time()-pmAll #

print(paste("finish block: ", block))

