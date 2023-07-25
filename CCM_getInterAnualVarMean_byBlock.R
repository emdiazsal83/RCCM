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

block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#print(block)
#q('no')
#block <- 14

pmAll <- proc.time()
print(paste("start block: ", block))

#server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"
# reposData <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/data/"
reposData <- paste("/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/data/world", sep="")
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



getMeanIAVar <- function(x, paramsFunc){
  seqsDays <- lapply(1:46, function(ini) seq(ini, by=46, length.out=11))
  varsIA <- sapply(seqsDays, function(seqsDay){
    res <- var(x[seqsDay], na.rm=T)
    return(res)
  })
  
  res <- mean(varsIA, na.rm=T)
  #if(is.nan(res)) res <- NA
  
  return(res)
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


pm <- proc.time()
fmeanIAvars_vals <- applyFuncStack(stck=miniStck, func=getMeanIAVar, paramsFunc=list(), mcCoresRow=5, mcCoresCol=5)
proc.time() - pm # 36 mins locally, 8.6 mins in erc
#as.array(fmeanIAvars_vals)
# plot(fmeanIAvars_vals)

fmeanIAvars_vals_nm <- paste("fmeanIAvars_vals", variable, min(years), max(years), block, sep="_")
assign(fmeanIAvars_vals_nm, fmeanIAvars_vals)


# save results

reposResults = "/home/emiliano/CCM/results_meanIAvar_byBlock/"
#reposResults <- paste("/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps_byBlock/",variable, "/", sep="")
#reposResults <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps/", sep="") 
fileName <- paste(paste(variable, min(years), max(years), block, "consec", "timeEmbed", sep="_"), "RData", sep=".")
save(list=c(fmeanIAvars_vals_nm), file=paste(reposResults, fileName, sep=""))

proc.time()-pmAll #

print(paste("finish block: ", block))

