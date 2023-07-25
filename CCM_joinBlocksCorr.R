# read in all .RData files containing maps of consecutive data points and time embedding for a certain variable
# and make a raster file of it
remove(list=ls())
library(raster)
library(ncdf4)

server <- "optimus.uv.es"
user <- "emiliano"

reposResults <- "/home/emiliano/Drives/erc/CCM/ccmCorrs_byBlock/"
#reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses_byBlock/", sep="")
#reposResults <- paste("timeEmbeddingMaps_byBlock/", paste(V1, "vs", V2, sep="_"), "/", sep="")

aux <- dir(reposResults)
sapply(aux, function(combo) length(dir(paste(reposResults, combo, sep=""))))
table(sapply(aux, function(combo) length(dir(paste(reposResults, combo, sep="")))))

variables <- c("soil_moisture","surface_moisture",
               "air_temperature_2m","gross_primary_productivity", 
               "radiation", "latent_energy", "precipitation")

variables <- c("surface_moisture","root_moisture")

combos <- as.data.frame(t(combn(variables, 2)))
colnames(combos) <- c("V1","V2")


for(i in 1:nrow(combos)){
  #i <- 1
  V1 <- combos$V1[i]
  V2 <- combos$V2[i]
  print("**************************")
  print(i)
  print(paste(V1, "vs", V2, sep="_"))

  reposResults <- "/home/emiliano/Drives/erc/CCM/ccmCorrs_byBlock/"
  #reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses_byBlock/", sep="")
  filesRData <- dir(paste(reposResults, paste(V1, "vs",V2,sep="_"), sep="/"))
  for(file in filesRData) load(paste(reposResults, paste(V1, "vs",V2,sep="_"), "/", file, sep=""))

  aux <- ls()
  indx_corr <- grep( paste("ccm_corr_vals", V1, V2, sep="_"), aux)
  length(indx_corr)
  

  aux2 <- strsplit(aux[indx_corr], split="_")
  years <- unique(sapply(aux2, function(el) paste(el[(length(el)-2):(length(el)-1)], collapse="_")))

  listRastsCorr <- lapply(indx_corr, function(i){
    res <- eval(parse(text=aux[i]))
    return(res)
  })
  
  
  rastCorr <- do.call(merge, listRastsCorr)

  names(rastCorr) <- c("pearson","spearman")
  plot(rastCorr)
  #plot(getValues(rastCause[[1]]),getValues(rastCause[[2]])); abline(a=0, b=1, col="red")
  #plot(stack(is.na(rastCause[[1]]),rastCause[[1]]>rastCause[[2]], rastCause[[3]]<rastCause[[4]]))

  #reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses/", sep="")
  #reposResults <- "/home/soulivanh/Documents/proyectos/CCM/ccmCauses/"
  reposResults <- "/home/emiliano/Drives/erc/CCM/ccmCorrs/"
  corr_nm <- paste("ccm_corr", paste(V1, "vs", V2, sep="_"), years, sep="_")


  png(file=paste(paste(reposResults, corr_nm, sep=""), "png", sep="."))
  plot(rastCorr)
  dev.off()

  #library(rgdal)
  
  reposResults <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/ccmCorrs/"
  #reposResults <- "/home/emiliano/Drives/erc/CCM/ccmCorrs/"
  #reposResults <- "/home/soulivanh/Documents/proyectos/CCM/ccmCorrs/"
  writeRaster(rastCorr, filename=paste(paste(reposResults, corr_nm, sep=""),"nc", sep="."), overwrite=T, format="CDF")
}

plot(rastCorr[["pearson"]], main = "pearson correlation")
plot(rastCorr[["spearman"]], main = "spearman correlation")
plot(rastCorr[["spearman"]]>0.7)

sum(as.array((rastCorr[["spearman"]]>0.7)==1), na.rm=T)/sum(!is.na(as.array(rastCorr[["spearman"]])))

reposResults <- "/home/emiliano/Drives/erc/CCM/ccmCorrs/"
#reposResults <- "/home/soulivanh/Documents/proyectos/CCM/ccmCauses/"
files <- dir(reposResults)
indx <- which(sapply(strsplit(files, "\\."), function(el) el[2])=="nc")
files <- files[indx]

rast <- stack( paste(reposResults, files, sep=""))
plot(rast)
dim(rast)
21*12

dim(rastCorr)
plot(rastCorr[[1]])
ext = drawExtent()

reposResults <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/ccmCorrs/"
png(file=paste(paste(reposResults, "corrSMs_pearson", sep=""), "png", sep="."), width=1200, height=720)
par(mar = c(0, 0, 0, 0.5))
plot(crop(rastCorr[["pearson"]],ext),axes = F, box = F, col = terrain.colors(100))
dev.off()

reposResults <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/ccmCorrs/"
png(file=paste(paste(reposResults, "corrSMs_spearman", sep=""), "png", sep="."), width=1200, height=720)
par(mar = c(0, 0, 0, 0.5))
plot(crop(rastCorr[["spearman"]],ext), axes = F, box = F, col = terrain.colors(100))
dev.off()
