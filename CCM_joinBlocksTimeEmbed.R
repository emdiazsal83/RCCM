# read in all .RData files containing maps of consecutive data points and time embedding for a certain variable
# and make a raster file of it
remove(list=ls())
library(raster)
library(ncdf4)

server <- "optimus.uv.es"
user <- "emiliano"
variable <- c("root_moisture","soil_moisture","surface_moisture",
              "air_temperature_2m","gross_primary_productivity", 
              "radiation", "latent_energy", "precipitation")[8]
reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps_byBlock/", variable, "/", sep="")
#reposResults <- paste("timeEmbeddingMaps_byBlock/", variable, "/", sep="")


filesRData <- dir(reposResults)
for(file in filesRData) load(paste(reposResults,file, sep=""))

aux <- ls()
indx_consec <- grep("fconsec_vals", aux)
indx_timeDim <- grep("ftimeDim_vals", aux)

aux2 <- strsplit(aux[indx_consec], split="_")
variable <- unique(sapply(aux2, function(el) paste(el[3:(length(el)-3)], collapse="_")))
years <- unique(sapply(aux2, function(el) paste(el[(length(el)-2):(length(el)-1)], collapse="_")))

listRastsConsec <- lapply(indx_consec, function(i){
  res <- eval(parse(text=aux[i]))
  return(res)
  })
rastConsec <- do.call(merge, listRastsConsec)

listRastsTimeDim <- lapply(indx_timeDim, function(i){
  res <- eval(parse(text=aux[i]))
  return(res)
})
rastTimeDim <- do.call(merge, listRastsTimeDim)

reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps/", sep="")
consec_nm <- paste("fconsec", variable, years, sep="_")
timeDim_nm <- paste("timeDim", variable, years, sep="_")

png(file=paste(paste(reposResults, consec_nm, sep=""), "png", sep="."))
plot(rastConsec)
dev.off()
png(file=paste(paste(reposResults, timeDim_nm, sep=""), "png", sep="."))
plot(rastTimeDim)
dev.off()

#library(rgdal)

reposResults <- "/home/soulivanh/Documents/proyectos/CCM/timeEmbeddingMaps/"
writeRaster(rastConsec, filename=paste(paste(reposResults, consec_nm, sep=""),"nc", sep="."), overwrite=T, format="CDF")
writeRaster(rastTimeDim, filename=paste(paste(reposResults, timeDim_nm, sep=""), "nc", sep="."), overwrite=TRUE, format="CDF")





