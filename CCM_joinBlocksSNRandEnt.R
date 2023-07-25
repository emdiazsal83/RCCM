# read in all .RData files containing maps of consecutive data points and time embedding for a certain variable
# and make a raster file of it
remove(list=ls())
library(raster)
library(ncdf4)

server <- "optimus.uv.es"
user <- "emiliano"
variable <- c("root_moisture","soil_moisture","surface_moisture",
              "air_temperature_2m","gross_primary_productivity", 
              "radiation", "latent_energy", "precipitation")[7]
reposResults <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/results_snr_ent_v2/"
#reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps_byBlock/", variable, "/", sep="")
#reposResults <- paste("timeEmbeddingMaps_byBlock/"variable, variable, "/", sep="")


filesRData <- dir(reposResults)
filesRData <- filesRData[grep(variable, filesRData)]
length(filesRData)
for(file in filesRData) load(paste(reposResults,file, sep=""))

aux <- ls()
indx_snrEnt <- grep("fsnrEnt_vals", aux)


aux2 <- strsplit(aux[indx_snrEnt], split="_")
variable <- unique(sapply(aux2, function(el) paste(el[3:(length(el)-1)], collapse="_")))

listRasts_snrEnt <- lapply(indx_snrEnt, function(i){
  res <- eval(parse(text=aux[i]))
  return(res)
  })
rast_snrEnt <- do.call(merge, listRasts_snrEnt)
plot(rast_snrEnt)


unlink(paste(reposResults,filesRData, sep=""))
snrEnt_nm <- paste("fsnrEnt", variable, sep="_")


png(file=paste(paste(reposResults, snrEnt_nm, sep=""), "png", sep="."))
plot(rast_snrEnt)
dev.off()

#library(rgdal)

writeRaster(rast_snrEnt, filename=paste(paste(reposResults, snrEnt_nm, sep=""),"nc", sep="."), overwrite=T, format="CDF")





