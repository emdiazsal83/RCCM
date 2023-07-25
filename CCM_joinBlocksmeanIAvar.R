# read in all .RData files containing maps of consecutive data points and time embedding for a certain variable
# and make a raster file of it
remove(list=ls())
library(raster)
library(ncdf4)

server <- "optimus.uv.es"
user <- "emiliano"
variable <- c("root_moisture","soil_moisture","surface_moisture",
              "air_temperature_2m","gross_primary_productivity", 
              "radiation", "latent_energy", "precipitation")[1]
reposResults <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/results_meanIAvar_byBlock/"
#reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/timeEmbeddingMaps_byBlock/", variable, "/", sep="")
#reposResults <- paste("timeEmbeddingMaps_byBlock/", variable, "/", sep="")


filesRData <- dir(reposResults)

indx <- which(regexpr(variable, filesRData, fixed=T)==1)
filesRData <- filesRData[indx]

for(file in filesRData) load(paste(reposResults,file, sep=""))

aux <- ls()
indx_meanIAvars <- grep("fmeanIAvars_vals", aux)


aux2 <- strsplit(aux[indx_meanIAvars], split="_")
variable <- unique(sapply(aux2, function(el) paste(el[3:(length(el)-3)], collapse="_")))
years <- unique(sapply(aux2, function(el) paste(el[(length(el)-2):(length(el)-1)], collapse="_")))

listRasts_meanIAvars <- lapply(indx_meanIAvars, function(i){
  res <- eval(parse(text=aux[i]))
  return(res)
  })
rastMeanIAvars <- do.call(merge, listRasts_meanIAvars)


reposResults <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/results_meanIAvar/"
meanIAvars_nm <- paste("fMeanIAvars", variable, years, sep="_")

save(rastMeanIAvars, file=paste(reposResults, meanIAvars_nm, ".RData", sep=""))


png(file=paste(paste(reposResults, meanIAvars_nm, sep=""), "png", sep="."))
plot(rastMeanIAvars)
dev.off()


#library(rgdal)

writeRaster(rastMeanIAvars, filename=paste(paste(reposResults, meanIAvars_nm, sep=""), "nc", sep="."), overwrite=TRUE, format="CDF")

##
# root moisture
load("/home/emiliano/Documents/ISP/proyectos/causality/CCM/results_meanIAvar/fMeanIAvars_root_moisture_2001_2011.RData")
root_moisture_rast <- rastMeanIAvars 
# surface moisture
load("/home/emiliano/Documents/ISP/proyectos/causality/CCM/results_meanIAvar/fMeanIAvars_surface_moisture_2001_2011.RData")
surface_moisture_rast <- rastMeanIAvars

hist(root_moisture_rast)
hist(log(root_moisture_rast))

plot(stack(root_moisture_rast, surface_moisture_rast))

plot(stack(clamp(root_moisture_rast, upper=0.01), clamp(surface_moisture_rast, upper=0.01)))

plot(stack(log(root_moisture_rast), log(surface_moisture_rast)))

rat_moisture_rast <- root_moisture_rast/surface_moisture_rast
hist(rat_moisture_rast)
plot(rat_moisture_rast)

plot(clamp(rat_moisture_rast, upper=1.1))
plot(rat_moisture_rast<1.1 & rat_moisture_rast>0.5)

plot(rat_moisture_rast < 0.5)

cutoff <- 0.0005
plot(stack(root_moisture_rast<cutoff, surface_moisture_rast<cutoff))
