# read in all .RData files containing maps of consecutive data points and time embedding for a certain variable
# and make a raster file of it
remove(list=ls())
library(raster)
library(ncdf4)

server <- "optimus.uv.es"
user <- "emiliano"


reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses_byBlock/", sep="")
#reposResults <- paste("timeEmbeddingMaps_byBlock/", paste(V1, "vs", V2, sep="_"), "/", sep="")

aux <- dir(reposResults)
sapply(aux, function(combo) length(dir(paste(reposResults, combo, sep=""))))
table(sapply(aux, function(combo) length(dir(paste(reposResults, combo, sep="")))))

variables <- c("soil_moisture","surface_moisture",
               "air_temperature_2m","gross_primary_productivity", 
               "radiation", "latent_energy", "precipitation")

combos <- as.data.frame(t(combn(variables, 2)))
colnames(combos) <- c("V1","V2")


for(i in 1:nrow(combos)){
  #i <- 1
  V1 <- combos$V1[i]
  V2 <- combos$V2[i]
  print("**************************")
  print(i)
  print(paste(V1, "vs", V2, sep="_"))

  reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses_byBlock/", sep="")
  filesRData <- dir(paste(reposResults, paste(V1, "vs",V2,sep="_"), sep="/"))
  for(file in filesRData) load(paste(reposResults, paste(V1, "vs",V2,sep="_"), "/", file, sep=""))

  aux <- ls()
  indx_cause <- grep( paste("ccm_cause_vals", V1, V2, sep="_"), aux)

  

  aux2 <- strsplit(aux[indx_cause], split="_")
  years <- unique(sapply(aux2, function(el) paste(el[(length(el)-2):(length(el)-1)], collapse="_")))

  listRastsCause <- lapply(indx_cause, function(i){
    res <- eval(parse(text=aux[i]))
    return(res)
  })
  rastCause <- do.call(merge, listRastsCause)

  plot(rastCause)
  #plot(getValues(rastCause[[1]]),getValues(rastCause[[2]])); abline(a=0, b=1, col="red")
  #plot(stack(is.na(rastCause[[1]]),rastCause[[1]]>rastCause[[2]], rastCause[[3]]<rastCause[[4]]))

  #reposResults <- paste("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/ODEs/CLEAN/code/ccm_rEDM/ccmCauses/", sep="")
  reposResults <- "/home/soulivanh/Documents/proyectos/CCM/ccmCauses/"
  cause_nm <- paste("ccm_cause", paste(V1, "vs", V2, sep="_"), years, sep="_")


  png(file=paste(paste(reposResults, cause_nm, sep=""), "png", sep="."))
  plot(rastCause)
  dev.off()

  #library(rgdal)

  reposResults <- "/home/soulivanh/Documents/proyectos/CCM/ccmCauses/"
  writeRaster(rastCause, filename=paste(paste(reposResults, cause_nm, sep=""),"nc", sep="."), overwrite=T, format="CDF")
}


reposResults <- "/home/soulivanh/Documents/proyectos/CCM/ccmCauses/"
files <- dir(reposResults)
indx <- which(sapply(strsplit(files, "\\."), function(el) el[2])=="nc")
files <- files[indx]

rast <- stack( paste(reposResults, files, sep=""))
plot(rast)
dim(rast)
21*12

