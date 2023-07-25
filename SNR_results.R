# look at SNR-var and SNR-ent for different CCM vars
# in map-world, map-south america and histogram form

remove(list=ls())
library(raster)
library(ncdf4)

server <- "optimus.uv.es"
user <- "emiliano"
variables <- c("root_moisture","soil_moisture","surface_moisture",
              "air_temperature_2m","gross_primary_productivity", 
              "radiation", "latent_energy", "precipitation")
reposResults <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/results_snr_ent/"

dir(reposResults)

# read in a raster

fileNm <- paste("fsnrEnt_", variables[1], ".nc", sep="")
stck <- stack(paste(reposResults, fileNm, sep=""))

names(stck) <- c("snr", "snr_sm", "h_ts_max", "h_ts_kdp", "h_ns_max", 
                 "h_ns_kdp", "h_ns_sm_max", "h_ns_sm_kdp")

plot(stck)

# get south america extent

plot(stck[["snr"]])
extentSA <- drawExtent()
plot(crop(stck[["snr"]], extentSA))

# read in all stacks, calculate snr-h with kdp and keep snr and snr-h only

stcks <- lapply(variables, function(var){
  # var <- variables[5]
  fileNm <- paste("fsnrEnt_", var, ".nc", sep="")
  stck <- stack(paste(reposResults, fileNm, sep=""))
  
  names(stck) <- c("snr", "snr_sm", "h_ts_max", "h_ts_kdp", "h_ns_max", 
                   "h_ns_kdp", "h_ns_sm_max", "h_ns_sm_kdp")
  snr_h <- 20*log(stck[["h_ts_kdp"]]/stck[["h_ns_kdp"]],10)
  
  plot(crop(stck[[c("h_ns_kdp","h_ns_max")]], extentSA))  
    
  #plot(crop(snr_h, extentSA)) 
  #res <- stack(stck[["snr"]], snr_h)
  #names(res) <- c("snr_v","snr_h")
  res <- stack(stck[["snr"]], stck[["h_ns_kdp"]])
  names(res) <- c("snr_v","h_ns")
  #plot(res)
  return(res)
})

names(stcks) <- variables

# plot SNR-var and SNR-ent for each variable separatley
# in a) map-world form, b)map-south america form c) histogram form

for(var in variables){
  # var <- variables[1]
  print("*************")
  print(var)
  # SNR-var
  plot(stcks[[var]][["snr_v"]], main=paste(var,"snrv"))
  plot(crop(stcks[[var]][["snr_v"]], extentSA), main=paste(var, "snrv"))
  hist(stcks[[var]][["snr_v"]], main=paste(var,"snrv"))
  
  # h_ns
  plot(stcks[[var]][["h_ns"]], main=paste(var,"h-ns"))
  plot(crop(stcks[[var]][["h_ns"]], extentSA), main=paste(var, "h-ns"))
  hist(stcks[[var]][["h_ns"]], main=paste(var,"h-ns"))
}


# get only snr in one stack

snrs <- stack(lapply(stcks, function(stck) stck[["snr_v"]]))

hs <- stack(lapply(stcks, function(stck) stck[["h_ns"]]))

plot(snrs[["precipitation"]]>-1)
plot(snrs[["root_moisture"]]>10 & snrs[["air_temperature_2m"]]>5 & 
       snrs[["gross_primary_productivity"]]>8 & snrs[["radiation"]]>3
     & snrs[["latent_energy"]]>10 & snrs[["precipitation"]]>-1)

highSNRmask <- snrs[["root_moisture"]]>10 & snrs[["air_temperature_2m"]]>5 & 
  snrs[["gross_primary_productivity"]]>8 & snrs[["radiation"]]>3 &
  snrs[["latent_energy"]]>10 & snrs[["precipitation"]]>-1

highSNRmask <- snrs[["root_moisture"]]>10 & snrs[["air_temperature_2m"]]>5 & 
  snrs[["gross_primary_productivity"]]>8  &
  snrs[["latent_energy"]]>10 

plot(highSNRmask)


extentAsia <- drawExtent()
extentAlaska <- drawExtent()
plot(highSNRmask)
plot(crop(highSNRmask, extentAsia))
plot(crop(highSNRmask, extentAlaska))
plot(crop(snrs[["gross_primary_productivity"]], extentAlaska))
plot(crop(snrs[["gross_primary_productivity"]], extentAsia))

plot(snrs[["gross_primary_productivity"]] < snrs[["radiation"]])
plot(crop(snrs[["gross_primary_productivity"]] < snrs[["radiation"]], extentSA))

plot(snrs[["precipitation"]] < snrs[["soil_moisture"]])
plot(crop(snrs[["precipitation"]] < snrs[["soil_moisture"]], extentSA))


reposCCM <- "/home/emiliano/Documents/ISP/proyectos/causality/CCM/ccmCorrs/"
V1 <- "precipitation"
V2  <- "surface_moisture"
fileCCM <- paste("ccm_corr_",V2, "_vs_",V1, "_2001_2011.nc", sep="") 
fileCCM %in% dir(reposCCM)
plot(stack(paste(reposCCM, fileCCM, sep="")))
plot(crop(stack(paste(reposCCM, fileCCM, sep="")), extentSA))

ccmStck <- stack(paste(reposCCM, fileCCM, sep=""))

plot(crop(ccmStck[[1]]<ccmStck[[2]], extentSA))

plot(crop(stack(snrs[["precipitation"]] < snrs[["surface_moisture"]], ccmStck[[1]]<ccmStck[[2]]), extentSA))


plot(crop(stack(hs[["gross_primary_productivity"]], 
                hs[["radiation"]]), extentSA))

plot(crop(stack(hs[["gross_primary_productivity"]]/hs[["radiation"]]), extentSA))
