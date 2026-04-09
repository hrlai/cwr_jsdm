library(dismo)
library(terra)
library(usdm)


# env
env_select <- vifstep(as(bioclim_soilgrids, "Raster"), th = 10)
env_select <- env_select@results$Variables
env_select <- stringr::str_replace(env_select, "0.5", "0-5")
env_final  <- bioclim_soilgrids[[names(bioclim_soilgrids) %in% env_select]]

# Presences
focals <- c("Oryza sativa", 
            "Zea mays", 
            "Triticum aestivum",
            "Solanum tuberosum",
            "Glycine max")
cwr_of_focals <- 
  bind_rows(
    read_csv("data/inventory/maize_cwr.csv"),
    read_csv("data/inventory/rice_cwr.csv"),
    read_csv("data/inventory/soybean_cwr.csv"),
    read_csv("data/inventory/potato_cwr.csv"),
    read_csv("data/inventory/wheat_cwr.csv")
  ) %>% 
  mutate(Species = paste(Genus, Specie)) %>% 
  pull(Species) %>% 
  unique
cwr_select <- names(which(table(cwr_clean$species) >= 10))
cwr_select <- setdiff(cwr_select, focals)
cwr_select <- intersect(cwr_of_focals, cwr_select)



fit_maxent <- function(species) {
  
  message("Fitting ", species)
  
  # subset data
  test <- cwr_clean[which(cwr_clean$species == species), ]
  test <- vect(test, geom = c("decimallongitude", "decimallatitude"))
  
  # creates a 4-decimal-degree buffer around the
  # occurrence data
  occ_buff <- buffer(test, 4)
  
  # crop study area to a manageable extent (rectangle shaped)
  studyArea <- crop(env_final, ext(occ_buff))  
  
  # the 'study area' created by extracting the buffer area from the raster stack
  studyArea <- mask(studyArea, occ_buff)
  # output will still be a raster stack, just of the study area
  
  # select background points from this buffered area; when the number provided 
  # to set.seed() function, the same random sample will be selected in the next line			
  # use this code before the sampleRandom function every time, if you want to get
  # the same "random samples"
  set.seed(1) 
  bg <- spatSample(x = studyArea,
                   size = 10000,
                   method = "random",
                   na.rm = T, #removes the 'Not Applicable' points  
                   as.points = TRUE) # return spatial points 

  # extracting env conditions for training occ from the raster
  # stack; a data frame is returned (i.e multiple columns)
  p <- terra::extract(env_final, test, ID = FALSE)
  # extracting env conditions for background
  a <- terra::extract(env_final, bg, ID = FALSE)
  
  # repeat the number 1 as many numbers as the number of rows
  # in p, and repeat 0 as the rows of background points
  pa <- c(rep(1, nrow(p)), rep(0, nrow(a)))
  
  # (rep(1,nrow(p)) creating the number of rows as the p data
  # set to have the number '1' as the indicator for presence;
  # rep(0,nrow(a)) creating the number of rows as the a data
  # set to have the number '0' as the indicator for absence;
  # the c combines these ones and zeros into a new vector that
  # can be added to the Maxent table data frame with the
  # environmental attributes of the presence and absence
  # locations
  pder <- as.data.frame(rbind(p, a))

  # train Maxent with tabular data
  mod <- maxent(x = pder, ## env conditions
                p = pa,   ## 1:presence or 0:absence
                path = paste0("out/maxent/", species), 
                args = prepPara(replicates = 5, ## 5 replicates
                                replicatetype = "crossvalidate"))
  pred <- lapply(mod@models, predict, x = env_final)
  pred_avg <- app(rast(pred), mean)
  
  write_rds(mod, paste0("out/maxent/", species, "/model.rds"))
  writeRaster(pred_avg, paste0("out/maxent/", species, "/average_prediction.tif"))
  
}

library(doMC)
library(parallel)
library(foreach)
registerDoMC(cores = 6)
foreach(i = cwr_select) %dopar% {
  fit_maxent(species = i)
}


# eval <- evaluate(crds(test), crds(bg), m@models[[1]], x = env_final)

mod_files <- 
  list.files(path = "out/maxent", 
             pattern = "model.rds", 
             full.names = TRUE, recursive = TRUE)
mod_list <- lapply(mod_files, read_rds)
pred_list <- lapply(mod_list, function(mod) {
  pred_tmp <- lapply(mod@models, predict, x = env_final)
  pred_avg <- app(rast(pred_tmp), mean)
})
names(pred_list) <- 
  list.files(path = "out/maxent", 
             pattern = "model.rds", 
             recursive = TRUE) %>% 
  str_remove("/model.rds")

avg_pred_files <- 
  list.files(path = "out/maxent", 
             pattern = "average_prediction.tif", 
             full.names = TRUE, recursive = TRUE)
avg_pred_list <- lapply(avg_pred_files, rast)
# pred_sum <- app(rast(avg_pred_list), sum)
# plot(pred_sum, col = map.pal("viridis", 10))
names(avg_pred_list) <-
  list.files(path = "out/maxent", 
             pattern = "average_prediction.tif", 
             recursive = TRUE) %>% 
  str_remove("/average_prediction.tif")

calc_richness <- function(inventory) {
  cwr_tmp <- 
    read_csv(inventory) %>% 
    mutate(Species = paste(Genus, Specie)) %>% 
    pull(Species) %>% 
    unique
  cwr_tmp <- intersect(names(avg_pred_list), cwr_tmp)
  message("There are ", 
          length(cwr_tmp), 
          " wild relatives:\n", paste(cwr_tmp, collapse = "\n"))
  richness <- app(rast(avg_pred_list[cwr_tmp]), sum)
  return(richness)
}

maize_cwr_richness <- calc_richness("data/inventory/maize_cwr.csv")
rice_cwr_richness <- calc_richness("data/inventory/rice_cwr.csv")
soybean_cwr_richness <- calc_richness("data/inventory/soybean_cwr.csv")
potato_cwr_richness <- calc_richness("data/inventory/potato_cwr.csv")
wheat_cwr_richness <- calc_richness("data/inventory/wheat_cwr.csv")

cwr_richness <- c(maize_cwr_richness,
                  rice_cwr_richness,
                  soybean_cwr_richness,
                  potato_cwr_richness,
                  wheat_cwr_richness)
names(cwr_richness) <- c("Maize", "Rice", "Soybean", "Potato", "Wheat")

plot(cwr_richness, col = map.pal("viridis", 10))

writeRaster(cwr_richness, 
            "out/maxent/cwr_richness.tif",
            overwrite = TRUE)















# A function that implements Maxent parameters using the general R manner
prepPara <- function(userfeatures=NULL, #NULL=autofeature, could be any combination of # c("L", "Q", "H", "H", "P")
                     responsecurves=TRUE,
                     jackknife=TRUE,      
                     outputformat="logistic",
                     outputfiletype="asc", 
                     projectionlayers=NULL,
                     randomseed=FALSE,
                     removeduplicates=TRUE,
                     betamultiplier=NULL,
                     biasfile=NULL,
                     testsamplesfile=NULL,
                     replicates=1,
                     replicatetype="crossvalidate",
                     writeplotdata=TRUE,
                     extrapolate=TRUE,
                     doclamp=TRUE,
                     beta_threshold=NULL,
                     beta_categorical=NULL,
                     beta_lqp=NULL,
                     beta_hinge=NULL,
                     applythresholdrule=NULL
){
  #20 & 29-33 features, default is autofeature
  if(is.null(userfeatures)){
    args_out <- c("autofeature")
  } else {
    args_out <- c("noautofeature")
    if(grepl("L",userfeatures)) args_out <- c(args_out,"linear") else args_out <- c(args_out,"nolinear")
    if(grepl("Q",userfeatures)) args_out <- c(args_out,"quadratic") else args_out <- c(args_out,"noquadratic")
    if(grepl("H",userfeatures)) args_out <- c(args_out,"hinge") else args_out <- c(args_out,"nohinge")
    if(grepl("P",userfeatures)) args_out <- c(args_out,"product") else args_out <- c(args_out,"noproduct")
    if(grepl("T",userfeatures)) args_out <- c(args_out,"threshold") else args_out <- c(args_out,"nothreshold")
  }
  
  #1 
  if(responsecurves) args_out <- c(args_out,"responsecurves") else args_out <- c(args_out,"noresponsecurves")
  #2
  #if(picture) args_out <- c(args_out,"pictures") else args_out <- c(args_out,"nopictures")
  #3
  if(jackknife) args_out <- c(args_out,"jackknife") else args_out <- c(args_out,"nojackknife")
  #4
  args_out <- c(args_out,paste0("outputformat=",outputformat))
  #5
  args_out <- c(args_out,paste0("outputfiletype=",outputfiletype))
  #7
  if(!is.null(projectionlayers))    args_out <- c(args_out,paste0("projectionlayers=",projectionlayers))
  #10
  if(randomseed) args_out <- c(args_out,"randomseed") else args_out <- c(args_out,"norandomseed")
  #16
  if(removeduplicates) args_out <- c(args_out,"removeduplicates") else args_out <- c(args_out,"noremoveduplicates")
  #20 & 53-56
  # check if negative
  betas <- c( betamultiplier,beta_threshold,beta_categorical,beta_lqp,beta_hinge)
  if(! is.null(betas) ){
    for(i in 1:length(betas)){
      if(betas[i] <0) stop("betamultiplier has to be positive")
    }
  }
  if (  !is.null(betamultiplier)  ){
    args_out <- c(args_out,paste0("betamultiplier=",betamultiplier))
  } else {
    if(!is.null(beta_threshold)) args_out <- c(args_out,paste0("beta_threshold=",beta_threshold))
    if(!is.null(beta_categorical)) args_out <- c(args_out,paste0("beta_categorical=",beta_categorical))
    if(!is.null(beta_lqp)) args_out <- c(args_out,paste0("beta_lqp=",beta_lqp))
    if(!is.null(beta_hinge)) args_out <- c(args_out,paste0("beta_hinge=",beta_hinge))
  }
  #22
  if(!is.null(biasfile))    args_out <- c(args_out,paste0("biasfile=",biasfile))
  #23
  if(!is.null(testsamplesfile))    args_out <- c(args_out,paste0("testsamplesfile=",testsamplesfile))
  #24&25
  replicates <- as.integer(replicates)
  if(replicates>1 ){
    args_out <- c(args_out,
                  paste0("replicates=",replicates),
                  paste0("replicatetype=",replicatetype) )
  }
  #37
  if(writeplotdata) args_out <- c(args_out,"writeplotdata") else args_out <- c(args_out,"nowriteplotdata")
  #39
  if(extrapolate) args_out <- c(args_out,"extrapolate") else args_out <- c(args_out,"noextrapolate")
  #42
  if(doclamp) args_out <- c(args_out,"doclamp") else args_out <- c(args_out,"nodoclamp")
  #60
  if(!is.null(applythresholdrule))    args_out <- c(args_out,paste0("applythresholdrule=",applythresholdrule))
  
  return(args_out)
}
