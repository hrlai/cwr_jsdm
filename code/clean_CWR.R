library(tidyverse)
library(readxl)
library(rgbif)
library(countrycode)
library(CoordinateCleaner)
library(sp)
library(terra)
library(geodata)

# crop wild relative inventory from GRIN
# https://npgsweb.ars-grin.gov/gringlobal/taxon/taxonomysearchcwr
cwr_list <- 
    read_xlsx("data/GBIF/CWR/Crop_Relatives_GRIN.xlsx",
              skip = 1)

cwr_names <- sort(unique(word(cwr_list$`CROP WILD RELATIVE`, 1, 2)))

# crop wild relative occurrences from GBIF
# https://doi.org/10.15468/jyrthk
cwr <- 
    read_tsv("data/GBIF/CWR/0284616-220831081235567.csv") %>% 
    #convert country code from ISO2c to ISO3c
    mutate(countryCode = countrycode(countryCode, 
                                     origin = "iso2c", 
                                     destination = "iso3c")) %>% 
    # filter CWR to those related to Savary's focal crops
    filter(species %in% cwr_names)

# flag problems
flags <- 
    cwr %>% 
    filter(!is.na(decimalLongitude),
           !is.na(decimalLatitude)) %>% 
    select(species, 
           countryCode, 
           decimallongitude = decimalLongitude, 
           decimallatitude = decimalLatitude) %>% 
    clean_coordinates(countries = "countryCode",
                           # species = "species",
                           tests = c("capitals", 
                                     "centroids", 
                                     "equal",
                                     "gbif", 
                                     "institutions",
                                     "zeros",
                                     "seas",
                                     "urban",
                                     "countries")) 


# visualise
wm <- borders("world", colour="gray50", fill="gray50")
ggplot() + 
    coord_fixed() + 
    wm +
    geom_point(data = flags, 
               aes(x = decimallongitude, 
                   y = decimallatitude,
                   colour = .summary), 
               size = 0.5)+
    theme_bw()

# output clean data
# we will trust clean_coordinates for now and use only the unflagged data points
cwr_clean <- 
  flags %>% 
  filter(.summary == TRUE) %>% 
  select(species, countryCode, starts_with("decimal")) %>% 
  as.data.frame()
# convert to spatial points
cwr_points <- 
  vect(cwr_clean, 
       geom = c("decimallongitude", "decimallatitude"),
       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# read the bioclim raster to use as reference grid
bioclim <- 
  worldclim_global(
    var = "bio", 
    res = 10, 
    path = "data/WorldClim")
# rasterise the point occurrences to raster counts
point_counts <- 
  rasterize(cwr_points, bioclim,
            field = "species",  
            by = "species",
            fun = "length",
            background = 0)
# obtain adjacency list
W <- adjacent(point_counts, 
              cells = 1:prod(dim(point_counts)[1:2]), 
              directions = "rook",  # or queen?
              pairs = TRUE)
# or turn the adjacency (edge) list into a sparse matrix
# NB: not 100% sure if this is correct
ww <- igraph::graph_from_edgelist(W)
ww <- igraph::get.adjacency(ww)
