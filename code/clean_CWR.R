library(tidyverse)
library(readxl)
library(rgbif)
library(countrycode)
library(CoordinateCleaner)
library(sp)

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
