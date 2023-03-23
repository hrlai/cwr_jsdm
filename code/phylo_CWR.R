# Script to build phylogeny for crop wild relatives

library(tidyverse)
library(readxl)
library(V.PhyloMaker2)
# library(Taxonstand)
library(lcvplants)

# crop wild relative inventory from GRIN
# https://npgsweb.ars-grin.gov/gringlobal/taxon/taxonomysearchcwr
cwr_list <- 
    read_xlsx("data/GBIF/CWR/Crop_Relatives_GRIN.xlsx",
              skip = 1)

# phylogeny
# harmonise nomenclature
phylo_in <- 
    data.frame(species = unique(word(cwr_list$`CROP WILD RELATIVE`, 1, 2))) %>% 
    filter(!str_detect(species, "×")) %>% 
    pull(species) %>% 
    lcvp_search()

# mine phylogeny
phylo <- phylo.maker(
    phylo_in %>% 
        select(species = Search,
               genus   = Input.Genus,
               family  = Family),
    tree = GBOTB.extended.LCVP,
    nodes = nodes.info.1.LCVP)
write.tree(phylo$scenario.3, "data/phylo/cwr.tre")

plot(phylo$scenario.3,
     "radial",
     font = 1, cex = 0.5)
