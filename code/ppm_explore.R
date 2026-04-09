library(inlabru)
library(INLA)
library(sf)

bru_safe_sp(force = TRUE)


# SpatialPolygons object for the boundary for the mesh

flat_env <- as.numeric(!is.na(bioclim_soilgrids[[26]]))
flat_env[flat_env == 0] <- NA
world <- 
  as.polygons(flat_env,
              values = FALSE, 
              extent = FALSE) %>% 
  as(., "Spatial")
plot(world)

world_mesh_segment <- as.inla.mesh.segment(world)
cwr_coord <- 
  st_as_sf(cwr_clean, coords = c("decimallongitude", "decimallatitude"))[,1] %>% 
  st_coordinates
max.edge <- diff(range(cwr_coord))/(3*5)

mesh <-
  inla.mesh.2d(
    loc = cwr_coord,
    boundary = world_mesh_segment,
    # cutoff = 0.2,
    max.edge = max.edge,
    # offset = c(0.1, 0.4)
  )

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

# Specify model

species <- cwr_select[1]
pres_coord <- 
  cwr_coord[which(cwr_clean$species == species), ] %>% 
  as.data.frame() %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  as(., "Spatial")


# Compile data
library(rnaturalearth)
map <- ne_countries(type = "countries", scale = "medium")
resolution <- 0.2
r <- raster(map, resolution = resolution)
r[] <- 0
tab <- table(cellFromXY(r, pres_coord))
r[as.numeric(names(tab))] <- tab
