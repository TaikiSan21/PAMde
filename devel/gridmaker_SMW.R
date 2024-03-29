# gridmaker.r --- creates uniform grid of decimal lat/longs based on 
#   input boundaries and pixel size (in degrees).
# Calculates area size of each pixel in square km

# Original code by Karin Forney and Elizabeth Becker
# Modified by Sam Woodman Nov 2018 so grid does not have to be rectangular.

# Note that there will be warning messages generated by some sf functions
#   regarding lat/long data


###############################################################################
# PREP
rm(list = ls())

# Run once to install needed packages
# install.packages(c("sf", "maps", "dplyr"))

library(sf)
library(maps)
library(dplyr)


###############################################################################
# VARIABLES SET BY USER

#------------------------------------------------------------------------------
### 1) Create map object for quick visualization of study area and grid

### EITHER ###
# Use this section if using longitude range [0, 360]
map.base <- st_as_sfc(maps::map('world2', plot = FALSE, fill = TRUE))

### OR ###
# Use this section if using longitude range [-180, 180]
# map.base <- st_as_sfc(maps::map('world', plot = FALSE, fill = TRUE))


#------------------------------------------------------------------------------
### 2) Define the study area (rectangular or otherwise)
# The centroids of the grid cells (NOT grid cell edges) will start 
#   at the bottom-left corner of the study area,
#   and thus proceed along the left and bottom edges of the study area.
# All grid cells whose centroids (NOT grid cell edges) are within 
#   the study area will be exported to the .csv file.
# Grid cells are not clipped at all by land

### EITHER ###
# Use this section to create a rectangular grid defined by min/max of lon/lat

# latmin <- 32.5
# latmax <- 42.0
# lonmin <- 360 - 125.0
# lonmax <- 360 - 118.0
# poly.df <- data.frame(
#   X1 = c(lonmax, lonmin, lonmin, lonmax, lonmax),
#   X2 = c(latmin, latmin, latmax, latmax, latmin)
# )

###  OR  ###
# Use this section to provide coordinates of a non-rectangular study area
# Coordinates must be entered as "c(lon, lat)" where all longitudes are in 
#   the same range (i.e. [0, 360] or [-180, 180]). 
# The first and the last set of coordinates must be the same.

list.vertices <- list(
  c(243, 32),
  c(238, 33),
  c(234, 38),
  c(234, 48.12),
  c(235.210364011, 48.5375173285),
  c(235.494102374, 48.3877105965),
  c(237, 48),
  c(237, 40),
  c(240, 36),
  c(243, 36),
  c(243, 32)
)
poly.df <- data.frame(do.call(rbind, list.vertices))


#------------------------------------------------------------------------------
### 3) Pixel size: the total length and width of each grid cell. 
# Thus, the shortest distance from a grid cell edge to 
#   the grid cell centroid will be (pixel / 2).

# pixel <- .225          # 25km
# pixel <- .090          # 10km
# pixel <- 0.10          # 0.1 degree for SeaGrant modeling project
# pixel <- 0.045         # 5km
pixel <- 0.027         # 3km
# pixel <- 0.018         # 2km


#------------------------------------------------------------------------------
### 4) Path and filename for output .csv file
# Output will be sorted by latitude and then by longitude

outfile <- paste0(
  "../whale-model-prep_data/Grid/", #Path to folder
  # paste0("Grid_Lat", latmin, "to", latmax, "_Lon", lonmin, "to", lonmax),
  "Grid_Nonrectangle", #Uncomment either this line or above line
  "_3km_WEAR.csv" # "_Step", pixel, "withArea.csv"
)


###############################################################################
# CODE RUN BY USER - NO CHANGES NEEDED

#------------------------------------------------------------------------------
# Create grid, etc.

### Create and visualize study area polygon; 4326 is WGS 84 coords
poly.bound <- st_sfc(st_polygon(list(as.matrix(poly.df))), crs = 4326)
plot(poly.bound, axes = TRUE, border = "red")
plot(map.base, add = TRUE, col = "tan")
# plot(poly.bound, add = TRUE, border = "red")

### Create grid (for areas) and get *centroids* of polygons.
grid <- st_make_grid(
  poly.bound, cellsize = pixel, 
  offset = st_bbox(poly.bound)[c("xmin", "ymin")] - pixel / 2
)
grid.cent <- st_set_precision(st_centroid(grid), 1e+10)

### Get the centroids that are within the study area
grid.cent.which <- unlist(st_intersects(poly.bound, grid.cent))

### Create data frame of centroid coordinates and area value
grid.cent.coords <- round(
  do.call(rbind, st_geometry(grid.cent)[grid.cent.which]), 
  10
)
lon <- grid.cent.coords[, 1]

grid.cent.df <- data.frame(
  lat = grid.cent.coords[, 2], 
  lon180 = ifelse(lon > 180, lon - 360, lon),
  lon360 = ifelse(lon < 0, lon + 360, lon),
  area_km = as.numeric(st_area(grid[grid.cent.which])) / 1e+06
) %>% 
  dplyr::arrange(lat, lon360) 


#------------------------------------------------------------------------------
# Save and visualize data

### Write data frame of centroid coordinates to .csv
write.csv(grid.cent.df, file = outfile, row.names = FALSE)

# ### Write grid to shapefile. Needs manual file output path
# grid.sa <- grid[grid.cent.which]
# st_write(grid.sa, "../whale-model-prep_data/shapefiles/grid_poly_027.shp")

# ### Visualize grid
# plot(grid[grid.cent.which], axes = TRUE)
# plot(poly.bound, add = TRUE, border = "red")
# plot(map.base, add = TRUE, col = "tan")

###############################################################################
# Exp
d <- eSDM::pts2poly_centroids(dplyr::select(grid.cent.df, lon180, lat, area_km), 0.027 / 2, crs = 4326) %>% 
  mutate(base_idx = 1:85869)
plot(st_geometry(d), border = "blue", add = TRUE)

x.land <- st_read("C:/SMW/Ensemble Case Study/GIS_Work/Shapefiles/World_countries_trunc.shp")
plot(x.land, add = TRUE, fill = "green")

d.e <- st_erase(d, x.land) %>% 
  mutate(area2_km = as.numeric(st_area(grid[grid.cent.which])) / 1e+06)
