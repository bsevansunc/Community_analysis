library(raster)

getwd()

setwd('/home/bsevans/spatial_data')

#---------------------------------------------------------------------------------------------*
# ---- Point data ----
#=============================================================================================*

# Get data:

pts = read.csv('spatial_data/nn_pts.csv')

# Set projection information:

proj1 = '+proj=utm +zone=18 +datum=WGS84 +no_defs +ellps=WGS84  +units=m'
proj2 = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '

# Convert point data to spatial points (data = site):

pts = SpatialPointsDataFrame(pts[,2:3], data = data.frame(pts[,1]), proj4string=CRS(proj1))
  names(pts) = 'site'

# Change projection to lat lon:

pts.ll = spTransform(pts, CRS(proj2))

#---------------------------------------------------------------------------------------------*
# ---- Raster data ----
#=============================================================================================*

# Create an extent object from the pts.ll file and expand by 1/10 degree:

e = extent(pts.ll)+c(-.1, .1, -.1, .1)

# Canopy cover:

can1 = raster('spatial_data/Hansen_GFC2013_treecover2000_40N_080W.tif')

# Crop to extent:

can = crop(can1, e)

# Impervious surface:

imp1 = raster('spatial_data/imp50')

# Change projection:

imp = projectRaster(imp1, can)

# Crop to study region:

imp = crop(imp, e)

#---------------------------------------------------------------------------------------------*
# ---- Extract LC data to points ----
#=============================================================================================*

# To ensure that distances are in meters, convert rasters to (or back to) UTM:

can.utm = projectRaster(can, imp1)

imp.utm = projectRaster(imp, imp1)

# I need to get rid of the water in can.utm (set as NA using the NA's in imp):

can.utm = mask(can.utm, imp.utm)

# Stack the rasters:

lc.stack = stack(can.utm,imp.utm)
  names(lc.stack) = c('can','imp')

# pts$imp = pts.imp[,2]

# Extract the raster stack to the points:

pts.lc = extract(lc.stack, pts, buffer = 100, fun = mean, na.rm = T)

# Bind with the points file:

pts.lc = cbind(data.frame(pts), pts.lc)

# Write to file:

write.csv(pts.lc, 'pts_lc100.csv', row.names = F)









