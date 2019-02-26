require(rlas)
require(lidR)
require(raster)
require(matlab) # for tic()/toc() function

#Note this will alter the original data

# MN lidar
#data <- '/Users/aaron/gdrive/projects/big_data_lidar_project/data/points.laz'
data <- '/Users/aaron/Desktop/thin_lidar_project/duncan_lake/1158-30-51.laz'

# Create a spatial index file (.lax) that sits next to the .las file during processing to speed things up
writelax(data)

# Read las; apply filters here if needed (see rlas:::lasfilterusage())
# Filtering out overlapping points (The MN LiDAR dataset has these Classified as 17)
las <- readLAS(data, filter="-drop_class 9 17")  # Water = 9 Overlapping = 17
las %>% grid_density %>% plot #Check the results

# Get basic info about las file
summary(las)
str(las)

# Normalize data (try reading las file in using only ground points readlas())
dtm <- grid_terrain(las, method = "knnidw", k = 10L)
lasnormalize(las=las, dtm=dtm)

# Check if normalization worked (vals should be in 0-30 range)
las@data$Z

# Look for outliers (e.g. birds)
bp <- boxplot(las@data$Z)
bp
bp$out # These are values beyond the extremes of the whiskers (the real outliers)

# TODO: Find out if I can use class 7 (Low Point (“low noise”)) and 8 (Model Keypoint (possibly high noise)) to help screen noise points
# https://community.esri.com/blogs/larry-zhang/2014/10/29/noise-removal-and-manual-classification-in-las-cloud-points
# Remove 0.15 > points < 40 (or, if you want DBH=1.37)
cleaned_las <- lasfilter(las, Z < 40 & Z >= 0.15)

# Generate a list of Grid Metrics
# https://github.com/Jean-Romain/lidR/wiki/grid_metrics

tic()
myMetrics = function(x,y,z, i)
{
  xmax = max(x)    # This is used to clean edge cells 
  ymax = max(y)    # https://github.com/Jean-Romain/lidR/wiki/grid_metrics-control
  xmin = min(x)
  ymin = min(y)
  
  lidrmetrics = stdmetrics_z(z)
  
  metrics = list(
    A       = (xmax - xmin)*(ymax - ymin),      # This is used to clean edge cells 
    zmean   = mean(z),
    imean   = mean(i),
    zwimean = sum(z*i)/sum(i),                  # Mean elevation weighted by intensities
    zimean  = mean(z*i),                        # Mean products of z by intensity
    zsqmean = sqrt(mean(z^2)),                  # Quadratic mean
    q01     = quantile(z, probs = c(0.01)),     # 1st percentile value for cell
    q05     = quantile(z, probs = c(0.05)),     # 5th percentile value for cell
    q10     = quantile(z, probs = c(0.1)),      # 10th percentile value for cell
    q20     = quantile(z, probs = c(0.2)),      # 20th percentile value for cell
    q30     = quantile(z, probs = c(0.3)),      # 30th percentile value for cell
    q40     = quantile(z, probs = c(0.4)),      # 40th percentile value for cell
    q50     = quantile(z, probs = c(0.5)),      # 50th percentile value for cell
    q60     = quantile(z, probs = c(0.6)),      # 60th percentile value for cell
    q70     = quantile(z, probs = c(0.7)),      # 70th percentile value for cell
    q75     = quantile(z, probs = c(0.75)),     # 75th percentile value for cell
    q80     = quantile(z, probs = c(0.8)),      # 80th percentile value for cell
    q90     = quantile(z, probs = c(0.9)),      # 90th percentile value for cell
    q95     = quantile(z, probs = c(0.95)),     # 95th percentile value for cell
    q99     = quantile(z, probs = c(0.99)),     # 99th percentile value for cell
    entropy = entropy(z),                       # TODO: Note the data stripe on bottom. Find a fix
    vci     = VCI(z, zmax = 40)
  )
  
  return(c(metrics, lidrmetrics))
}

# Run the grid metrics
metrics <- grid_metrics(cleaned_las, myMetrics(X,Y,Z, Intensity), 10) #spatial resolution = 10 here
toc()
plot(metrics, "vci") # Just checking...

# Clean the metrics so that there are no edge effects
# https://github.com/Jean-Romain/lidR/wiki/grid_metrics-control
# TODO: Need to figure out why there are NoData pixels being assigned throughout raster
#cleaned_metrics <- metrics[A > 100 * 0.90]
#plot(cleaned_metrics, "vci") # Now see if the cleaning worked

# Convert to raster (raster stack if there area several layers)
raster <- as.raster(metrics, z = "vci")
#cleaned_raster <- as.raster(cleaned_metrics, z = "vci")

# Convert NA to 0's
raster[is.na(raster[])] <- 0  #TODO: does it matter if the collar is 0 too?

# Assign spatial reference 
crs(raster) <- las@crs@projargs # Pulling the CRS from the las header

# Write to 8bit unsigned integer (Note may want to use floating point for some metrics???)
writeRaster(raster, file.path('~/Desktop/', "vcitest"), datatype = "FLT4S", format = "HFA") # format = "HFA" for .img
writeRaster(cleaned_raster, file.path('~/Desktop/', "vci_cleaned"), datatype = "FLT4S", format = "HFA") # format = "HFA" for .img




