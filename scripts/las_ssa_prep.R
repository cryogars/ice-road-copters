library(raster)
library(lidR)
library(dplyr)
library(readr)
library(terra)


# WRAP THIS ALL INTO A FUNCTION
compute_rfl_and_angles <- function(las_path, road_cal_factor, crs){

    # PATH TO SURFACE RETURNS
    las <- readLAS(las_path)
    st_crs(las) <- crs

    # MERGE INFORMATION
    # all the paths to the rasters to be made... TODO: 
    x <- raster("./data/fl_230209_225728/ground_data/n_i.tif")
    y <- raster("./data/fl_230209_225728/ground_data/n_j.tif")
    z <- raster("./data/fl_230209_225728/ground_data/n_k.tif")
    slope <- raster("./data/fl_230209_225728/ground_data/slope.tif")
    las <- merge_spatial(las, x, attribute = "n_i")
    las <- merge_spatial(las, y, attribute = "n_j")
    las <- merge_spatial(las, z, attribute = "n_k")
    las <- merge_spatial(las, slope, attribute = "slope")

    # Payload - Convert to dataframe and clean bad data
    df <- payload(las)
    df<- filter(df, NumberOfReturns == 1)
    df<- filter(df, ReturnNumber == 1)
    df<- filter(df, n_i != "NA")

    # PLACEHOLDER FOR HELI IMU DATA
    # I WILL ACTUALLY SYNC WITH GPS TIME HERE
    df$X_h <-df$X+1000
    df$Y_h <-df$Y+1000 
    df$Z_h <-df$Z + 15000 

    # Compute cosi and rfl for all the points
    df <- df %>%
        mutate(
          cosi = ((X_h-X)*n_i + (Y_h-Y)*n_j + (Z_h-Z)*n_k) / (sqrt( (X_h-X)^2 + (Y_h-Y)^2 + (Z_h-Z)^2) * sqrt(n_i^2 + n_j^2 +n_k^2)) , 
          rfl = 10^(Reflectance/10) * road_cal_factor
            )

    # CONVERT TO RASTER /  GRID HERE

}



####
# NEXT FUNCTION
# Uses the buffer shape to clip only road points and where abs(ScanAngleRank) < 5
# Compute median value as the calibration factor. 
# Return csv of all points for computing histogram later
#las = readLAS("file.las")
#spdf <- readOGR(dsn = "...", layer = "...")
#clipped_las = lasclip(ctg, spdf)



####
# NEXT FUNCTION
# Load in rasters and compute Slope , Aspect, and normal vector
# Save to disk

