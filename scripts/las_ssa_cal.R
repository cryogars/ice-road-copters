library(raster)
library(lidR)
library(dplyr)
library(readr)
library(terra)
library(sf)


args <- commandArgs(trailingOnly = TRUE)
cal_las <- args[1]
crs <- args[2]
shp_fp_rfl <- args[3]
n_e_d_shift <- args[4]
output_csv <- args[5]

# Uses the buffer shape to clip only road points and where abs(ScanAngleRank) < 5
# Return csv of all points for computing histogram and median

# MY SHAPE: /Users/brent/Code/ice-road-copters/transform_area/Eagle/eagle_res_buffered.shp
# layer = eagle_res_buffered
# dsn = /Users/brent/Code/ice-road-copters/transform_area/Eagle/

# Assumes this pandas df is reasonable size and can be saved
crs <- as.numeric(crs)
las <- readLAS(cal_las)
st_crs(las) <- crs

n_e_d_shift <- strsplit(n_e_d_shift, ",")
n_shift <- as.numeric(n_e_d_shift[[1]][1])
e_shift <- as.numeric(n_e_d_shift[[1]][2])
d_shift <- as.numeric(n_e_d_shift[[1]][3])

# SHIFT X, Y, and Z based on ASP PC_ALIGN
las@data$X <- las@data$X  + e_shift
las@data$Y <- las@data$Y  + n_shift
las@data$Z <- las@data$Z  - d_shift

spdf <- st_read(file.path(shp_fp_rfl))

las <- clip_roi(las, spdf)

df <- payload(las)

df <- filter(df, abs(ScanAngleRank) <= 5)

# COMPUTE RFL
df <- df %>%
    mutate( 
      rfl = 10^(Reflectance/10)
        )

write.csv(df, output_csv)


