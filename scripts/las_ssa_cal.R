library(raster)
library(lidR)
library(rlas)
library(dplyr)
library(readr)
library(terra)
library(sf)
library(data.table)

#LOAD ARGS
args <- commandArgs(trailingOnly = TRUE)
cal_las <- args[1]
crs <- args[2]
shp_fp_rfl <- args[3]
n_e_d_shift <- args[4]
output_csv <- args[5]
imu_data <- args[6]

# LOAD CAL LAS
crs <- as.numeric(crs)
las <- readLAS(cal_las)
st_crs(las) <- crs


# PREP SHIFT
n_e_d_shift <- strsplit(n_e_d_shift, ",")
n_shift <- as.numeric(n_e_d_shift[[1]][1])
e_shift <- as.numeric(n_e_d_shift[[1]][2])
d_shift <- as.numeric(n_e_d_shift[[1]][3])

# SHIFT X, Y, and Z based on ASP PC_ALIGN
las@data$X <- las@data$X  + e_shift
las@data$Y <- las@data$Y  + n_shift
las@data$Z <- las@data$Z  - d_shift

# CLIP TO TARGET
spdf <- st_read(file.path(shp_fp_rfl))
las <- clip_roi(las, spdf)

# COMPUTE SURFACE NORMALS AND INCLUDE IN LAS
dtm <- rasterize_terrain(las, algorithm = knnidw(k = 10L, p = 2))
slope <- terrain(dtm, v = c("slope"), unit = "radians")
aspect <- terrain(dtm, v = c("aspect"), unit = "radians")
x <- sin(aspect) * sin(slope)
y <- cos(aspect) * sin(slope) 
z <- cos(slope)
las <- merge_spatial(las, x, attribute = "n_i")
las <- merge_spatial(las, y, attribute = "n_j")
las <- merge_spatial(las, z, attribute = "n_k")

# GET AS DATATABLE
df <- payload(las)

# PLACEHOLDER FOR HELI IMU DATA
# I WILL ACTUALLY SYNC WITH GPS TIME HERE
# imu <- read_csv(imu_data)
# names(imu)[names(imu) == 'X?'] <- 'X_h'
# names(imu)[names(imu) == 'Y?'] <- 'Y_h'
# names(imu)[names(imu) == 'Z?'] <- 'Z_h'
# joined_df <- merge(df, imu, by = 'gpstime')
df$X_h<-df$X+1000 #INCLUDE SHIFT HERE
df$Y_h<-df$Y+1000 #INCLUDE SHIFT HERE
df$Z_h<-df$Z+15000 #INCLUDE SHIFT HERE


# COMPUTE RFL (NORMALIZED BY INC ANGLE)
df <- df %>%
    mutate( 
      cosi = ((X_h-X)*n_i + (Y_h-Y)*n_j + (Z_h-Z)*n_k) / (sqrt( (X_h-X)^2 + (Y_h-Y)^2 + (Z_h-Z)^2) * sqrt(n_i^2 + n_j^2 +n_k^2)),
      rfl = 10^(Reflectance/10) / cosi
)

# WRITE CSV
write.csv(df, output_csv)


