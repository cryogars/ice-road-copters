library(raster)
library(lidR)
library(rlas)
library(dplyr)
library(readr)
library(terra)
library(sf)
library(data.table)
options(digits = 22)

#LOAD ARGS
args <- commandArgs(trailingOnly = TRUE)
cal_las <- args[1]
crs <- args[2]
shp_fp_rfl <- args[3]
n_e_d_shift <- args[4]
output_csv <- args[5]
imu_data <- args[6]
pix_size <- args[7]
alpha <- args[8]

# LOAD CAL LAS
crs <- as.numeric(crs)
pix_size <- as.numeric(pix_size)
alpha <- as.numeric(alpha)

las <- readLAS(cal_las)
st_crs(las) <- crs

# PREP SHIFT
n_e_d_shift <- strsplit(n_e_d_shift, ",")
n_shift <- as.numeric(n_e_d_shift[[1]][1])
e_shift <- as.numeric(n_e_d_shift[[1]][2])
d_shift <- as.numeric(n_e_d_shift[[1]][3])

# SHIFT X, Y, and Z based on ASP PC_ALIGN
# NOTE: TRANSLATION IS APPLIED FIRST HERE IN ORDER TO CLIP TO ROI IN REAL WORLD COORD.
las@data$X <- las@data$X  + e_shift
las@data$Y <- las@data$Y  + n_shift
las@data$Z <- las@data$Z  - d_shift

# CLIP TO TARGET
spdf <- st_read(file.path(shp_fp_rfl))
las <- clip_roi(las, spdf)

# HARD SET TO TELL RLIDR THESE ARE GROUND POINTS
las@data$Classification <- 2

# COMPUTE SURFACE NORMALS AND INCLUDE IN LAS
dtm <- rasterize_terrain(las, algorithm = knnidw(k = 10L, p = 2), res=pix_size)
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

# HELI IMU DATA
imu <- read_csv(imu_data,show_col_types = FALSE)
imu <- as.data.table(imu, TRUE) 

#REDUCE THE NUMBER OF ROWS FOR JOIN
max_time <- max(df$gpstime, na.rm = TRUE)
min_time <- min(df$gpstime, na.rm = TRUE)

# SET COL NAMES FOR JOIN
names(imu)[names(imu) == 'Easting[m]'] <- 'X_h'
names(imu)[names(imu) == 'Northing[m]'] <- 'Y_h'
names(imu)[names(imu) == 'Height[m]'] <- 'Z_h'
names(imu)[names(imu) == 'Time[s]'] <- 'gpstime'

# APPLY FILTER TO INCREASE SPEED
imu <- filter(imu, gpstime >= min_time)
imu <- filter(imu, gpstime <= max_time)

# SET KEYS FOR JOIN
setkey(df,gpstime)
setkey(imu,gpstime)

# JOIN BY NEAREST
df <- imu[df,roll = "nearest"]

# APPLY ASP PC-ALIGN TRANSLATION
df$X_h <- df$X + e_shift 
df$Y_h <- df$Y + n_shift 
df$Z_h <- df$Z - d_shift 

# COMPUTE RFL (NORMALIZED BY INC ANGLE)
df <- df %>%
    mutate( 
      range = sqrt((X_h-X)^2 + (Y_h-Y)^2 + (Z_h-Z)^2),
      trans = exp(alpha * range),
      cosi = ((X_h-X)*n_i + (Y_h-Y)*n_j + (Z_h-Z)*n_k) / (range * sqrt(n_i^2 + n_j^2 +n_k^2)),
      rfl = (10^(Reflectance/10) / cosi) * (1 / trans**2)
)

# WRITE CSV
write.csv(df, output_csv)


