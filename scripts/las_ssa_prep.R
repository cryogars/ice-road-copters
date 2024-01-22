library(raster)
library(lidR)
library(rlas)
library(dplyr)
library(readr)
library(terra)
library(sf)
library(data.table)
options(digits = 22)

# LOAD ARGS
args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
crs <- args[2]
ni_fp <- args[3]
nj_fp <- args[4]
nk_fp <- args[5]
rfl_fp <- args[6]
n_e_d_shift <- args[7]
road_cal_factor <- args[8]
imu_data <- args[9]
alpha <- args[10]

# PATH TO SURFACE RETURNS AND SET CRS
crs <- as.numeric(crs)
road_cal_factor <- as.numeric(road_cal_factor)
alpha <- as.numeric(alpha)

las <- readLAS(f)
st_crs(las) <- crs

# PREP SHIFT X, Y (ignore Z) ARGS BASED ON ASP PC_ALIGN
n_e_d_shift <- strsplit(n_e_d_shift, ",")
n_shift <- as.numeric(n_e_d_shift[[1]][1])
e_shift <- as.numeric(n_e_d_shift[[1]][2])
d_shift <- as.numeric(n_e_d_shift[[1]][3])

# APPLY SHIFT
las@data$X <- las@data$X + e_shift # ASP shift applied here (E)
las@data$Y <- las@data$Y + n_shift # ASP shift applied here (N)
las@data$Z <- las@data$Z - d_shift # ASP shift applied here (D)

# MERGE INFORMATION INTO LAS
x <- raster(ni_fp)
y <- raster(nj_fp)
z <- raster(nk_fp)
las <- merge_spatial(las, x, attribute = "n_i")
las <- merge_spatial(las, y, attribute = "n_j")
las <- merge_spatial(las, z, attribute = "n_k")

# CONVERT TO DF. FILTER OUT NA AND SELECT ONLY FIRST RETURNS.
df <- payload(las)
df<- filter(df, NumberOfReturns == 1)
df<- filter(df, ReturnNumber == 1)
df<- filter(df, n_i != "NA")

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

# APPLY ASP SHIFT TO HELI POSITION
df$X_h <- df$X_h + e_shift # ASP shift applied here (E)
df$Y_h <- df$Y_h + n_shift # ASP shift applied here (N)
df$Z_h <- df$Z_h - d_shift # ASP shift applied here (D)

# COMPUTE RFL (NORMALIZED BY INC ANGLE)
df <- df %>%
    mutate( 
      range = sqrt((X_h-X)^2 + (Y_h-Y)^2 + (Z_h-Z)^2),
      trans = exp(alpha * range),
      cosi = ((X_h-X)*n_i + (Y_h-Y)*n_j + (Z_h-Z)*n_k) / (range * sqrt(n_i^2 + n_j^2 +n_k^2)),
      rfl = (10^(Reflectance/10) / cosi) * (road_cal_factor  / trans**2)
)

# REMOVE COSI < 0.50 (EMPIRICAL/SENS)
df <- filter(df, cosi>=0.50)

# REMOVE RFL> 1 or less than 0
df <- filter(df, rfl>=0.0)
df <- filter(df, rfl<=1.0)

# MAKE NEW TEMP LAS FILE THAT HOLD RFL AS THE Z VALUE
lasheader = header_create(las)
df$Z <- df$rfl
write.las(rfl_fp, lasheader, df)