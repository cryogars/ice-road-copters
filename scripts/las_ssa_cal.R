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

# GET AS DATATABLE
df <- payload(las)

# REMOVE LARGE SCAN ANGLES
df <- filter(df, abs(ScanAngleRank) <= 5)

# COMPUTE RFL
df <- df %>%
    mutate( 
      rfl = 10^(Reflectance/10)
)

# WRITE CSV
write.csv(df, output_csv)


