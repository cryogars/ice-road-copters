library(raster)
library(lidR)
library(rlas)
library(dplyr)
library(readr)
library(terra)
library(sf)
library(data.table)

# LOAD ARGS
args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
crs <- args[2]
ni_fp <- args[3]
nj_fp <- args[4]
nk_fp <- args[5]
rfl_fp <- args[6]
n_e_d_shift <- args[7]
cosi_fp <- args[8]
road_cal_factor <- args[9]

# PATH TO SURFACE RETURNS AND SET CRS
crs <- as.numeric(crs)
road_cal_factor <- as.numeric(road_cal_factor)
pix_size <- as.numeric(pix_size)
las <- readLAS(f)
st_crs(las) <- crs

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

# PLACEHOLDER FOR HELI IMU DATA
# I WILL ACTUALLY SYNC WITH GPS TIME HERE
df$X_h <-df$X+1000
df$Y_h <-df$Y+1000 
df$Z_h <-df$Z + 15000 


# COMPUTE ALL RFL AND INCIDENCE ANGLES
df <- df %>%
    mutate(
      cosi = ((X_h-X)*n_i + (Y_h-Y)*n_j + (Z_h-Z)*n_k) / (sqrt( (X_h-X)^2 + (Y_h-Y)^2 + (Z_h-Z)^2) * sqrt(n_i^2 + n_j^2 +n_k^2)) , 
      rfl = 10^(Reflectance/10) * road_cal_factor
        )

# MAKE NEW TEMP LAS FILES THAT WILL HOLD COSI AND RFL AS THE Z VALUE
lasheader = header_create(las)
df$Z <- df$cosi # put the data into acceptable headers
write.las(cosi_fp, lasheader, df)

df$Z <- df$rfl
write.las(rfl_fp, lasheader, df)

# PREP SHIFT X, Y (ignore Z) ARGS BASED ON ASP PC_ALIGN
n_e_d_shift <- strsplit(n_e_d_shift, ",")
n_shift <- as.numeric(n_e_d_shift[[1]][1])
e_shift <- as.numeric(n_e_d_shift[[1]][2])

# APPLY SHIFT TO COSI
las <- readLAS(cosi_fp)
st_crs(las) <- crs
las@data$X <- las@data$X + e_shift # ASP shift applied here (E)
las@data$Y <- las@data$Y + n_shift # ASP shift applied here (N)
writeLAS(las, cosi_fp)

# APPLY SHIFT TO RFL
las <- readLAS(rfl_fp)
st_crs(las) <- crs
las@data$X <- las@data$X + e_shift # ASP shift applied here (E)
las@data$Y <- las@data$Y + n_shift # ASP shift applied here (N)
writeLAS(las, rfl_fp)
