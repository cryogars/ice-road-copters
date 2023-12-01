install.packages('raster')
install.packages('lidR')
install.packages('dplyr')
install.packages('readr')
install.packages('terra')
library(raster)
library(lidR)
library(dplyr)
library(readr)
library(terra)


# PATH TO SURFACE RETURNS AND SET CRS
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

# SHIFT X, Y, and Z based on ASP PC_ALIGN
n_e_d_shift <- strsplit(n_e_d_shift, ",")
n_shift <- as.numeric(n_e_d_shift[[1]])
e_shift <- as.numeric(n_e_d_shift[[2]])
d_shift <- as.numeric(n_e_d_shift[[3]])
df$X <- df$X + e_shift
df$Y <- df$Y + n_shift
df$Z <- df$Z - d_shift


# CONVERT TO RASTER /  GRID HERE
# file paths:
# rfl_fp and cosi_fp
