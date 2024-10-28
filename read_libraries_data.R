# script to load the commonly used packages and data

library(deSolve); library(DescTools);
library(tidyverse); library(zipfR); library(GoFKernel); library(ggnewscale)
library(cowplot); library(suncalc); library(chillR); library(ggpmisc)
library(lubridate); library(zoo); library(rio); library(RColorBrewer);library(foreach); library(doParallel)
library(readxl); library(rstan); library(dfoptim); library(viridis); library(pROC);
library(doBy); library(shinystan); library(truncnorm);

source(file = "utils/data_functions.R")

###########
### EIP ###
###########
# getting the EIP values
fit <- readRDS(file = "data/fit_mSOS_temp_only_f2_f3.rds")
fit_ <- rstan::extract(fit)
# index 1 is the lowest temperature - 17
# index 11 is the highest temperature - 30
# single k value
temps <- c(17, 18, 19, 20, 21, 23, 25, 27, 28, 29, 30)
n_t <- length(temps)

# for all temperatures
temps_all <- seq(17, 30, 0.01)
mean_temp <- mean_temp_g <-  23.06936
sd_temp <- sd_temp_g <- 4.361642
scaled_temps_all <- (temps_all - mean_temp) / sd_temp # scaling so on same scale as parameter fits
scaled_temps <- (temps - mean_temp) / sd_temp
# getting the EIP value
params_temp <- rstan::extract(fit)

# for scaling the parameters to the relevant data scaling
mean_temp_s <- 27.9032
sd_temp_s <- 3.471223

# for all temperatures
temps_all <- seq(17, 30, 0.01)
scaled_temps_all <- (temps_all - mean_temp_g) / sd_temp_g # scaling so on same scale as parameter fits

# getting the EIP value

EIP_index <- get_EIP(params_temp, scaled_temps, 10000)
EIP_index_all <- get_EIP(params_temp, scaled_temps_all, 10000)

# calculating the mean EIP

# index <- seq(1, length(temps_all))
# mean_EIP <- as.data.frame(t(sapply(index, calc_mean_EIP, EIP_index = EIP_index_all)))
# mean_EIP$temp <- temps_all
#saveRDS(mean_EIP, file = "data/mean_EIP.rds")
mean_EIP <- readRDS(file = "data/mean_EIP.rds")

####################
### EIP function ###
####################
# assumes a linear interpolation for missing values
# separate EIP functions for the posterior quantiles - different mean values
EIP_fun <- vector(mode = "list", length = 3)
EIP_fun[[1]] <- approxfun(mean_EIP$temp, mean_EIP$`2.5%`, yleft = max(mean_EIP$`2.5%`), yright = min(mean_EIP$`2.5%`))
EIP_fun[[2]] <- approxfun(mean_EIP$temp, mean_EIP$`50%`, yleft = max(mean_EIP$`50%`), yright = min(mean_EIP$`50%`)) # extrapolating beyond
EIP_fun[[3]] <- approxfun(mean_EIP$temp, mean_EIP$`97.5%`, yleft = max(mean_EIP$`97.5%`), yright = min(mean_EIP$`97.5%`))