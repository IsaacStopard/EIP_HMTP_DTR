# script to wrangle and interpolate all the temperature data
# Author: Isaac J Stopard
# Version: 1.0
# Notes: 

rm(list = ls())

library(deSolve); library(tidyverse); library(zipfR); library(GoFKernel); library(ggnewscale)
library(cowplot); library(suncalc); library(chillR); library(ggpmisc)
library(lubridate); library(zoo); library(rio); library(RColorBrewer);library(foreach); library(doParallel)
library(readxl)

########################
### temperature data ###
########################
f <- c("data/FW_EHT_data_analysis/19_09_to_26_10_20_Tengrela.xlsx",
       "data/FW_EHT_data_analysis/16_08_to_18_08_20_Tiefora.xlsx",
       "data/FW_EHT_data_analysis/16_08_to_18_08_20_Tengrela.xlsx", # max and min temp - actually NA
       "data/FW_EHT_data_analysis/19_09_to_26_10_20_Tiefora.xlsx",
       "data/FW_EHT_data_analysis/31_10_to_1_12_20_Tiefora.xlsx",
       "data/FW_EHT_data_analysis/31_10_to_09_11_20_Tengrela.xlsx")

lapply(f, function(file_){
  data <- rio::import_list(file = file_, setclass = "tbl_df", rbind = TRUE)
  print(c(paste0("Date NA:", sum(is.na(data$Date))),
          paste0("Temp NA:", sum(is.na(data$Temperature))),
          paste0("Max temp NA:", sum(is.na(data$`Maximum Temperature`))),
          paste0("Min temp NA:", sum(is.na(data$`Minimum Temperature`)))))
})

temp_data <- rbind(rio::import_list(file = "data/FW_EHT_data_analysis/16_08_to_18_08_20_Tiefora.xlsx", setclass = "tbl_df", rbind = TRUE),
                   rio::import_list(file = "data/FW_EHT_data_analysis/16_08_to_18_08_20_Tengrela.xlsx", setclass = "tbl_df", rbind = TRUE),
                   rio::import_list(file = "data/FW_EHT_data_analysis/19_09_to_26_10_20_Tengrela.xlsx", setclass = "tbl_df", rbind = TRUE),
                   rio::import_list(file = "data/FW_EHT_data_analysis/19_09_to_26_10_20_Tiefora.xlsx", setclass = "tbl_df", rbind = TRUE),
                   rio::import_list(file = "data/FW_EHT_data_analysis/31_10_to_1_12_20_Tiefora.xlsx", setclass = "tbl_df", rbind = TRUE),
                   rio::import_list(file = "data/FW_EHT_data_analysis/31_10_to_09_11_20_Tengrela.xlsx", setclass = "tbl_df", rbind = TRUE)) %>% 
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), # ignore daylight savings
         year = year(Date),
         month = month(Date),
         day = day(Date),
         hour = hour(Date),
         date = as.Date(Date),
         quarter = as.yearqtr(date, format =  "%Y-%m-%d %H:%M:%S"),
         f_date = floor_date(Date, unit = "hour"))

temp_data[which(temp_data$Location == "Intdoor"), "Location"] <- "Indoor"

# getting the start of the day and the start of night
temp_data[,"night_start_Ti"] <- suncalc::getSunlightTimes(date = temp_data$date, lat = 10.633333, lon = -4.55, keep = c("night"))$night
temp_data[,"night_start_Te"] <- suncalc::getSunlightTimes(date = temp_data$date, lat = 10.673024326627694, lon = -4.817461017395222, keep = c("night"))$night
temp_data[,"night_end_Ti"] <- suncalc::getSunlightTimes(date = temp_data$date, lat = 10.633333, lon = -4.55, keep = c("nightEnd"))$nightEnd
temp_data[,"night_end_Te"] <- suncalc::getSunlightTimes(date = temp_data$date, lat = 10.673024326627694, lon = -4.817461017395222, keep = c("nightEnd"))$nightEnd

temp_data <- temp_data %>% rowwise() %>% mutate(day_time = ifelse(Site == "Tiefora",
                                                                  ifelse(Date <= night_start_Ti & Date > night_end_Ti, "day", "night"),
                                                                  ifelse(Date <= night_start_Te & Date > night_end_Te, "day", "night")
))

hourly_temp_data <- rbind(temp_data %>% group_by(f_date, Site, Location, Hut) %>% summarise(temp = mean(Temperature),
                                                                                            max_temp = mean(`Maximum Temperature`),
                                                                                            min_temp = mean(`Minimum Temperature`)) %>% 
                            mutate(date = as.Date(f_date),
                                   hour = hour(f_date)), 
                          subset(temp_data, Location == "Indoor") %>% group_by(f_date, Site, Location) %>% summarise(temp = mean(Temperature),
                                                                                                                     max_temp = mean(`Maximum Temperature`),
                                                                                                                     min_temp = mean(`Minimum Temperature`)) %>% 
                            mutate(date = as.Date(f_date),
                                   hour = hour(f_date)) %>% mutate(Hut = "mean"), # calculated the hourly temperature data ignoring the house
                          subset(temp_data, Location == "Outdoor") %>% group_by(f_date, Site, Location) %>% summarise(temp = mean(Temperature),
                                                                                                                      max_temp = mean(`Maximum Temperature`),
                                                                                                                      min_temp = mean(`Minimum Temperature`)) %>% 
                            mutate(date = as.Date(f_date),
                                   hour = hour(f_date)) %>% mutate(Hut = "mean")) 

# ERA5 data
BF_data <- as.data.frame(read_csv(file = "data/ERA5/t2m_2020_BF.csv") %>% #filter(time >= min(hourly_temp_data$f_date) & time <= max(hourly_temp_data$f_date)) %>% 
                           mutate(temp = t2m -  273.5,
                                  min_temp = NA,
                                  max_temp = NA,
                                  f_date = time,
                                  Site = "Tiefora",
                                  Location = "Indoor",
                                  Hut = "ERA5",
                                  date = as.Date(time),
                                  hour = hour(time)))

hourly_temp_data <- rbind(as.data.frame(hourly_temp_data), BF_data[, colnames(hourly_temp_data)])

##############################################
### interpolating the missing temperatures ###
##############################################

u_t <- as.data.frame(unique(hourly_temp_data[,c("Site", "Location", "Hut")]))
u_t[,"index"] <- seq(1, nrow(u_t))

n_na <- rep(NA, nrow(u_t))
temp_data <- vector(mode = "list", length = nrow(u_t))
temp_data_fun <- vector(mode = "list", length = nrow(u_t))

for(i in 1:nrow(u_t)){
  
  data <- subset(hourly_temp_data, 
                 Location == u_t[i, "Location"] &
                   Site == u_t[i, "Site"] & 
                   Hut == u_t[i, "Hut"]) %>% 
    mutate(Year = year(f_date),
           Month = month(f_date),
           Day = day(f_date),
           Hour = hour,
           Temp = temp)
  
  ht <- as.data.frame(data[,c("Year", "Month", "Day", "Hour", "Temp")])
  
  c <- interpolate_gaps_hourly(hourtemps = ht, latitude = 10.633333)
  
  temp_data[[i]] <- as.data.frame(c$weather) %>% mutate(f_date = ymd_h(paste(Year, Month, Day, Hour, sep= ' ')),
                                                        Site = u_t[i, "Site"],
                                                        Location = u_t[i, "Location"],
                                                        Hut = u_t[i, "Hut"],
                                                        house = ifelse(Hut == "ERA5" | Hut == "mean", Hut, str_sub(u_t[i, "Hut"], start = -1)),
                                                        interpolate = ifelse(is.na(Temp_measured) == 1, "yes", "no"),
                                                        index = i)
  # ordering by the date
  temp_data[[i]] <- temp_data[[i]][order(temp_data[[i]]$f_date),]
  temp_data_fun[[i]] <- approxfun(x = seq(0, nrow(temp_data[[i]])-1)/24, y = temp_data[[i]]$Temp, yright = NA, yleft = NA)
  
  n_na[i] <- round(sum(is.na(temp_data[[i]]$Temp_measured))/nrow(temp_data[[i]])*100, digits = 2)
  
  rm(list = c("data", "c"))
}

saveRDS(list("temp_data" = temp_data, "n_na" = n_na, "temp_data_fun" = temp_data_fun, "hourly_temp_data" = hourly_temp_data, "u_t" = u_t), file = "data/temp_data.rds")
