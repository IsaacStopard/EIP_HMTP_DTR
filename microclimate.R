rm(list = ls())

source(file = "utils/data_functions.R"); source(file = "utils/model_functions.R"); source(file = "read_libraries_data.R")

######################################################
##### mosquito biting time from the Burkina Faso #####
######################################################
source(file = "read_bt_data.R")

png(file = "results/figures/bite_time_dens.png", height = 400, width = 700)
ggplot(data = rbind(bt_density_outdoor %>% mutate(location = "Outdoor"),
                    bt_density_indoor %>% mutate(location = "Indoor")), 
       aes(x = Hour, y = d)) +
  ylab("Probability mass") +
  geom_bar(stat = "identity", alpha = 0.75, position = "dodge", fill = "grey70",
           col = "grey25") + 
  theme_bw() + theme(text = element_text(size = 18),
                     axis.text.x = element_text(angle = -90)) +
  facet_wrap(~location) +
  scale_y_continuous(limits = c(0,0.25))
dev.off()

##### temperature data #####
temp_data_rds <- readRDS(file = "data/temp_data.rds")
n_na <- temp_data_rds$n_na
temp_fun <- temp_data_rds$temp_data_fun
temp_data <- temp_data_rds$temp_data
hourly_temp_data <- temp_data_rds$hourly_temp_data
u_t <- temp_data_rds$u_t

temp_data_plot <- bind_rows(temp_data) #

temp_data_plot <- rbind(temp_data_plot, subset(temp_data_plot, house == "ERA5") %>% mutate(Location = "Outdoor"))

DTR_data <- temp_data_plot %>% mutate(s_date = as.Date(f_date)) %>% group_by(s_date, Site, Location, Hut, house, Month) %>% 
  summarise(DTR = max(Temp) - min(Temp)) %>% mutate(Location = ifelse(house == "ERA5", "ERA5", Location))

DTR_data$Location <- factor(DTR_data$Location, levels = c("ERA5", "Outdoor", "Indoor"))

png(file = "results/figures/DTR_micro_plot.png", height = 450, width = 700)
ggplot(data = DTR_data %>% subset(Site == "Tiefora" & Month > 8),
       aes(x = house, y = DTR, fill = Location, group = interaction(house, Location))) +
  geom_boxplot(alpha = 0.75) +
  theme_bw() + theme(text = element_text(size = 18)) +
  xlab("House") +
  ylab("Diurnal temperature range (DTR) (°C)") +
  scale_fill_manual(values = c("#56B4E9","#CC79A7",  "#E69F00"), name = "Temperature\ndata") +
  scale_y_continuous(breaks = seq(0, 35, 5))
dev.off()  

t_plot <- 
  plot_grid(
    ggplot(data = subset(temp_data_plot, Site == "Tiefora" & Day == 1 & Month > 8) %>% 
             mutate(date_plot = paste0("Date: ",format(as.Date(f_date, format = "%d/%m/%Y"), "%d-%b-%Y"))),
           aes(x = Hour, y = Temp, col = factor(house))) +
      geom_line(linewidth = 1.5) +
      geom_point(data = subset(temp_data_plot, Site == "Tiefora" & interpolate == "yes" & Day == 1), 
                 aes(x = Hour, y = Temp), col = "#E69F00", size = 0.25) +
      facet_grid( vars(Location), vars(date_plot)) +
      theme_bw() + theme(text = element_text(size = 18)) +
      scale_colour_viridis_d(name = "house") +
      xlab("Time (hour)") +
      ylab("Temperature (°C)"),
    
    ggplot(data = subset(temp_data_plot, Site == "Tiefora"& Month > 8) %>% 
             mutate(Location = ifelse(house=="ERA5","ERA5",Location)) ,
           aes(x = f_date, y = Temp)) +
      geom_line(linewidth = 0.5, alpha = 0.75, col = "#56B4E9") +
      geom_point(data = subset(temp_data_plot, Site == "Tiefora" & interpolate == "yes"), 
                 aes(x = f_date, y = Temp), col =  "#E69F00", size = 0.1, alpha = 0.75) +
      facet_wrap(~Location + house) +
      theme_bw() + theme(text = element_text(size = 18)) +
      #geom_smooth(formula = y ~ x, se = FALSE, method = "loess", col = "grey50", alpha = 0.25, size = 0.5, linetype = 2) +
      xlab("Date and time (hour - day - month") +
      ylab("Temperature (°C)") +
      scale_x_datetime(date_labels = "%d %b",
                   date_breaks = "2 month"),
    nrow = 2, labels = c("A", "B"), rel_heights = c(1, 1.5)
  )

png(file = "results/figures/BF_temp_data.png", height = 850, width = 1100)
t_plot
dev.off()

temp_data_plot_tiefora <- subset(temp_data_plot, Site == "Tiefora" & house == "mean" | Site == "Tiefora" & house == "ERA5" & Location == "Outdoor") %>% 
         mutate(Location = ifelse(house == "ERA5", "ERA5", Location))

temp_data_plot_tiefora$Location <- factor(temp_data_plot_tiefora$Location, levels = c("ERA5", "Outdoor", "Indoor"))

temp_data_plot_tiefora_daily <- temp_data_plot_tiefora %>% mutate(s_date = as.Date(f_date)) %>% 
  group_by(Site, Location, Hut, house, index, s_date) %>% summarise(r_temp = mean(Temp)) %>% 
  mutate(f_date = as.POSIXct(paste(s_date, "00:00:00"), format="%Y-%m-%d %H:%M:%S"))

temp_data_plot_tiefora <- temp_data_plot_tiefora %>% dplyr::mutate(r_temp = rollmean(temp_data_plot_tiefora$Temp,
                                                                              24 * 3,
                                                                              align = "center",
                                                                              fill = NA),
                                                                   s_date = as.Date(f_date))

# mean monthly temp
monthly_temp <- temp_data_plot_tiefora %>%
  group_by(Year, Month, Site, Location) %>% dplyr::summarise(m_temp = mean(Temp)) %>% 
  dplyr::left_join(temp_data_plot_tiefora, by = c("Year", "Month", "Site", "Location"))

# mean daily temp
daily_temp <- temp_data_plot_tiefora %>%
  group_by(Year, Month, s_date, Site, Location) %>% dplyr::summarise(d_temp = mean(Temp)) %>% 
  dplyr::left_join(temp_data_plot_tiefora, by = c("Year", "s_date", "Month", "Site", "Location"))

temp_plot_df <- rbind(temp_data_plot_tiefora[,c("f_date", "Location", "Temp")] %>% mutate(name = "Hourly"),
                      #daily_temp[,c("f_date", "Location", "d_temp")] %>% rename(Temp = d_temp) %>% mutate(name = "Daily"),
                      monthly_temp[,c("f_date", "Location", "m_temp")] %>% rename(Temp = m_temp) %>% mutate(name = "Monthly"))

temp_plot_df$name <- factor(temp_plot_df$name, levels = c("Hourly", "Daily", "Monthly"))

temp_plot <- ggplot(data = temp_plot_df,
       aes(x = f_date, y = Temp, col = Location, group = interaction(Location, name))) +
  geom_line(aes(linewidth = name, linetype = name), alpha = 1) +
  scale_linetype_manual(values = c(1, 2)) +
  # geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp"), se = FALSE, 
  #             linewidth = 1.3) +
  theme_bw() + theme(text = element_text(size = 18)) +
  scale_linewidth_manual(values = c(0.18, 1)) +
  #geom_smooth(formula = y ~ x, se = FALSE, method = "gam", alpha = 0.25, size = 0.75) +
  xlab("Date") +
  ylab("Temperature (°C)") +
  #facet_wrap(~Location, nrow = 3) +
  scale_y_continuous(breaks = seq(15, 50, 5), limits = c(15, 50)) +
  scale_colour_manual(values = c("#56B4E9","#CC79A7",  "#E69F00"), name = "Temperature\ndata") +
  #scale_colour_viridis_d(option = "mako") + 
  scale_x_datetime(date_labels = "%d %b",
               date_breaks = "1 month")

#############################################
### estimating the EIP given the DTR data ###
#############################################

ht_data_subset <- as.data.frame(temp_data_plot %>% mutate(date = as.Date(f_date)) %>%
                                  group_by(Site, Location, Hut) %>% 
                                  filter(date != min(date) & date < max(date) - days(20)) %>% ungroup())

# only selecting
u_f_start <- unique(ht_data_subset[,c("date", "index")])

p_u_f <- subset(u_f_start, index %in% u_t[which(u_t$Hut == "mean" & u_t$Site == "Tiefora" | u_t$Hut == "ERA5" & u_t$Site == "Tiefora"),"index"])

s_times <- seq(0, 23, 1)

u_f <- bind_rows(lapply(seq(1, length(s_times)), function(i, df, s_times){
  return(df %>% mutate(s_time = s_times[i]))},
  df = p_u_f,
  s_times = s_times
  ))

u_f$f_date <- ymd_h(paste(u_f$date, u_f$s_time))

#s_time = 0 # biting time
params <- c(mu = 0.1, shape = 47, v_EIP = TRUE)
state <- c(U = 0, E = c(1, rep(0, params[["shape"]]-1)),  I = 0)

out_m <- vector(mode = "list", length = nrow(u_f))
# out_l <- vector(mode = "list", length = nrow(u_f))
# out_u <- vector(mode = "list", length = nrow(u_f))

for(i in 1:nrow(u_f)){ 
  
  print(i)
  
  ind <- u_f[i, "index"]
  
  start_hour <- which(temp_data[[ind]]$f_date == u_f[i, "f_date"]) 
  start_time <- (start_hour  - 1)/24
  
  times <- seq(0, 480)/24 # must be in days
  
  # checking that the function is correct
  if(round(temp_data[[ind]][start_hour,"Temp"], digits = 4) == round(temp_fun[[ind]](start_time), digits = 4)){
    
    out <- NULL
    attempt <- 0
    while(is.null(out) && attempt <= 5){
      attempt <- attempt + 1
      try(
        out <- as.data.frame(ode(y = c(state), times = times, func = model_4_microclimate, parms = c(params, index = ind, EIP_CI = 2,
                                                                                            start_time_in = start_time))) %>% 
          mutate(start_hour = start_hour,
                 s_date = u_f[i, "date"][[1]],
                 s_time = u_f[i, "s_time"][[1]],
                 Location = u_t[u_f[i, "index"][[1]], "Location"],
                 Hut = u_t[u_f[i, "index"][[1]], "Hut"],
                 Site = u_t[u_f[i, "index"][[1]], "Site"]) %>%
          rowwise() %>%
          mutate(M = U + sum(c_across("E1":paste0("E",params[["shape"]]))) + I,
                 s_prev = I / M)
      )
    }
    
    if(is.null(out)){
      out_m[[i]] <- NA
    } else{
        out_m[[i]] <- out
      }
    
    rm(list = c("start_hour", "ind", "start_time", "times", "out"))
  } else{
    out_m[[i]] <- NA
    rm(list = c("start_hour", "ind", "start_time", "times"))
  }
}

#saveRDS(out_m, file = "results/DTR_season_model_outputs_m_r.rds")

out_m <- readRDS(file = "results/DTR_season_model_outputs_m_r.rds")

na_inds <- which(is.na(out_m))

# checking that EIP 90 has been reached or whether more time is required
all_EIP90_c <- lapply(seq(1, length(out_m)),
                   function(i, out){
                     s_prev <- out[[i]]$s_prev
                     return(max(s_prev) > 0.9)},
                   out = out_m)

sum(unlist(all_EIP90_c)) == length(out_m)

times <- seq(0, 480)/24

adjust_out <- function(out){
  out <- discard(out, ~ all(is.na(.x)))
  out <- bind_rows(out)
  s_out <- out %>% group_by(Location, Site, Hut, s_date, s_time) %>% 
    summarise(i_EIP_10 = which.min(abs(s_prev - 0.1)),
              d_EIP_10 = min(abs(s_prev - 0.1), na.rm = T),
              s_prev_10 = s_prev[i_EIP_10],
              EIP_10 = time[i_EIP_10],
              i_EIP_50 = which.min(abs(s_prev - 0.5)),
              d_EIP_50 = min(abs(s_prev - 0.5), na.rm = T),
              s_prev_50 = s_prev[i_EIP_50],
              EIP_50 = time[i_EIP_50],
              i_EIP_90 = which.min(abs(s_prev - 0.9)),
              d_EIP_90 = min(abs(s_prev - 0.9), na.rm = T),
              s_prev_90 = s_prev[i_EIP_90],
              EIP_90 = time[i_EIP_90])
  s_out$house <- ifelse(s_out$Hut == "ERA5" | s_out$Hut == "mean", s_out$Hut, str_sub(s_out$Hut, -1))
  return(s_out)
}

s_out_m <- adjust_out(out_m) %>% mutate(p = "median")

s_out_m$f_date <- as.POSIXct(ymd_h(paste(s_out_m$s_date, s_out_m$s_time)))

s_out_m_bt <- s_out_m
s_out_m <- subset(s_out_m, s_time %in% unique(bt_density_all$s_time))

s_out_m <- rbind(s_out_m, subset(s_out_m, Hut == "ERA5") %>% mutate(Location = "Outdoor"))

s_b_inds <- match(interaction(s_out_m$s_time, s_out_m$Location), interaction(bt_density_all$s_time, bt_density_all$Location))
s_out_m$scaled_EIP_50 <- s_out_m$EIP_50 * bt_density_all[s_b_inds, "d"]
s_out_m$scaled_EIP_10 <- s_out_m$EIP_10 * bt_density_all[s_b_inds, "d"]
s_out_m$scaled_EIP_90 <- s_out_m$EIP_90 * bt_density_all[s_b_inds, "d"]

sum_out_m <- s_out_m %>% group_by(Location, Site, Hut, s_date) %>% summarise(EIP_50 = sum(scaled_EIP_50),
                                                                             EIP_10 = sum(scaled_EIP_10),
                                                                             EIP_90 = sum(scaled_EIP_90))

sum_out_m <- sum_out_m %>% subset(Hut=="ERA5" & Location!="Indoor" | Hut == "mean")
sum_out_m <- sum_out_m %>% mutate(Location = ifelse(Hut == "ERA5", "ERA5", Location),
                              month = lubridate::month(s_date),
                              year = lubridate::year(s_date))



subset(sum_out_m, Location == "ERA5")[subset(sum_out_m, Location == "ERA5")$EIP_50 == min(subset(sum_out_m, Location == "ERA5")$EIP_50),]
subset(sum_out_m, Location == "ERA5")[subset(sum_out_m, Location == "ERA5")$EIP_50 == max(subset(sum_out_m, Location == "ERA5")$EIP_50),]

sum_out_m$Location <- factor(sum_out_m$Location, levels = c("ERA5", "Outdoor", "Indoor"))

# mean monthly EIP
# estimates from the hourly values
sum_out_m_m <- sum_out_m %>% group_by(Location, Site, Hut, month, year) %>% 
  summarise(m_EIP_50 = mean(EIP_50)) %>% dplyr::left_join(sum_out_m, by = c("Location", "Site", "Hut", "month", "year"))

# estimating the thermal performance curve values at each mean temperature
EIP_fun[[2]] <- Vectorize(EIP_fun[[2]])

monthly_temp$EIP <- EIP_fun[[2]](monthly_temp$m_temp)
monthly_temp$s_date <- as.Date(monthly_temp$f_date, format = "%Y-%m-%d")

daily_temp$EIP <- EIP_fun[[2]](daily_temp$d_temp)
daily_temp$s_date <- as.Date(daily_temp$f_date, format = "%Y-%m-%d")

EIP_plot_df <- rbind(sum_out_m_m[,c("s_date", "Location", "EIP_50")] %>% rename(EIP = EIP_50) %>% mutate(model = "DTR dependent: daily"),
                     #daily_temp[,c("s_date", "Location", "EIP")] %>% mutate(model = "DTR independent: daily"),
                     monthly_temp[,c("s_date", "Location", "EIP")] %>% mutate(model = "DTR independent: monthly")
                     )

EIP_plot <- ggplot(data = EIP_plot_df,
                   aes(x = s_date,
                       fill = Location, col =  Location, y = EIP, linetype = model))+
  geom_line(linewidth = 1, alpha = 1) +
  xlab("Date") + 
  scale_linetype_manual(values = c(1, 2)) +
  ylab(expression(paste("Predicted EIP (days)"))) +  
  theme_bw() + theme(text = element_text(size = 18)) +
  scale_y_continuous(breaks = seq(8, 16, 2), limits = c(8, 16)) +
  scale_colour_manual(values = c("#56B4E9","#CC79A7",  "#E69F00"), name = "Temperature\ndata") +
  scale_x_date(date_labels = "%d %b",
               date_breaks = "2 month",
               limits = as.Date(c("01/01/2020", "31/12/2020"), format = "%d/%m/%Y"))

###################################################
##### differences in the EIP with biting time #####
###################################################

s_out_m_plot <- s_out_m_bt %>% subset(Hut=="ERA5" & Location=="Indoor" | Hut == "mean") %>% mutate(Location = ifelse(Hut == "ERA5", "ERA5", Location))

png(file = "results/figures/EIP_bite_time.png", height = 570, width = 1040)
ggplot(data = s_out_m_plot %>% gather(key = "percentile", value = "EIP", EIP_10, EIP_50, EIP_90) %>% 
                           mutate(percentile = ifelse(percentile == "EIP_10", "EIP[10]",
                                                      ifelse(percentile == "EIP_50", "EIP[50]", "EIP[90]"))), 
                        aes(x = factor(s_time), y = EIP, group = interaction(s_time, percentile))) +
  geom_boxplot(alpha = 0.25, fill = "grey70") +
  facet_grid(vars(percentile), vars(Location), scales = "free", labeller = label_parsed) +
  theme_bw() + theme(text = element_text(size = 18)) +
  xlab("Hour of infection") + ylab("Predicted EIP (days)") +
  scale_y_continuous(breaks = seq(8, 18, 2), limits = c(7, 18))+
  scale_x_discrete(breaks = seq(0, 22, 2))
dev.off()

###########################################################
### estimating delta given the DTR data and biting time ###
###########################################################

s_times_b <- seq(0, 23.9, 0.1)

u_f_b <- bind_rows(lapply(seq(1, length(s_times_b)), function(i, df, s_times){
  return(df %>% mutate(s_time = s_times[i]))},
  df = p_u_f,
  s_times = s_times_b
))

u_f_b$f_date <- ymd_h(paste(u_f_b$date, floor(u_f_b$s_time)))

pred_DTR <- lapply(seq(1, nrow(u_f_b)), gen_delta_DTR, u_f_b = u_f_b, temp_fun = temp_fun, temp_data = temp_data)

#saveRDS(pred_DTR, file = "results/pred_DTR.rds")
pred_DTR <- readRDS(file = "results/pred_DTR.rds")
sum(is.na(pred_DTR))

pred_DTR <- bind_rows(pred_DTR)

u_f_b <- cbind(u_f_b, pred_DTR) 

u_f_b$Location <- u_t[u_f_b$index, "Location"]
u_f_b$Site = u_t[u_f_b$index, "Site"]
u_f_b$Hut = u_t[u_f_b$index, "Hut"]

u_f_b$s_time_f <- floor(u_f_b$s_time)

hourly_u_f_b <- u_f_b %>% group_by(date, index, f_date, Location, Site, Hut, s_time_f) %>% summarise(mean = mean(`50%`))

u_f_b_plot <- u_f_b %>% group_by(date, index, f_date, Location, Site, Hut, s_time_f) %>% summarise(mean = mean(`50%`)) %>% 
  mutate(Location = ifelse(Hut == "ERA5", "ERA5", Location),
                               day_m = lubridate::mday(date),
                               month = lubridate::month(date)) %>% 
         subset(month %in% seq(2, 11, 3)) %>% 
         mutate(month = lubridate::month(date, label = TRUE))

u_f_b_plot$Location <- factor(u_f_b_plot$Location, levels = c("ERA5", "Outdoor", "Indoor"))

png(file = "results/figures/pred_HMTP_s_time.png", height = 700, width = 950)
ggplot(data = u_f_b_plot, 
       aes(x = s_time_f, y = mean, group = interaction(Location, Site, Hut, month, s_time_f),
           fill = Location, col = Location)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5, col = "grey40") +
  #geom_bar(stat = "identity", alpha = 0.25, col = "black", position = "dodge") +
  #geom_line() +
  facet_wrap(~month) +
  theme_bw() + theme(text = element_text(size = 18)) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 22, 2)) +
  ylab("Predicted HMTP for each hourl and day of the month") +
  xlab("Hour of infection") +
  scale_fill_manual(values = c("#56B4E9","#CC79A7",  "#E69F00"), name = "Temperature\ndata") +
  scale_colour_manual(values = c("#56B4E9","#CC79A7",  "#E69F00"), name = "Temperature\ndata")
dev.off()


# assuming mosquitoes only bite during the recorded times
hourly_u_f_b <- hourly_u_f_b %>% subset(s_time_f %in% bt_density_all$s_time) # calculating the daily expected value due to differences in the biting times

hourly_u_f_b$bt_density <- bt_density_all[match(interaction(hourly_u_f_b$Location, hourly_u_f_b$s_time_f), interaction(bt_density_all$Location, bt_density_all$s_time)),"d"]
hourly_u_f_b <- hourly_u_f_b %>% mutate(scaled_mean = mean * bt_density)

daily_delta <- hourly_u_f_b %>% group_by(date, index, Location, Site, Hut) %>% 
  summarise(delta = sum(scaled_mean),
            tot_dens = sum(bt_density)) %>% mutate(Location = ifelse(Hut == "ERA5", "ERA5", Location),
                                                 week = week(date),
                                                 month = month(date),
                                                 year = year(date))

unique(daily_delta$tot_dens) # should equal 1

daily_delta <- daily_delta %>% group_by(index, Location, Site, Hut, month, year) %>% summarise(m_delta = mean(delta)) %>% 
  left_join(daily_delta, by = c("index", "Location", "Site", "Hut", "month", "year"))

monthly_temp <- monthly_temp %>% rowwise() %>% mutate(s_temp = (m_temp - 23.06936) / 4.361642,
                        delta_m = gen_delta(fit = fit, temp = s_temp)[2])

daily_temp <- daily_temp %>% rowwise() %>% mutate(s_temp = (d_temp - 23.06936) / 4.361642,
                                                      delta_d = gen_delta(fit = fit, temp = s_temp)[2])

delta_plot_df <- rbind(daily_delta[,c("date", "Location", "Site", "delta")] %>% mutate(model = "DTR dependent: daily"),
                                  #daily_temp[,c("s_date", "Location", "Site", "delta_d")] %>% rename(date = s_date, delta = delta_d) %>% mutate(model = "DTR independent: daily"),
                                  monthly_temp[,c("s_date", "Location", "Site", "delta_m")] %>% rename(date = s_date, delta = delta_m) %>% mutate(model = "DTR independent: monthly"))

delta_plot <- ggplot(data = delta_plot_df, 
                     aes(x = date, y = delta, col = Location, group = interaction(Location, model), linetype = model)) +
  geom_line(size = 1, alpha = 0.8) +
  scale_linetype_manual(values = c(1, 2)) +
  ylab("Predicted HMTP") +
  xlab("Date") + theme_bw() + theme(text = element_text(size = 18)) +
  scale_colour_manual(values = c("#56B4E9","#CC79A7",  "#E69F00"), name = "Temperature\ndata") + #"#00AFBB", "grey70",  "#FFDB6D"
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_date(date_labels = "%d %b",
               date_breaks = "2 month",
               limits = as.Date(c("01/01/2020", "31/12/2020"), format = "%d/%m/%Y"))


png(file = "results/figures/pred_EIP_HMTP.png", height = 1000, width = 1100)
plot_grid(
  plot_grid(
    temp_plot + theme(legend.position = "none"),
    EIP_plot + theme(legend.position = "none"),
    delta_plot + theme(legend.position = "none"),
    nrow = 3, 
    labels = c("A", "B", "C")),
  get_legend(EIP_plot),
  nrow = 1, 
  rel_widths = c(0.9, 0.28)
  )
dev.off()
