# Author: Isaac J Stopard
# Version: 1.0.0
# Last updated: February 2025
# Notes: script to fit the Ross-Macdonald type malaria transmission model to the sporozoite prevalence data from Tiefora

rm(list = ls())

source(file = "utils/model_functions.R") 
source(file = "utils/read_libraries_data.R")
source(file = "utils/read_bt_data.R")

###########################################
##### reading in the temperature data #####
###########################################

read_ERA5 <- function(file_){
  return(
    as.data.frame(read_csv(file = file_) %>%
                           mutate(temp = t2m -  273.5,
                                  f_date = time,
                                  Site = "Tiefora",
                                  Location = "ERA5",
                                  Hut = "ERA5",
                                  date = as.Date(time),
                                  hour = hour(time)))
    )
}

temp_data <- rbind(read_ERA5(file_ = "data/ERA5/t2m_2015_BF.csv"),
      read_ERA5(file_ = "data/ERA5/t2m_2016_BF.csv"),
      read_ERA5(file_ = "data/ERA5/t2m_2017_BF.csv"),
      read_ERA5(file_ = "data/ERA5/t2m_2018_BF.csv"),
      read_ERA5(file_ = "data/ERA5/t2m_2019_BF.csv")) %>% mutate(month = month(date),
                                     year = year(date))

dates <- unique(temp_data$date)

temp_data <- left_join(temp_data, temp_data %>% group_by(month, year) %>% summarise(monthly_temp = mean(temp)), by = c("month", "year")) 
temp_data <- left_join(temp_data, temp_data %>% group_by(date) %>% summarise(daily_temp = mean(temp)), by = c("date")) 

# x is the start of the hour 
temp_fun <- approxfun(x = seq(0, nrow(temp_data)-1)/24, y = temp_data$temp, yright = NA, yleft = NA, method = "linear") # hourly temperature fluctuations

temp_fun_day <- approxfun(x = seq(0, nrow(temp_data)-1)/24, y = temp_data$daily_temp, yright = NA, yleft = NA, method = "constant") # mean daily temperature fluctuations

temp_fun_month <- approxfun(x = seq(0, nrow(temp_data)-1)/24, y = temp_data$monthly_temp, yright = NA, yleft = NA, method = "constant") # mean monthly temperature fluctuations

############################################
### EIP thermal performance curve values ###
############################################

# assumes a linear interpolation for values between simulated values
EIP_fun <- approxfun(mean_EIP$temp, mean_EIP$`50%`, yleft = max(mean_EIP$`50%`), yright = min(mean_EIP$`50%`)) #vector(mode = "list", length = 3)

######################################################################
### calculating the daily HMTP given the temp and biting time data ###
######################################################################

### DTR-dependent model
delta_df <- expand.grid("date" = dates, "hour" = seq(0, 23.9, 0.1)) %>%
  mutate(days = difftime(date, "2015-01-01", units = "days")) %>% subset(date < "2019-01-01")

gen_delta_DTR_all <- function(i, delta_df, temp_fun, temp_data, fit, mean_temp = 23.06936, sd_temp = 4.361642){
  if(i %% 1000 == 0){print(i)}
  
  start_hour <- delta_df[i,"hour"]
  
  start_day <- delta_df[i,"days"][[1]]
  
  start_time <- start_day + start_hour/24 
  
  temp_in <- max(temp_fun(seq(start_time, start_time + (10/24), 0.0001)))
  temp_in <- (temp_in - mean_temp) / sd_temp # scaling so on same scale as parameter fits
  
  if(round(subset(temp_data, date == delta_df[i, "date"] & hour == floor(delta_df[i,"hour"]))$temp, digits = 4) == round(temp_fun(start_day + floor(start_hour)/24), digits = 4)){
    return(gen_delta(fit = fit, temp = temp_in))} else{
      return(NA)
    }
}

pred_DTR <- lapply(seq(1, nrow(delta_df)), 
                   gen_delta_DTR_all, 
                   delta_df = delta_df, 
                   temp_fun = temp_fun, 
                   temp_data = temp_data,
                   fit = fit)

saveRDS(pred_DTR, file = "results/pred_DTR_EIR_fit.rds")
pred_DTR <- readRDS(file = "results/pred_DTR_EIR_fit.rds")

pred_DTR <- bind_rows(pred_DTR)

# calculating the mean hourly HMTP given the DTR-dependent model
hourly_delta <- cbind(delta_df, pred_DTR) %>% mutate(hour = floor(hour)) %>% group_by(date, hour, days) %>%
  summarise(delta = mean(`50%`))

bt_density$hour <- c(seq(19, 23), seq(0, 5)) 

hourly_delta[,"d"] <- bt_density[match(hourly_delta$hour, bt_density$hour), "d"]

hourly_delta$d[is.na(hourly_delta$d)==1] <- 0

# calculating the mean daily HMTP given the DTR-dependent model and differences in mosquito biting times
daily_delta <- hourly_delta %>% mutate(scaled_delta = delta * d) %>% 
  group_by(date, days) %>% summarise(delta = sum(scaled_delta),
                                     d = sum(d))

daily_delta <- daily_delta[order(daily_delta$date),]

delta_fun_DTR_dependent <- approxfun(x = daily_delta$days, y = daily_delta$delta, yleft = NA, yright = NA, method = "constant") # function

### DTR independent model HMTP

temp_data$DTR_ind_monthly_delta <- t(sapply((temp_data$monthly_temp - mean_temp)/sd_temp, gen_delta, fit = fit))[,2]

delta_fun_DTR_ind_monthly <- approxfun(x = seq(0, nrow(temp_data)-1)/24, y = temp_data$DTR_ind_monthly_delta, yright = NA, yleft = NA, method = "constant")

#########################################
##### daily mosquito mortality rate #####
#########################################

# Marten's model 2 from: https://link.springer.com/article/10.1186/1756-3305-6-20

martens_fun <- Vectorize(function(temp){
  if(temp < 4 | temp > 39){
    return(NA)
  } else{
    surv_prob <- exp(-(1/(- 4.4 + 1.31 * temp - 0.03 * temp^2))) # daily survival probability
    mu <- -log(surv_prob)
  return(mu)
  }
})

# DTR independent daily model
temp_data <- temp_data %>% mutate(daily_mu = martens_fun(daily_temp),
                                  monthly_mu = martens_fun(monthly_temp))
mu_fun_DTR_ind_daily <- approxfun(x = seq(0, nrow(temp_data)-1)/24, y = temp_data$daily_mu, yright = NA, yleft = NA, method = "constant")

# DTR independent daily model
mu_fun_DTR_ind_monthly <- approxfun(x = seq(0, nrow(temp_data)-1)/24, y = temp_data$monthly_mu, yright = NA, yleft = NA, method = "constant")

# constant model
c_mu <- mean(temp_data$daily_mu)

##############################################
##### estimating the mosquito birth rate #####
##############################################

# fitted sporozoite prevalence, HBR and EIR data
spz_data <- readRDS(file = "data/format_spz_data.rds")
spz_data$raw$day_plot <- difftime(spz_data$raw$Date, "2015-01-01", units = c("days"))
spz_data$raw_m$day_plot <- difftime(spz_data$raw_m$Date, "2015-01-01", units = c("days"))

spline_ma <- subset(spz_data$pred_m, Location == "IN") %>% mutate(day = day_y)

# replicating the same mosquito seasonality for each year
spline_ma_pad <- spline_ma[rep(seq(1, nrow(spline_ma)),5),] %>% mutate(day_o = day,
                                                                               day = seq(0, (n()-1)))

# adult mosquito emergence rates
beta_c <- diff(log(spline_ma_pad$p)) + c_mu # mosquito birth rates assuming a constant mortality rate

beta_daily <- diff(log(spline_ma_pad$p)) + subset(temp_data, hour == 0)$daily_mu[1:(nrow(spline_ma_pad)-1)]
beta_monthly <- diff(log(spline_ma_pad$p)) + subset(temp_data, hour == 0)$monthly_mu[1:(nrow(spline_ma_pad)-1)]

beta_fun_c <- approxfun(x = seq(0, (nrow(spline_ma_pad)-2)), y = beta_c, yright = NA, yleft = NA, method = "constant")

beta_fun_DTR_ind_daily <- approxfun(x = seq(0, (nrow(spline_ma_pad)-2)), y = beta_daily, yright = NA, yleft = NA, method = "constant")

beta_fun_DTR_ind_monthly <- approxfun(x = seq(0, (nrow(spline_ma_pad)-2)), y = beta_monthly, yright = NA, yleft = NA, method = "constant")

# model
# Adapted Aron & May version of the classic malaria transmission model
# S: susceptible mosquitoes (absolute numbers)
# E: exposed mosquitoes  (absolute numbers)
# I: infectious mosquitoes  (absolute numbers)
# z: infectious mosquitoes  (prevalence)
# Ih: infected people (absolute numbers)
# pop: total human population
# x: infected people (prevalence)
# a: per-mosquito biting rate per day
# beta: per-mosquito birth rate per day
# b: mosquito-to-human transmission probability
# c: human-to-mosquito transmission probability
# r: per-person recovery rate per day
# mu: per-capita mosquito mortality rate per day

classic_malaria_model <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    
    if(model == 1){#"DTR_dependent"
      temp <- temp_fun(t) # temperature at time t
      c <- delta_fun_DTR_dependent(t) 
      
      mu <- mu_fun_DTR_ind_daily(t)
      beta <- beta_fun_DTR_ind_daily(t)
      
    } else if(model == 2){#"DTR_independent_monthly"
      temp <- temp_fun_month(t)
      c <- delta_fun_DTR_ind_monthly(t)
      mu <- mu_fun_DTR_ind_monthly(t)
      beta <- beta_fun_DTR_ind_monthly(t)
      
    } else if(model == 3){#"constant"
      temp <- temp_fun(t)
      c <- m_HMTP
      mu <- m_mu
      beta <- beta_fun_c(t)
    }
    
    if(model != 3){"not constant"
      rate <- shape/EIP_fun(temp) # PDR at time t
    } else{
      rate <- m_rate
    }
    
    c <- c * scale_HMTP # HMTP at time t
    
    # calculating the total numbers of exposed mosquitoes
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    z <- I/M
    
    dI_h <- a * b * M * z * (pop - I_h) - r * I_h
    
    x <- I_h / pop
    
    dS <- beta * M - (mu * S) - a * c * S * x
    
    dE1 <- a * c * S * x - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    EIR <- M * a * z
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI, dI_h), temp = temp, rate = rate, M = M, c = c, beta = beta, x = x, z = z, mu = mu, EIR = EIR)
  })
}

# parameters
c_mu_in <- c_mu
c_HMTP_in <- mean(daily_delta$delta)
c_rate_in <- mean(47/EIP_fun(temp_data$temp))
spline_ma_pad_in <- spline_ma_pad

run_malaria_model <- function(scale_m, 
                              scale_HMTP,
                              r,
                              a,
                              model,
                              spline_ma_pad = spline_ma_pad_in,
                              c_HMTP = c_HMTP_in,
                              c_rate = c_rate_in,
                              c_mu = c_mu_in
                              ){
  
  b <- 0.5 # references
  
  #r <- 1/100
  shape <- 47
  
  params <- c(model = model,
              m_mu = c_mu,
              b = b,
              r = r,
              a = a,
              m_HMTP = c_HMTP,
              shape = shape,
              m_rate = c_rate,
              pop = 1000,
              scale_HMTP = scale_HMTP)
  
  M1 <- spline_ma_pad[1, "p"] * scale_m
  
  state <- c(S = M1, E = c(0, rep(0, shape-1)), I = 0, I_h = 1)
  
  times <- seq(0, 1460, 0.5)
  
  out <- NULL
  attempt <- 0
  while(is.null(out) && attempt <= 20){
    attempt <- attempt + 1
    try(
      out <- ode(y = c(state), times = times, func = classic_malaria_model, parms = c(params)) %>% 
        as.data.frame() %>%
        #rowwise() %>%
        mutate(I = ifelse(I < 0 & abs(I) < 1E-05, 0, I),
          #check_M = S + sum(c_across("E1":paste0("E",params[["shape"]]))) + I,
               s_prev = I / M)
    )
  }
  
  out$date <- lubridate::ymd_hms("2015 01 01 00 00 00") + days(floor(out$time)) + hours(out$time %%1 * 24)
  out$year <- lubridate::year(out$date)
  out$day_y <- lubridate::yday(out$date)
  
return(out %>% as.data.frame())
}

#ov_stan <- median(rstan::extract(spz_data$t_EIR$stanfit, "reciprocal_dispersion")$reciprocal_dispersion)

# function to calculate the likelihoods
calc_ll <- function(par, # overdispersion parameter
                    model,
                    a,
                    r,
                    v_M,
                    #f_EIR,
                    #v_ov,
                    #ov_in = ov_stan,
                    data = spz_data$raw,
                    min_year = 2015, # only fit to data from 2017
                    max_year = 2018
                    ){
  
  if(v_M == TRUE){
      scale_m <- par[1]
      scale_HMTP <- 1
  } else if(v_M == "both"){
    scale_m <- par[3]
    scale_HMTP <- par[1]
  } else{
    scale_m <- 1/a # multiplied by 3 to account for the biting rate
    scale_HMTP <- par[1]
  }
  
  ov <- par[2]
  
  data <- subset(data, year > min_year & year < max_year) %>% as.data.frame()
  
  m_sim <- run_malaria_model(scale_m = scale_m,
                             scale_HMTP = scale_HMTP,
                             model = model,
                             r = r,
                             a = a
                             )
  
  data <- data %>% mutate(pred_z = as.vector(m_sim[match(data$day_plot, m_sim$time),"z"]),
                            log_likelihood_spz = dbinom(x = tot_p, size = tot, prob = pred_z, log = TRUE))
    
    ll_z <- -sum(data$log_likelihood_spz)
    
  return(ll_z) 
}

combs <- expand.grid(model = c(1, 2, 3),
                     r = c(1/50, 1/100),
                     a = c(1/2, 1/3)) 

combs_in_v_HMTP <- combs %>% dplyr::mutate(v_M = FALSE, upper = case_when(model == 1 ~ 1/max(delta_fun_DTR_dependent(daily_delta$days)),
                                                                          model == 2 ~ 1/max(delta_fun_DTR_ind_monthly(daily_delta$days)),
                                                                          model == 3 ~ 1/c_HMTP_in)) 


v_list <- c("temp_fun", "temp_fun_month",
            "c_rate_in",
            "delta_fun_DTR_dependent", "delta_fun_DTR_ind_monthly", "c_HMTP_in", 
            "beta_fun_DTR_ind_daily", "beta_fun_DTR_ind_monthly", "beta_fun_c",
            "mu_fun_DTR_ind_daily", "mu_fun_DTR_ind_monthly", "c_mu_in", "c_mu", 
            "spline_ma_pad_in", 
            "EIP_fun",
            "run_malaria_model",
            "calc_ll", 
            "spz_data", 
            "classic_malaria_model", 
            #"ov_stan", 
            "combs_in_v_HMTP")

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = v_list)

pred_ll_v_M_t <- foreach(i=1:nrow(combs_in_v_HMTP),
                       .packages = (.packages())
) %dopar% {
  #dfoptim::nmkb
  optimize(#par = c(1), 
           f = calc_ll,
           lower = c(0),
           upper = c(combs_in_v_HMTP[i, "upper"]), # upper limit prevent the HMTP being greater than one
           model = combs_in_v_HMTP[i, "model"],
           r = combs_in_v_HMTP[i, "r"],
           a = combs_in_v_HMTP[i, "a"],
           #v_ov = FALSE,
           v_M = FALSE,
           #f_EIR = FALSE,
           max_year = 2017,
           data = spz_data$raw,
           #control = list(tol = 1e-04)
           tol = 1e-04
           )
}
saveRDS(pred_ll_v_M_t, 
        file = "results/EIR_fit_scale_hmtp_uni_day_n.rds")
stopCluster(cl)

pred_ll_v_M_t <- readRDS(file = "results/EIR_fit_scale_hmtp_uni_day_n.rds")

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c(v_list, "pred_ll_v_M_t"))

pred_vals <- foreach(i=1:nrow(combs_in_v_HMTP),
                     .packages = (.packages())
) %dopar% {
  run_malaria_model(scale_m = 1/combs_in_v_HMTP[i, "a"],
                    scale_HMTP = pred_ll_v_M_t[[i]]$minimum,
                    model = combs_in_v_HMTP[i, "model"],
                    r = combs_in_v_HMTP[i, "r"],
                    a = combs_in_v_HMTP[i, "a"]
                    
  )
}
saveRDS(pred_vals, 
        file = "results/EIR_fit_scale_hmtp_pred_day_n.rds")
stopCluster(cl)

pred_vals <- readRDS(file = "results/EIR_fit_scale_hmtp_pred_day_n.rds")

for(i in 1:length(pred_vals)){
  pred_vals[[i]] <- subset(pred_vals[[i]], year != 2015 & year != 2019 & time %% 1 == 0)
  pred_vals[[i]] <- pred_vals[[i]] %>% mutate(hour = lubridate::hour(date),
                                              date_p = as.Date("2016-12-31", format = "%Y-%m-%d") + day_y)
  pred_vals[[i]]$s_prev <- pred_vals[[i]]$s_prev * 100
  pred_vals[[i]] <- pred_vals[[i]][-nrow(pred_vals[[i]]),] # removing the last row
}


spz_data$raw$z <- (spz_data$raw$tot_p / spz_data$raw$tot) * 100
tot_hlc_m <- spz_data$raw_m %>% subset(Location == "IN") %>% group_by(Date, day_plot) %>% summarise(m_tot_hlc = mean(tot_hlc))
spz_data$raw$tot_hlc <- unlist(unname(as.vector(tot_hlc_m[match(spz_data$raw$Date, tot_hlc_m$Date), "m_tot_hlc"])))
spz_data$raw$EIR <- (spz_data$raw$tot_p / spz_data$raw$tot) * spz_data$raw$tot_hlc

calc_MAE <- function(data = spz_data$raw, pred, v1 = "EIR", v2 = "EIR", year_in = 2017){
  
  data <- subset(data, year > year_in) %>% as.data.frame()
  inds <- match(data$Date, pred$date)
  
  return(mean(abs(data[, v1] - pred[inds, v2])))
}

calc_auc <- function(data = spz_data$raw, year_in = 2017, pred){
  data <- subset(data, year > year_in) %>% as.data.frame()
  inds <- match(data$Date, pred$date)
  
  data$pred_z <- pred[inds, "z"]
  p <- c()
  pred_z <- c()
  for(i in 1:nrow(data)){
    p <- c(p, rep(1, data[i, "tot_p"]))
    pred_z <- c(pred_z, rep(data[i, "pred_z"], data[i, "tot_p"]))
    p <- c(p, rep(0, (data[i, "tot"] - data[i, "tot_p"])))
    pred_z <- c(pred_z, rep(data[i, "pred_z"], (data[i, "tot"] - data[i, "tot_p"])))
  }
  
  return(auc(roc(p, pred_z)))
  
}

for(i in 1:nrow(combs)){
  combs[i, "s_HMTP"] <- round(pred_ll_v_M_t[[i]]$minimum, digits = 2)
  combs[i, "log-likelihood"] <- - round(pred_ll_v_M_t[[i]]$objective, digits = 2)
  combs[i, "auc_z"] <- round(calc_auc(data = spz_data$raw, pred = pred_vals[[i]], year_in = 2017), digits = 2)
  combs[i, "MAE_z"] <- round(calc_MAE(data = spz_data$raw, pred = pred_vals[[i]], v1 = "z", v2 = "s_prev", year_in = 2017), digits = 2)
}

combs <- combs %>% mutate(model = case_when(model == 1 ~ "DTR-dependent",
                                   model == 2 ~ "DTR-independent",
                                   model == 3 ~ "Constant"))

write.csv(combs, file = "results/EIR_fit.csv")

#################
##### plots #####
#################

adjust_date <- function(df){
  
  df <- df[rep(seq(1, nrow(df)),5),] %>% mutate(day_o = day_y,
                                                day = seq(0, (n()-1)),
                                                Date = lubridate::ymd_hms(paste(2015, 01, 01, 00, 00, 00, sep = " ")) + lubridate::days(day),
                      date_p = as.Date("2016-12-31", format = "%Y-%m-%d") + day,
                      year = lubridate::year(Date)) %>% 
    filter(year %in% c(2016, 2017, 2018))
  return(df)
}

spz_data$pred_s_day <- adjust_date(spz_data$pred_s)

spz_data$pred_m_day <- rbind(adjust_date(subset(spz_data$pred_m, Location == "IN")) %>% mutate(Location = "Inside"),
                             adjust_date(subset(spz_data$pred_m, Location == "OUT")) %>% mutate(Location = "Outside"))

spz_data$raw <- spz_data$raw %>% mutate(date_p = as.Date("2016-12-31", format = "%Y-%m-%d") + day_y,
                                        train = ifelse(year == 2018, "testing", "training"),
                                        lower_p = DescTools::BinomCI(x = tot_p, n = tot, conf.level = 0.95,
                                                                     method = c("clopper-pearson"), sides = c("two.sided"))[,"lwr.ci"],
                                        upper_p = DescTools::BinomCI(x = tot_p, n = tot, conf.level = 0.95,
                                                                     method = c("clopper-pearson"), sides = c("two.sided"))[,"upr.ci"])

spz_data$raw_m <- spz_data$raw_ %>% mutate(date_p = as.Date("2016-12-31", format = "%Y-%m-%d") + day_y,
                                        train = ifelse(year == 2018, "testing", "training"))


mos_plot <- ggplot() +
  xlab("Date") +
  ylab("HBR per person per day") +
  geom_line(data = spz_data$pred_m_day, aes(x = Date, y = p, col = Location, group = Location), linewidth = 1.25) +
  geom_ribbon(data = spz_data$pred_m_day, aes(x = Date, ymin = lower, ymax = upper, fill = Location), alpha = 0.1) +
  theme_bw() + theme(text = element_text(size = 18),
                     legend.text = element_text(size = 10),
                     legend.title = element_text(size = 12)) +
  scale_colour_manual(name = "Location", values = c("#000000", "grey50")) +
  scale_y_continuous(limits = c(0, 450)) +
  scale_fill_manual(name = "Location", values = c("#000000", "grey50")) +
  coord_cartesian(xlim = as.POSIXct(c("01/01/2016", "01/01/2019"), format = "%d/%m/%Y")) +
  geom_point(data = spz_data$raw_m %>% mutate(Location = ifelse(Location == "IN", "Inside", "Outside")), 
             aes(x = Date, y = tot_hlc, group = factor(year), col = Location), size = 4, alpha = 0.8) +
  theme(legend.position = c(0.8, 0.9), legend.box = "horizontal")
  
spz_plot <- ggplot() +
  xlab("Date") +
  ylab("Sporozoite prevalence") +
  #geom_ribbon(data = spz_data$pred_s_day, aes(x = Date, ymin = lower, ymax = upper), alpha = 0.1) +
  geom_pointrange(data = spz_data$raw, aes(x = Date, y = tot_p/tot, ymin = lower_p, ymax = upper_p, group = factor(year), shape = train), 
                  size = 0.95, alpha = 0.85) +
  #geom_line(data = spz_data$pred_s_day, aes(x = Date, y = p, col = "GAM"), linewidth = 1.25) +
  theme_bw() + theme(text = element_text(size = 18),
                     legend.text = element_text(size = 10),
                     legend.title = element_text(size = 12)) +
  scale_colour_manual(name = "Model", values = c("#0072B2", "#009E73", "#E69F00")) +
  scale_shape_manual(name = "Data", values = c(16, 15)) +
  #scale_fill_manual(name = "year", values = c("#000000", "#0072B2", "#009E73")) +
  coord_cartesian(xlim = as.POSIXct(c("01/01/2016", "01/01/2019"), format = "%d/%m/%Y")) +
  scale_y_sqrt(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  geom_line(data = pred_vals[[10]], aes(x = date, y = z, col = "DTR-dependent"), linewidth = 1.25) +
  geom_line(data = pred_vals[[11]], aes(x = date, y = z, col = "DTR-independent"), linewidth = 1.25) +
  geom_line(data = pred_vals[[12]], aes(x = date, y = z, col = "Constant"), linewidth = 1.25) +
  theme(legend.position = c(0.15, 0.85), legend.box = "horizontal")

ggsave(file = "results/figures/Figure_5.pdf",
plot_grid(
    mos_plot,
    spz_plot,
  ncol = 1,
  labels = c("A", "B")
), 
width = 800/30, height = 800/30, units = "cm", device = "pdf")

