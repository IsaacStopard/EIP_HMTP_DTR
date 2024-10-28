#################
### functions ###
#################

# function to calculate the mean EIP from the get_EIP function
calc_mean_EIP <- function(EIP_index, i){
  p <- data.frame("a" = EIP_index$shape_total_S[,i],
                  "b" = EIP_index$rate_total_S[,i],
                  "mu" = EIP_index$mu[,i],
                  "k" = EIP_index$k[,i]) %>% rowwise() %>% mutate(mean = tryCatch({integrate(mean_func, a = a, b = b,
                                                                                             mu = mu, k = k, lower = 0, upper = Inf)[[1]]},
                                                                                  error = function(e){
                                                                                    NA
                                                                                  }))
  return(c("mean" = mean(na.omit(p$mean)), quantile(na.omit(p$mean), c(0.025, 0.5, 0.975)),
           "n_na" = sum(is.na(p$mean))))
}

# new method for extracting the EIP
generate_EIP <- function(fit_all, l_n = 11){
  fit_ <- rstan::extract(fit_all)
  EIP_df <- data.frame("a" = fit_$shape_total_S[,1], "b" = fit_$rate_total_S[,1], "mu" = fit_$mu[,1], "k" = fit_$k[,1], "index_temp" = rep(1, (iterations - warmup)*chains))
  EIP_df[,"EIP[10]"] <- Inv_EIP_CDF(EIP_df$a, EIP_df$b, EIP_df$mu, EIP_df$k, 0.1)
  EIP_df[,"EIP[50]"] <- Inv_EIP_CDF(EIP_df$a, EIP_df$b, EIP_df$mu, EIP_df$k, 0.5)
  EIP_df[,"EIP[90]"] <- Inv_EIP_CDF(EIP_df$a, EIP_df$b, EIP_df$mu, EIP_df$k, 0.9)
  
  EIP_quantiles <-cbind(data.frame("index_temp" = rep(1, 3), "EIP" = c("EIP[10]", "EIP[50]", "EIP[90]"), "mean" = c(mean(EIP_df[,"EIP[10]"]),
                                                                                                                    mean(EIP_df[,"EIP[50]"]),
                                                                                                                    mean(EIP_df[,"EIP[90]"]))),
                        as.data.frame(rbind(quantile(EIP_df[,"EIP[10]"], c(0.025, 0.5, 0.975)),
                                            quantile(EIP_df[,"EIP[50]"], c(0.025, 0.5, 0.975)),
                                            quantile(EIP_df[,"EIP[90]"], c(0.025, 0.5, 0.975)))))
  colnames(EIP_quantiles) <- c("index_temp", "EIP", "mean", "lower", "median", "upper")
  
  for(i in 2:l_n){
    placeholder_df <- data.frame("a" = fit_$shape_total_S[,i], "b" = fit_$rate_total_S[,i], "mu" = fit_$mu[,i], "k" = fit_$k[,i], "index_temp" = rep(i, (iterations - warmup)*chains))
    placeholder_df[,"EIP[10]"] <- Inv_EIP_CDF(placeholder_df$a, placeholder_df$b, placeholder_df$mu, placeholder_df$k, 0.1)
    placeholder_df[,"EIP[50]"] <- Inv_EIP_CDF(placeholder_df$a, placeholder_df$b, placeholder_df$mu, placeholder_df$k, 0.5)
    placeholder_df[,"EIP[90]"] <- Inv_EIP_CDF(placeholder_df$a, placeholder_df$b, placeholder_df$mu, placeholder_df$k, 0.9)
    EIP_df <- rbind(EIP_df, placeholder_df)
    
    placeholder_quantiles <-cbind(data.frame("index_temp" = rep(i, 3), "EIP" = c("EIP[10]", "EIP[50]", "EIP[90]"), "mean" = c(mean(placeholder_df[,"EIP[10]"]),
                                                                                                                              mean(placeholder_df[,"EIP[50]"]),
                                                                                                                              mean(placeholder_df[,"EIP[90]"]))),
                                  as.data.frame(rbind(quantile(placeholder_df[,"EIP[10]"], c(0.025, 0.5, 0.975)),
                                                      quantile(placeholder_df[,"EIP[50]"], c(0.025, 0.5, 0.975)),
                                                      quantile(placeholder_df[,"EIP[90]"], c(0.025, 0.5, 0.975)))))
    colnames(placeholder_quantiles) <- c("index_temp", "EIP", "mean", "lower", "median", "upper")
    EIP_quantiles <- rbind(EIP_quantiles, placeholder_quantiles)
    rm(list = c("placeholder_df", "placeholder_quantiles"))
  }
  
  EIP_quantiles$temp <- all_data$unique_temp[EIP_quantiles$index_temp]
  return(list("EIP_df" = EIP_df, "EIP_quantiles" = EIP_quantiles))
}

EIP_PDF <- function(t, a, b, mu, k){
  return(-(((1/b)^-a) * exp(-b * t) * (t^(-1+a)) * mu * ((k / (k + mu))^(1 + k)) * (((k + mu * Rgamma(a, b*t)) / (k + mu)) ^ (-1 - k))
           / ((-1 + (k / (k + mu))^k) * gamma(a))))
}

EIP_CDF <- function(a, b, mu, k, t){
  return(
    (1 - k^k * (k + mu * Rgamma(a, b * t, lower = TRUE))^-k) / 
      (1 - (k / (k + mu))^k)
  )
}

# function to calculate the mean of the EIP PDF 
# to get the mean this function must be integrated
mean_func <- function(t, a, b, mu, k){
  return(t * EIP_PDF(t, a, b, mu, k))
}

int_mean_func <- function(a, b, mu, k){
  return(tryCatch({integrate(mean_func, a = a, b = b,
                             mu = mu, k = k, lower = 0, upper = Inf)[[1]]},
                  error = function(cond){
                    return(NA)
                  }))
}

v.int_mean_func <- Vectorize(int_mean_func) 

# function to calculate the variance of the EIP PDF 
# to get the variance this function must be integrated
var_func <- function(t, a, b, mu, k){
  return(t^2 * EIP_PDF(t, a, b, mu, k))
}

int_var_func <- function(a, b, mu, k){
  return(tryCatch({integrate(var_func, a = a, b = b,
                             mu = mu, k = k, lower = 0, upper = Inf)[[1]] -
      int_mean_func(a = a, b = b, mu = mu, k = k)^2
  },
  error = function(cond){
    return(NA)
  }))
}

v.int_var_func <- Vectorize(int_var_func)

# function that calculates the EIP_X given the percentile X
Inv_EIP_CDF <- function(a, b, mu, k, p){
  one <- ((-(p * (1 - (k / (k + mu))^k) - 1)/(k^k))^(-1/k) - k) / mu
  return((Rgamma.inv(a, one, lower = TRUE)/b))
}

# function required to calculate the continuous 
c_eq_I_m <- function(t, a, b, mu, k, m){
  return(exp(-(t * m)) * EIP_PDF(t, a, b, mu, k))
}

# function that calculate the EIP from mSOS params
get_EIP <- function(params, temps, iter){
  shape_S <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  rate_S <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  mu <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  shape_O <- as.data.frame(replicate(length(temps), params$shape_O))
  rate_O <- as.data.frame(replicate(length(temps), params$rate_O))
  
  k <- as.data.frame(replicate(length(temps), params$k))
  
  for(i in 1:length(temps)){
    shape_S[,i] <- temps[i]^2 * params$a_shape_S + temps[i] * params$b_shape_S + params$c_shape_S
    rate_S[,i] <- temps[i] * params$m_rate_S + params$c_rate_S
    mu[,i] <- (1 / (1 + exp(-(temps[i]^2 * params$a_mu + temps[i] * params$b_mu + params$c_mu)))) * params$g_mu
    #mu[,i] <- params$a_mu * exp(-((temps[i] - params$b_mu)^2 / (2 * params$c_mu)^2))
  }
  
  mu_total_S <- (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S)
  sigma_sq_S  <- (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2)
  shape_total_S <- mu_total_S^2 / sigma_sq_S;
  rate_total_S <- mu_total_S / sigma_sq_S;
  mean_total_S <- shape_total_S / rate_total_S
  
  EIP_10 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  EIP_50 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  EIP_90 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  for(i in 1:length(temps)){
    EIP_10[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.1)
    EIP_50[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.5)
    EIP_90[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.9)
  }
  
  return(list("shape_O" = shape_O, "rate_O" = rate_O, "shape_S" = shape_S, "rate_S" = rate_S, "shape_total_S" = shape_total_S, 
              "rate_total_S" = rate_total_S, "mean_total_S" = mean_total_S, "EIP_10" = EIP_10, "EIP_50" = EIP_50, "EIP_90" = EIP_90,
              "mu" = mu, "k" = k))
}

get_EIP_stephensi <- function(params, temps, iter){
  #mu <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  rate_O <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  shape_S <- as.data.frame(replicate(length(temps), params$shape_sporozoite))
  rate_S <- as.data.frame(replicate(length(temps), params$rate_sporozoite))
  shape_O <- as.data.frame(replicate(length(temps), params$shape_oocyst))
  k <- as.data.frame(replicate(length(temps), params$k_NB))
  mu <- as.data.frame(replicate(length(temps), params$mu_NB))
  
  for(i in 1:length(temps)){
    rate_O[,i] <- temps[i] * params$m_rate_oocyst + params$c_rate_oocyst
  }
  
  mu_total_S <- (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S)
  sigma_sq_S  <- (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2)
  shape_total_S <- mu_total_S^2 / sigma_sq_S
  rate_total_S <- mu_total_S / sigma_sq_S
  mean_total_S <- shape_total_S / rate_total_S
  
  EIP_10 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  EIP_50 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  EIP_90 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  for(i in 1:length(temps)){
    EIP_10[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.1)
    EIP_50[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.5)
    EIP_90[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.9)
  }
  
  return(list("shape_O" = shape_O, "rate_O" = rate_O, "shape_S" = shape_S, "rate_S" = rate_S, "shape_total_S" = shape_total_S, 
              "rate_total_S" = rate_total_S, "mean_total_S" = mean_total_S, "EIP_10" = EIP_10, "EIP_50" = EIP_50, "EIP_90" = EIP_90,
              "mu" = mu, "k" = k))
}

# function that calculates the PDF for each MCMC iteration
get_PDF_multi <- function(fit_, gt = c(1, 11), temp = c(17, 30)){
  times <- seq(0,170,0.1)
  n_t <- length(times)
  n_iter <- length(fit_$shape_total_S[,1])
  
  EIP_MCMC_list <- list()
  
  PDF_quantiles_out <- data.frame(matrix(data = NA, nrow = 0, ncol = 6))
  colnames(PDF_quantiles_out) <- c("time", "mean", "lower", "median", "upper", "index_gt")
  
  n_gt <- length(gt)
  
  for(j in 1:n_gt){
    i <- gt[j]
    PDF_quantiles <- as.data.frame(matrix(data = NA, nrow = n_t, ncol = 6))
    colnames(PDF_quantiles) <- c("time", "mean", "lower", "median", "upper", "index_gt")
    PDF_quantiles[,1] <- times
    PDF_quantiles[,6] <- rep(i, n_t)
    
    shape <- matrix(rep(fit_$shape_total_S[,i], n_t), ncol = n_t)
    rate <- matrix(rep(fit_$rate_total_S[,i], n_t), ncol = n_t)
    mu <- matrix(rep(fit_$mu[,i], n_t), ncol = n_t)
    k <- matrix(rep(fit_$k, n_t), ncol = n_t)
    times_in <- matrix(rep(times, n_iter), nrow = n_iter, byrow = TRUE)
    PDF <- EIP_PDF(times_in, a = shape, 
                   b = rate,
                   mu = mu,
                   k = k)
    EIP_MCMC_list[[j]] <- PDF
    PDF_t <- t(PDF)
    PDF_quantiles[1:n_t,2] <- apply(PDF_t, 1, mean)
    qs <- (apply(PDF_t, 1, quantile, probs = c(0.025, 0.5, 0.975)))
    PDF_quantiles[1:n_t, 3] <- qs[1,]
    PDF_quantiles[1:n_t, 4] <- qs[2,]
    PDF_quantiles[1:n_t, 5] <- qs[3,]
    PDF_quantiles$temp <- rep(temp[j], nrow(PDF_quantiles))
    PDF_quantiles_out <- rbind(PDF_quantiles_out, PDF_quantiles)
    rm(list = c("PDF_quantiles"))
  }
  return(list("MCMC_all" = EIP_MCMC_list, "PDF_quantiles" = PDF_quantiles_out))
}


gen_quantiles <- function(df, temps){
  out <- as.data.frame(matrix(NA, nrow = length(temps), ncol = 5))
  
  for(i in 1:length(temps)){
    out[i,] <- c(temps[i], quantile(df[,i], c(0.025, 0.5, 0.975)), mean(df[,i]))
  }
  colnames(out) <- c("temp", "lower", "median", "upper", "mean")
  return(out)
}

#####################
### PDF functions ###
#####################

# get the parameters to calculate the EIP PDF
get_EIP_params <- function(temp, species, params){
  if(species == "stephensi"){
    out <- data.frame(
      "rate_O" = temp * params$m_rate_oocyst + params$c_rate_oocyst,
      "mu" = params$mu_NB,
      "k" = params$k_NB,
      "shape_S" = params$shape_sporozoite,
      "rate_S" = params$rate_sporozoite,
      "shape_O" = params$shape_oocyst
    )
  } else if(species == "gambiae"){
    out <- data.frame(
      "shape_S" = temp^2 * params$a_shape_S + temp * params$b_shape_S + params$c_shape_S,
      "rate_S" = temp * params$m_rate_S + params$c_rate_S,
      "mu" = (1/(1+exp(-(params$a_mu * temp^2 + 
                           params$b_mu * temp + params$c_mu)))) * params$g_mu,
      "k" = params$k,
      "shape_O" = params$shape_O,
      "rate_O" = params$rate_O
    )
  }
  
  out <- out %>% mutate(
    mu_total_S = (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S),
    sigma_sq_S  = (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2),
    shape_total_S = mu_total_S^2 / sigma_sq_S,
    rate_total_S = mu_total_S / sigma_sq_S,
    mean_total_S = shape_total_S / rate_total_S,
  )
  return(out)
}

# function to calculate the mean EIP from the get_EIP function
calc_mean_EIP <- function(EIP_index, i){
  p <- data.frame("a" = EIP_index$shape_total_S[,i],
                  "b" = EIP_index$rate_total_S[,i],
                  "mu" = EIP_index$mu[,i],
                  "k" = EIP_index$k[,i]) %>% rowwise() %>% mutate(mean = tryCatch({integrate(mean_func, a = a, b = b,
                                                                                             mu = mu, k = k, lower = 0, upper = Inf)[[1]]},
                                                                                  error = function(e){
                                                                                    NA
                                                                                  }))
  return(c("mean" = mean(na.omit(p$mean)), quantile(na.omit(p$mean), c(0.025, 0.5, 0.975)),
           "n_na" = sum(is.na(p$mean))))
}

calc_mean_EIP_nls <- function(params, temp, temps){
  i <- which(temps == temp)
  
  a <- params$shape_total_S[,i]
  b <- params$rate_total_S[,i]
  mu <- params$mu[,i]
  k <- params$k 
  EIP_mean <- v.int_mean_func(a = a, b = b, mu = mu, k = k)
  return(EIP_mean)
}

calc_mean_EIP_single_temp <- function(EIP_index, i){
  p <- data.frame("a" = EIP_index$shape_total_S[,i],
                  "b" = EIP_index$rate_total_S[,i],
                  "mu" = EIP_index$mu[,i],
                  "k" = EIP_index$k[,i]) %>% mutate(mean = v.int_mean_func(a = a, b = b, mu = mu, k = k))
  
  return(c("mean" = mean(na.omit(p$mean)), quantile(na.omit(p$mean), c(0.025, 0.5, 0.975)),
           "n_na" = sum(is.na(p$mean))))
}

v.calc_mean_EIP_single_temp <- Vectorize(calc_mean_EIP_single_temp)

calc_var_EIP <- function(params, temp, temps){
  i <- which(temps == temp)
  #a <- temp^2 * params$a_shape_S + temp * params$b_shape_S + params$c_shape_S
  #b <- temp * params$m_rate_S + params$c_rate_S
  #mu <- (1 / (1 + exp(-(temp^2 * params$a_mu + temp * params$b_mu + params$c_mu)))) * params$g_mu
  a <- params$shape_total_S[,i]
  b <- params$rate_total_S[,i]
  mu <- params$mu[,i]
  k <- params$k 
  EIP_var <- v.int_var_func(a = a, b = b, mu = mu, k = k)
  return(EIP_var)
}

get_EIP_mean <- function(params, temps, iter){
  shape_S <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  rate_S <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  mu <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  shape_O <- as.data.frame(replicate(length(temps), params$shape_O))
  rate_O <- as.data.frame(replicate(length(temps), params$rate_O))
  
  k <- as.data.frame(replicate(length(temps), params$k))
  
  for(i in 1:length(temps)){
    shape_S[,i] <- temps[i]^2 * params$a_shape_S + temps[i] * params$b_shape_S + params$c_shape_S
    rate_S[,i] <- temps[i] * params$m_rate_S + params$c_rate_S
    mu[,i] <- (1 / (1 + exp(-(temps[i]^2 * params$a_mu + temps[i] * params$b_mu + params$c_mu)))) * params$g_mu
    #mu[,i] <- params$a_mu * exp(-((temps[i] - params$b_mu)^2 / (2 * params$c_mu)^2))
  }
  
  mu_total_S <- (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S)
  sigma_sq_S  <- (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2)
  shape_total_S <- mu_total_S^2 / sigma_sq_S;
  rate_total_S <- mu_total_S / sigma_sq_S;
  mean_total_S <- shape_total_S / rate_total_S
  
  EIP_mean <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  for(i in 1:length(temps)){
    EIP_mean[,i] <- v.int_mean_func(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i])
  }
  
  return(list("EIP_mean" = EIP_mean))
}



###########################################################################################################
##### calculates the probability one EIP PDF is greater than another for a single posterior iteration #####
###########################################################################################################

calc_pr <- function(o_g, o_s, iter_g, iter_s){
  df <- data.frame("g" = Inv_EIP_CDF(a = o_g[iter_g,"shape_total_S"], 
                                     b = o_g[iter_g, "rate_total_S"], 
                                     mu = o_g[iter_g, "mu"], 
                                     k = o_g[iter_g, "k"], 
                                     p = runif(10000)),
                   "s" = Inv_EIP_CDF(a = o_s[iter_s,"shape_total_S"], 
                                     b = o_s[iter_s, "rate_total_S"], 
                                     mu = o_s[iter_s, "mu"], 
                                     k = o_s[iter_s, "k"], 
                                     p = runif(10000))) %>% 
    mutate(p = ifelse(g > s, 1, 0))
  out <- sum(df$p)/nrow(df)
  return(out)
}

# calculate the probability the difference is greater
calc_pr_val <- function(o_g, o_s, iter_g, iter_s){
  df <- data.frame("g" = Inv_EIP_CDF(a = o_g[iter_g,"shape_total_S"], 
                                     b = o_g[iter_g, "rate_total_S"], 
                                     mu = o_g[iter_g, "mu"], 
                                     k = o_g[iter_g, "k"], 
                                     p = runif(10000)),
                   "g_2" = Inv_EIP_CDF(a = o_g[iter_g,"shape_total_S"], 
                                       b = o_g[iter_g, "rate_total_S"], 
                                       mu = o_g[iter_g, "mu"], 
                                       k = o_g[iter_g, "k"], 
                                       p = runif(10000)),
                   "s" = Inv_EIP_CDF(a = o_s[iter_s,"shape_total_S"], 
                                     b = o_s[iter_s, "rate_total_S"], 
                                     mu = o_s[iter_s, "mu"], 
                                     k = o_s[iter_s, "k"], 
                                     p = runif(10000))) %>% 
    mutate(d_g = abs(g-g_2),
           d_s = abs(g - s),
           p = ifelse(d_g >= d_s, 1, 0))
  out <- sum(df$p)/nrow(df)
  return(out)
}

# calculate the probability one is greater than the other by sampling
run_pr <- function(temp, p_s, p_g, scaled_temps_s, temps_s, scaled_temps_g, temps_g, s_one = "stephensi",
                   s_two = "gambiae"){
  temp <- round(temp, digits = 3)
  temps_s <- round(temps_s, digits = 3)
  temps_g <- round(temps_g, digits = 3)
  print(temp)
  # set the seed 
  set.seed(12345)
  out_s <- get_EIP_params(scaled_temps_s[which(temps_s == temp)], 
                          s_one,
                          p_s)
  
  out_g <- get_EIP_params(scaled_temps_g[which(temps_g == temp)], 
                          s_two,
                          p_g)
  
  all <- data.frame(
    "i_g" = sample(seq(1, nrow(out_g), 1), 10000, replace = FALSE), # 10000 samples
    "i_s" = sample(seq(1, nrow(out_s), 1), 10000, replace = FALSE)
  ) %>% rowwise() %>% 
    mutate(
      p = calc_pr(o_g = out_g, o_s = out_s, iter_g = i_g, iter_s = i_s)
    )
  return(data.frame("mean" = mean(all$p),
                    "median" = median(all$p),
                    "lower" = quantile(all$p, c(0.025))[[1]],
                    "upper" = quantile(all$p, c(0.975))[[1]]))
}

# calculate the probability the difference is greater by sampling
run_pr_val <- function(temp, p_s, p_g, scaled_temps_s, temps_s, scaled_temps_g, temps_g, s_one = "stephensi",
                   s_two = "gambiae"){
  
  temp <- round(temp, digits = 3)
  temps_s <- round(temps_s, digits = 3)
  temps_g <- round(temps_g, digits = 3)
  
  print(temp)
  # set the seed 
  set.seed(12345)
  out_s <- get_EIP_params(scaled_temps_s[which(temps_s == temp)], 
                          s_one,
                          p_s)
  
  out_g <- get_EIP_params(scaled_temps_g[which(temps_g == temp)], 
                          s_two,
                          p_g)
  
  all <- data.frame(
    "i_g" = sample(seq(1, nrow(out_g), 1), 10000, replace = FALSE), # 10000 samples # should this be replaced = TRUE - no because this is uniform distribution?
    "i_s" = sample(seq(1, nrow(out_s), 1), 10000, replace = FALSE)
  ) %>% rowwise() %>% 
    mutate(
      p = calc_pr_val(o_g = out_g, o_s = out_s, iter_g = i_g, iter_s = i_s)
    )
  return(data.frame("mean" = mean(all$p),
                    "median" = median(all$p),
                    "lower" = quantile(all$p, c(0.025))[[1]],
                    "upper" = quantile(all$p, c(0.975))[[1]]))
}

calc_PDF <- function(p_x, t_ = t){
  
  g_PDF <- as.data.frame(t(mapply(EIP_PDF, 
                                  a = as.list(p_x[, "shape_total_S"]),
                                  b = as.list(p_x[, "rate_total_S"]),
                                  mu = as.list(p_x[, "mu"]),
                                  k = as.list(p_x[, "k"]),
                                  MoreArgs = list(t = t_)))) # %>% mutate(t = seq(0, 30, 0.1)) 
  
  return(data.frame("t" = t_, 
                    "mean" = sapply(g_PDF, mean),
                    "median" = sapply(g_PDF, median),
                    "lower" = sapply(g_PDF, quantile, probs = c(0.025)),
                    "upper" = sapply(g_PDF, quantile, probs = c(0.975))))
}

##############################################################################################
##### function that takes the EIP PDF parameters (from the MCMC) and discretizes the PDF #####
##############################################################################################

discretize_EIP_multi <- function(nc, shape, rate, mu, k, lo = 0.00001, up = 0.99999){
  
  n_iter <- length(shape)
  
  ps <- seq(lo, up, length.out = (nc + 1))
  shape_ <- matrix(rep(shape, nc + 1), ncol = nc + 1)
  rate_ <- matrix(rep(rate, nc + 1), ncol = nc + 1)
  mu_ <- matrix(rep(mu, nc + 1), ncol = nc + 1)
  k_ <- matrix(rep(k, nc + 1), ncol = nc + 1)
  
  ps_ <- matrix(rep(ps, nrow(shape_)), ncol = (nc + 1), byrow = TRUE)
  t_ps <- Inv_EIP_CDF(a = shape_, b = rate_, mu = mu_, k = k_, ps_)
  
  EIP_vals <- as.data.frame(array(NA, dim = c(n_iter * nc, 7)))
  colnames(EIP_vals) <- c("n", "p_", "lower", "upper", "p", "comp", "iteration")
  EIP_vals$lower <- as.vector(t(t_ps[,-(nc + 1)])) # lower
  EIP_vals$upper <- as.vector(t(t_ps[,-1])) # upper
  EIP_vals$comp <- rep(seq(1, nc, 1), n_iter)
  EIP_vals$iteration <- rep(seq(1, n_iter, 1), each = nc)
  EIP_vals$p <- rep(1/nc, n_iter * nc)
  EIP_vals$shape <- rep(shape, each = nc)
  EIP_vals$rate <- rep(rate, each = nc)
  EIP_vals$mu <- rep(mu, each = nc)
  EIP_vals$k <- rep(k, each = nc)
  
  EIP_vals <- EIP_vals %>% rowwise() %>% mutate(p_ = integrate(EIP_PDF, 
                                                               lower = lower, 
                                                               upper = upper,
                                                               a = shape, b = rate,
                                                               mu = mu, k = k)[[1]]) %>% 
    mutate(n = integrate(mean_func, lower = lower, upper = upper,
                         a = shape, b = rate,
                         mu = mu, k = k)[[1]]/p_)
  
  # gives the value as the mean
  EIP_vals_sum <- EIP_vals %>% group_by(comp) %>% summarise(n = mean(n),
                                                            lower_n = quantile(n, 0.025)[[1]],
                                                            upper_n = quantile(n, 0.975)[[1]],
                                                            lower = mean(lower),
                                                            lower_lower = quantile(lower, 0.025)[[1]],
                                                            upper_lower = quantile(lower, 0.975)[[1]],
                                                            upper = mean(upper),
                                                            lower_upper = quantile(upper, 0.025)[[1]],
                                                            upper_upper = quantile(upper, 0.975)[[1]])
  
  EIP_vals_sum$p = rep(1/nc, nrow(EIP_vals_sum))
  
  return(list("EIP_vals" = as.data.frame(EIP_vals),
              "EIP_vals_sum" = as.data.frame(EIP_vals_sum)))
}

# change to probability of EIP being different days
discretize_EIP_multi_t <- function(diff = 1, shape, rate, mu, k, lo = 0.00001, up = 0.99999){
  
  n_iter <- length(shape)
  
  # limits
  lo_ps <- min(floor(Inv_EIP_CDF(a = shape, b = rate, mu = mu, k = k, lo)))
  up_ps <- max(ceiling(Inv_EIP_CDF(a = shape, b = rate, mu = mu, k = k, up)))
  
  ps <- seq(lo_ps, up_ps, diff)
  
  nc <- length(ps)
  
  EIP_vals <- as.data.frame(array(NA, dim = c(n_iter * (nc-1), 7)))
  
  colnames(EIP_vals) <- c("n", "p_", "lower", "upper", "p", "comp", "iteration")
  
  EIP_vals$lower <- rep(ps[-nc], n_iter) # lower
  EIP_vals$upper <- rep(ps[-1]) # upper
  EIP_vals$comp <- rep(seq(1, (nc-1), 1), n_iter)
  EIP_vals$iteration <- rep(seq(1, n_iter, 1), each = (nc-1))
  
  #EIP_vals$p <- rep(1/nc, n_iter * nc)
  EIP_vals$shape <- rep(shape, each = (nc-1))
  EIP_vals$rate <- rep(rate, each = (nc-1))
  EIP_vals$mu <- rep(mu, each = (nc-1))
  EIP_vals$k <- rep(k, each = (nc-1))
  
  EIP_vals <- EIP_vals %>% rowwise() %>% mutate(p = integrate(EIP_PDF, 
                                                              lower = lower, 
                                                              upper = upper,
                                                              a = shape, b = rate,
                                                              mu = mu, k = k)[[1]]) %>% 
    mutate(n = integrate(mean_func, lower = lower, upper = upper,
                         a = shape, b = rate,
                         mu = mu, k = k)[[1]]/p)
  
  # gives the value as the mean
  EIP_vals_sum <- EIP_vals %>% group_by(comp, lower, upper) %>% summarise(n = mean(n),
                                                                          lower_n = quantile(n, 0.025)[[1]],
                                                                          upper_n = quantile(n, 0.975)[[1]],
                                                                          p = mean(p),
                                                                          lower_p = quantile(p, 0.025)[[1]],
                                                                          upper_p = quantile(p, 0.975)[[1]]) %>% 
    mutate(n_comp = max(comp))
  
  return(list("EIP_vals" = as.data.frame(EIP_vals),
              "EIP_vals_sum" = as.data.frame(EIP_vals_sum)))
}

##########################################
### helper functions to run the models ###
##########################################
# calculates the sporozoite prevalence for each compartment for the discrete model at equilibrium
# where the EIP is not stratified by day
# human prevalence included as parameter
run_ind_comp <- function(EIP_vals_sum, mu, temp_){
  df <- cbind(subset(EIP_vals_sum, temp == temp_)[,c("comp", "n", "p", "temp", "height", "comp_total")], 
              data.frame("mu" = rep(mu, EIP_vals_sum[1, "comp_total"]))) %>% 
    rowwise() %>% 
    mutate(
      s_disc = single_discrete_eq_I_m(a = 1/3,
                                      c = 0.5,
                                      I_h = 0.3,
                                      mu = 1/mu,
                                      n = n,
                                      p = p)
    )
  return(df)
}

# calculates the probability of surviving the EIP (EIP not stratified by day)
run_ind_comp_exp <- function(EIP_vals_sum, mu, temp_){
  df <- cbind(subset(EIP_vals_sum, temp == temp_)[,c("comp", "n", "p", "temp", "height", "comp_total")], 
              data.frame("mu" = rep(mu, EIP_vals_sum[1, "comp_total"]))) %>% 
    rowwise() %>% 
    mutate(
      s_disc = exp(-n*mu)*p
    )
  return(df)
}

# calculates the total equilibrium I given the discrete model
# human prevalence included as parameter
run_ind_comp_all <- function(EIP_vals_sum, mu, temp_){
  EIP <- subset(EIP_vals_sum, temp == temp_)[,c("comp", "n", "p", "temp", "height", "comp_total")]
  return(discrete_eq_I_m(a = 1/3,
                         c = 0.5,
                         I_h = 0.3,
                         mu = 1/mu,
                         n = EIP[,"n"],
                         p = EIP[,"p"]
  ))
}

# calculates the probability of surviving the EIP (EIP stratified by day)
run_ind_comp_exp_day <- function(EIP_vals_sum, mu, temp_){
  df <- subset(EIP_vals_sum_t, temp == temp_)[,c("comp", "n", "p", "temp", "comp_total", "lower", "upper")] %>%
    mutate(mu = mu) %>% 
    rowwise() %>% 
    mutate(
      s_disc = exp(-n*mu)*p
    )
  return(df)
}

# calculates the probability of surviving the EIP for each EIP (EIP stratified by day)
# human prevalence not included as a parameter
run_ind_comp_day <- function(EIP_vals_sum, mu_in, temp_){
  df <- subset(EIP_vals_sum_t, temp == temp_)[,c("comp", "n", "p", "temp", "comp_total", "lower", "upper")] 
  n_all <- df[,"n"]
  p_all <- df[,"p"]
  df <- df %>%
    mutate(mu = mu_in) %>%
    rowwise() %>% 
    mutate(
      s_disc = single_discrete_eq_I_m_full(a = 1/3,
                                           b = 0.35, 
                                           c = 0.35, 
                                           m = 20, 
                                           r = 1/20,
                                           mu = 1/mu,
                                           n = n_all,
                                           p = p_all,
                                           n_i = n,
                                           p_i = p)
    )
  return(df)
}

#############################################
### functions for DTR SMFA data wrangling ###
#############################################

process_prevalence_data <- function(data){
  out <- data.frame("Experiment" = integer(),
                    "Cup" = integer(),
                    "DPI" = double(),
                    "presence" = integer(),
                    "temp" = double(),
                    "DTR" = double(),
                    "gametocytemia" = double())
  
  for(i in 1:nrow(data)){
    n_inf <- data[i, "Infected"]
    n_un <- data[i, "Dissected"] - data[i, "Infected"]
    if(n_inf > 0){
      place_inf <- data.frame("Experiment" = rep(data[i, "Experiment"], n_inf),
                              "Cup" = rep(data[i, "Cup.ID"], n_inf),
                              "DPI" = rep(data[i, "DPI"], n_inf),
                              "presence" = rep(1, n_inf),
                              "temp" = rep(data[i, "temp"], n_inf),
                              "DTR" = rep(data[i, "DTR"], n_inf),
                              "gametocytemia" = rep(data[i, "gametocytemia"], n_inf))
      out <- rbind(out, place_inf)
    }
    if(n_un > 0){
      place_un <- data.frame("Experiment" = rep(data[i, "Experiment"], n_un),
                             "Cup" = rep(data[i, "Cup.ID"], n_un),
                             "DPI" = rep(data[i, "DPI"], n_un),
                             "presence" = rep(0, n_un),
                             "temp" = rep(data[i, "temp"], n_un),
                             "DTR" = rep(data[i, "DTR"], n_un),
                             "gametocytemia" = rep(data[i, "gametocytemia"], n_un))
      out <- rbind(out, place_un)
    }
  }
  
  # checking the numbers match up
  if(nrow(out) == sum(data$Dissected) &  sum(data$Infected) == sum(out$presence)){
    return(out)
  } else{
    return(print("numbers do not match up"))
  }
}

generate_prevalence_DTR <- function(data){
  totals <- unique(data[,c("DPI","index_temp", "temp", "DTR", "bt", "study")])
  
  for(i in 1:nrow(totals)){
    totals[i,"sample"] <- length(which(data[,"DPI"] == totals[i, "DPI"]
                                       & data[,"index_temp"] == totals[i, "index_temp"] &
                                         data[,"DTR"] == totals[i, "DTR"] &
                                         data[,"bt"] == totals[i, "bt"]  &
                                         data[,"study"] == totals[i, "study"]))
    
    totals[i,"positive"] <- length(which(data[,"DPI"] == totals[i, "DPI"] 
                                         & data[,"presence"] > 0 & data[,"index_temp"] == totals[i, "index_temp"] &
                                           data[,"DTR"] == totals[i, "DTR"] &
                                           data[,"bt"] == totals[i, "bt"] &
                                           data[,"study"] == totals[i, "study"]))
  }
  
  rm(i)
  
  totals <- mutate(totals, prevalence = positive / sample) # prevalence
  totals <- mutate(totals, lower = prevalence - (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 5% CI
  totals <- mutate(totals, upper = prevalence + (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 95% CI
  totals[which(totals$lower < 0), "lower"] <- 0 # preventing the lower confidence interval being below 0
  
  return(totals)
}

index_fun <- function(.df, variable, v_out){
  u_g <- unique(.df[, variable])
  for(i in 1:length(u_g)){
    .df[which(.df[, variable] == u_g[i]), v_out] <- i
  }
  return(.df)
}

# Parton Logan temperature model
day_temp <- function(t, Tmin, Tmax, D, p){
  return(Tmin + (Tmax - Tmin) * sin(pi * (t - 12 + D/2)/(D + 2 * p))) # t - 12 + D/2 = t - 6, so start at 6
}

night_temp <- function(t, Tmin, Tset, N, NTC, tset){
  return((Tmin - Tset * exp(-N/NTC) + (Tset - Tmin) * exp(-(t - tset)/NTC))/(1 - exp(-N/NTC)))
}

# starts at the start of feeding (18pm) and generates the temperature each hour post feeding for 24 hours
# adapted so biting time can be added
gen_temp <- function(mid, DTR, n_days, bt, d = 0.1){
  
  one_day <- c(night_temp(seq(18, 18 + 12 - d, d), mid - DTR/2, 
                          day_temp(18, mid - DTR/2, mid + DTR/2, D = 12, p = 1.5),
                          N = 12, NTC = 4, tset = 18),
               day_temp(seq(6, 18, d), mid - DTR/2, mid + DTR/2, D = 12, p = 1.5)
  )
  
  t <- c(seq(12, 23.9, d), seq(0, 12, d))
  
  n <- length(one_day)
  
  if(bt == 12){
    return(
      c(rep(one_day[-n], n_days - 1), one_day)
    )
  } else{
    one_day <- one_day[-n]
    i_bt <- which(t == bt)[1]
    one_day_s <- one_day[seq(1,i_bt)]
    one_day_f <- one_day[seq(i_bt,length(one_day))]
    one_day <- c(one_day_f, one_day_s)
    return(
      c(rep(one_day[-n], n_days - 1), one_day)
    )
  }
}

v.gen_temp <- Vectorize(gen_temp)

