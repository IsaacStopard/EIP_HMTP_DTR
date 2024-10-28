# Model functions for the sensitivity analysis of the impact of EIP (mean and variance) 
# and seasonality on sporozoite prevalence
# Author: Isaac J Stopard
# Version: 1.0 
# Last updated: 14/06/2021
# Notes: 

####################################
### classic model at equilibrium ###
####################################
### equilibrium of classic model
eq_I_m <- function(a, c, I_h, mu, n){
  return((a * c * I_h * exp(- n * mu))/ (a * c * I_h + mu)) 
}

v.eq_I_m <- Vectorize(eq_I_m)

eq_I_m_full <- function(a, b, c, mu, m, n, r){
  return((-r * mu + a^2 * b * c * m * exp(-n * mu)) / (a * b * m * (a * c + mu)))
}

v.eq_I_m_full <- Vectorize(eq_I_m_full)

eq_I_h <- function(a, b, c, mu, m, n, r){
  p <- exp(- n * mu)
  return((-r * mu + a^2 * b * c * m * p)/(a * c * (r + a * b * m * p)))
}

v.eq_I_h <- Vectorize(eq_I_h)

R0_M1 <- function(m, a, b, c, mu, n, r){
  return((m * a^2 * b * c * exp(-mu * n))/(r * mu))
}

v.R0_M1 <- Vectorize(R0_M1)

### discretized EIP PDF
discrete_eq_I_m <- function(a, c, I_h, mu, n, p){
  return(sum(exp(-n * mu) * p) * a * c * I_h /
           (a * c * I_h + mu))
}

discrete_eq_I_m_full <- function(a, b, c, mu, m, r, n, p){
  return((-r * mu + a^2 * b * c * m * sum(exp(-n * mu) * p)) / (a * b * m * (a * c + mu)))
}

single_discrete_eq_I_m <- function(a, c, I_h, mu, n, p){
  return(exp(-n * mu) * p * a * c * I_h /
           (a * c * I_h + mu))
}

single_discrete_eq_I_m_full <- function(a, b, c, mu, m, r, n, p, n_i, p_i){
  return((exp(-n_i * mu)*p_i*(-r * mu + a^2 * b * c * m * sum(exp(-n * mu) * p))) / (a * b * m * (a * c + mu) * sum(exp(-n * mu) * p))) 
}

v.single_discrete_eq_I_m_full <- Vectorize(single_discrete_eq_I_m_full)

# continuous PDF model
prop_surv <- function(t, mu, shape, rate, mu_PL, k){
  return(exp(-mu * t) * EIP_PDF(t, shape, rate, mu_PL, k))
}

continuous_eq_I_m <- function(a, c, I_h, mu, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return(a * c * I_h * p / (a * c * I_h + mu))
}

continuous_eq_I_m_full <- function(a, b, c, mu, m, r, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return((-r * mu + a^2 * b * c * m * p)/ (a * b * m * (a * c + mu)))
}

v.continuous_eq_I_m_full <- Vectorize(continuous_eq_I_m_full)

prop_surv_int <- function(mu, shape, rate, mu_PL, k){
  return(integrate(prop_surv, lower = 0, upper = Inf, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]])
}

v.prop_surv_int <- Vectorize(prop_surv_int)

continuous_eq_I_h_full <- function(a, b, c, mu, m, r, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return((-r * mu + a^2 * b * c * m * p)/(a * c * (r + a * b * m * p)))
}

v.continuous_eq_I_h_full <- Vectorize(continuous_eq_I_h_full)

R0 <- function(a, b, c, mu, m, r, n){
  return((m * a^2 * b * c * exp(-mu * n))/(r * mu))
}

v.R0 <- Vectorize(R0)

continuous_R0 <- function(a, b, c, mu, m, r, n, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return((m * a^2 * b * c * p)/(r * mu))
}
  
v.continuous_R0 <- Vectorize(continuous_R0)

discrete_R0 <- function(a, b, c, mu, m, r, n, p){
  return((m * a^2 * b * c * sum(exp(-n * mu) * p))/(r * mu))
}

v.discrete_R0 <- Vectorize(discrete_R0)
######################
### dynamic models ###
######################

# desolve models

# DDE model with a single EIP value
# no births
single_model_nb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    # EIP delay
    EIP_lag <- t - n_1
    if(EIP_lag < 0){ # this needs to be the time of the first
      state_lag <- rep(0,length(state)) # the initial values
    } else{
      state_lag <- lagvalue(EIP_lag)
    }
    
    ############
    ### DDEs ###
    ############
    
    M <- S + E + I
    
    dS <- -(a * c * I_h * S)  - (mu * S)
    
    dE <- (a * c * I_h * S) - (a * c * I_h * state_lag[1] * exp(-mu * n_1)) - (mu * E)
    
    dI <- (a * c * I_h * state_lag[1] * exp(-mu * n_1)) - (mu * I)
    
    list(c(dS, dE, dI))
  })
}

# DDE model with single EIP value
# constant mosquito births
single_model_b <- function(t, state, params, init_vals = state_single){
  with(as.list(c(state, params)),{
    EIP_lag <- t - n_1
    if(EIP_lag < 0){
      state_lag <- rep(0, length(state_vals)) # the initial values
    } else{
      state_lag <- lagvalue(EIP_lag)
    }
    
    dE_m <- (a * c * I_h * (1 - E_m - I_m)) - 
      (a * c * I_h * (1 - state_lag[1] - state_lag[2]) * exp(-mu * n_1)) -
      (mu * E_m)
    
    dI_m <- (a * c * I_h * (1 - state_lag[1] - state_lag[2]) * exp(-mu * n_1)) - (mu * I_m)
    
    list(c(dE_m, dI_m))
  })
}

# ODE model
# Erlang distribution model
# no mosquito births
multi_comp_model_nb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    dS <- - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}

# Erlang distribution EIP model with constant births
multi_comp_model_cb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    #M <- S + E_ + I
    S <- 1 - E_ - I
      
    dI_h <- a * b * m * I * (1 - I_h) - r * I_h
    
    #dS <- M * (beta + mu) - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(E_out, dI, dI_h))
  })
}

# model 3 - Erlang distribution model with time-varying EIP and constant births
model_3 <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    rate <- shape/EIP_fun_t(t)
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    S <- 1 - E_ - I
    #M <- S + E_ + I
    
    #dS <- M * (beta + mu) - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    dI_h <- a * b * m * I * (1 - I_h) - r * I_h
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(E_out, dI, dI_h))
  })
}

multi_comp_model_b <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    beta <- beta_fun(t)
    
    dS <- M * (beta + mu) - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}

multi_comp_model_eb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    beta <- beta_fun(t)
    
    dS <- beta - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}

# model 4
# SMFA model: no births
# needs to be on the log scale so the numbers aren't below zero
model_4 <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    if(v_EIP == TRUE){
      rate <- shape/EIP_fun[[EIP_CI]](temp_fun[[index]](t))
    }
    
    #rate <- shape/EIP_fun(temp_fun[[i]](t), a, b, c, t0)
    dU <- - mu * U
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    dE1 <- - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dU, E_out, dI))
  })
}

model_4_microclimate <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    if(v_EIP == TRUE){
      temp <- temp_fun[[index]](t + start_time_in)
      rate <- shape/EIP_fun[[EIP_CI]](temp)
    }
    
    #rate <- shape/EIP_fun(temp_fun[[i]](t), a, b, c, t0)
    dU <- - mu * U
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    dE1 <- - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dU, E_out, dI), temp = temp, rate = rate)
  })
}


#
model_4_simple <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    if(v_EIP == TRUE){
      rate <- shape/EIP_fun[[EIP_CI]](temp_fun(t))
    }
    
    #rate <- shape/EIP_fun(temp_fun[[i]](t), a, b, c, t0)
    dU <- - mu * U
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    dE1 <- - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dU, E_out, dI))
  })
}

# model with varying 
model_4_ht <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    
    rate <- shape/EIP_fun_l[[index_EIP]](temp_fun_l[[index_temp]](t))
    #rate <- shape/EIP_fun(temp_fun[[i]](t), a, b, c, t0)
    dU <- - mu * U
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    dE1 <- - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dU, E_out, dI))
  })
}


multi_comp_model <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    dS <- - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}



# extending to the multiple compartment version
multiple_model <- function(t, state, params){
  with(as.list(c(state, params)),{
    
    EIP_lag <- rep(NA, length_n)
    state_lag <- matrix(NA, nrow = length_n + 1, ncol = length_n + 1)
    for(i in 1:length_n){
      EIP_lag[i] <- t - get(paste0("n",i))
      if(EIP_lag[i] < 0){
        state_lag[i+1,] <- init_vals
      } else{
        state_lag[i+1,] <- lagvalue(EIP_lag[i]) # values of the state variables at previous times
      }
    }
    if(t > 0){
      state_lag[1,] <- lagvalue(t) # the current time
    } else{
      state_lag[1,] <- init_vals
    }
    
    for(i in 1:length_n){
      assign(paste0("dE_m",i),
             (a * c * I_h * (1 - sum(state_lag[1,c(1:(length_n + 1))])) * get(paste0("p",i))) -
               (a * c * I_h * (1 - sum(state_lag[i+1, c(1:(length_n + 1))])) * get(paste0("p",i)) * exp(-mu * get(paste0("n",i)))) -
               (mu * get(paste0("E_m",i))))
    }
    
    dI_m <- (a * c * I_h * (1 - sum(state_lag[2, c(1:(length_n + 1))])) * get(paste0("p",1)) * exp(-mu * get(paste0("n",1))))
    for(i in 2:length_n){
      dI_m <- dI_m + 
        (a * c * I_h * (1 - sum(state_lag[i+1, c(1:(length_n + 1))])) * get(paste0("p",i)) * exp(-mu * get(paste0("n",i))))
    }
    dI_m <- dI_m - (mu * I_m)
    
    # creating the list to return
    out <- rep(NA, length_n + 1)
    for(i in 1:length_n){
      out[i] <- get(paste0("dE_m",i))
    }
    out[length_n + 1] <- dI_m
    list(out)
  })
}

#######################
##### SMFA models #####
#######################

run_SMFA_model <- function(temp, 
                           DTR, 
                           bt,
                           p_ind = 2,
                           c, # if constant
                           t = seq(0, max_time, 0.1),
                           unique_t,
                           EIP_p = NULL,
                           v_EIP = TRUE,
                           rate = NULL){
  
  
  index <- which(unique_t$temp == temp & unique_t$DTR == DTR & unique_t$bt == bt & unique_t$p_ind == p_ind & unique_t$c == c)
  
  #params <- c(mu = 0.1, shape = 47, i = index, a = EIP_p[1], b = EIP_p[2], c = EIP_p[3], t0 = EIP_p[4])
  if(v_EIP == TRUE){
    params <- c(mu = 0.1, shape = 47, index = index, EIP_CI = p_ind, v_EIP = v_EIP)
  } else{
    params <- c(mu = 0.1, shape = 47, index = index, EIP_CI = p_ind, v_EIP = v_EIP, rate = rate)
  }
  
  delta <- unique_t[index,"delta"]
  #delta <- gen_delta(fit_temp, (mean(temp_fun[[index]](seq(0, delta_time, 0.1))) - 23.06936)/4.361642)
  state <- c(U = (1 - delta)*1, E = c(delta*1, rep(0, params[["shape"]]-1)),  I = 0)
  
  out_mc <- as.data.frame(ode(y = state, times = t, func = model_4,
                              parms = params)) %>% rowwise() %>%
    mutate(M = U + sum(c_across("E1":paste0("E",params[["shape"]]))) + I,
           s_prev = I / M,
           temp = unique_t[index, "temp"],
           DTR = unique_t[index, "DTR"],
           bt = unique_t[index, "bt"],
           p_ind = unique_t[index, "p_ind"])
  return(out_mc)
}

gen_delta <- function(fit, temp){
  placeholder <- (1/(1+exp(-(rstan::extract(fit, "a_delta")[[1]] * temp^2 + 
                               rstan::extract(fit, "b_delta")[[1]] * temp + 
                               rstan::extract(fit, "c_delta")[[1]])))) *
    (1/(1+exp(-(rstan::extract(fit, "a_delta_S")[[1]] * temp^2 + 
                  rstan::extract(fit, "b_delta_S")[[1]] * temp + 
                  rstan::extract(fit, "c_delta_S")[[1]]))))
  out <- c(mean(placeholder), quantile(placeholder, c(0.5, 0.025, 0.975)))
  names(out)[[1]] <- "mean"
  return(out)
}

# calculating the delta values over different time-periods
# could be the mean over different time-frames - so calculate the mean between a range?
run_diff_delta <- function(h,
                           u_l,
                           delta_vt = TRUE, # if delta is determined by variable temperature = TRUE
                           EIP_vt = TRUE, # if EIP is determined by the variable temperature = TRUE, previously ll_calc
                           delta_fun = "max",
                           #ll_calc = TRUE
                           unique_t_DTR,
                           fit,
                           mean_temp,
                           sd_temp,
                           v_EIP = TRUE,
                           d = 0.1){
  
  # rather than hours should this be the time until the same % of the sporogony is complete?
  # delta estimation
  if(delta_vt == TRUE){
    n_ <- c(rep(seq(1,u_l,1), 2)) # so that the delta values are the same for each replicate
  } else{
    n_ <- c(rep(seq((u_l+1),nrow(unique_t_DTR),1), 2))
  }
  
  unique_t_DTR$n <- n_
  
  if(delta_fun == "max"){
    unique_t_DTR <- unique_t_DTR %>% rowwise() %>% 
      mutate(m_temp = max(temp_fun[[n]](seq(0, h, d)/24)))
    unique_t_DTR <- as.data.frame(unique_t_DTR %>% mutate(delta = gen_delta(fit = fit, temp = (m_temp - mean_temp)/sd_temp)['50%'][[1]]))
  } else if(delta_fun == "min"){
    unique_t_DTR <- unique_t_DTR %>% rowwise() %>% 
      mutate(m_temp = min(temp_fun[[n]](seq(0, h, d)/24)))
    unique_t_DTR <- as.data.frame(unique_t_DTR %>% mutate(delta = gen_delta(fit = fit, temp = (m_temp - mean_temp)/sd_temp)['50%'][[1]]))
  } else if(delta_fun == "mean"){
    unique_t_DTR <- unique_t_DTR %>% rowwise() %>% 
      mutate(m_temp = mean(temp_fun[[n]](seq(0, h, d)/24)))
    unique_t_DTR <- as.data.frame(unique_t_DTR %>% mutate(delta = gen_delta(fit = fit, temp = (m_temp - mean_temp)/sd_temp)['50%'][[1]]))
  } else if(delta_fun == "mean_delta"){
    for(i in 1:nrow(unique_t_DTR)){
      p <- data.frame("m_temp" = temp_fun[[unique_t_DTR[i, "n"]]](seq(0, h, d)/24)) %>% rowwise() %>% mutate(t_ = (m_temp - mean_temp)/sd_temp,
                                                                                                               delta = gen_delta(fit = fit, temp = t_)['50%'][[1]])
      unique_t_DTR[i, "delta"] <- mean(p$delta)
    }
    
  } else if(delta_fun == "min_delta"){
    for(i in 1:nrow(unique_t_DTR)){
      p <- data.frame("m_temp" = temp_fun[[unique_t_DTR[i, "n"]]](seq(0, h, d)/24)) %>% rowwise() %>% mutate(t_ = (m_temp - mean_temp)/sd_temp,
                                                                                                               delta = gen_delta(fit = fit, temp = t_)['50%'][[1]])
      unique_t_DTR[i, "delta"] <- min(p$delta)
    }
  } else if(delta_fun == "max_delta"){
    for(i in 1:nrow(unique_t_DTR)){
      p <- data.frame("m_temp" = temp_fun[[unique_t_DTR[i, "n"]]](seq(0, h, d)/24)) %>% rowwise() %>% mutate(t_ = (m_temp - mean_temp)/sd_temp,
                                                                                                               delta = gen_delta(fit = fit, temp = t_)['50%'][[1]])
      unique_t_DTR[i, "delta"] <- max(p$delta)
    }
  }  else if(delta_fun == "geom_m_delta"){
    for(i in 1:nrow(unique_t_DTR)){
      p <- data.frame("m_temp" = temp_fun[[unique_t_DTR[i, "n"]]](seq(0, h, d)/24)) %>% rowwise() %>% mutate(t_ = (m_temp - mean_temp)/sd_temp,
                                                                                                               delta = gen_delta(fit = fit, temp = t_)['50%'][[1]])
      unique_t_DTR[i, "delta"] <- exp(mean(log(p$delta)))
    }
  }
  
  if(EIP_vt == TRUE){
    u_t_in <- unique_t_DTR[1:u_l,]
  } else{
    u_t_in <- unique_t_DTR[(u_l+1):nrow(unique_t_DTR),]
  }
  
  if(v_EIP == TRUE){
    out_mc <- as.data.frame(bind_rows(mapply(run_SMFA_model, 
                                             temp = u_t_in[,"temp"], 
                                             DTR = u_t_in[,"DTR"], 
                                             bt = u_t_in[,"bt"],
                                             p_ind = u_t_in[,"p_ind"], 
                                             c = u_t_in[,"c"], SIMPLIFY = FALSE,
                                             MoreArgs = list(unique_t = unique_t_DTR, v_EIP = TRUE)))) # p_ind = 2 so this model is with the median mean EIP 
  } else{
    out_mc <- as.data.frame(bind_rows(mapply(run_SMFA_model, 
                                             temp = u_t_in[,"temp"], 
                                             DTR = u_t_in[,"DTR"], 
                                             bt = u_t_in[,"bt"],
                                             p_ind = u_t_in[,"p_ind"], 
                                             c = u_t_in[,"c"], 
                                             rate = u_t_in[,"rate"],
                                             SIMPLIFY = FALSE,
                                             MoreArgs = list(unique_t = unique_t_DTR, v_EIP = FALSE))))
  }
  
  out_mc$temp <- ceiling(out_mc$temp)
  
  out_mc <- out_mc %>% mutate(DPI = time,
                              delta_vt = delta_vt,
                              EIP_vt = EIP_vt)
  
  unique_t_DTR <- unique_t_DTR[,!names(unique_t_DTR) %in% c("n", "m_temp", "delta")]
  rm(list = c("u_t_in"))
  return(out_mc)
}

# functions to calculate the maximum likelihoods
calc_ll_DTR <- function(s_totals_ = s_totals_l, out){
  index <- match(interaction(round(s_totals_$DPI, digits = 1), s_totals_$temp, s_totals_$bt, s_totals_$DTR), interaction(round(out$DPI, digits = 1), out$temp, out$bt, out$DTR))
  p <- out[index, "s_prev"]
  ll <- dbinom(x = s_totals_[,"positive"], size = s_totals_[,"sample"], prob = p, log = TRUE)
  return(ll)
}

calc_ll_all <- function(delta_vt_in, EIP_vt_in, delta_fun_in, min_h = 1, max_h = 24, s_h = 1, v_EIP_in = TRUE){
  likelihoods <- bind_rows(lapply(seq(min_h, max_h, s_h), function(x){
    out_fc <- NULL
    attempt <- 0
    while(is.null(out_fc) && attempt <= 5){
      attempt <- attempt + 1
      try(
        out_fc <- run_diff_delta(h = x,
                                 u_l = u_l,
                                 delta_vt = delta_vt_in,
                                 EIP_vt = EIP_vt_in,
                                 delta_fun = delta_fun_in,
                                 unique_t_DTR = unique_t_DTR,
                                 fit = fit,
                                 mean_temp = mean_temp,
                                 sd_temp = sd_temp,
                                 v_EIP = v_EIP_in)
      )
    }
    l <- data.frame("h" = x,
                    "ll" = sum(calc_ll_DTR(out = out_fc)),
                    "attempt" = attempt)
    rm(list = c("attempt", "out_fc"))
    return(l)
  }))
  return(likelihoods)
}

calc_ll_ml <- function(x, delta_vt_in, EIP_vt_in, delta_fun_in, v_EIP_in = TRUE){
  
  out_fc <- NULL
  attempt <- 0
  while(is.null(out_fc) && attempt <= 5){
    attempt <- attempt + 1
    try(
      out_fc <- run_diff_delta(h = x, 
                               u_l = u_l, 
                               delta_vt = delta_vt_in, 
                               EIP_vt = EIP_vt_in, 
                               delta_fun = delta_fun_in, 
                               unique_t_DTR = unique_t_DTR,  
                               fit = fit, 
                               mean_temp = mean_temp, 
                               sd_temp = sd_temp,
                               v_EIP = v_EIP_in)
    )
  }
  l <- sum(calc_ll_DTR(out = out_fc))
  rm(list = c("attempt", "out_fc"))
  return(l * -1) # because optim function
}



