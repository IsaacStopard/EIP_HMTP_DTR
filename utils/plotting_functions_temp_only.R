get_EIP_params_temp <- function(temp, params){
  out <- data.frame(
    "shape_S" = temp^2 * params$a_shape_S + temp * params$b_shape_S + params$c_shape_S,
    "rate_S" = temp * params$m_rate_S + params$c_rate_S,
    "mu" = (1/(1+exp(-(params$a_mu * temp^2 + 
                         params$b_mu * temp + params$c_mu)))) * params$g_mu,
    "k" = params$k,
    "shape_O" = params$shape_O,
    "rate_O" = params$rate_O
  ) %>% mutate(
    mu_total_S = (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S),
    sigma_sq_S  = (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2),
    shape_total_S = mu_total_S^2 / sigma_sq_S,
    rate_total_S = mu_total_S / sigma_sq_S,
    mean_total_S = shape_total_S / rate_total_S,
  )
  return(out)
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
