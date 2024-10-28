# Comparison of the Erlang distribution variance to the EIP PDF variances (with the same means) across a range of temperatures
# Author: Isaac J Stopard
# Version: 1.0
# Last updated: 15/07/2024
# Notes: 

rm(list = ls())

source(file = "utils/functions_temp_only.R"); source(file = "utils/model_functions.R"); 
source(file = "utils/vector_comp_functions.R");
source(file = "read_libraries_data.R")

#################
### functions ###
#################

# calculates the variance of an Erlang distribution with a different mean value (in this case the same as EIP50)

get_var <- function(comp, temp, params = params_temp){
  rate <- comp / gen_quantiles(na.omit(get_EIP_mean(params, temps = (temp - mean_temp_g) / sd_temp_g, 10000)$EIP_mean), temps = temp)$median
  return(comp / rate^2)
}

#################################
### temperature dependent EIP ###
#################################
e_temps <- c(17, 18, 19,
               20, 21, 23, 25,
               27, 28, 29, 30)

scaled_e_temps <- (e_temps - mean_temp_g) / sd_temp_g

### finding the number of compartments that give a similar variance to the observed data
# for each temperature and iteration this calculates mean and variance of the EIP PDF
# then calculates the mean mean and mean variance of all the MCMC iterations

m_df <- as.data.frame(sapply(e_temps, calc_mean_EIP_nls, params = params_temp, temps = e_temps))
colnames(m_df) <- e_temps

v_df <- as.data.frame(sapply(e_temps, calc_var_EIP, params = params_temp, temps = e_temps))
colnames(v_df) <- e_temps

var_df <- cbind(m_df %>% gather(key = temp, value = mean) %>% na.omit() %>% group_by(temp) %>% summarise(m = mean(mean)),
                v_df %>% gather(key = temp, value = var_) %>% na.omit() %>% group_by(temp) %>% summarise(v = mean(var_)))[,-c(1,3)] %>%
  mutate(temp = e_temps)
 
data_plot_v <- v_df %>% gather(key = temp, value = var_) %>% na.omit() %>% group_by(temp) %>% summarise(mean = mean(var_),
                                                                                                      median = median(var_),
                                                                                                      lower = quantile(var_, probs = c(0.1))[[1]],
                                                                                                      upper = quantile(var_, probs = c(0.9))[[1]]) %>%
  mutate(temp = e_temps)
  
### Non-linear squares to find compartment number for different EIP_X values
nls(formula = v ~ get_var(comp, temp), data = var_df, start = list(comp = 20), control = list(tol = 1e-03), trace = TRUE)

plot_v_df <- data.frame("temp" = seq(17, 30, 0.5)) %>% #
  mutate("47" = get_var(47, temp))
colnames(plot_v_df)[2] <- "v"

######################################
##### running for the SMFA model #####
######################################

##### delta values

delta_df <- vector(mode = "list", length = length(scaled_e_temps))

for(i in 1:length(delta_df)){
  delta_df[[i]] <- cbind(data.frame("temp" = e_temps[i], "scaled_temp" = scaled_e_temps[i]),
                         as.data.frame(t(gen_delta(fit = fit, temp = scaled_e_temps[i]))))
}

delta_df <- bind_rows(delta_df)
colnames(delta_df) <- c("temp", "scaled_temp", "mean", "median", "lower", "upper")

# loading in the data from the constant temperature SMFAs
s_data_C <- read.csv(file = "data/processed/ES_new_constant_temp_spz_processed.csv") %>% 
  mutate(gametocytemia = round(gametocytemia, digits = 5), bt = 12, study = 1) %>% index_fun(variable = "gametocytemia", v_out = "index_g") %>% filter(index_g!=1)

s_totals_c <- generate_prevalence_DTR(s_data_C[,-c(1)] %>% index_fun(variable = "temp", v_out = "index_temp"))

# specifying the credible interval HMTP estimates that correspond to the credible interval mean EIP estimates
unique_t <- replicate(3,unique(s_totals_c[,c("temp", "DTR")]),  simplify = FALSE)
for(i in 1:3){
  unique_t[[i]]$p_ind <- i
  if(i == 1){
    unique_t[[i]]$delta <- delta_df[match(unique_t[[i]][, "temp"], delta_df$temp), "lower"]
  } else if(i == 2){
    unique_t[[i]]$delta <- delta_df[match(unique_t[[i]][, "temp"], delta_df$temp), "median"]
  } else{
    unique_t[[i]]$delta <- delta_df[match(unique_t[[i]][, "temp"], delta_df$temp), "upper"]
  }
}
unique_t <- bind_rows(unique_t)
unique_t$bt = 12
unique_t$c = 1

max_time <- 50
l <- nrow(unique_t)
temp_fun <- vector(mode = "list", length = l)
# creating functions that describe the changes in temperature with time - should be no changes
for(i in 1:l){
  temp_fun[[i]] <- approxfun(seq(0, 24 * max_time, 0.1)/24, 
                             gen_temp(unique_t[i, "temp"], unique_t[i, "DTR"], max_time, bt = 12))
}

# running the Erlang distribution model for each constant temperature
# unique_t is subset to only run for the median posterior values
out_mc <- as.data.frame(bind_rows(mapply(run_SMFA_model, temp = unique_t[,"temp"], DTR = unique_t[,"DTR"], p_ind = unique_t[,"p_ind"], 
                                         MoreArgs = list(bt = 12, unique_t = unique_t, c = 1), SIMPLIFY = FALSE)))

# combining the plots 
plot_df <- out_mc[,c("time", "s_prev", "temp", "DTR", "bt", "p_ind")] %>% spread(key = p_ind, value = s_prev)
colnames(plot_df) <- c("DPI", "temp", "DTR", "bt", "lower", "median", "upper")

### getting the sporozoite prevalence values predicted by mSOS
# posterior predictive distribution times
PPD_times_g <- readRDS(file = "data/PPD_times_g")
length_ppd_times <- length(PPD_times_g)

iterations <- 5500
warmup <- 3000
chains <- 4

g_fits <- run_prop_ppd_df(fit, "pooled", "S_prevalence_ppd", length_ppd_times, PPD_times_g, unique_temp = e_temps)

# subsetting the predicted values so they are only 2.5 days after the maximum collected day
m_DPI <- as.data.frame(s_totals_c %>% group_by(temp) %>% summarise(m = max(DPI)))

g_fits <- bind_rows(lapply(seq(1, nrow(m_DPI)),
       function(i, m_DPI, fit){
         subset(fit, temp == m_DPI[i,"temp"] & DPI < (m_DPI[i, "m"]+2.5))
       }, 
       m_DPI = m_DPI, 
       fit = g_fits))

plot_df <- bind_rows(lapply(seq(1, nrow(m_DPI)),
                            function(i, m_DPI, fit){
                              subset(fit, temp == m_DPI[i,"temp"] & DPI < (m_DPI[i, "m"]+2.5))
                            }, 
                            m_DPI = m_DPI, 
                            fit = plot_df))


### predictive accuracy
s_totals_c[,"m3"] <- plot_df[match(interaction(s_totals_c$DPI, s_totals_c$temp), interaction(plot_df$DPI, plot_df$temp)),"median"]
s_totals_c[,"mSOS"] <- g_fits[match(interaction(s_totals_c$DPI, s_totals_c$temp), interaction(g_fits$DPI, g_fits$temp)),"median"]

s_data_C[,"m3"] <- plot_df[match(interaction(s_data_C$DPI, s_data_C$temp), interaction(plot_df$DPI, plot_df$temp)),"median"]
s_data_C[,"mSOS"] <- g_fits[match(interaction(s_data_C$DPI, s_data_C$temp), interaction(g_fits$DPI, g_fits$temp)),"median"]

scaled_brier <- function(obs, pre){
  brier <- sum((obs - pre)^2) / length(obs)
  #brier_max <- mean(obs) * (1-mean(obs))
  brier_max <- sum((obs - mean(obs))^2) / length(obs) 
  return(1-brier/brier_max)
}

round(scaled_brier(s_data_C$presence, 
             s_data_C$m3), digits = 3)

round(scaled_brier(s_data_C$presence,
         s_data_C$mSOS), digits = 3)

# ROC curves
get_roc <- function(plot_df, i, dates, species){
  b <- get_pred_act(plot_df = plot_df, i = i, dates = dates, species = species)
  return(roc(b$y, b$predicted))
  
}

mSOS_roc <- roc(s_data_C$presence, s_data_C$mSOS)
m3_roc <- roc(s_data_C$presence, s_data_C$m3)

auc(mSOS_roc)
auc(m3_roc)

roc_curves <- ggroc(list("m3" = m3_roc,
                         "mSOS" = mSOS_roc),
               size = 1.75, alpha = 0.8) +
    theme_bw() +
    scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "", labels = c("model 3 (SMFA)", "mSOS")) +
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), size = 1, linetype = 2, col = "grey75", inherit.aes = FALSE) +
    ylab("Sensitivity") + xlab("Specificity") + theme(text = element_text(size = 18))
  
actual_fitted_plot <- ggplot(data = s_totals_c %>% gather(key = "model", value = "value", m3, mSOS), aes(x = prevalence, y = value, colour = model)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = value, col = model), size = 0.5, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, size = 1) +
  geom_point(size = 2) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = 1.5) +
  theme_bw() + theme(text = element_text(size = 15)) +
  xlab("Actual sporozoite prevalence") +
  ylab("Predicted sporozoite prevalence") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "", labels = c("model 3 (SMFA)", "mSOS")) +


png(file = "results/figures/var_plot.png", height = 725, width = 725)
plot_grid(
  plot_grid(ggplot() + 
  geom_pointrange(data = data_plot_v, aes(x = temp, y = mean, ymin = lower, ymax = upper), size = 0.85) + 
    geom_line(data = plot_v_df,
                     aes(x = temp, y = v), size = 1.5, alpha = 0.9, col = "#56B4E9") + theme_bw() +
  theme(text = element_text(size = 15)) + ylab("Variance in the EIP") + xlab("Temperature (°C)") + 
  scale_y_continuous(trans = "sqrt") +
  scale_x_continuous(breaks = seq(18, 30, 2)),
  
  roc_curves + theme(legend.position = "none"),
  ncol = 2, rel_widths = c(1, 1), labels = c("A", "B")),
  
  ggplot() +
  geom_pointrange(data = s_totals_c %>% mutate(temp = paste0(temp, "°C")), #temp %in% c(19, 23, 27)
                  aes(x = DPI, y = prevalence, ymin = lower, ymax = upper), 
                  alpha = 0.75, shape = 21, fill = "grey30", size = 0.75) + 
  facet_wrap(~temp, scales = "free_x") + 
  theme_bw() +
  geom_line(data =  rbind(g_fits[,c("DPI", "median", "temp")] %>% mutate(model = "mSOS",
                                                                         temp = paste0(temp, "°C")),
      plot_df[,c("DPI", "median", "temp")] %>% mutate(model = "model M3 (SMFA)",
                                                      temp = paste0(temp, "°C"))), 
      aes(x = DPI, y = median, col = model, linetype = model), 
      size = 1.8, alpha = 0.9) +
  theme(text = element_text(size = 15)) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "") +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  xlab("Days post infection") + 
  ylab("Sporozoite prevalence") +
  scale_y_continuous(labels = scales::percent),
  ncol = 1,
  labels = c("","C"), rel_heights = c(1, 1.35)
)
dev.off()



