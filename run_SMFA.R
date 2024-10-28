# DTR impact within standard membrane feeding assays 
# Author: Isaac J Stopard
# Version: 1.0
# Last updated: 06/08/2021
# Notes: 
rm(list = ls())

source(file = "utils/model_functions.R")
source(file = "read_libraries_data.R")

###########
### EIP ###
###########

# getting the EIP value
EIP_10 <- gen_quantiles(EIP_index_all$EIP_10, temps_all)
EIP_50 <- gen_quantiles(EIP_index_all$EIP_50, temps_all)
EIP_90 <- gen_quantiles(EIP_index_all$EIP_90, temps_all)

########################
### fitting the data ###
########################

# loading the data
# constant temperature data
s_data_C <- read.csv(file = "data/processed/ES_new_constant_temp_spz_processed.csv") %>% 
  mutate(gametocytemia = round(gametocytemia, digits = 5), bt = 12, study = 1) %>% index_fun(variable = "gametocytemia", v_out = "index_g") %>% filter(index_g!=1)

s_data_DTR <- process_prevalence_data(
  read.csv(file = "data/DTR/DTR_spz_data.csv")
) %>% 
  mutate(gametocytemia = round(gametocytemia, digits = 5),
         bt = 12, study = 1) %>% index_fun(variable = "gametocytemia", v_out = "index_g") %>% filter(index_g!=1)

# data from the other studies

w_data <- readRDS(file = "data/DTR/waite_low.rds") %>% mutate(data_source = "Waite et al 2019")

e_data <- readRDS(file = "data/DTR/Suh_NEE_data_fig_1.rds") %>% mutate(data_source = "Suh et al (2020)")

m_data <- readRDS(file = "data/DTR/murdock_a_g_spz.rds") %>% mutate(data_source = "Murdock et al (2016)")

s_data <- rbind(s_data_DTR %>% mutate(data_source = "novel"), w_data, e_data, m_data) %>% index_fun(variable = "temp", v_out = "index_temp")

s_totals_test <- generate_prevalence_DTR(s_data)

########################################
### generating the SMFA temperatures ###
########################################

u_t <- unique(s_totals_test[,c("temp", "DTR", "bt")])

# DTR temperature plots
# biting time at 18:00 or ZT=12
n_days <- 1
t <- seq(0, 24 * n_days, 0.1)
n_t <- length(t)

temps_DTR <- v.gen_temp(mid = u_t[,"temp"], DTR = u_t[,"DTR"], n_days, bt = u_t[,"bt"])
m_temps_DTR <- colMeans(temps_DTR)

t_plot_df <- vector(mode = "list", length = nrow(u_t))

for(i in 1:nrow(u_t)){
  t_plot_df[[i]] <- data.frame("mid" = rep(u_t[i, "temp"], n_t * 2),
                        "DTR" = c(rep(u_t[i, "DTR"], n_t * 2)),
                        "temp" = c(temps_DTR[,i],
                                   gen_temp(m_temps_DTR[i], 0, n_days, bt = 12)),
                        "DPI" = rep(t, 2),
                        "bt" = rep(u_t[i, "bt"], 2),
                        "mean" = c(rep("FALSE", n_t), rep("TRUE", n_t))
)
}

t_plot_df <- bind_rows(t_plot_df)

########################
##### delta values #####
########################

# HMTP

delta_df <- vector(mode = "list", length = length(scaled_temps))

for(i in 1:length(delta_df)){
  delta_df[[i]] <- cbind(data.frame("temp" = temps[i], "scaled_temp" = scaled_temps[i]),
                         as.data.frame(t(gen_delta(fit = fit, temp = scaled_temps[i]))))
}

delta_df <- bind_rows(delta_df)
colnames(delta_df) <- c("temp", "scaled_temp", "mean", "median", "lower", "upper")

###########################################################################
### running the model with time-varying temperature for all experiments ###
###########################################################################
max_time <- 70

# delta values - (1) mean temperature, (2) mean, min or max temperature over a range of times and (3) mean, min or max delta over a range of times
# EIP values - (1) mean temperature, (2) fluctuating temperature and (3) mean parasite development rate

# run four model combinations and calculate the likelihood for each
# 1: delta - 1 and EIP - 1
# 2: delta - 1 and EIP - 2
# 3: delta - 2 and EIP - 1
# 4: delta - 2 and EIP - 2
# 5: delta - 3 and EIP - 1
# 6: delta - 3 and EIP - 2
# 7: delta - 1 and EIP - 3
# 8: delta - 2 and EIP - 3
# 9: delta - 3 and EIP - 3

# when delta is determined by the time-dependent changes in the temperature this is repeated for the mean, min and max temperatures over are range of hours post infection

u_l <- nrow(u_t) # actual number of unique combinations

unique_t_DTR <- rbind(as.data.frame(u_t %>% mutate(p_ind = 2,
                                                   DTR_in = DTR,
                                                   c = 1)),
                      as.data.frame(u_t %>% mutate(p_ind = 2,
                               temp = m_temps_DTR[row_number()],
                               DTR_in = 0,
                               c = 2)))

# both the EIP models are run at the same time
# temperature functions
# uses a linear interpolation
l <- nrow(unique_t_DTR)
temp_fun <- vector(mode = "list", length = l)
for(i in 1:l){
  temp_fun[[i]] <- approxfun(seq(0, 24 * max_time, 0.1)/24, 
                             gen_temp(mid = unique_t_DTR[i, "temp"], 
                                      DTR = unique_t_DTR[i, "DTR_in"],
                                      max_time, 
                                      bt = unique_t_DTR[i, "bt"],
                                      d = 0.1))
}

# calculating the mean rates over the first 24 hours
rates <- rep(NA, u_l)
for(i in 1:u_l){
  rates[i] <- mean(47/EIP_fun[[2]](temp_fun[[i]](seq(0, 24, 0.1)/24)))
}

unique_t_DTR$rate <- c(rates, rates)

######################################
##### calculating the likelihood #####
######################################

# testing dataset
s_totals_l <- generate_prevalence_DTR(s_data %>% mutate(study = 1))

# models 1 and 2 do not have a variable delta so the time should not affect the likelihood - this is checked
# for each model (ll_1 and ll_2) the max, min and mean shouldn't effect the results
# model 1
ll_1 <- calc_ll_all(delta_vt_in = FALSE, EIP_vt_in = FALSE, delta_fun_in = "mean", min_h = 1, max_h = 21, s_h = 10)

# model 2
ll_2 <- calc_ll_all(delta_vt_in = FALSE, EIP_vt_in = TRUE, delta_fun_in = "mean", min_h = 1, max_h = 21, s_h = 10)

#saveRDS(list("ll_1" = ll_1, "ll_2" = ll_2), file = "data/ml_DTR_simulations_mean.rds")

ll_list <- readRDS(file = "data/ml_DTR_simulations_mean.rds")

# models 3 and 4 have variable delta values, so we estimate the time-frame over which this occurs using maximum likelihood
# model 3
# doing in parallel
mod_sims <- expand.grid(EIP_vt_in = c(FALSE, TRUE), delta_fun_in = c("max", "mean", "min", 
                                                                     "max_delta", "mean_delta", "min_delta")) %>% mutate(delta_vt_in = TRUE,
                                                                                                       v_EIP_in = TRUE)

cl <- makeCluster(4)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("mod_sims", "calc_ll_ml", "calc_ll_DTR", "model_4", "run_diff_delta", "run_SMFA_model", "EIP_fun", "temp_fun"))

ml_all <- foreach(i=1:nrow(mod_sims),
                               .packages = (.packages())
) %dopar% {
  tryCatch({optim(par = 8, calc_ll_ml, 
                  delta_vt_in = TRUE, 
                  EIP_vt_in = mod_sims[i, "EIP_vt_in"], 
                  delta_fun_in = mod_sims[i, "delta_fun_in"],
                  v_EIP_in = TRUE,
                  lower = 0.01, upper = 24, 
                  method = "Brent")},
  error = function(cond){
    return(NA)
  })
}
stopCluster(cl)

saveRDS(ml_all,      
        file = "data/ml_DTR_simulations_parallel.rds")

ml_all <- readRDS(file = "data/ml_DTR_simulations_parallel.rds")

# extracting the parameter values
for(i in 1:nrow(mod_sims)){
  mod_sims[i, "likelihood"] <- round(ml_all[[i]]$value * -1, digits = 2)
  mod_sims[i, "h_in"] <- floor(ml_all[[i]]$par/0.1)*0.1
  mod_sims[i, "h"] <- floor(ml_all[[i]]$par/0.1)*0.1
}

ml_in <- rbind(data.frame(EIP_vt_in = c(FALSE, TRUE),
                          delta_fun_in = c("mean", "mean"),
                          delta_vt_in = c(FALSE, FALSE),
                          v_EIP_in = c(TRUE, TRUE),
                          likelihood = c(round(ll_list$ll_1$ll[1], digits = 2), round(ll_list$ll_2$ll[1], digits = 2)),
                          h = rep("h*", 2),
                          h_in = rep(0, 2)
                          ),
               mod_sims) %>% 
  mutate(delta = case_when(delta_vt_in == FALSE ~ "a (mean temperature during first 24 hours)",
                           delta_vt_in == TRUE & delta_fun_in == "mean" ~ "b (mean)",
                           delta_vt_in == TRUE & delta_fun_in == "max" ~ "b (maximum)",
                           delta_vt_in == TRUE & delta_fun_in == "min" ~ "b (minimum)",
                           delta_vt_in == TRUE & delta_fun_in == "mean_delta" ~ "c (mean)",
                           delta_vt_in == TRUE & delta_fun_in == "max_delta" ~ "c (maximum)",
                           delta_vt_in == TRUE & delta_fun_in == "min_delta" ~ "c (minimum)"),
         EIP = case_when(EIP_vt_in == FALSE ~ "a (mean temperature)",
                         EIP_vt_in == TRUE ~ "b (fluctuating temperature)"))

# generating the model predictions
pred_out <- lapply(seq(1, nrow(ml_in)), function(i){
  
  out_fc <- NULL
  attempt <- 0
  while(is.null(out_fc) && attempt <= 5){
    attempt <- attempt + 1
    try(
      out_fc <- run_diff_delta(h = ml_in[i, "h_in"], 
                               u_l = u_l, 
                               delta_vt = ml_in[i, "delta_vt_in"], 
                               EIP_vt = ml_in[i, "EIP_vt_in"], 
                               delta_fun = ml_in[i, "delta_fun_in"], 
                               unique_t_DTR = unique_t_DTR,  
                               fit = fit,
                               mean_temp = mean_temp, 
                               sd_temp = sd_temp, 
                               v_EIP = ml_in[i, "v_EIP_in"])
    )
  }
  
  return(out_fc)
  
})

#saveRDS(pred_out, file = "data/pred_out_novel_data.rds")
pred_out <- readRDS(file = "data/pred_out.rds")

max_DPI <- left_join(unique(s_totals_l[,c("temp", "DTR", "bt")]), 
                     as.data.frame(s_totals_l %>% group_by(temp) %>% summarise(mDPI = max(DPI))), by = "temp")

for(j in 1:length(pred_out)){
  pred_out[[j]] <- bind_rows(lapply(seq(1, nrow(max_DPI)), function(i){
  df <- subset(pred_out[[j]], temp == max_DPI[i, "temp"] &
                 DTR == max_DPI[i, "DTR"] &
                 bt == max_DPI[i, "bt"])
  df <- df[-which(df$DPI > max_DPI[i, "mDPI"] + 5),]
  return(df)}
  ))
  }

scaled_brier <- function(obs, pre){
  brier <- sum((obs - pre)^2) / length(obs)
  #brier_max <- mean(obs) * (1-mean(obs))
  brier_max <- sum((obs - mean(obs))^2) / length(obs) 
  return(1-brier/brier_max)
}

# ROC curves
get_roc <- function(plot_df, i, dates, species){
  b <- get_pred_act(plot_df = plot_df, i = i, dates = dates, species = species)
  return(roc(b$y, b$predicted))
  
}

roc_curves <- vector(mode = "list", length = nrow(ml_in))

for(i in 1:nrow(ml_in)){
  pred <- pred_out[[i]]
  pred_vals <- pred[match(interaction(s_data$temp, s_data$DTR, s_data$bt, s_data$DPI), interaction(pred$temp, pred$DTR, pred$bt, pred$DPI)),"s_prev"]
  ml_in[i, "Brier"] <- round(scaled_brier(s_data$presence, pred_vals), digits = 2)
  roc_curves[[i]] <- roc(s_data$presence, pred_vals)
  ml_in[i, "AUC"] <- round(auc(roc_curves[[i]]), digits = 2)
}

ml_in$likelihood_r <- round(ml_in$likelihood, digits = 2)

#ml_in$EIP <- gsub("a (mean temperature)", "mean", ml_in$EIP)
#ml_in$EIP <- gsub("b (fluctuating temperature)", "fluctuating", ml_in$EIP)
write.csv(ml_in[,c("EIP", "delta", "h", "likelihood_r", "Brier", "AUC")], file = "results/likelihood_table.csv")

png(file = "results/figures/DTR_temperature_data.png", height = 700, width = 850)
ggplot(data = subset(t_plot_df %>% mutate(mid = paste0("Temperature: ", mid, "°C"),
                                            bt = paste0("Biting time: ZT",bt)), DPI <= 24 & mean == "FALSE"), aes(x = DPI, y = temp, col = factor(DTR))) + 
  geom_line(linewidth = 1.5, alpha = 0.95) +
    theme_bw() + scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 4)) +
    xlab("Hours post infection") + ylab("Temperature (°C)") +
    scale_colour_manual(values = c("#CC79A7", "#009E73", "#E69F00", "#56B4E9"), name = "DTR") +
    facet_wrap(~mid + bt) +
    scale_y_continuous(limits = c(12.5, 37.5), breaks = seq(15, 35, 5)) +
    theme(text = element_text(size = 15))
dev.off()

s_totals_test <- s_totals_test %>% mutate(data_source = case_when(study == 1 ~ "Novel",
                                                 study == 2 ~ "Waite et al (2019)",
                                                 study == 3 ~ "Suh et al (2020)",
                                                 study == 4 ~ "Murdock et al (2016)"))

s_totals_test$data_source <- factor(s_totals_test$data_source, levels = c("Novel", "Murdock et al (2016)", "Suh et al (2020)", "Waite et al (2019)"))

smfa_plot <- function(i, n = 2){
  ggplot() + 
    geom_pointrange(data = s_totals_test %>% mutate(temp = paste0("Temperature: ", temp, "°C"),
                                                 bt = paste0("Biting time: ZT",bt)), 
                    aes(x = DPI, y = prevalence, ymin = lower, ymax = upper, size = sample, 
                        fill = factor(DTR), shape = data_source), col = "black") +
    geom_line(data = pred_out[[i]] %>% mutate(temp = paste0("Temperature: ", temp, "°C"),
                                              bt = paste0("Biting time: ZT",bt)), 
              aes(x = DPI, y = s_prev, col = factor(DTR)), linewidth = 1.075) + 
    facet_wrap(~temp + bt, scales = "free_x", nrow = n) +
    scale_size_continuous(range = c(0.1, 1.0), name = "sample\nsize") + 
    theme_bw() + theme(text = element_text(size = 20)) +
    scale_colour_manual(values = c("#CC79A7", "#009E73", "#E69F00",  "#56B4E9"), name = "DTR") +
    scale_fill_manual(values = c("#CC79A7", "#009E73", "#E69F00", "#56B4E9"), name = "DTR") +
    scale_shape_manual(values = c(22, 21, 23, 24), name = "study") +
    ylab("Sporozoite prevalence") + xlab("Days post infection") +
    scale_y_continuous(labels = scales::percent)
}

bf_plot <- smfa_plot(4)

png(file = "results/figures/DTR_simulations_testing_data.png", height = 550, width = 1000)
bf_plot
dev.off()

png(file = "results/figures/DTR_simulations_all.png", height = 1400, width = 1550)
plot_grid(
  plot_grid(
  smfa_plot(2, n = 2) + theme(legend.position = "none", plot.title = element_text(face = "bold")) +
    ggtitle("model 2")+
    theme(title = element_text(size = 14)),
  
  smfa_plot(6, n = 2) + theme(legend.position = "none", plot.title = element_text(face = "bold")) +
    ggtitle("model 6") +
    theme(title = element_text(size = 14)),
  
  smfa_plot(8, n = 2) + theme(legend.position = "none", plot.title = element_text(face = "bold")) + 
    ggtitle("model 8") +
    theme(title = element_text(size = 14)),
  
  smfa_plot(10, n = 2) + theme(legend.position = "none", plot.title = element_text(face = "bold")) +
    ggtitle("model 10") +
    theme(title = element_text(size = 14)),
  
  smfa_plot(12, n = 2) + theme(legend.position = "none", plot.title = element_text(face = "bold")) +
    ggtitle("model 12") +
    theme(title = element_text(size = 14)),
  
  smfa_plot(14, n = 2) + theme(legend.position = "none", plot.title = element_text(face = "bold")) +
    ggtitle("model 14") +
    theme(title = element_text(size = 14)),
  nrow = 3, 
  labels = c("A", "B", "C", "D", "E", "F")
  ),
  get_legend(bf_plot), ncol = 2, rel_widths = c(1, 0.2)
)
dev.off()

########################
### explanatory plot ###
########################

# hours <- seq(0, 24, 0.1)
# 
# df <- data.frame(hour = hours, 
#                  m_temp = c(temp_fun[[1]](hours/24)),
#                  t = rep(19, length(hours))) %>% rowwise() %>% mutate(t_ = (m_temp - mean_temp)/sd_temp,
#                                                                               delta = gen_delta(fit = fit, temp = t_)['50%'][[1]],
#                                                                               PDR = 47/EIP_fun[[2]](m_temp))
# 
# df_delta <- bind_rows(lapply(seq(1, nrow(df)),function(i, df){
#   return(data.frame("max_delta" = max(df[seq(1,i,1),"delta"]),
#              "min_delta" = min(df[seq(1,i,1),"delta"]),
#              "mean_delta" = mean(df[seq(1,i,1),"delta"]),
#              "max_temp" = max(df[seq(1,i,1),"m_temp"]),
#              "min_temp" = min(df[seq(1,i,1),"m_temp"]),
#              "mean_temp_in" = mean(df[seq(1,i,1),"m_temp"])))
# }, df = as.data.frame(df))) %>% rowwise() %>% mutate(max_temp_delta = gen_delta(fit = fit, temp = (max_temp - mean_temp)/sd_temp)['50%'][[1]],
#                                   min_temp_delta = gen_delta(fit = fit, temp = (min_temp - mean_temp)/sd_temp)['50%'][[1]],
#                                   mean_temp_delta = gen_delta(fit = fit, temp = (mean_temp_in - mean_temp)/sd_temp)['50%'][[1]])
# 
# df <- cbind(df, df_delta)
# 
# h_df <- data.frame("y_" = seq(11, 14, 1), "yend_" = seq(11, 14, 1), "xend_" = c(24, 16, 8, 0))
# 
# h_df_d <- data.frame("y_" = rep(0, 4),
#                      "yend_" = c(df[which(df[,"hour"] == 0),"mean_delta"],
#                                  df[which(df[,"hour"] == 8),"mean_delta"],
#                                  df[which(df[,"hour"] == 16),"mean_delta"],
#                                  df[which(df[,"hour"] == 24),"mean_delta"]), 
#                      "xend_" = c(0, 8, 16, 24), 
#                      "x" = c(0, 8, 16, 24),
#                      "label" = c("h = 0", "h = 8", "h = 16", "h = 24"))
# 
# out_h <- bind_rows(lapply(seq(0, 24, 8), function(h_){
#   run_diff_delta(h = h_, 
#                         u_l = u_l, 
#                         delta_vt = TRUE, EIP_vt = TRUE, delta_fun = "mean_delta", 
#                         unique_t_DTR = unique_t_DTR,  fit = fit, mean_temp = mean_temp, sd_temp = sd_temp, v_EIP = FALSE) %>% mutate(h = h_)
# }))
# 
# 
# 
# 
# png(file = "report/h_plot.png", height = 1000, width = 800)
# plot_grid(
#   plot_grid(ggplot() + 
#     geom_line(data = df, aes(x = hour, y = m_temp, linetype = factor(t)), size = 1.5, alpha = 0.75) +
#   theme_bw() + theme(text = element_text(size = 15)) + 
#   geom_segment(data = h_df, aes(x = 0, y =  y_, xend = xend_, yend = yend_, col = factor(xend_)),
#                  arrow = arrow(length = unit(0.3, "cm")), size = 1.5) +
#     scale_colour_viridis_d(name = "h") +
#     scale_linetype_manual(name = "diurnal\ntemperature (°C)", labels = c("19±5"), values = c(1)) +
#     ylab("Temperature (°C)") + xlab("Hours post infection") +
#     scale_x_continuous(breaks = seq(0,24,4)), NULL, rel_widths = c(1, 0.05)),
#   
#   plot_grid(ggplot(data = df %>% gather(key = "delta_in", value = "delta_plot", max_delta, min_delta, mean_delta, max_temp_delta, min_temp_delta, mean_temp_delta)) + 
#     geom_line(data = df, aes(x = hour, y = delta, linetype = delta_in), size = 1.7, alpha = 0.7, linetype = 2) +
#               geom_line(aes(x = hour, y = delta_plot, colour = delta_in), size = 1.7, alpha = 0.7) +
#     theme_bw() + theme(text = element_text(size = 15)) + 
#   scale_colour_brewer(name = "Method to\nestimate the HMTP", 
#                       palette = "BrBG",
#                 #         values = scales::seq_gradient_pal("#E69F00", "#56B4E9", "Lab")(seq(0,1,length.out=7)),
#                          labels = c(expression(paste("maximum ", delta)),
#                                     expression(paste("maximum temperature ", delta)),
#                                     expression(paste("mean ", delta)),
#                                     expression(paste("mean temperature ", delta)),
#                                     expression(paste("minimum ", delta)),
#                                     expression(paste("minimum temperature ", delta)))
#                       ) +
#     ylab("Estimated HMTP") + xlab("The hours post infection (h) over which the HMTP is determined") +
#   geom_segment(data = h_df_d, aes(x = x, y = yend_ - 0.05, xend = xend_, yend = yend_),
#                arrow = arrow(length = unit(0.3, "cm")), size = 1.08) +
#   geom_text(data = h_df_d, aes(x = x, y = yend_ + 0.05, label = label), size = 5) +
#     scale_x_continuous(breaks = seq(0,24,4)),NULL, rel_widths = c(1, 0)),
#   
#   plot_grid(ggplot(data = na.omit(subset(out_h, temp %in% c(19) & DTR == 10 & bt == 12 & DPI < 50)),
#          aes(x = DPI, y = s_prev, col = factor(h))) + geom_line(size = 1.5) +
#     scale_colour_viridis_d(name = "h") +
#     xlab("Days post infection") + ylab("Sporozoite prevalence") +
#     scale_y_continuous(labels = scales::percent) + 
#     theme_bw() + theme(text = element_text(size = 15)) +
#       scale_x_continuous(limits = c(0, 50), breaks = seq(0,50,5)), 
#     NULL, rel_widths = c(1, 0.1725)),
#   
#   nrow = 3, labels = c("A", "B", "C"), rel_widths = c(0.25, 1, 0.25), scale = 0.97
# )
# dev.off()

# #################################################
# ### running the model for higher temperatures ###
# #################################################
# 
# # EIP models with different cutoff
# c_temps <- seq(35, 38)
# m_EIP_all <- bind_rows(sapply(c(Inf, c_temps), 
#                  function(c_temp){
#                    t <- seq(15, 45, 0.01)
#                    EIP <- data.frame("temp" = t,
#                                      "EIP" = EIP_fun[[2]](t), # to get the mean EIP
#                                      "m_EIP" = rep(c_temp, length(t)))
#                    EIP[which(t > c_temp), "EIP"] <- Inf
#                    return(EIP)
#                    }, 
#                  simplify = FALSE))
# png(file = "report/cutoff_sensitivity_PDF.png", height = 500, width = 850)
# ggplot(data = m_EIP_all, aes(x = temp, y = 47/EIP)) +
#   geom_line(size = 1.25) + theme_bw() + facet_wrap(~m_EIP) +
#   theme(text = element_text(size = 15)) +
#   xlab("Temperature (°C)") + ylab("Parasite development rate (shape / EIP)")
# dev.off()
# 
# # EIP functions
# EIP_fun_l <- vector(mode = "list", length = length(c_temps) + 1)
# EIP_fun_l[[1]] <- EIP_fun[[2]]
# p <- data.frame("temp" = seq(15, 45, 0.01),
#                 "EIP" = EIP_fun[[2]](seq(15, 45, 0.01)))
# 
# for(i in 1:length(c_temps)){
#   place <- p[which(p$temp <= c_temps[i]),]
#   EIP_fun_l[[i+1]] <- approxfun(place$temp, place$EIP, yleft = max(mean_EIP$`50%`), yright = Inf, method = "linear")
# }
# 
# # temperature functions
# g_temps <- c(32.5)
# temp_fun_l <- lapply(g_temps, function(i){
#   approxfun(seq(0, 24 * max_time, 0.1)/24, 
#                              gen_temp(i, 10, max_time, bt = 12))
# })
# 
# # delta value of 1 - so just the infected mosquitoes
# 
# index_all <- expand.grid(seq(1,length(g_temps)), seq(1, (length(c_temps)+1)))
# index_temp_ <- index_all[, 1]
# index_EIP_ <- index_all[, 2]
#  
# out <- bind_rows(mapply(function(index_temp, 
#                 index_EIP, 
#                 state = c(U = (1 - 1)*1, E = c(1*1, rep(0, 47-1)),  I = 0),
#                 t = seq(0, 70, 0.1)){
#    params <- c(mu = 0.1, shape = 47, index_EIP = index_EIP, index_temp = index_temp)
#     return(out_mc <- as.data.frame(ode(y = state, times = t, func = model_4_ht,
#                                parms = params)) %>% rowwise() %>%
#      mutate(M = U + sum(c_across("E1":paste0("E",params[["shape"]]))) + I,
#             s_prev = I / M,
#             index_temp = index_temp,
#             index_EIP = index_EIP))
# }, index_temp_, index_EIP_, SIMPLIFY = FALSE))
# 
# c_temps_i <- c(Inf, c_temps)
# out$m_EIP <- c_temps_i[out$index_EIP]
# out$temp <- g_temps[out$index_temp]
# 
# png(file = "report/DTR_simulations_sensitivity.png", height = 750, width = 500)
# plot_grid(
#   ggplot(data = m_EIP_all, aes(x = temp, y = 47/EIP)) +
#   geom_line(size = 1.25) + theme_bw() + facet_wrap(~m_EIP, nrow = 3) +
#   theme(text = element_text(size = 15)) +
#   xlab("Temperature (°C)") + ylab("Parasite development rate (shape / EIP)"),
#   
#   plot_grid(NULL,
#   ggplot(data = out %>% mutate(m_EIP = ifelse(is.infinite(m_EIP), NA, m_EIP)), 
#        aes(x = time, y = s_prev, col = m_EIP, group = factor(m_EIP))) +
#   theme_bw() + theme(text = element_text(size = 15)) +
#   geom_line(size = 0.75) +
#   coord_cartesian(xlim = c(0, 30)) +
#   ylab("Sporozoite prevalence") +
#   xlab("Days post infection") + scale_y_continuous(labels = scales::percent) +
#   scale_colour_gradient(low = "#56B4E9", high = "#E69F00", na.value = "black",
#                         name = ""), NULL, rel_heights = c(0.1, 0.75, 0.1), nrow = 3),
#  nrow = 2, labels = c("A", "B")
# )
# dev.off()
