# Description: R script to compare the Anopheles stephensi and Anopheles gambiae EIP PDFs
# Version: 1.0
# Date: 15/07/2024

rm(list = ls())

#################
### functions ###
#################

source(file = "utils/functions_temp_only.R")
source(file = "utils/vector_comp_functions.R")
source(file = "read_libraries_data.R")

############
### fits ###
############

# number of iterations used to fit the gambiae data
iterations <- 5500
warmup <- 3000
chains <- 4

###########
### EIP ###
###########

# temperature to simulate the EIP PDF
temps_s <- seq(21, 34, 0.1)
scaled_temps_s <- (temps_s - mean_temp_s) / sd_temp_s
temps_g <- seq(17, 30, 0.1)
scaled_temps_g <- (temps_g - mean_temp_g) / sd_temp_g

# calculating the EIP
EIP_index_g <- get_EIP(rstan::extract(fit_gambiae), scaled_temps_g, 10000)
EIP_index_s <- get_EIP_stephensi(rstan::extract(fit_stephensi), scaled_temps_s, 12000)

EIP_q_df <- rbind(cbind(gen_quantiles(EIP_index_g$EIP_10, temps_g), data.frame("EIP" = rep("EIP[10]", length(temps_g)), "vector" = rep("A. gambiae", length(temps_g)))),
                 cbind(gen_quantiles(EIP_index_g$EIP_50, temps_g), data.frame("EIP" = rep("EIP[50]", length(temps_g)), "vector" = rep("A. gambiae", length(temps_g)))),
                 cbind(gen_quantiles(EIP_index_g$EIP_90, temps_g), data.frame("EIP" = rep("EIP[90]", length(temps_g)), "vector" = rep("A. gambiae", length(temps_g)))),
                 
                 cbind(gen_quantiles(EIP_index_s$EIP_10, temps_s), data.frame("EIP" = rep("EIP[10]", length(temps_s)), "vector" = rep("A. stephensi", length(temps_s)))),
                 cbind(gen_quantiles(EIP_index_s$EIP_50, temps_s), data.frame("EIP" = rep("EIP[50]", length(temps_s)), "vector" = rep("A. stephensi", length(temps_s)))),
                 cbind(gen_quantiles(EIP_index_s$EIP_90, temps_s), data.frame("EIP" = rep("EIP[90]", length(temps_s)), "vector" = rep("A. stephensi", length(temps_s)))))

# function to calculate the EIP whilst accounting for uncertainty in the model fits
calc_EIP_v <- function(EIP_index, temps, n_iter){
  
  set.seed(12345)
  x <- runif(n_iter)
  n_temp <- length(temps)
  
  shape_total_S <- EIP_index$shape_total_S
  rate_total_S <- EIP_index$rate_total_S
  mu <- EIP_index$mu
  k <- EIP_index$k
  
  out <- bind_rows(
    lapply(seq(1, n_temp),
         function(j, n_iter, temp, shape_total_S, rate_total_S, mu, k, x){
           print(temp[j])
           EIP <- lapply(seq(1, n_iter), 
                         function(i, j, temp, shape_total_S, rate_total_S, mu, k, x){
                           Inv_EIP_CDF(shape_total_S[i, j], rate_total_S[i, j], mu[i, j], k[i, j], x)
                           }, 
                         j = j,
                         temp = temp,
                         shape_total_S = shape_total_S, 
                         rate_total_S = rate_total_S, 
                         mu = mu, 
                         k = k, 
                         x = x)
    EIP <- unlist(EIP)
    out_i <- data.frame("EIP_50" = quantile(EIP, probs = c(0.5))[[1]],
             "EIP_10" = quantile(EIP, probs = c(0.1))[[1]],
             "EIP_90" = quantile(EIP, probs = c(0.9))[[1]],
             "mean" = mean(EIP),
             "temp" = temps[j])
    rm(list = c("EIP"))
    return(out_i)},
    
    n_iter = n_iter,
    temp = temps,
    shape_total_S = shape_total_S, 
    rate_total_S = rate_total_S, 
    mu = mu, 
    k = k, 
    x = x
  ))
  
 return(out)   
}

EIP_plot_g <- calc_EIP_v(EIP_index = EIP_index_g, temps = temps_g, n_iter = 10000)
EIP_plot_s <- calc_EIP_v(EIP_index = EIP_index_s, temps = temps_s, n_iter = 12000)
saveRDS(EIP_plot_g, file = "EIP_gambiae.rds")      
saveRDS(EIP_plot_s, file = "EIP_stephensi.rds")      

EIP_plot_g <- readRDS(file = "EIP_gambiae.rds")
EIP_plot_s <- readRDS(file = "EIP_stephensi.rds")

EIP_plot <- ggplot() +
  geom_ribbon(data = rbind(EIP_plot_g %>% mutate(vector_species = "gambiae"),
                    EIP_plot_s %>% mutate(vector_species = "stephensi")),
       aes(x = temp, ymin = EIP_10, ymax = EIP_90, fill = vector_species), 
       alpha = 0.325) +
  geom_line(data = rbind(EIP_plot_g %>% mutate(vector_species = "gambiae"),
                         EIP_plot_s %>% mutate(vector_species = "stephensi")),
            aes(x = temp, y = EIP_50, colour = vector_species), linewidth = 1.25) +
  geom_line(data = EIP_plot_g, aes(x = temp, y = mean), col = "#56B4E9", linetype = 2, linewidth = 1.25) +
  geom_line(data = EIP_plot_s, aes(x = temp, y = mean), col = "#E69F00", linetype = 2, linewidth = 1.25) +
  theme_bw() + theme(text = element_text(size = 15)) +
  ylab("EIP") +
  xlab("Temperature (°C)") +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_y_continuous(breaks = seq(0, 80, 10)) +
  scale_x_continuous(breaks = seq(18, 34, 2)) 

################
##### HMTP #####
################

params_g <- rstan::extract(fit_gambiae)
params_s <- rstan::extract(fit_stephensi)

HMTP_s <- matrix(NA, nrow = length(params_s$m_delta), ncol = length(scaled_temps_s))

for(i in 1:length(scaled_temps_s)){
  print(i)
  HMTP_s[,i] <- 1/(1 + exp(-(params_s$m_delta * scaled_temps_s[i] + params_s$c_delta)))
}

HMTP_s <- t(apply(HMTP_s, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame() %>% mutate(temp = temps_s,
                                                                                                     species = "stephensi")

colnames(HMTP_s) <- c("lower", "median", "upper", "temp", "species")

HMTP_g <- matrix(NA, nrow = length(params_g$c_delta), ncol = length(scaled_temps_g))

for(i in 1:length(scaled_temps_g)){
  print(i)
  HMTP_g[,i] <- (1/(1 + exp(-(params_g$a_delta * scaled_temps_g[i]^2 + params_g$b_delta * scaled_temps_g[i] + params_g$c_delta)))) * 
    (1/(1 + exp(-(params_g$a_delta_S * scaled_temps_g[i]^2 + params_g$b_delta_S * scaled_temps_g[i] + params_g$c_delta_S))))
}

HMTP_g <- t(apply(HMTP_g, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame() %>% mutate(temp = temps_g,
                                                                                                     species = "gambiae")

colnames(HMTP_g) <- c("lower", "median", "upper", "temp", "species")

HMTP_df <- rbind(HMTP_s, HMTP_g)

HMTP_plot <- ggplot(data = HMTP_df, aes(x = temp, y = median, ymin = lower, ymax = upper, fill = species)) +
  geom_ribbon(alpha = 0.325) +
  geom_line(aes(colour = species), linewidth = 1.25) +
  theme_bw() + theme(text = element_text(size = 15)) +
  ylab("HMTP") +
  xlab("Temperature (°C)") +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(18, 34, 2)) 

############################
### visualising the fits ###
############################

# gambiae data
g_s_data <- read.csv(file = "data/processed/ES_new_constant_temp_spz_processed.csv") %>% mutate(vector_species = "gambiae",
                                                                                                        gametocytemia = round(gametocytemia, digits = 5),
                                                                                                        ref = "Eunho_Suh")
g_s_data <- g_s_data[,c(2:ncol(g_s_data))]

g_s_data$gametocytemia <- round(g_s_data$gametocytemia, digits = 5)
g_s_data <- subset(g_s_data, gametocytemia!=0.00024)

g_s_data <- g_s_data %>% mutate(index_temp = match(g_s_data$temp, sort(unique(g_s_data$temp))))

g_s_data_in <- generate_prevalence_temp(g_s_data)

# posterior predictive distribution times
PPD_times_g <- readRDS(file = "data/PPD_times_g")
length_ppd_times <- length(PPD_times_g)

g_fits <- run_prop_ppd_df(fit_gambiae, "pooled", "S_prevalence_ppd", length_ppd_times, PPD_times_g, unique_temp = sort(unique(g_s_data_in$temp)))

# stephensi data
stephensi_data <- read.csv("data/processed/mSOS_all_parasite_data.csv") %>% subset(DTR == 0) %>% 
  mutate(gametocytemia = ifelse(ref == "Shapiro_L_Thomas" | ref == "Murdock_Thomas" | ref == "Shapiro_Thomas", 0.08, NA))

stephensi_s_data <- subset(stephensi_data, lifestage == "gland-spz")
stephensi_s_data <- subset(stephensi_s_data, vector_species!="gambiae")

s_fits <- prop_ppd_function(as.data.frame(fit_stephensi), 
                            n_unq_gt = length(sort(unique(stephensi_data$temp))), 
                            length_ppd_times = length(seq(0,30,0.1)),
                            PPD_times = seq(0, 30, 0.1),
                            iterations = 4500, 
                            warmup = 1500, 
                            chains = 4,
                            Stan_data_name = "sporozoite_prevalence_ppd")

s_fits <- s_fits %>% mutate(temp = sort(unique(stephensi_data$temp))[index_gt])

fits_df <- rbind(g_fits %>% mutate("vector_species" = "gambiae"),
                 s_fits %>% mutate("model" = "pooled",
                   "vector_species" = "stephensi"))

stephensi_s_data <- stephensi_s_data[,c("day_post_inf", "presence", "temp", "gametocytemia", "vector_species", "ref")]
colnames(stephensi_s_data)[1] <- "DPI"

s_data_in <- stephensi_s_data %>% group_by(DPI, temp, vector_species) %>% summarise(positive = sum(presence),
                                                                                         sample = n()) %>% 
  ungroup() %>% 
  mutate(prevalence = positive/sample,
         lower = prevalence - (1.96 * sqrt(prevalence * (1 - prevalence) / sample)),
         upper = prevalence + (1.96 * sqrt(prevalence * (1 - prevalence) / sample)))

s_data <- rbind(g_s_data_in[,-which(colnames(g_s_data_in)=="index_temp")] %>% mutate(vector_species = "gambiae"), s_data_in)

m_days <- as.data.frame(s_data %>% group_by(temp) %>% summarise(m_dpi = max(DPI)))

fits_plot_df <- bind_rows(lapply(seq(1, nrow(m_days)), function(i, fits_df, m_days){
   subset(fits_df, DPI < (m_days[i, "m_dpi"] + 5) & temp == m_days[i, "temp"])
},
fits_df = fits_df,
m_days = m_days))

fit_plot_lim <- ggplot(data = fits_plot_df %>% subset(temp %in% c(21, 27, 30)) %>% 
                     mutate(temp_label = paste0(temp, "°C")), aes(x = DPI, y = median, ymin = lower, ymax = upper, group = vector_species)) +
  geom_pointrange(data = s_data %>% subset(temp %in% c(21, 27, 30)) %>% 
                    mutate(temp_label = paste0(temp, "°C")), 
                  aes(x = DPI, y = prevalence, ymin = lower, ymax = upper, 
                      fill = vector_species), alpha = 0.75, shape = 21, col = "grey30", size = 0.8) +
  geom_ribbon(alpha = 0.325, aes(fill = vector_species)) +
  geom_line(aes(col = vector_species), linewidth = 1.25) +
  theme_bw() + theme(text = element_text(size = 15)) +
  facet_wrap(~temp_label) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  ylab("Sporozoite prevalence") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.2)) +
  xlab("Days post infection")

png(file = "results/figures/vector_plot.png", height = 450, width = 700)
plot_grid(
 fit_plot_lim,
 plot_grid(
   EIP_plot, HMTP_plot, nrow = 1, labels = c("B", "C")
 ),
 nrow = 2, labels = c("A",""),
 rel_heights = c(0.35, 0.5)
)

dev.off()

