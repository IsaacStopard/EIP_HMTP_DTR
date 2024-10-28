# gambiae functions

generate_prevalence_temp <- function(data){
  totals <- unique(data[,c("DPI","index_temp", "temp")])
  
  for(i in 1:nrow(totals)){
    totals[i,"sample"] <- length(which(data[,"DPI"] == totals[i, "DPI"]
                                       & data[,"index_temp"] == totals[i, "index_temp"]))
    
    totals[i,"positive"] <- length(which(data[,"DPI"] == totals[i, "DPI"] 
                                         & data[,"presence"] > 0 & data[,"index_temp"] == totals[i, "index_temp"]))
  }
  
  rm(i)
  
  totals <- mutate(totals, prevalence = positive / sample) # prevalence
  totals <- mutate(totals, lower = prevalence - (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 5% CI
  totals <- mutate(totals, upper = prevalence + (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 95% CI
  totals[which(totals$lower < 0), "lower"] <- 0 # preventing the lower confidence interval being below 0
  
  return(totals)
}

prop_ppd_function <- function(fit_df, n_unq_gt, length_ppd_times, PPD_times, iterations, warmup, chains, Stan_data_name){
  prop_ppd <- array(NaN, c(length_ppd_times, ((iterations - warmup) * chains), n_unq_gt))
  for(i in 1:length_ppd_times){
    for(j in 1:n_unq_gt){
      prop_ppd[i,,j] <- fit_df[,paste0(Stan_data_name,"[",i,",",j,"]")]
    }
  }
  prop_ppd_df <- data.frame()
  labs_gt_ind <- c()
  for(i in 1:n_unq_gt){
    place <- as.data.frame(prop_ppd[,,i])
    prop_ppd_df <- rbind(prop_ppd_df, place)
    labs_gt_ind <-  append(labs_gt_ind, rep(i, length_ppd_times))
  }
  
  day_post_inf <- rep(PPD_times, n_unq_gt)
  prop_ppd_df[,"DPI"] <- day_post_inf
  prop_ppd_df[,"index_gt"] <- labs_gt_ind
  prop_ppd_df$median <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, median)
  prop_ppd_df$lower <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.025))
  prop_ppd_df$upper <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.975))
  prop_ppd_df$mean <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, mean)
  
  #x <- ((iterations - warmup) * chains) 
  #prop_ppd_df <- prop_ppd_df %>% gather("iteration", "value", 1:x)
  #prop_ppd_df$percentile <- rep("iteration", nrow(prop_ppd_df))
  
  #prop_quantile_ppd_df <- subset(prop_ppd_df, iteration == "V1")
  #prop_quantile_ppd_df <- prop_quantile_ppd_df[,c("DPI", "index_gt", "median", "lower", "upper", "mean")]
  
  return(prop_ppd_df[,c("DPI", "index_gt", "median", "lower", "upper", "mean")])
}

run_prop_ppd_df <- function(fit, model, Stan_data_name, length_ppd_times, PPD_times, unique_temp){
  df <- prop_ppd_function(as.data.frame(fit),length(unique_temp), length_ppd_times, PPD_times, iterations, warmup, chains, Stan_data_name)
  df$model <- rep(model, nrow(df))
  df$temp <- unique_temp[df$index_gt]
  return(df)
}



