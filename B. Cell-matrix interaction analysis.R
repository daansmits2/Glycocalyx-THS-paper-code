## install necessary packages------
library(ggplot2)
library(gridExtra)
library(tidyr)
library(matrixStats)
library(plyr)
library(dplyr)
library(ggpubr)
library(Rmisc)
library(svDialogs)
library(cowplot)
library(ggh4x)
library(ggpmisc)
library(nortest)
library(moments)
library(stats)
library(tseries)
library(multcomp)


## chosen parameters ----

pxsize = 0.042
Xthreshold = 3

xlim = 15

Xthr= 10

Clustersize_threshold = 0.3


bleb_length_segregation_threshold = 1/3
retraction_fiber_length_segregation_threshold = 4/5
protrusion_length_segregation_threshold = 0.25

## load necessary functions ------
zscore <- function(data){
  if(class(data)=="numeric" | class(data)=="integer"){
    zscore <- (data-mean(data, na.rm=TRUE))/sd(data, na.rm=TRUE)
  }
  if(class(data)=="data.frame"){
    zscore <- sapply(data, function(data) (data-mean(data, na.rm=TRUE))/sd(data, na.rm=TRUE))
  }
  return(zscore)
}

check_bad_values <- function(data, x_var, y_var) {
  # Convert column names to actual column references
  x_col <- data[[x_var]]
  y_col <- data[[y_var]]
  
  # Calculate counts of bad values
  x_na_count <- sum(is.na(x_col))
  x_nan_count <- sum(is.nan(x_col))
  x_inf_count <- sum(is.infinite(x_col))
  
  y_na_count <- sum(is.na(y_col))
  y_nan_count <- sum(is.nan(y_col))
  y_inf_count <- sum(is.infinite(y_col))
  
  # Create a list to return the results
  results <- list(
    x_var = list(
      NA_count = x_na_count,
      NaN_count = x_nan_count,
      Inf_count = x_inf_count
    ),
    y_var = list(
      NA_count = y_na_count,
      NaN_count = y_nan_count,
      Inf_count = y_inf_count
    )
  )
  
  return(results)
}

filter_bad_values <- function(data, x_var, y_var){
  
  library(dplyr)
  library(rlang)
  x_var_sym <- sym(x_var)
  y_var_sym <- sym(y_var)
  
  data %>% 
    mutate(
      !!x_var_sym := ifelse(!!x_var_sym == 0, 0.00001, !!x_var_sym),
      !!y_var_sym := ifelse(!!y_var_sym == 0, 0.00001, !!y_var_sym)
    ) %>%
    filter(!is.na(!!x_var_sym) & !is.infinite(!!x_var_sym) & !is.nan(!!x_var_sym) &
             !is.na(!!y_var_sym) & !is.infinite(!!y_var_sym) & !is.nan(!!y_var_sym))
}

transform_smooth <- function(value) {
  log_smooth_min <- log(smooth_min)
  log_smooth_max <- log(smooth_max)
  log_smooth <- log_smooth_min + (log_smooth_max - log_smooth_min) * (value - data_min) / (data_max - data_min)
  log_smooth <- exp(log_smooth)
  return(log_smooth)
}

glht.table <- function(x) {
  library(multcomp)
  pq <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  pname <- switch(x$alternativ, 
                  less = paste("Pr(<", ifelse(x$df == 0, "z", "t"), ")", sep = ""), 
                  greater = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""), 
                  two.sided = paste("Pr(>|", ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
  colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df == 0, "z value", "t value"), pname)
  
  # Add a "Significance" column based on p-values
  significance <- ifelse(mtests[, "Pr(>|t|)"] < 0.001, "***",
                         ifelse(mtests[, "Pr(>|t|)"] < 0.01, "**",
                                ifelse(mtests[, "Pr(>|t|)"] < 0.05, "*", "")))
  mtests <- cbind(mtests, Significance = significance)
  
  return(mtests)
}

delete_entry <- function(file, B1int_GCX_ratio, number){
  B1int_GCX_ratio <- B1int_GCX_ratio[B1int_GCX_ratio$number!=number,]
  file <- file[file$number!=number,]
}

bleb_length_segregation <- function(B1int_GCX_ratio, bleb_length_segregation_threshold){
  
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%  group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>% 
    mutate(X_max = max (X),
           Type = ifelse(X<=(X_max*bleb_length_segregation_threshold), "Bleb", "Cell body"))
  
  return(B1int_GCX_ratio)
}


protrusion_length_segregation <- function(B1int_GCX_ratio, protrusion_length_segregation_threshold){
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>% 
    arrange(number, X) %>% 
    group_by(number) %>%
    mutate(X_relative = X_micron / max(X_micron))
  
  sub_protrusion <- subset(B1int_GCX_ratio, Orientation == "Leading edge protrusion")
  sub_protrusion <- sub_protrusion %>%
    arrange(number, X) %>%
    group_by(number) %>%
    mutate(Length_label = ifelse(X_relative < (protrusion_length_segregation_threshold), "Protrusion tip", "Intermediate zone"),
           Length_label = ifelse(X_relative > (1-protrusion_length_segregation_threshold), "Protrusion base", Length_label)
    )
  
  
  df <- subset(B1int_GCX_ratio, Orientation!="Leading edge protrusion")
  
  df <- rbind(df, sub_protrusion)
  return(df)
}

define_focalized_b1int_leadingedge <- function(B1int_GCX_ratio, Xdiff, Zdiff){
  library(pracma)
  leadingedge_data <- B1int_GCX_ratio %>%  filter(Orientation=="Leading edge protrusion")
  

  
  maxima <- leadingedge_data  %>%  group_by(Orientation, number, RoiNR, CellNR, 
                                            Signal_type, Position_relative_to_fiber) %>% 
    summarise(
      mean_B1int = mean(B1int, na.rm=TRUE),
      sd_B1int = sd(B1int, na.rm=TRUE),
      maxima_values = list(findpeaks(B1int, nups=Xdiff, ndowns=Xdiff, minpeakheight=mean_B1int+Zdiff*sd_B1int))
    ) %>%
    dplyr::select(-c(mean_B1int, sd_B1int)) %>% 
    ungroup()
  

  # Explode the list of maxima into separate rows
  maxima <- maxima %>%
    tidyr::unnest(maxima_values) 
  
  
  maxima_values <- maxima$maxima_values
  maxima <- maxima[,-ncol(maxima)]
  maxima_values <- maxima_values[,1]
  
  maxima <- cbind(maxima, maxima_values)
  
  maxima <- maxima %>%
    rename(Peaks = maxima_values)
  
  maxima <- maxima %>%  mutate(Type="Focalized")
  
  leadingedge_data <- leadingedge_data %>%
    mutate(Type = ifelse(B1int %in% maxima$Peaks, "Focalized", "Unspecified"))
  
  
  
  
  # Find minima using findpeaks with appropriate options
  minima <- leadingedge_data %>%
    group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>%
    summarise(
      minima_values = list(findpeaks(-B1int, nups=Xdiff, ndowns=Xdiff))
    ) %>%
    ungroup()
  
  # Explode the list of minima into separate rows
  minima <- minima %>%
    tidyr::unnest(minima_values)

  
  
  minima_values <- minima$minima_values
  minima <- minima[,-ncol(minima)]
  minima_values <- minima_values[,1]
  
  minima <- cbind(minima, minima_values)
  
  minima <- minima %>%
    rename(Peaks = minima_values)
  
  minima <- minima %>%  mutate(Type="Non-focalized")
  
  
  
  # Label the rows in leadingedge_data based on the minima
  leadingedge_data <- leadingedge_data %>%
    mutate(Type = ifelse(-B1int %in% minima$Peaks, "Non-focalized", Type))
  
  nonleadingedge_data <- subset(B1int_GCX_ratio, Orientation!="Leading edge protrusion")
  barry <- rbind(leadingedge_data, nonleadingedge_data)
  return(barry)
}



define_focalized_b1int_leadingedge_v2 <- function(B1int_GCX_ratio, smooth=TRUE, smooth_span=0.15){
  library(pracma)
  library(dplyr)
  library(tidyr)
  
  leadingedge_data <- B1int_GCX_ratio
  
  if(smooth==TRUE){
    print("Mean smoothening applied for peak detection")
    leadingedge_data <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>% 
      mutate(B1int_smooth = ifelse(Signal_type=="foreground",
             loess(B1int ~ seq_along(B1int), span=smooth_span, family="gaussian")$fitted, 
             B1int)
      )
  }else{
    print("Peak detection without smoothening")
    leadingedge_data <- leadingedge_data %>% 
      mutate(B1int_smooth  = B1int)
  }
  
  maxima <- leadingedge_data %>%
    group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>%
    summarise(
      maxima_values = list(pracma::findpeaks(B1int_smooth))  # Assuming default values for minima
    ) %>%
    ungroup() %>%
    tidyr::unnest(maxima_values)
  
  new_Peak_ID <- rep(1:nrow(maxima))
  maxima <- cbind(maxima$number, maxima$Signal_type, new_Peak_ID, as.data.frame(maxima$maxima_values))
  colnames(maxima) <- c("number", "Signal_type", "new_Peak_ID", "B1int_smooth", 
                        "Peaklocation_X", "Peak_left_X", "Peak_right_X")  
  
  
  selection <- maxima %>%  dplyr::select(Signal_type, number, Peaklocation_X) %>%  rename(X=Peaklocation_X)
  
  merged_data <- leadingedge_data %>%
    group_by(Signal_type, number) %>%
    left_join(selection, keep=TRUE) %>%
    dplyr::rename(X = X.x,
                  Peaklocation_X = X.y,
                  Signal_type = Signal_type.x,
                  number = number.x) %>% 
    dplyr::select(-c("Signal_type.y", "number.y"))
  
  selection <- maxima %>%  dplyr::select(Signal_type, number, Peak_left_X) %>%  rename(X=Peak_left_X)
  
  merged_data <- merged_data %>%
    group_by(Signal_type, number) %>%
    left_join(selection, keep=TRUE) %>%
    dplyr::rename(X = X.x,
                  Peak_left_X = X.y,
                  Signal_type = Signal_type.x,
                  number = number.x) %>% 
    dplyr::select(-c("Signal_type.y", "number.y"))
  
  
  selection <- maxima %>%  dplyr::select(Signal_type, number, Peak_right_X) %>%  rename(X=Peak_right_X)
  
  merged_data <- merged_data %>%
    group_by(Signal_type, number) %>%
    left_join(selection, keep=TRUE) %>%
    dplyr::rename(X = X.x,
                  Peak_right_X = X.y,
                  Signal_type = Signal_type.x,
                  number = number.x) %>% 
    dplyr::select(-c("Signal_type.y", "number.y"))
  
  
  merged_data <- merged_data %>% 
    group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>%
    arrange(X) %>% 
    mutate(new_Peak_ID_left = cumsum(!is.na(Peak_left_X)),
           new_Peak_ID_right= cumsum(!is.na(Peak_right_X))+1,
           new_Peak_ID = ifelse(new_Peak_ID_left == new_Peak_ID_right, new_Peak_ID_left, NA),
           new_Peak_ID = as.character(new_Peak_ID)
    ) %>% 
    dplyr::select(-c(new_Peak_ID_left, new_Peak_ID_right))
  
  return(merged_data)
}

define_focalized_b1int_leadingedge_v2 <- function(B1int_GCX_ratio, smooth=TRUE, smooth_span=0.15){
  library(pracma)
  library(dplyr)
  library(tidyr)
  
  grouped_coordinates = 48
  
  leadingedge_data <- B1int_GCX_ratio
  
  if(smooth==TRUE){
    print("Mean smoothening applied for peak detection")
    leadingedge_data <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>% 
      mutate(X_group = (X-1) %/% grouped_coordinates + 1,
             X_group = ifelse(X_group == max(X_group, na.rm=TRUE), X_group-1, X_group)
             ) %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber, X_group) %>% 
      mutate(B1int_smooth = ifelse(Signal_type=="foreground",
                                   loess(B1int ~ seq_along(B1int), span=smooth_span, family="gaussian")$fitted, 
                                   B1int)
      )
  }else{
    print("Peak detection without smoothening")
    leadingedge_data <- leadingedge_data %>% 
      mutate(B1int_smooth  = B1int)
  }
  
  maxima <- leadingedge_data %>%
    group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>%
    summarise(
      maxima_values = list(pracma::findpeaks(B1int_smooth))  # Assuming default values for minima
    ) %>%
    ungroup() %>%
    tidyr::unnest(maxima_values)
  
  new_Peak_ID <- rep(1:nrow(maxima))
  maxima <- cbind(maxima$number, maxima$Signal_type, new_Peak_ID, as.data.frame(maxima$maxima_values))
  colnames(maxima) <- c("number", "Signal_type", "new_Peak_ID", "B1int_smooth", 
                        "Peaklocation_X", "Peak_left_X", "Peak_right_X")  
  
  
  selection <- maxima %>%  dplyr::select(Signal_type, number, Peaklocation_X) %>%  rename(X=Peaklocation_X)
  
  merged_data <- leadingedge_data %>%
    group_by(Signal_type, number) %>%
    left_join(selection, keep=TRUE) %>%
    dplyr::rename(X = X.x,
                  Peaklocation_X = X.y,
                  Signal_type = Signal_type.x,
                  number = number.x) %>% 
    dplyr::select(-c("Signal_type.y", "number.y"))
  
  selection <- maxima %>%  dplyr::select(Signal_type, number, Peak_left_X) %>%  rename(X=Peak_left_X)
  
  merged_data <- merged_data %>%
    group_by(Signal_type, number) %>%
    left_join(selection, keep=TRUE) %>%
    dplyr::rename(X = X.x,
                  Peak_left_X = X.y,
                  Signal_type = Signal_type.x,
                  number = number.x) %>% 
    dplyr::select(-c("Signal_type.y", "number.y"))
  
  
  selection <- maxima %>%  dplyr::select(Signal_type, number, Peak_right_X) %>%  rename(X=Peak_right_X)
  
  merged_data <- merged_data %>%
    group_by(Signal_type, number) %>%
    left_join(selection, keep=TRUE) %>%
    dplyr::rename(X = X.x,
                  Peak_right_X = X.y,
                  Signal_type = Signal_type.x,
                  number = number.x) %>% 
    dplyr::select(-c("Signal_type.y", "number.y"))
  
  
  merged_data <- merged_data %>% 
    group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>%
    arrange(X) %>% 
    mutate(new_Peak_ID_left = cumsum(!is.na(Peak_left_X)),
           new_Peak_ID_right= cumsum(!is.na(Peak_right_X))+1,
           new_Peak_ID = ifelse(new_Peak_ID_left == new_Peak_ID_right, new_Peak_ID_left, NA),
           new_Peak_ID = as.character(new_Peak_ID)
    ) %>% 
    dplyr::select(-c(new_Peak_ID_left, new_Peak_ID_right))
  
  return(merged_data)
}



define_focalized_b1int_all_v2 <- function(B1int_GCX_ratio, Zdiff=1, smooth=TRUE, smooth_span=0.15){
  library(pracma)
  library(dplyr)
  library(tidyr)
  
  leadingedge_data <- B1int_GCX_ratio %>%  filter(Signal_type=="foreground")
  
  if(smooth==TRUE){
    print("Mean smoothening applied for peak detection")
    leadingedge_data <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>% 
      mutate(B1int_smooth = loess(B1int ~ seq_along(B1int), span=smooth_span, family="gaussian")$fitted)
  }else{
    print("Peak detection without smoothening")
    leadingedge_data <- leadingedge_data %>% 
      mutate(B1int_smooth  = B1int)
  }
  
  # Initialize a dataframe to store peaks and their max Xdiff
  peaks_xdiff <- data.frame(Peak = numeric(), MaxXdiff = numeric(), Orientation = character(), number = character())
  
  
  for (xdiff in 1:70) {
    print(xdiff)
    
    temp_maxima <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>%
      summarise(
        mean_B1int_smooth = mean(B1int_smooth, na.rm = TRUE),
        sd_B1int_smooth = sd(B1int_smooth, na.rm = TRUE),
        B1int_norm.profile = B1int/ mean(B1int, na.rm=TRUE),
        maxima_values = list(pracma::findpeaks(B1int_smooth, nups = xdiff, ndowns = xdiff))
      ) %>%
      ungroup() 
    
    temp_maxima <- unnest(temp_maxima, maxima_values) 
    
    temp_maxima <- temp_maxima %>% 
      dplyr::select(Peak = maxima_values)
    
    print(temp_maxima)
    
    # Update peaks_xdiff with the maximum Xdiff, Orientation, and number for each peak
    for (peak in temp_maxima$Peak) {
      if (!peak %in% peaks_xdiff$Peak) {
        new_row <- data.frame(Peak = peak, MaxXdiff = xdiff, 
                              Orientation = temp_maxima$Orientation[1],  # Assuming all orientations are the same for a given peak
                              number = temp_maxima$number[1])  # Assuming all numbers are the same for a given peak
        peaks_xdiff <- dplyr::bind_rows(peaks_xdiff, new_row)
      } else {
        peaks_xdiff$MaxXdiff[peaks_xdiff$Peak == peak] <- max(peaks_xdiff$MaxXdiff[peaks_xdiff$Peak == peak], xdiff)
        peaks_xdiff$Orientation[peaks_xdiff$Peak == peak] <- temp_maxima$Orientation[1]  # Assuming all orientations are the same for a given peak
        peaks_xdiff$number[peaks_xdiff$Peak == peak] <- temp_maxima$number[1]  # Assuming all numbers are the same for a given peak
      }
    }
  }
  
  
  
  # Merge the peaks_xdiff with leadingedge_data
  leadingedge_data <- leadingedge_data %>%
    left_join(peaks_xdiff, by = c("B1int_smooth" = "Peak")) %>%
    mutate(Type = ifelse(!is.na(MaxXdiff), "Focalized", "Unspecified"))
  
  bg <- subset(B1int_GCX_ratio, Signal_type=="background")
  combined_data <- rbind(leadingedge_data, bg)
  
  
  return(combined_data)
}

define_focalized_b1int_all_v2 <- function(B1int_GCX_ratio, Zdiff=1, smooth=TRUE, smooth_span=0.15){
  library(pracma)
  library(dplyr)
  library(tidyr)
  
  leadingedge_data <- B1int_GCX_ratio %>%  filter(Signal_type=="foreground")
  
  if(smooth==TRUE){
    print("Mean smoothening applied for peak detection")
    leadingedge_data <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>% 
      mutate(B1int_smooth = loess(B1int ~ X, span=smooth_span, family="gaussian")$fitted)
  }else{
    print("Peak detection without smoothening")
    leadingedge_data <- leadingedge_data %>% 
      mutate(B1int_smooth  = B1int)
  }
  
  # Initialize a dataframe to store peaks and their max Xdiff
  peaks_xdiff <- data.frame(Peak = numeric(), MaxXdiff = numeric(), Orientation = character(), number = character())
  
  
  for (xdiff in 1:70) {
    print(xdiff)
    
    temp_maxima <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>%
      summarise(
        mean_B1int_smooth = mean(B1int_smooth, na.rm = TRUE),
        sd_B1int_smooth = sd(B1int_smooth, na.rm = TRUE),
        B1int_norm.profile = B1int/ mean(B1int, na.rm=TRUE),
        maxima_values = list(pracma::findpeaks(B1int_smooth, nups = xdiff, ndowns = xdiff))
      ) %>%
      ungroup() 
    
    temp_maxima <- unnest(temp_maxima, maxima_values) 
    
    temp_maxima <- temp_maxima %>% 
      dplyr::select(Peak = maxima_values)
    
    print(temp_maxima)
    
    # Update peaks_xdiff with the maximum Xdiff, Orientation, and number for each peak
    for (peak in temp_maxima$Peak) {
      if (!peak %in% peaks_xdiff$Peak) {
        new_row <- data.frame(Peak = peak, MaxXdiff = xdiff, 
                              Orientation = temp_maxima$Orientation[1],  # Assuming all orientations are the same for a given peak
                              number = temp_maxima$number[1])  # Assuming all numbers are the same for a given peak
        peaks_xdiff <- dplyr::bind_rows(peaks_xdiff, new_row)
      } else {
        peaks_xdiff$MaxXdiff[peaks_xdiff$Peak == peak] <- max(peaks_xdiff$MaxXdiff[peaks_xdiff$Peak == peak], xdiff)
        peaks_xdiff$Orientation[peaks_xdiff$Peak == peak] <- temp_maxima$Orientation[1]  # Assuming all orientations are the same for a given peak
        peaks_xdiff$number[peaks_xdiff$Peak == peak] <- temp_maxima$number[1]  # Assuming all numbers are the same for a given peak
      }
    }
  }
  
  
  
  # Merge the peaks_xdiff with leadingedge_data
  leadingedge_data <- leadingedge_data %>%
    left_join(peaks_xdiff, by = c("B1int_smooth" = "Peak")) %>%
    mutate(Type = ifelse(!is.na(MaxXdiff), "Focalized", "Unspecified"))
  
  bg <- subset(B1int_GCX_ratio, Signal_type=="background")
  combined_data <- rbind(leadingedge_data, bg)
  
  
  return(combined_data)
}

define_focalized_b1int_cellbody_v2 <- function(B1int_GCX_ratio, Zdiff=1, smooth=TRUE, smooth_span=0.15){
  library(pracma)
  library(dplyr)
  library(tidyr)
  
  leadingedge_data <- B1int_GCX_ratio %>% filter(Orientation == "Cellbody")
  
  leadingedge_data$Signal_type <- as.factor(leadingedge_data$Signal_type)
  leadingedge_data$Signal_type <- factor(leadingedge_data$Signal_type, levels=c("foreground", "background"))
  
  if(smooth==TRUE){
    print("Mean smoothening applied for peak detection")
    leadingedge_data <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>% 
      mutate(B1int_smooth = loess(B1int ~ seq_along(B1int), span=smooth_span, family="gaussian")$fitted)
  }else{
    print("Peak detection without smoothening")
    leadingedge_data <- leadingedge_data %>% 
      mutate(B1int_smooth  = B1int)
  }
  
  # Initialize a dataframe to store peaks and their max Xdiff
  peaks_xdiff <- data.frame(Peak = numeric(), MaxXdiff = numeric())
  
  for (xdiff in 1:70) {
    temp_maxima <- leadingedge_data %>% 
      group_by(Orientation, number, LocationNR, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>%
      summarise(
        mean_B1int_smooth = mean(B1int_smooth, na.rm = TRUE),
        sd_B1int_smooth = sd(B1int_smooth, na.rm = TRUE),
        B1int_norm.profile = B1int/ mean(B1int, na.rm=TRUE),
        maxima_values = list(pracma::findpeaks(B1int_smooth, nups = xdiff, ndowns = xdiff))
      ) %>%
      ungroup() 
    
    temp_maxima <- unnest(temp_maxima, maxima_values) 
    
    temp_maxima <- temp_maxima %>% 
      dplyr::select(Peak = maxima_values)
    
    # Update peaks_xdiff with the maximum Xdiff for each peak
    for (peak in temp_maxima$Peak) {
      if (!peak %in% peaks_xdiff$Peak) {
        peaks_xdiff <- rbind(peaks_xdiff, data.frame(Peak = peak, MaxXdiff = xdiff))
      } else {
        peaks_xdiff$MaxXdiff[peaks_xdiff$Peak == peak] <- max(peaks_xdiff$MaxXdiff[peaks_xdiff$Peak == peak], xdiff)
      }
    }
  }
  
  
  
  # Merge the peaks_xdiff with leadingedge_data
  leadingedge_data <- leadingedge_data %>%
    left_join(peaks_xdiff, by = c("B1int_smooth" = "Peak")) %>%
    mutate(Type = ifelse(!is.na(MaxXdiff), "Focalized", "Unspecified"))
  
  # Find minima using findpeaks with appropriate options
  minima <- leadingedge_data %>%
    group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>%
    summarise(
      minima_values = list(pracma::findpeaks(-B1int, nups = 1, ndowns = 1))  # Assuming default values for minima
    ) %>%
    ungroup() %>%
    tidyr::unnest(minima_values)
  
  # Process the identified minima
  minima <- minima %>%
    mutate(Type = "Non-focalized")
  
  # Combine the processed leading edge data with non-leading edge data
  nonleadingedge_data <- subset(B1int_GCX_ratio, Orientation != "Cellbody")
  combined_data <- rbind(leadingedge_data, minima, nonleadingedge_data)
  
  
  return(combined_data)
}


additional_peak_filter <- function(protrusion, ratiothr=0, sd_thr=0){
  protrusion <- protrusion %>% group_by(number, Signal_type, PeakID_perprofile) %>% 
    mutate(B1int_smooth_edge = mean(B1int_smooth[Peak_position=="Edge"]),
           B1int_smooth_max = max(B1int_smooth[Peak_position=="Inside"]),
           B1int_max = max(B1int[Peak_position=="Inside"], na.rm=TRUE),
           B1int_diff = B1int_max / B1int_smooth_edge,
           B1int_sd = sd(B1int, na.rm=TRUE),
           PeakID_perprofile = ifelse(B1int_diff < ratiothr, NA, PeakID_perprofile),
           PeakID_perprofile = ifelse(B1int_sd < sd_thr, NA, PeakID_perprofile)
           ) %>% 
    dplyr::select(-c(B1int_smooth_edge, B1int_smooth_max))
  
  return(protrusion)
}

calc_local_segregation_leadingedge <- function(B1int_GCX_ratio, Xdiff) {
  protrusion_data <- B1int_GCX_ratio %>%  filter(Orientation=="Leading edge protrusion")
  
  
  # First, identify the X positions of the focalized peaks
  protrusion_data <- protrusion_data %>% group_by(number, Signal_type) %>% arrange(X) %>% 
    mutate(lag_Focalized = Type== "Focalized") 
  protrusion_data <- protrusion_data %>%  
    mutate(PartOfPeak = lag(lag_Focalized, Xdiff) | lead(lag_Focalized, Xdiff)) %>%  
    mutate(PartOfPeak = ifelse(Type=="Focalized", TRUE, PartOfPeak)) %>% 
    dplyr::select(-lag_Focalized)
  
  
  
  protrusion_data <- protrusion_data %>%
    group_by(number,CellNR, LocationNR, RoiNR, Type) %>% 
    mutate(
      Peak_ID = if_else(Type == "Focalized", cumsum(Type == "Focalized"), NA_integer_)
    )
  
  protrusion_data <- protrusion_data %>%
    group_by(Peak_ID, number) %>%
    mutate(Peak_ID = ifelse(is.na(Peak_ID), NA, cur_group_id())) %>%
    ungroup()
  
  protrusion_data <- protrusion_data %>% 
    group_by(Peak_ID, Signal_type) %>% 
    mutate(Peak_position = ifelse(Type=="Focalized" & PartOfPeak==TRUE, "Inside FA", NA)) %>% 
    mutate(Peak_position = ifelse(Type!="Focalized" & PartOfPeak==TRUE, "Outside FA", Peak_position))
  
  
  protrusion_data <-  protrusion_data %>% 
    group_by(Peak_ID,Peak_position, Signal_type) %>%  
    mutate(ratio_perpeak = ifelse(PartOfPeak==TRUE, B1int/Gcx, NA)) %>% 
    mutate(ratio_perpeak = ifelse(is.infinite(ratio_perpeak), NA, ratio_perpeak)) %>% 
    ungroup() %>% group_by(Signal_type) %>% 
    mutate(ratio_perpeak = ratio_perpeak/ mean(ratio_perpeak[Peak_position=="Outside FA"], na.rm=TRUE))
  
  
  
  return(protrusion_data)
}

calc_local_segregation_leadingedge <- function(B1int_GCX_ratio, Xdiff) {
  protrusion_data <- B1int_GCX_ratio 
  
  
  # First, identify the X positions of the focalized peaks
  protrusion_data <- protrusion_data %>% group_by(number, Signal_type) %>% arrange(X) %>% 
    mutate(lag_Focalized = Type== "Focalized") 
  protrusion_data <- protrusion_data %>%  
    mutate(PartOfPeak = lag(lag_Focalized, Xdiff) | lead(lag_Focalized, Xdiff)) %>%  
    mutate(PartOfPeak = ifelse(Type=="Focalized", TRUE, PartOfPeak)) %>% 
    dplyr::select(-lag_Focalized)
  
  
  
  protrusion_data <- protrusion_data %>%
    group_by(number,CellNR, LocationNR, RoiNR, Type) %>% 
    mutate(
      Peak_ID = if_else(Type == "Focalized", cumsum(Type == "Focalized"), NA_integer_)
    )
  
  protrusion_data <- protrusion_data %>%
    group_by(Peak_ID, number) %>%
    mutate(Peak_ID = ifelse(is.na(Peak_ID), NA, cur_group_id())) %>%
    ungroup()
  
  protrusion_data <- protrusion_data %>% 
    group_by(Peak_ID, Signal_type) %>% 
    mutate(Peak_position = ifelse(Type=="Focalized" & PartOfPeak==TRUE, "Inside FA", NA)) %>% 
    mutate(Peak_position = ifelse(Type!="Focalized" & PartOfPeak==TRUE, "Outside FA", Peak_position))
  
  
  protrusion_data <-  protrusion_data %>% 
    group_by(Peak_ID,Peak_position, Signal_type) %>%  
    mutate(ratio_perpeak = ifelse(PartOfPeak==TRUE, B1int/Gcx, NA)) %>% 
    mutate(ratio_perpeak = ifelse(is.infinite(ratio_perpeak), NA, ratio_perpeak)) %>% 
    ungroup() %>% group_by(Signal_type) %>% 
    mutate(ratio_perpeak = ratio_perpeak/ mean(ratio_perpeak[Peak_position=="Outside FA"], na.rm=TRUE))
  
  
  
  return(protrusion_data)
}

label_peak_id_grouped <- function(df) {
  # Initialize a vector to store the updated Peak_IDs
  new_Peak_ID <- df$Peak_ID
  
  for (i in which(df$Type == "Focalized")) {
    # Calculate range
    start_index = max(1, i - df$MaxXdiff[i])
    end_index = min(nrow(df), i + df$MaxXdiff[i])
    
    start_index = start_index + 1
    end_index = end_index -1
    
    # Assign the same Peak_ID to the range
    new_Peak_ID[start_index:end_index] <- df$Peak_ID[i]
  }
  
  return(new_Peak_ID)
}

define_focalized_b1int_retractionfiber <- function(B1int_GCX_ratio, Xdiff){
  library(pracma)
  leadingedge_data <- B1int_GCX_ratio %>%  filter(Orientation=="Retraction fiber")
  
  
  
  maxima <- leadingedge_data  %>%  group_by(Orientation, number, RoiNR, CellNR, Signal_type, Position_relative_to_fiber) %>% 
    summarise(
      maxima_values = list(findpeaks(B1int, nups=Xdiff, ndowns=Xdiff))
    ) %>%
    ungroup()
  
  
  # Explode the list of maxima into separate rows
  maxima <- maxima %>%
    tidyr::unnest(maxima_values) 
  
  
  maxima_values <- maxima$maxima_values
  maxima <- maxima[,-ncol(maxima)]
  maxima_values <- maxima_values[,1]
  
  maxima <- cbind(maxima, maxima_values)
  
  maxima <- maxima %>%
    rename(Peaks = maxima_values)
  
  maxima <- maxima %>%  mutate(Type="Focalized")
  
  leadingedge_data <- leadingedge_data %>%
    mutate(Type = ifelse(B1int %in% maxima$Peaks, "Focalized", "Unspecified"))
  
  
  
  
  # Find minima using findpeaks with appropriate options
  minima <- leadingedge_data %>%
    group_by(Orientation, number, RoiNR, CellNR, Signal_type) %>%
    summarise(
      minima_values = list(findpeaks(-B1int, nups=Xdiff, ndowns=Xdiff))
    ) %>%
    ungroup()
  
  # Explode the list of minima into separate rows
  minima <- minima %>%
    tidyr::unnest(minima_values)
  
  
  
  minima_values <- minima$minima_values
  minima <- minima[,-ncol(minima)]
  minima_values <- minima_values[,1]
  
  minima <- cbind(minima, minima_values)
  
  minima <- minima %>%
    rename(Peaks = minima_values)
  
  minima <- minima %>%  mutate(Type="Non-focalized")
  
  
  
  # Label the rows in leadingedge_data based on the minima
  leadingedge_data <- leadingedge_data %>%
    mutate(Type = ifelse(-B1int %in% minima$Peaks, "Non-focalized", Type))
  
  nonleadingedge_data <- subset(B1int_GCX_ratio, Orientation!="Retraction fiber")
  barry <- rbind(leadingedge_data, nonleadingedge_data)
  return(barry)
}

retraction_fiber_length_segregation <- function(B1int_GCX_ratio, retraction_fiber_length_segregation_threshold) {
  B1int_GCX_ratio <- B1int_GCX_ratio %>% 
    arrange(number, X) %>% 
    group_by(number) %>%
    mutate(X_relative = row_number() / n())
  
  sub_retraction_fibers <- subset(B1int_GCX_ratio, Orientation == "Retraction fiber")
  sub_retraction_fibers <- sub_retraction_fibers %>%
    arrange(number, X) %>%
    group_by(number) %>%
    mutate(Fiber_type = ifelse(X_relative < retraction_fiber_length_segregation_threshold, "Retraction fiber", "Cell body"))
 
  
  df <- subset(B1int_GCX_ratio, Orientation!="Retraction fiber")
  
  df <- rbind(df, sub_retraction_fibers)
  return(df)
}

ratio_normalization_blebs <- function(B1int_GCX_ratio){

  B1int_GCX_ratio <- B1int_GCX_ratio %>% group_by(number, Signal_type) %>%  mutate(ratio_normalized_bleb = 
                                                                        ifelse(Orientation=="Bleb",
                                                                             ratio/ mean(ratio[Type=="Cell body"]),
                                                                             NA))
  
  return(B1int_GCX_ratio)
  
  
}

ratio_normalization_retractionfiber <- function(B1int_GCX_ratio){
  B1int_GCX_ratio <- B1int_GCX_ratio %>% group_by(number, Signal_type) %>%  mutate(ratio_normalized_retractionfiber = 
                                                                        ifelse(Orientation=="Retraction fiber",
                                                                               ratio/ mean(ratio[Type=="Cell body"]),
                                                                               NA))
  
  return(B1int_GCX_ratio)
}

ratio_normalization_leadingedgeprotrusion <- function(B1int_GCX_ratio, Xthr) {
  sub <- subset(B1int_GCX_ratio, Orientation == "Leading edge protrusion")
  sub <- sub %>%
    group_by(number, Position_relative_to_fiber, Signal_type) %>%
    mutate(ratio_normalized_leadingedgeprotrusion = ratio / mean(ratio[X>Xthr & X<20]))
  
  summary <- sub %>%
    group_by(number, Position_relative_to_fiber, Signal_type) %>%
    dplyr::summarise(ratio_normalized_leadingedgeprotrusion = mean(ratio_normalized_leadingedgeprotrusion, na.rm=TRUE)) %>% 
    filter(Signal_type=="foreground")
  
  print(summary)
  
  
  df <- subset(B1int_GCX_ratio, Orientation != "Leading edge protrusion")
  
  df <- df  %>%
    mutate(ratio_normalized_leadingedgeprotrusion = NA)
  
  df <- rbind(df, sub)
  return(df)
}

label_function <- function(labels) {
  paste0("Nested: ", labels)
}

calculate_residuals_and_tests <- function(x, y) {
  # Linear relationship
  linear_formula <- "y ~ x"
  linear_model <- lm(linear_formula, data = data.frame(x = x, y = y))
  linear_residuals <- residuals(linear_model)
  linear_jarque <- jarque.test(linear_residuals)
  linear_sw <- shapiro.test(linear_residuals)
  linear_r_squared <- summary(linear_model)$r.squared
  linear_adj_r_squared <- summary(linear_model)$adj.r.squared
  linear_aic <- AIC(linear_model)
  linear_bic <- BIC(linear_model)
  linear_coeffs <- coef(linear_model)
  linear_intercept <- linear_coeffs[1]
  linear_slope <- linear_coeffs[2]
  linear_intercept_p_value <- summary(linear_model)$coefficients[1, "Pr(>|t|)"]
  linear_slope_p_value <- summary(linear_model)$coefficients[2, "Pr(>|t|)"]
  linear_intercept_sig <- ifelse(linear_intercept_p_value < 0.05, "*", "n.s")
  linear_slope_sig <- ifelse(linear_slope_p_value < 0.05, "*", "n.s")
  
  # Exponential relationship
  exponential_formula <- "y ~ exp(x)"
  exponential_model <- lm(exponential_formula, data = data.frame(x = x, y = y))
  exponential_residuals <- residuals(exponential_model)
  exponential_jarque <- jarque.test(exponential_residuals)
  exponential_sw <- shapiro.test(exponential_residuals)
  exponential_r_squared <- summary(exponential_model)$r.squared
  exponential_adj_r_squared <- summary(exponential_model)$adj.r.squared
  exponential_aic <- AIC(exponential_model)
  exponential_bic <- BIC(exponential_model)
  exponential_coeffs <- coef(exponential_model)
  exponential_intercept <- exponential_coeffs[1]
  exponential_slope <- exponential_coeffs[2]
  exponential_intercept_p_value <- summary(exponential_model)$coefficients[1, "Pr(>|t|)"]
  exponential_slope_p_value <- summary(exponential_model)$coefficients[2, "Pr(>|t|)"]
  exponential_intercept_sig <- ifelse(exponential_intercept_p_value < 0.05, "*", "n.s")
  exponential_slope_sig <- ifelse(exponential_slope_p_value < 0.05, "*", "n.s")
  
  # Logarithmic relationship
  logarithmic_formula <- "y ~ log(x)"
  logarithmic_model <- lm(logarithmic_formula, data = data.frame(x = x, y = y))
  logarithmic_residuals <- residuals(logarithmic_model)
  logarithmic_jarque <- jarque.test(logarithmic_residuals)
  logarithmic_sw <- shapiro.test(logarithmic_residuals)
  logarithmic_r_squared <- summary(logarithmic_model)$r.squared
  logarithmic_adj_r_squared <- summary(logarithmic_model)$adj.r.squared
  logarithmic_aic <- AIC(logarithmic_model)
  logarithmic_bic <- BIC(logarithmic_model)
  logarithmic_coeffs <- coef(logarithmic_model)
  logarithmic_intercept <- logarithmic_coeffs[1]
  logarithmic_slope <- logarithmic_coeffs[2]
  logarithmic_intercept_p_value <- summary(logarithmic_model)$coefficients[1, "Pr(>|t|)"]
  logarithmic_slope_p_value <- summary(logarithmic_model)$coefficients[2, "Pr(>|t|)"]
  logarithmic_intercept_sig <- ifelse(logarithmic_intercept_p_value < 0.05, "*", "n.s")
  logarithmic_slope_sig <- ifelse(logarithmic_slope_p_value < 0.05, "*", "n.s")
  
  # Create a data frame to store the results
  results_df <- data.frame(
    Formula = c(linear_formula, exponential_formula, logarithmic_formula),
    JB_p_value = c(linear_jarque$p.value, exponential_jarque$p.value, logarithmic_jarque$p.value),
    SW_p_value = c(linear_sw$p.value, exponential_sw$p.value, logarithmic_sw$p.value),
    R_squared = c(linear_r_squared, exponential_r_squared, logarithmic_r_squared),
    Adj_R_squared = c(linear_adj_r_squared, exponential_adj_r_squared, logarithmic_adj_r_squared),
    AIC = c(linear_aic, exponential_aic, logarithmic_aic),
    BIC = c(linear_bic, exponential_bic, logarithmic_bic),
    Coefficient = c(linear_slope, exponential_slope, logarithmic_slope),
    Intercept = c(linear_intercept, exponential_intercept, logarithmic_intercept),
    Coefficient_Significance = c(linear_slope_sig, exponential_slope_sig, logarithmic_slope_sig),
    Intercept_Significance = c(linear_intercept_sig, exponential_intercept_sig, logarithmic_intercept_sig),
    Coefficient_p_value = c(linear_slope_p_value, exponential_slope_p_value, logarithmic_slope_p_value),
    Intercept_p_value = c(linear_intercept_p_value, exponential_intercept_p_value, logarithmic_intercept_p_value),
    stringsAsFactors = FALSE
  )
  
  return(results_df)
}

calculate_residuals_and_tests_limitcorrected <- function(x, y) {
  library(tseries)
  
  # Linear relationship
  linear_formula <- "y ~ x"
  linear_model <- lm(linear_formula, data = data.frame(x = x, y = y))
  linear_residuals <- residuals(linear_model)
  linear_jarque <- jarque.bera.test(linear_residuals)
  linear_sw <- shapiro.test(linear_residuals)
  linear_r_squared <- summary(linear_model)$r.squared
  linear_adj_r_squared <- summary(linear_model)$adj.r.squared
  linear_aic <- AIC(linear_model)
  linear_bic <- BIC(linear_model)
  linear_coeffs <- coef(linear_model)
  linear_intercept <- linear_coeffs[1]
  linear_slope <- linear_coeffs[2]
  linear_intercept_p_value <- summary(linear_model)$coefficients[1, "Pr(>|t|)"]
  linear_slope_p_value <- summary(linear_model)$coefficients[2, "Pr(>|t|)"]
  linear_intercept_sig <- ifelse(linear_intercept_p_value < 0.05, "*", "n.s")
  linear_slope_sig <- ifelse(linear_slope_p_value < 0.05, "*", "n.s")
  
  # Exponential relationship
  exponential_formula <- "y ~ exp(x)"
  exponential_model <- NULL
  exponential_residuals <- NULL
  exponential_jarque <- NULL
  exponential_sw <- NULL
  exponential_r_squared <- NULL
  exponential_adj_r_squared <- NULL
  exponential_aic <- NULL
  exponential_bic <- NULL
  exponential_coeffs <- NULL
  exponential_intercept <- NULL
  exponential_slope <- NULL
  exponential_intercept_p_value <- NULL
  exponential_slope_p_value <- NULL
  exponential_intercept_sig <- NULL
  exponential_slope_sig <- NULL
  
  if (any(!is.na(x) & is.finite(x) & x > log(.Machine$integer.max))) {
    warning("Exponential regression skipped due to high values.")
  } else {
    exponential_formula <- "y ~ exp(x)"
    exponential_model <- lm(exponential_formula, data = data.frame(x = x, y = y))
    exponential_residuals <- residuals(exponential_model)
    exponential_jarque <- jarque.bera.test(exponential_residuals)
    exponential_sw <- shapiro.test(exponential_residuals)
    exponential_r_squared <- summary(exponential_model)$r.squared
    exponential_adj_r_squared <- summary(exponential_model)$adj.r.squared
    exponential_aic <- AIC(exponential_model)
    exponential_bic <- BIC(exponential_model)
    exponential_coeffs <- coef(exponential_model)
    exponential_intercept <- exponential_coeffs[1]
    exponential_slope <- exponential_coeffs[2]
    exponential_intercept_p_value <- summary(exponential_model)$coefficients[1, "Pr(>|t|)"]
    exponential_slope_p_value <- summary(exponential_model)$coefficients[2, "Pr(>|t|)"]
    exponential_intercept_sig <- ifelse(exponential_intercept_p_value < 0.05, "*", "n.s")
    exponential_slope_sig <- ifelse(exponential_slope_p_value < 0.05, "*", "n.s")
  }
  
  # Logarithmic relationship
  logarithmic_formula <- "y ~ log(x)"
  logarithmic_model <- lm(logarithmic_formula, data = data.frame(x = x, y = y))
  logarithmic_residuals <- residuals(logarithmic_model)
  logarithmic_jarque <- jarque.bera.test(logarithmic_residuals)
  logarithmic_sw <- shapiro.test(logarithmic_residuals)
  logarithmic_r_squared <- summary(logarithmic_model)$r.squared
  logarithmic_adj_r_squared <- summary(logarithmic_model)$adj.r.squared
  logarithmic_aic <- AIC(logarithmic_model)
  logarithmic_bic <- BIC(logarithmic_model)
  logarithmic_coeffs <- coef(logarithmic_model)
  logarithmic_intercept <- logarithmic_coeffs[1]
  logarithmic_slope <- logarithmic_coeffs[2]
  logarithmic_intercept_p_value <- summary(logarithmic_model)$coefficients[1, "Pr(>|t|)"]
  logarithmic_slope_p_value <- summary(logarithmic_model)$coefficients[2, "Pr(>|t|)"]
  logarithmic_intercept_sig <- ifelse(logarithmic_intercept_p_value < 0.05, "*", "n.s")
  logarithmic_slope_sig <- ifelse(logarithmic_slope_p_value < 0.05, "*", "n.s")
  
  # Square root relationship (y ~ sqrt(x))
  sqrt_formula <- "y ~ sqrt(x)"
  sqrt_model <- lm(sqrt_formula, data = data.frame(x = x, y = y))
  sqrt_residuals <- residuals(sqrt_model)
  sqrt_jarque <- jarque.bera.test(sqrt_residuals)
  sqrt_sw <- shapiro.test(sqrt_residuals)
  sqrt_r_squared <- summary(sqrt_model)$r.squared
  sqrt_adj_r_squared <- summary(sqrt_model)$adj.r.squared
  sqrt_aic <- AIC(sqrt_model)
  sqrt_bic <- BIC(sqrt_model)
  sqrt_coeffs <- coef(sqrt_model)
  sqrt_intercept <- sqrt_coeffs[1]
  sqrt_slope <- sqrt_coeffs[2]
  sqrt_intercept_p_value <- summary(sqrt_model)$coefficients[1, "Pr(>|t|)"]
  sqrt_slope_p_value <- summary(sqrt_model)$coefficients[2, "Pr(>|t|)"]
  sqrt_intercept_sig <- ifelse(sqrt_intercept_p_value < 0.05, "*", "n.s")
  sqrt_slope_sig <- ifelse(sqrt_slope_p_value < 0.05, "*", "n.s")
  
  # x/sqrt(x) relationship (y ~ x/sqrt(x))
  x_div_sqrtx_formula <- "y ~ x/sqrt(x)"
  x_div_sqrtx_model <- lm(x_div_sqrtx_formula, data = data.frame(x = x, y = y))
  x_div_sqrtx_residuals <- residuals(x_div_sqrtx_model)
  x_div_sqrtx_jarque <- jarque.bera.test(x_div_sqrtx_residuals)
  x_div_sqrtx_sw <- shapiro.test(x_div_sqrtx_residuals)
  x_div_sqrtx_r_squared <- summary(x_div_sqrtx_model)$r.squared
  x_div_sqrtx_adj_r_squared <- summary(x_div_sqrtx_model)$adj.r.squared
  x_div_sqrtx_aic <- AIC(x_div_sqrtx_model)
  x_div_sqrtx_bic <- BIC(x_div_sqrtx_model)
  x_div_sqrtx_coeffs <- coef(x_div_sqrtx_model)
  x_div_sqrtx_intercept <- x_div_sqrtx_coeffs[1]
  x_div_sqrtx_slope <- x_div_sqrtx_coeffs[2]
  x_div_sqrtx_intercept_p_value <- summary(x_div_sqrtx_model)$coefficients[1, "Pr(>|t|)"]
  x_div_sqrtx_slope_p_value <- summary(x_div_sqrtx_model)$coefficients[2, "Pr(>|t|)"]
  x_div_sqrtx_intercept_sig <- ifelse(x_div_sqrtx_intercept_p_value < 0.05, "*", "n.s")
  x_div_sqrtx_slope_sig <- ifelse(x_div_sqrtx_slope_p_value < 0.05, "*", "n.s")
  
  # Create a data frame to store the results
  results_df <- data.frame(
    Formula = c(linear_formula, exponential_formula, logarithmic_formula, sqrt_formula, x_div_sqrtx_formula),
    JB_p_value = c(linear_jarque$p.value, ifelse(is.null(exponential_jarque), NA, exponential_jarque$p.value), logarithmic_jarque$p.value, sqrt_jarque$p.value, x_div_sqrtx_jarque$p.value),
    SW_p_value = c(linear_sw$p.value, ifelse(is.null(exponential_sw), NA, exponential_sw$p.value), logarithmic_sw$p.value, sqrt_sw$p.value, x_div_sqrtx_sw$p.value),
    R_squared = c(linear_r_squared, ifelse(is.null(exponential_r_squared), NA, exponential_r_squared), logarithmic_r_squared, sqrt_r_squared, x_div_sqrtx_r_squared),
    Adj_R_squared = c(linear_adj_r_squared, ifelse(is.null(exponential_adj_r_squared), NA, exponential_adj_r_squared), logarithmic_adj_r_squared, sqrt_adj_r_squared, x_div_sqrtx_adj_r_squared),
    AIC = c(linear_aic, ifelse(is.null(exponential_aic), NA, exponential_aic), logarithmic_aic, sqrt_aic, x_div_sqrtx_aic),
    BIC = c(linear_bic, ifelse(is.null(exponential_bic), NA, exponential_bic), logarithmic_bic, sqrt_bic, x_div_sqrtx_bic),
    Coefficient = c(linear_slope, ifelse(is.null(exponential_slope), NA, exponential_slope), logarithmic_slope, sqrt_slope, x_div_sqrtx_slope),
    Intercept = c(linear_intercept, ifelse(is.null(exponential_intercept), NA, exponential_intercept), logarithmic_intercept, sqrt_intercept, x_div_sqrtx_intercept),
    Coefficient_Significance = c(linear_slope_sig, ifelse(is.null(exponential_slope_sig), NA, exponential_slope_sig), logarithmic_slope_sig, sqrt_slope_sig, x_div_sqrtx_slope_sig),
    Intercept_Significance = c(linear_intercept_sig, ifelse(is.null(exponential_intercept_sig), NA, exponential_intercept_sig), logarithmic_intercept_sig, sqrt_intercept_sig, x_div_sqrtx_intercept_sig),
    Coefficient_p_value = c(linear_slope_p_value, ifelse(is.null(exponential_slope_p_value), NA, exponential_slope_p_value), logarithmic_slope_p_value, sqrt_slope_p_value, x_div_sqrtx_slope_p_value),
    Intercept_p_value = c(linear_intercept_p_value, ifelse(is.null(exponential_intercept_p_value), NA, exponential_intercept_p_value), logarithmic_intercept_p_value, sqrt_intercept_p_value, x_div_sqrtx_intercept_p_value),
    stringsAsFactors = FALSE
  )
  
  return(results_df)
}

calculate_residuals_and_tests_limitcorrected_noexponential <- function(x, y) {
  library(tseries)
  
  # Linear relationship
  linear_formula <- "y ~ x"
  linear_model <- lm(linear_formula, data = data.frame(x = x, y = y))
  linear_residuals <- residuals(linear_model)
  linear_jarque <- jarque.bera.test(linear_residuals)
  linear_sw <- shapiro.test(linear_residuals)
  linear_r_squared <- summary(linear_model)$r.squared
  linear_adj_r_squared <- summary(linear_model)$adj.r.squared
  linear_aic <- AIC(linear_model)
  linear_bic <- BIC(linear_model)
  linear_coeffs <- coef(linear_model)
  linear_intercept <- linear_coeffs[1]
  linear_slope <- linear_coeffs[2]
  linear_intercept_p_value <- summary(linear_model)$coefficients[1, "Pr(>|t|)"]
  linear_slope_p_value <- summary(linear_model)$coefficients[2, "Pr(>|t|)"]
  linear_intercept_sig <- ifelse(linear_intercept_p_value < 0.05, "*", "n.s")
  linear_slope_sig <- ifelse(linear_slope_p_value < 0.05, "*", "n.s")
  
  # Logarithmic relationship
  logarithmic_formula <- "y ~ log(x)"
  logarithmic_model <- lm(logarithmic_formula, data = data.frame(x = x+1e-10, y = y))
  logarithmic_residuals <- residuals(logarithmic_model)
  logarithmic_jarque <- jarque.bera.test(logarithmic_residuals)
  logarithmic_sw <- shapiro.test(logarithmic_residuals)
  logarithmic_r_squared <- summary(logarithmic_model)$r.squared
  logarithmic_adj_r_squared <- summary(logarithmic_model)$adj.r.squared
  logarithmic_aic <- AIC(logarithmic_model)
  logarithmic_bic <- BIC(logarithmic_model)
  logarithmic_coeffs <- coef(logarithmic_model)
  logarithmic_intercept <- logarithmic_coeffs[1]
  logarithmic_slope <- logarithmic_coeffs[2]
  logarithmic_intercept_p_value <- summary(logarithmic_model)$coefficients[1, "Pr(>|t|)"]
  logarithmic_slope_p_value <- summary(logarithmic_model)$coefficients[2, "Pr(>|t|)"]
  logarithmic_intercept_sig <- ifelse(logarithmic_intercept_p_value < 0.05, "*", "n.s")
  logarithmic_slope_sig <- ifelse(logarithmic_slope_p_value < 0.05, "*", "n.s")
  
  # Square root relationship (y ~ sqrt(x))
  sqrt_formula <- "y ~ sqrt(x)"
  sqrt_model <- lm(sqrt_formula, data = data.frame(x = x, y = y))
  sqrt_residuals <- residuals(sqrt_model)
  sqrt_jarque <- jarque.bera.test(sqrt_residuals)
  sqrt_sw <- shapiro.test(sqrt_residuals)
  sqrt_r_squared <- summary(sqrt_model)$r.squared
  sqrt_adj_r_squared <- summary(sqrt_model)$adj.r.squared
  sqrt_aic <- AIC(sqrt_model)
  sqrt_bic <- BIC(sqrt_model)
  sqrt_coeffs <- coef(sqrt_model)
  sqrt_intercept <- sqrt_coeffs[1]
  sqrt_slope <- sqrt_coeffs[2]
  sqrt_intercept_p_value <- summary(sqrt_model)$coefficients[1, "Pr(>|t|)"]
  sqrt_slope_p_value <- summary(sqrt_model)$coefficients[2, "Pr(>|t|)"]
  sqrt_intercept_sig <- ifelse(sqrt_intercept_p_value < 0.05, "*", "n.s")
  sqrt_slope_sig <- ifelse(sqrt_slope_p_value < 0.05, "*", "n.s")
  
  # x/sqrt(x) relationship (y ~ x/sqrt(x))
  x_div_sqrtx_formula <- "y ~ x/sqrt(x)"
  x_div_sqrtx_model <- lm(x_div_sqrtx_formula, data = data.frame(x = x, y = y))
  x_div_sqrtx_residuals <- residuals(x_div_sqrtx_model)
  x_div_sqrtx_jarque <- jarque.bera.test(x_div_sqrtx_residuals)
  x_div_sqrtx_sw <- shapiro.test(x_div_sqrtx_residuals)
  x_div_sqrtx_r_squared <- summary(x_div_sqrtx_model)$r.squared
  x_div_sqrtx_adj_r_squared <- summary(x_div_sqrtx_model)$adj.r.squared
  x_div_sqrtx_aic <- AIC(x_div_sqrtx_model)
  x_div_sqrtx_bic <- BIC(x_div_sqrtx_model)
  x_div_sqrtx_coeffs <- coef(x_div_sqrtx_model)
  x_div_sqrtx_intercept <- x_div_sqrtx_coeffs[1]
  x_div_sqrtx_slope <- x_div_sqrtx_coeffs[2]
  x_div_sqrtx_intercept_p_value <- summary(x_div_sqrtx_model)$coefficients[1, "Pr(>|t|)"]
  x_div_sqrtx_slope_p_value <- summary(x_div_sqrtx_model)$coefficients[2, "Pr(>|t|)"]
  x_div_sqrtx_intercept_sig <- ifelse(x_div_sqrtx_intercept_p_value < 0.05, "*", "n.s")
  x_div_sqrtx_slope_sig <- ifelse(x_div_sqrtx_slope_p_value < 0.05, "*", "n.s")
  
  # Create a data frame to store the results
  results_df <- data.frame(
    Formula = c(linear_formula, logarithmic_formula, sqrt_formula, x_div_sqrtx_formula),
    JB_p_value = c(linear_jarque$p.value, logarithmic_jarque$p.value, sqrt_jarque$p.value, x_div_sqrtx_jarque$p.value),
    SW_p_value = c(linear_sw$p.value, logarithmic_sw$p.value, sqrt_sw$p.value, x_div_sqrtx_sw$p.value),
    R_squared = c(linear_r_squared, logarithmic_r_squared, sqrt_r_squared, x_div_sqrtx_r_squared),
    Adj_R_squared = c(linear_adj_r_squared, logarithmic_adj_r_squared, sqrt_adj_r_squared, x_div_sqrtx_adj_r_squared),
    AIC = c(linear_aic, logarithmic_aic, sqrt_aic, x_div_sqrtx_aic),
    BIC = c(linear_bic, logarithmic_bic, sqrt_bic, x_div_sqrtx_bic),
    Coefficient = c(linear_slope, logarithmic_slope, sqrt_slope, x_div_sqrtx_slope),
    Intercept = c(linear_intercept, logarithmic_intercept, sqrt_intercept, x_div_sqrtx_intercept),
    Coefficient_Significance = c(linear_slope_sig, logarithmic_slope_sig, sqrt_slope_sig, x_div_sqrtx_slope_sig),
    Intercept_Significance = c(linear_intercept_sig, logarithmic_intercept_sig, sqrt_intercept_sig, x_div_sqrtx_intercept_sig),
    Coefficient_p_value = c(linear_slope_p_value, logarithmic_slope_p_value, sqrt_slope_p_value, x_div_sqrtx_slope_p_value),
    Intercept_p_value = c(linear_intercept_p_value, logarithmic_intercept_p_value, sqrt_intercept_p_value, x_div_sqrtx_intercept_p_value),
    stringsAsFactors = FALSE
  )
  
  return(results_df)
}

custom_lineplot <- function(data, x, y, Orientation, Type, Fiber_type){
  
  distvsratio_protrusion <- ggplot(data=protrusion_data,
                                   aes_string(x=x, y=y)) +
    geom_point(size=1, stroke=0, alpha=0.5) +
    geom_smooth(method="lm", formula=y~log(x), fullrange=TRUE, se=TRUE, color="black", linewidth=0.5)+
    labs(x="Distance from most distal point -> cell body (micron)") +
    scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
    scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
    theme(legend.position = c(0.7,0.7)) +
    #geom_vline(xintercept=Xthreshold, linetype=2) +
    
    theme(axis.text.y = element_text(size = 8),
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
          axis.text = element_text(size = 8),
          legend.position = c(0.85, 0.55),
          legend.justification = c(1, 0),
          legend.box.just = "right",
          legend.title=element_blank())
  
  distvsratio_protrusion <- distvsratio_protrusion +
    facet_grid(Fiber_type ~Orientation+Type)
  
  
  
  
  print(distvsratio_protrusion)
  
  return(distvsratio_protrusion)
}

custom_lineplot_with_rsquared <- function(data, x_var, y_var, Orientation, Type, Fiber_type) {
  library(ggplot2)
  library(dplyr)
  # Your existing ggplot code
  distvsratio_protrusion <- ggplot(data, aes(x = {{x_var}}, y = {{y_var}})) +
    geom_point(size = 1, stroke = 0, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ log(x), fullrange = TRUE, se = TRUE, color = "black", linewidth = 0.5) +
    labs(x = "Distance from most distal point -> cell body (micron)") +
    scale_colour_viridis_d(option = "turbo", direction = -1, begin = 0.2, end = 0.8) +
    scale_fill_viridis_d(option = "turbo", direction = -1, begin = 0.2, end = 0.8) +
    theme(legend.position = c(0.7, 0.7)) +
    theme(
      axis.text.y = element_text(size = 8),
      text = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(c(2, 1), "mm"),
      axis.text = element_text(size = 8),
      legend.position = c(0.85, 0.55),
      legend.justification = c(1, 0),
      legend.box.just = "right",
      legend.title = element_blank()
    ) +
    facet_grid(Fiber_type ~ Orientation + Type)
  
  # Calculate R-squared values per facet
  rsquared_p_values <- data %>%
    group_by(Fiber_type, Orientation, Type) %>%
    summarize(
      rsquared = summary(lm({{y_var}} ~ log({{x_var}})))$r.squared,
      adj_r_squared = 1 - (1 - rsquared) * (length(unique({{x_var}})) - 1) / (length({{x_var}}) - length(coef(lm({{y_var}} ~ log({{x_var}}))))),
      p_value = coef(summary(lm({{y_var}} ~ log({{x_var}}))))[, "Pr(>|t|)"]
    ) %>%
    filter(p_value < 0.05) %>%
    ungroup()
  
  # Add R-squared annotations to the original plot
  final_plot <- distvsratio_protrusion +
    geom_text(data = rsquared_p_values, aes(label = sprintf("R-squared = %.2f", adj_r_squared)), 
              x = Inf, y = -Inf, hjust = 1, vjust = 0, size=3)
  
  print(final_plot)
  
  return(final_plot)
  
}

changePlotSize <- function(plot) {
  
  plotname <- readline(prompt = "Enter plot name: ")
  
  # Open a new graphics device
  dev.new()
  
  # Set the size of the plot
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
  
  # Display the plot
  print(plot)
  
  # Interactively adjust the plot size in the "Plots" pane
  tryCatch({
    locator(1)
  }, error = function(e) {
    # Handle the error when the plot window is closed
    warning("Plot window closed. Operation aborted.")
    return(invisible())
  })
  
  
  # Save the modified plot as a PDF
  dev.copy2pdf(file = paste0(plotname,".pdf"))
  
  # Close the graphics device
  print(plot)
  saveRDS(plot, file=paste0(plotname, ".RData"))
  dev.off()
}

performTwoWayANOVA <- function(data, responseVar, factor1, factor2) {
  anova_model <- aov(reformulate(paste(factor1, factor2, sep = "*"), responseVar), data = data)
  print(summary(anova_model))
  
  posthoc <- TukeyHSD(anova_model)
  print(posthoc)
  
  return(list(anova = anova_model, posthoc = posthoc))
}

custom_plot <- function(bleb_ratio, x, y, Orientation, B1int_status){
  distvsratio_bleb_example <-  ggplot(subset(bleb_ratio, B1int_status!="Unchanged"), aes_string(x=x, y=y))+
    #geom_point(size=2, stroke=0, alpha=0.5,aes(color=Type))+
    geom_line(size=0.5, aes(group=number, color=Type))+
    
    #geom_smooth(method="lm", formula = y~(x), fullrange=TRUE, se=TRUE) +
    
    #labs(y="B1 integrin / Glycocalyx ratio") +
    
    labs(x="Relative distance (fraction)")+
    scale_x_continuous(breaks=c(0,0.5,1))+
    
    
    scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    theme(axis.text.y = element_text(size = 8),
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
          axis.text = element_text(size = 8),
          legend.position = c(0.95, 0.8),
          legend.title=element_blank())
  
  
  distvsratio_bleb_example <- distvsratio_bleb_example +
    facet_nested(. ~ Orientation+B1int_status) +
    geom_vline(xintercept=0.33, linetype=2, size=0.5)+custom_theme
  
  
  print(distvsratio_bleb_example)
  return(distvsratio_bleb_example)
}


# Define new transform
my_transform <- function(threshold = 25, squeeze_factor = 10) {
  force(threshold)
  force(squeeze_factor)
  my_transform <- trans_new(
    name = "trans_squeeze",
    transform = function(x) {
      ifelse(x > threshold, 
             ((x - threshold) * (1 / squeeze_factor)) + threshold, 
             x) 
    },
    inverse = function(x) {
      ifelse(x > threshold, 
             ((x - threshold) * squeeze_factor) + threshold, 
             x)
    }
  )
  return(my_transform)
}


custom_lineplot_with_rsquared <- function(data, x_var, y_var, 
                                          regression_type = "linear",
                                          color_var = NULL,
                                          group_vars = NULL,
                                          alpha = 0.5, size = 1) {
  library(ggplot2)
  library(dplyr)
  print("Function started")
  
  # Convert x_var and y_var to symbols
  x_var_sym <- rlang::ensym(x_var)
  y_var_sym <- rlang::ensym(y_var)
  
  # Use aes_string instead of aes
  distvsratio_protrusion <- ggplot(data, aes_string(x = rlang::as_string(x_var_sym), y = rlang::as_string(y_var_sym)))
  
  if (!is.null(color_var)) {
    distvsratio_protrusion <- distvsratio_protrusion +
      geom_point(shape=16, size = size, stroke = 0, alpha = alpha, aes_string(color = color_var))+
      scale_colour_gradientn(colours = viridis(50, direction=-1, begin=0.05, end=0.95))+
      scale_fill_gradientn(colours = viridis(50, direction=-1, begin=0.05, end=0.95))
  } else {
    distvsratio_protrusion <- distvsratio_protrusion +
      geom_point(size = size, stroke = 0, alpha = alpha)
  }
  
  distvsratio_protrusion <- distvsratio_protrusion +
    theme(legend.position = c(0.7, 0.7)) +
    theme(
      axis.text.y = element_text(size = 8),
      text = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(c(2, 1), "mm"),
      axis.text = element_text(size = 8),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank()
    )
  print(distvsratio_protrusion)
  
  # Construct the facet_nested statement
  if (!is.null(group_vars)) {
    facet_nested_expr <- paste(".~", paste(group_vars, collapse = "+"))
    distvsratio_protrusion <- distvsratio_protrusion + facet_nested(as.formula(facet_nested_expr))
  }
  
  # Calculate R-squared values per facet
  x_var_str <- deparse(substitute(x_var))
  y_var_str <- deparse(substitute(y_var))
  
  lm_formula <- switch(
    regression_type,
    "linear" = as.formula(paste(y_var_str, "~", x_var_str)),
    "log" = as.formula(paste(y_var_str, "~", "log", "(", x_var_str, ")")),
    "exp" = as.formula(paste(y_var_str, "~", "exp(", x_var_str, ")"))
  )
  
  print(lm_formula)
  rsquared_p_values <- data %>%
    filter(!is.na(.data[[x_var_str]]) & !is.infinite(.data[[x_var_str]]) & !is.nan(.data[[x_var_str]]) &
             !is.na(.data[[y_var_str]]) & !is.infinite(.data[[y_var_str]]) & !is.nan(.data[[y_var_str]]) &
             if (regression_type == "exp") .data[[x_var_str]] < 700 else TRUE) %>%
    group_by_at(group_vars) %>%
    summarize(
      rsquared = summary(lm(lm_formula, data = .))$r.squared,
      adj_r_squared = 1 - (1 - rsquared) * (length(unique(.data[[x_var_str]])) - 1) / (length(.data[[x_var_str]]) - length(coef(lm(lm_formula, data = .)))),
      p_values = coef(summary(lm(lm_formula, data = .)))[, "Pr(>|t|)"],
      values = coef(summary(lm(lm_formula, data = .)))[, "Estimate"],
      names = c("Intercept", "Coefficient")
    ) %>%
    ungroup()
  
  rsquared_p_values$regression_type <- rep(regression_type, nrow(rsquared_p_values))
  
  print(rsquared_p_values)
  
  # Adjust p-values using Bonferroni correction
  num_comparisons <- nrow(rsquared_p_values)  # Assuming each row corresponds to a comparison
  rsquared_p_values$p_values_adj <- p.adjust(rsquared_p_values$p_values, method = "bonferroni", n = num_comparisons)
  
  categorize_significance <- function(p_value_adj) {
    if (is.na(p_value_adj) || !is.numeric(p_value_adj)) {
      return("NA")
    } else if (p_value_adj < 0.001) {
      return("***")
    } else if (p_value_adj < 0.01) {
      return("**")
    } else if (p_value_adj < 0.05) {
      return("*")
    } else {
      return("N.S")
    }
  }
  
  # Apply the function to the rsquared_p_values dataframe
  rsquared_p_values <- rsquared_p_values %>%
    mutate(
      significance_intercept = categorize_significance(p_values_adj[1]),
      significance_coefficient = categorize_significance(p_values_adj[2])
    )
  
  print(rsquared_p_values)
  
  # Prepare a dataframe for annotations
  annotations_df <- rsquared_p_values %>%
    mutate(annotation_text = paste(
      sprintf("R^2= %.2f", adj_r_squared[2]),
      #sprintf("p.adj.coef= %s", significance_coefficient[2]),
      sprintf("P.adj= %.2e", p_values_adj[2]),
      sep = "\n"),
      x_coord = Inf, y_coord = Inf)
  
  print(annotations_df)  
  
  # Create a dataframe from the rsquared_p_values
  rsquared_values_df <- as.data.frame(rsquared_p_values)
  
  basic_formula <- switch(
    regression_type,
    "linear" = as.formula(paste("y ~ x")),
    "log" = as.formula(paste("y ~ log(x)")),
    "exp" = as.formula(paste("y ~ exp(x)"))
  )
  print(basic_formula)
  
  # Add linear regression line if significant
  if (regression_type %in% c("linear", "log", "exp")) {
    if (regression_type == "exp") {
      final_plot <- distvsratio_protrusion +
        geom_smooth(data = filter(data, .data[[x_var_str]] < 700),
                    method = "lm", formula = basic_formula, fullrange = TRUE, se = TRUE, color = "black", linewidth = 0.5, span = 0.1)
    } else {
      final_plot <- distvsratio_protrusion +
        geom_smooth(data = data,
                    method = "lm", formula = basic_formula, fullrange = TRUE, se = TRUE, color = "black", linewidth = 0.5, span = 0.1)
    }
  }
  
  # Add annotations
  final_plot <- final_plot +
    geom_text(data = annotations_df, aes(x = x_coord, y = y_coord, label = annotation_text),
              hjust = 1.2, vjust = 1.2, size = 3)
  
  # Return both the plot and the R-squared values dataframe
  return(list(plot = final_plot, rsquared_values = rsquared_values_df))
}

custom_lineplot_with_rsquared <- function(data, x_var, y_var, 
                                          regression_type = "linear",
                                          color_var = NULL,
                                          group_vars = NULL,
                                          alpha = 0.5, size = 1,
                                          add.extra.plot = NULL,
                                          annotate_signif=TRUE) {
  library(ggplot2)
  library(dplyr)
  print("Function started")
  
  # Convert x_var and y_var to symbols
  x_var_sym <- rlang::ensym(x_var)
  y_var_sym <- rlang::ensym(y_var)
  
  # Use aes_string instead of aes
  distvsratio_protrusion <- ggplot(data, aes_string(x = rlang::as_string(x_var_sym), y = rlang::as_string(y_var_sym)))+
    theme_classic()
  
  if(!is.null(add.extra.plot)){
    distvsratio_protrusion <- distvsratio_protrusion + add.extra.plot
  }
  
  if (!is.null(color_var)) {
    distvsratio_protrusion <- distvsratio_protrusion +
      geom_point(shape=16, size = size, stroke = 0, alpha = alpha, aes_string(color = color_var))+
      scale_colour_gradientn(colours = viridis(50, direction=-1, begin=0.05, end=0.95))+
      scale_fill_gradientn(colours = viridis(50, direction=-1, begin=0.05, end=0.95))
  } else {
    distvsratio_protrusion <- distvsratio_protrusion +
      geom_point(size = size, stroke = 0, alpha = alpha)
  }
  
  distvsratio_protrusion <- distvsratio_protrusion +
    theme(legend.position = c(0.7, 0.7)) +
    theme(
      axis.text.y = element_text(size = 8),
      text = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(c(2, 1), "mm"),
      axis.text = element_text(size = 8),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank()
    )
  print(distvsratio_protrusion)
  
  # Construct the facet_nested statement
  if (!is.null(group_vars)) {
    facet_nested_expr <- paste(".~", paste(group_vars, collapse = "+"))
    distvsratio_protrusion <- distvsratio_protrusion + facet_nested(as.formula(facet_nested_expr))
  }
  
  # Calculate R-squared values per facet
  x_var_str <- deparse(substitute(x_var))
  y_var_str <- deparse(substitute(y_var))
  
  lm_formula <- switch(
    regression_type,
    "linear" = as.formula(paste(y_var_str, "~", x_var_str)),
    "log" = as.formula(paste(y_var_str, "~", "log", "(", x_var_str, ")")),
    "exp" = as.formula(paste(y_var_str, "~", "exp(", x_var_str, ")"))
  )
  
  print(lm_formula)
  rsquared_p_values <- data %>%
    filter(!is.na(.data[[x_var_str]]) & !is.infinite(.data[[x_var_str]]) & !is.nan(.data[[x_var_str]]) &
             !is.na(.data[[y_var_str]]) & !is.infinite(.data[[y_var_str]]) & !is.nan(.data[[y_var_str]]) &
             if (regression_type == "exp") .data[[x_var_str]] < 700 else TRUE) %>%
    group_by_at(group_vars) %>%
    summarize(
      rsquared = summary(lm(lm_formula, data = .))$r.squared,
      adj_r_squared = 1 - (1 - rsquared) * (length(unique(.data[[x_var_str]])) - 1) / (length(.data[[x_var_str]]) - length(coef(lm(lm_formula, data = .)))),
      p_values = coef(summary(lm(lm_formula, data = .)))[, "Pr(>|t|)"],
      values = coef(summary(lm(lm_formula, data = .)))[, "Estimate"],
      names = c("Intercept", "Coefficient")
    ) %>%
    ungroup()
  
  rsquared_p_values$regression_type <- rep(regression_type, nrow(rsquared_p_values))
  
  print(rsquared_p_values)
  
  # Adjust p-values using Bonferroni correction
  num_comparisons <- nrow(rsquared_p_values)  # Assuming each row corresponds to a comparison
  rsquared_p_values$p_values_adj <- p.adjust(rsquared_p_values$p_values, method = "bonferroni", n = num_comparisons)
  
  categorize_significance <- function(p_value_adj) {
    if (is.na(p_value_adj) || !is.numeric(p_value_adj)) {
      return("NA")
    } else if (p_value_adj < 0.001) {
      return("***")
    } else if (p_value_adj < 0.01) {
      return("**")
    } else if (p_value_adj < 0.05) {
      return("*")
    } else {
      return("N.S")
    }
  }
  
  # Apply the function to the rsquared_p_values dataframe
  rsquared_p_values <- rsquared_p_values %>%
    mutate(
      significance_intercept = categorize_significance(p_values_adj[1]),
      significance_coefficient = categorize_significance(p_values_adj[2])
    )
  
  print(rsquared_p_values)
  
  # Prepare a dataframe for annotations
  annotations_df <- rsquared_p_values %>%
    mutate(annotation_text = paste(
      sprintf("R^2= %.2f", adj_r_squared[2]),
      #sprintf("p.adj.coef= %s", significance_coefficient[2]),
      sprintf("P.adj= %.2e", p_values_adj[2]),
      sep = "\n"),
      x_coord = Inf ,
      y_coord = Inf
    )
  
  print(annotations_df)  
  
  # Create a dataframe from the rsquared_p_values
  rsquared_values_df <- as.data.frame(rsquared_p_values)
  
  basic_formula <- switch(
    regression_type,
    "linear" = as.formula(paste("y ~ x")),
    "log" = as.formula(paste("y ~ log(x)")),
    "exp" = as.formula(paste("y ~ exp(x)"))
  )
  print(basic_formula)
  
  # Add linear regression line if significant
  if (regression_type %in% c("linear", "log", "exp")) {
    if (regression_type == "exp") {
      final_plot <- distvsratio_protrusion +
        geom_smooth(data = filter(data, .data[[x_var_str]] < 700),
                    method = "lm", formula = basic_formula, fullrange = FALSE, se = TRUE, color = "black", linewidth = 0.5, span = 0.1)
    } else {
      final_plot <- distvsratio_protrusion +
        geom_smooth(data = data,
                    method = "lm", formula = basic_formula, fullrange = FALSE, se = TRUE, color = "black", linewidth = 0.5, span = 0.1)
    }
  }
  
  if(annotate_signif==TRUE){
  # Add annotations
  final_plot <- final_plot +
    geom_text(data = annotations_df, aes(x = x_coord, y = y_coord, label = annotation_text),
              hjust = 1.2, vjust = 1.2, size = 3)
  }
  # Return both the plot and the R-squared values dataframe
  return(list(plot = final_plot, rsquared_values = rsquared_values_df))
}

custom_lineplot_with_rsquared <- function(data, x_var, y_var, 
                                          regression_type = "linear",
                                          color_var = NULL,
                                          group_vars = NULL,
                                          alpha = 0.5, size = 1,
                                          add.extra.plot = NULL,
                                          annotate_signif=TRUE,
                                          shape_var=NULL) {
  library(ggplot2)
  library(dplyr)
  print("Function started")
  
  # Convert x_var and y_var to symbols
  x_var_sym <- rlang::ensym(x_var)
  y_var_sym <- rlang::ensym(y_var)
  
  # Use aes_string instead of aes
  distvsratio_protrusion <- ggplot(data, aes_string(x = rlang::as_string(x_var_sym), y = rlang::as_string(y_var_sym)))+
    theme_classic()
  
  if(!is.null(add.extra.plot)){
    distvsratio_protrusion <- distvsratio_protrusion + add.extra.plot
  }
  
  
  if (!is.null(color_var)) {
    distvsratio_protrusion <- distvsratio_protrusion +
      geom_point(shape=16, size = size, stroke = 0, alpha = alpha, aes_string(color = color_var))+
      scale_colour_gradientn(colours = viridis(50, direction=-1, begin=0.05, end=0.95))+
      scale_fill_gradientn(colours = viridis(50, direction=-1, begin=0.05, end=0.95))
  } else {
    distvsratio_protrusion <- distvsratio_protrusion +
      geom_point(size = size, stroke = 0, alpha = alpha, aes_string(shape=shape_var))
  }
  
  distvsratio_protrusion <- distvsratio_protrusion +
    theme(legend.position = c(0.7, 0.7)) +
    theme(
      axis.text.y = element_text(size = 8),
      text = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(c(2, 1), "mm"),
      axis.text = element_text(size = 8),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank()
    )
  print(distvsratio_protrusion)
  
  # Construct the facet_nested statement
  if (!is.null(group_vars)) {
    facet_nested_expr <- paste(".~", paste(group_vars, collapse = "+"))
    distvsratio_protrusion <- distvsratio_protrusion + facet_nested(as.formula(facet_nested_expr))
  }
  
  # Calculate R-squared values per facet
  x_var_str <- deparse(substitute(x_var))
  y_var_str <- deparse(substitute(y_var))
  
  lm_formula <- switch(
    regression_type,
    "linear" = as.formula(paste(y_var_str, "~", x_var_str)),
    "log" = as.formula(paste(y_var_str, "~", "log", "(", x_var_str, ")")),
    "exp" = as.formula(paste(y_var_str, "~", "exp(", x_var_str, ")"))
  )
  
  print(lm_formula)
  rsquared_p_values <- data %>%
    filter(!is.na(.data[[x_var_str]]) & !is.infinite(.data[[x_var_str]]) & !is.nan(.data[[x_var_str]]) &
             !is.na(.data[[y_var_str]]) & !is.infinite(.data[[y_var_str]]) & !is.nan(.data[[y_var_str]]) &
             if (regression_type == "exp") .data[[x_var_str]] < 700 else TRUE) %>%
    group_by_at(group_vars) %>%
    summarize(
      rsquared = summary(lm(lm_formula, data = .))$r.squared,
      adj_r_squared = 1 - (1 - rsquared) * (length(unique(.data[[x_var_str]])) - 1) / (length(.data[[x_var_str]]) - length(coef(lm(lm_formula, data = .)))),
      p_values = coef(summary(lm(lm_formula, data = .)))[, "Pr(>|t|)"],
      values = coef(summary(lm(lm_formula, data = .)))[, "Estimate"],
      names = c("Intercept", "Coefficient")
    ) %>%
    ungroup()
  
  rsquared_p_values$regression_type <- rep(regression_type, nrow(rsquared_p_values))
  
  print(rsquared_p_values)
  
  # Adjust p-values using Bonferroni correction
  num_comparisons <- nrow(rsquared_p_values)  # Assuming each row corresponds to a comparison
  rsquared_p_values$p_values_adj <- p.adjust(rsquared_p_values$p_values, method = "bonferroni", n = num_comparisons)
  
  categorize_significance <- function(p_value_adj) {
    if (is.na(p_value_adj) || !is.numeric(p_value_adj)) {
      return("NA")
    } else if (p_value_adj < 0.001) {
      return("***")
    } else if (p_value_adj < 0.01) {
      return("**")
    } else if (p_value_adj < 0.05) {
      return("*")
    } else {
      return("N.S")
    }
  }
  
  # Apply the function to the rsquared_p_values dataframe
  rsquared_p_values <- rsquared_p_values %>%
    mutate(
      significance_intercept = categorize_significance(p_values_adj[1]),
      significance_coefficient = categorize_significance(p_values_adj[2])
    )
  
  print(rsquared_p_values)
  
  # Prepare a dataframe for annotations
  annotations_df <- rsquared_p_values %>%
    mutate(annotation_text = paste(
      sprintf("R^2= %.2f", adj_r_squared[2]),
      #sprintf("p.adj.coef= %s", significance_coefficient[2]),
      sprintf("P.adj= %.2e", p_values_adj[2]),
      sep = "\n"),
      x_coord = Inf ,
      y_coord = Inf
    )
  
  print(annotations_df)  
  
  # Create a dataframe from the rsquared_p_values
  rsquared_values_df <- as.data.frame(rsquared_p_values)
  
  basic_formula <- switch(
    regression_type,
    "linear" = as.formula(paste("y ~ x")),
    "log" = as.formula(paste("y ~ log(x)")),
    "exp" = as.formula(paste("y ~ exp(x)"))
  )
  print(basic_formula)
  
  # Add linear regression line if significant
  if (regression_type %in% c("linear", "log", "exp")) {
    if (regression_type == "exp") {
      final_plot <- distvsratio_protrusion +
        geom_smooth(data = filter(data, .data[[x_var_str]] < 700),
                    method = "lm", formula = basic_formula, fullrange = FALSE, se = TRUE, color = "black", linewidth = 0.5, span = 0.1)
    } else {
      final_plot <- distvsratio_protrusion +
        geom_smooth(data = data,
                    method = "lm", formula = basic_formula, fullrange = FALSE, se = TRUE, color = "black", linewidth = 0.5, span = 0.1)
    }
  }
  
  if(annotate_signif==TRUE){
    # Add annotations
    final_plot <- final_plot +
      geom_text(data = annotations_df, aes(x = x_coord, y = y_coord, label = annotation_text),
                hjust = 1.2, vjust = 1.2, size = 3)
  }
  # Return both the plot and the R-squared values dataframe
  return(list(plot = final_plot, rsquared_values = rsquared_values_df))
}


perform_tukey_hsd <- function(data, group_vars, x_var, y_var) {
  # Ensure necessary libraries are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
  
  library(dplyr)
  library(rstatix)
  
  # Convert group_vars to symbols
  group_syms <- rlang::syms(group_vars)
  
  # Create a formula for the analysis
  formula <- as.formula(paste0(y_var, " ~ ", x_var))
  
  # Perform the analysis
  result <- data %>%
    group_by(!!!group_syms) %>%
    rstatix::tukey_hsd(formula) %>%
    add_significance() %>%
    add_xy_position()
  
  return(result)
}

custom_violinplot_with_anova_tukey <- function(data, group_vars, x_var, y_var, log=FALSE, jitter=FALSE, 
                                               y_adjust_min=NULL, y_adjust_max=NULL) {
  
  
  # Ensure necessary libraries are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
  
  library(dplyr)
  library(rstatix)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  library(ggpubr)
  library(scales)
  library(viridis)
  
  
  # Construct the facet_nested statement
  facet_nested_expr <- paste(".~", paste(group_vars, collapse = "+"))
  
  print(facet_nested_expr)
  
  violin <- ggplot(data, aes_string(x=x_var, y=y_var))
  
  if(jitter==TRUE){
    violin <- violin + geom_jitter(size=1, stroke=0)
  }
  
  violin <- violin +
    geom_violin(scale = "width", aes_string(fill=x_var)) +
    geom_boxplot(outlier.shape=NA, width = 0.2, position=position_dodge(.9), aes_string(fill=x_var))
  
  pretty_breaks <- pretty_breaks(n=3)
  violin <- violin + scale_y_continuous(breaks=pretty_breaks, expand = expansion(mult = c(0.05, 0.1))) 
  
  # Convert strings to symbols
  group_syms <- rlang::syms(group_vars)
  x_sym <- rlang::sym(x_var)
  y_sym <- rlang::sym(y_var)
  
  
  # Create formula
  formula <- as.formula(paste0(y_var, " ~ ", x_var))
  
  print(paste("tukey performed with grouping:", group_syms))
  
  print(paste("tukey performed with formula:", y_var," ~", x_var))
  
  
  # Perform the analysis
  result <- data %>%
    group_by(!!!group_syms) %>%
    rstatix::tukey_hsd(formula) %>%
    add_significance() %>%
    add_xy_position()
  
  print("tukey worked")
  print(result)  
  
  result <- result %>%
    mutate(y.position = 1.1*(y.position))
  
  
  if (log) {
    result <- result %>%
      mutate(y.position = 1.2*(log10(y.position)))
    violin <- violin+ scale_y_continuous(trans="log10", n.breaks=4, 
                                         expand = expansion(mult = c(0.05, 0.1)))
  }

  
  if (!is.null(y_adjust_min)& !is.null(y_adjust_max)) {
    violin <- violin + scale_y_continuous(limits = c(y_adjust_min, y_adjust_max*1.2))
    
    result <- result %>%
      mutate(y.position = y.position*(y_adjust_max/y.position*1.1))
  }
  
  violin <- violin+
    stat_pvalue_manual(
      result, label = "p.adj.signif"
    ) +
    scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
    scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)
  
  
  
  if(!is.null(group_vars)){
    violin <- violin + 
      facet_nested(as.formula(facet_nested_expr))
  }
  
  
  
  
  print(violin)
  
  
  # Create a dataframe the tukey values
  result <- as.data.frame(result)
  
  print(result)
  
  # Return both the plot and the tukey values
  return(list(plot = violin, tukey_values = result))
  
  
}




custom_violinplot_with_statistics <- function(data, group_vars, x_var, y_var, log=FALSE, jitter=FALSE, fill_var=NULL,
                                              y_adjust_min=NULL, y_adjust_max=NULL, annotate_signif=TRUE, annotate_only_signif=TRUE, plot_mean=FALSE) {
  
  # Ensure necessary libraries are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
  
  library(dplyr)
  library(rstatix)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  library(ggpubr)
  library(scales)
  library(viridis)
  
  # Construct the facet_nested statement
  facet_nested_expr <- paste(".~", paste(group_vars, collapse = "+"))
  
  print(facet_nested_expr)
  
  violin <- ggplot(data, aes_string(x=x_var, y=y_var))
  
  if(jitter==TRUE){
    violin <- violin + geom_jitter(size=1, stroke=0)
  }
  
  if (plot_mean) {
    violin <- violin +
      geom_violin(scale = "width", aes_string(fill=fill_var)) +
      stat_summary(fun=mean, geom="crossbar", aes_string(group=fill_var), size=0.5, width=0.5, position=position_dodge(width=0.9))
  } else {
    violin <- violin +
      geom_violin(scale = "width", aes_string(fill=fill_var)) +
      geom_boxplot(outlier.shape=NA, width = 0.2, position=position_dodge(.9), aes_string(fill=fill_var))
  }
  
  pretty_breaks <- pretty_breaks(n=3)
  violin <- violin + scale_y_continuous(breaks=pretty_breaks, expand = expansion(mult = c(0.05, 0.1))) 
  
  # Convert strings to symbols
  group_syms <- rlang::syms(group_vars)
  x_sym <- rlang::sym(x_var)
  y_sym <- rlang::sym(y_var)
  
  # Create formula
  formula <- as.formula(paste0(y_var, " ~ ", x_var))
  
  print(paste("Grouping:", group_syms))
  
  print(paste("Statistical test performed with formula:", y_var," ~", x_var))
  
  
  # Check the number of levels in x_var
  num_levels <- length(unique(data[[x_var]]))
  
  if(is.null(group_vars)){
    
    if(num_levels > 2){
      # Perform ANOVA
      result <- data %>%
        group_by(!!!group_syms) %>%
        rstatix::tukey_hsd(formula) %>%
        add_significance() %>%
        add_xy_position() %>%
        mutate(test_type = "ANOVA-Tukey")
      
      print("More than two levels in x_var - ANOVA performed")
    } else {
      # Perform T-test
      result <- data %>%
        group_by(!!!group_syms) %>% 
        rstatix::t_test(formula) %>%
        add_significance() %>%
        add_xy_position() %>%
        mutate(test_type = "T-test")
      
      print("Only two levels in x_var - T test performed")
    }
  }
  
  if(!is.null(group_vars)){
    # Perform the analysis
    result <- data %>%
      group_by(!!!group_syms) %>%
      rstatix::tukey_hsd(formula) %>%
      add_significance() %>%
      add_xy_position() %>%
      mutate(test_type = "ANOVA-Tukey")
    
    print("More than two groups compared - ANOVA with tukey post-hoc test performed")
  }
  
  print(result)
  
  if (log) {
    result <- result %>%
      mutate(y.position = 1.2*(log10(y.position)))
    violin <- violin+ scale_y_continuous(trans="log10", n.breaks=4, 
                                         expand = expansion(mult = c(0.05, 0.1)))
  }
  
  if (!is.null(y_adjust_min)& !is.null(y_adjust_max)) {
    violin <- violin + scale_y_continuous(limits = c(y_adjust_min, y_adjust_max*1.1))
    
    result <- result %>%
      mutate(y.position = y.position*(y_adjust_max/y.position*1.1))
  }
  
  if(annotate_signif==TRUE){
    
    if(!is.null(group_vars) | (num_levels > 2)){
      if(annotate_only_signif==TRUE){
        result <- result %>%  filter(p.adj.signif != "ns")
      }
      violin <- violin +
        stat_pvalue_manual(
          result, label = "p.adj.signif"
        )
    } else {
      if(annotate_only_signif==TRUE){
        result <- result %>%  filter(p.signif != "ns")
      }
      violin <- violin +
        stat_pvalue_manual(
          result, label = "p.signif"
        )
    }
  }
  
  if(!is.null(group_vars)){
    violin <- violin + 
      facet_nested(as.formula(facet_nested_expr))
  }
  
  print(violin)
  
  # Create a dataframe the tukey values
  result <- as.data.frame(result)
  
  print(result)
  
  # Return both the plot and the tukey values
  return(list(plot = violin, statistical_test = result))
  
}




custom_plot_with_statistics <- function(data, group_vars, x_var, y_var, log=FALSE, jitter=FALSE, fill_var=NULL, size=2, alpha=1,
                                              y_adjust_min=NULL, y_adjust_max=NULL, annotate_signif=TRUE, annotate_only_signif=TRUE, plot_mean=FALSE, violinplot=FALSE) {
  
  # Ensure necessary libraries are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
  
  library(dplyr)
  library(rstatix)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  library(ggpubr)
  library(scales)
  library(viridis)
  
  # Construct the facet_nested statement
  facet_nested_expr <- paste(".~", paste(group_vars, collapse = "+"))
  
  print(facet_nested_expr)
  
  violin <- ggplot(data, aes_string(x=x_var, y=y_var))
  

  if(violinplot==TRUE){
    violin <- violin +
      geom_violin(scale = "width", aes_string(fill=fill_var))
      
  }
  
  if(jitter==TRUE){
    violin <- violin + geom_jitter(size=1, stroke=0)
  }
  
  
  if (plot_mean) {
    violin <- violin +
      stat_summary(fun=mean, geom="crossbar", aes_string(group=fill_var), size=0.5, width=0.5, position=position_dodge(width=0.9))
  } else {
    violin <- violin +
      geom_boxplot(outlier.size=size, outlier.stroke=0, outlier.alpha=alpha, outlier.shape = NA, width = 0.5, position=position_dodge(.9), aes_string(fill=fill_var))
  }
  
  pretty_breaks <- pretty_breaks(n=3)
  violin <- violin + scale_y_continuous(breaks=pretty_breaks, expand = expansion(mult = c(0.05, 0.1))) 
  
  # Convert strings to symbols
  group_syms <- rlang::syms(group_vars)
  x_sym <- rlang::sym(x_var)
  y_sym <- rlang::sym(y_var)
  

  # Create formula
  formula <- as.formula(paste0(y_var, " ~ ", x_var))
  
  print(paste("Grouping:", group_syms))
  
  print(paste("Statistical test performed with formula:", y_var," ~", x_var))
  
  # Check the number of levels in x_var
  num_levels <- length(unique(data[[x_var]]))
  print(paste("Number of levels in", x_var, ":", num_levels))
  
  # Check the number of levels in group_vars
  if(!is.null(group_vars)){
  group_levels <- length(unique(data[[group_vars]]))
  print(paste("Number of levels in", group_vars, ":", group_levels))
  print(levels(group_vars))
  }
  
  
  if(is.null(group_vars) | !is.null(group_vars) & group_levels==1){
    
    if(num_levels > 2){
      print("More than two levels in x_var - ANOVA with tukey post-hoc performed")
      
      # Perform ANOVA
      result <- data %>%
        group_by(!!!group_syms) %>%
        rstatix::tukey_hsd(formula) %>%
        add_significance() %>%
        add_xy_position() %>%
        mutate(test_type = "ANOVA-Tukey"
        )
      
    } else {
      # Perform T-test
      print("Only two levels in x_var - T test performed")
      
      result <- data %>%
        group_by(!!!group_syms) %>% 
        rstatix::t_test(formula) %>%
        add_significance() %>%
        add_xy_position() %>%
        mutate(test_type = "T-test",
               p.adj.signif = p.signif
        )
    }
  }
  
  if(!is.null(group_vars) & group_levels>1){
    # Perform the analysis
    print("More than two groups compared - ANOVA with tukey post-hoc test performed")
    
    result <- data %>%
      group_by(!!!group_syms) %>%
      rstatix::tukey_hsd(formula) %>%
      add_significance() %>%
      add_xy_position() %>%
      mutate(test_type = "ANOVA-Tukey"
      )
    
  }
  
  print(result)
  
  
  
  
  if (log) {
    result <- result %>%
      mutate(y.position = 1.2*(log10(y.position)))
    violin <- violin+ scale_y_continuous(trans="log10", n.breaks=4, 
                                         expand = expansion(mult = c(0.05, 0.1)))
  }

  
  if (!is.null(y_adjust_min)& !is.null(y_adjust_max)) {
    violin <- violin + scale_y_continuous(limits = c(y_adjust_min, y_adjust_max*1.1))
    
    result <- result %>%
      mutate(y.position = y.position*(y_adjust_max/y.position*1.1))
  }
  
  result <- result %>%
    mutate(y.position = y.position*1.1)
  
  if(annotate_signif==TRUE){
    
    if(!is.null(group_vars) | (num_levels > 2)){
      if(annotate_only_signif==TRUE){
        result <- result %>%  filter(p.adj.signif != "ns")
      }
      violin <- violin +
        stat_pvalue_manual(
          result, label = "p.adj.signif"
        )
    } else {
      if(annotate_only_signif==TRUE){
        result <- result %>%  filter(p.signif != "ns")
      }
      violin <- violin +
        stat_pvalue_manual(
          result, label = "p.signif",
          step.increase = 0.5
        )
    }
  }
  
  if(!is.null(group_vars)){
    violin <- violin + 
      facet_nested(as.formula(facet_nested_expr))
  }
  
  print(violin)
  
  # Create a dataframe the tukey values
  result <- as.data.frame(result)
  
  print(result)
  
  # Return both the plot and the tukey values
  return(list(plot = violin, statistical_test = result))
  
}






create_histogram <- function(data, x_var, facet_var, num_facets = NULL, num_bins = 30, factor = 1) {
  if (!is.null(num_facets) && num_facets < length(unique(data[[facet_var]]))) {
    set.seed(123)
    selected_facets <- sample(unique(data[[facet_var]]), num_facets)
    data <- data[data[[facet_var]] %in% selected_facets, ]
  }
  
  stats <- data %>%
    group_by(!!rlang::sym(facet_var)) %>%
    summarise(mean = mean(!!rlang::sym(x_var)),
              sd = sd(!!rlang::sym(x_var)))
  
  plot_data <- stats %>%
    mutate(upper_threshold = mean + factor * sd,
           lower_threshold = mean - factor * sd) %>%
    dplyr::select(-sd) %>%
    tidyr::pivot_longer(cols = c("mean", "upper_threshold", "lower_threshold"), names_to = "statistic", values_to = "value")
  
  histograms <- ggplot2::ggplot(data, ggplot2::aes(x=!!rlang::sym(x_var))) +
    ggplot2::geom_histogram(ggplot2::aes(y=after_stat(ndensity)), bins=num_bins, position="identity", fill="grey") +
    ggplot2::geom_density(ggplot2::aes(y=after_stat(scaled)), adjust=1, color="black") +
    ggplot2::facet_wrap(as.formula(paste0('.~', facet_var)), scales = "fixed", ncol = 4) +  # Set ncol to 4
    ggplot2::geom_vline(data=plot_data, ggplot2::aes(xintercept=value, color=statistic), linetype="dashed") +
    ggplot2::scale_color_manual(values=c("mean"="blue", "upper_threshold"="red", "lower_threshold"="red")) +
    scale_linetype_manual(values=c("dashed", "solid", "dashed"))+
    ylim(0,1) +
    ggplot2::theme_classic()
  
  print(histograms)
}



## import data------
workdir <- choose.dir(getwd(), "Choose folder containing subdirectories with .csv files")
parentdirname <- basename(workdir)

## if available, import existing file--
existingfilename <- basename(list.files(workdir, pattern="Data_filtered.csv"))
existingfilename_alldata <- basename(list.files(workdir, pattern="Unprocessed_data.csv"))

if(length(existingfilename==1)){
  existingfile <- read.csv2(paste0(workdir,"/",existingfilename), header=T, sep = ",")
  analyzedrois <- unique(existingfile$name)
  existingfile <- existingfile[2:ncol(existingfile)]
  
  
  #existingfile_alldata <- read.csv2(paste0(workdir,"/",existingfilename_alldata), header=T, sep = ",")
  #existingfile_alldata <- existingfile_alldata[2:ncol(existingfile_alldata)]
  
  analyzedprofiles <- length(unique(existingfile$number))
}else{
  analyzedprofiles = 0
}


protrusiontypes <- list.dirs(path=workdir, recursive=FALSE)

for(p in 1:length(protrusiontypes)){
  
  protrusiontype <- protrusiontypes[p]
  protrusionbasename <- basename(protrusiontype)
  
  dirlist <- list.dirs(path=protrusiontype, recursive=FALSE)
  ndirs <- length(dirlist)
  
  for(d in 1:ndirs){
    print(d)
    dir <- dirlist[d]
    setwd(dir)
    Filebasename <- basename(list.files(dir, pattern=".csv"))
    for(n in 1:length(Filebasename)){
      filename <- Filebasename[n]
      filepath <- paste0(dir,"/", filename)
      filepath <- normalizePath(filepath)
      file <- read.csv2(filepath, header=T, sep = ",", quote="")
      filename <- substring(filename, 1, (nchar(filename)-4))
      file$name <- filename
      

      if(n==1){
        
      dirfiles <- file
      
      }else{
        dirfiles2 <- file
        
        dirfiles <-rbind(dirfiles,dirfiles2)
      }
    }
    setwd(workdir)
    if(d==1){
    file1 <- dirfiles
   
    
    }else{
    file2 <- dirfiles
    
    file1 <- rbind(file1,file2)
    }
  }
  
  if(p==1){
  print(protrusionbasename)
  data <- cbind(file1, rep(protrusionbasename, nrow(file1)))
  colnames(data)[ncol(data)] <- "Protrusiontype"
  }else{
    
    print(protrusionbasename)
    data2 <- cbind(file1, rep(protrusionbasename, nrow(file1)))
    colnames(data2)[ncol(data2)] <- "Protrusiontype"
    
    data <- rbind(data,data2)
  }
}

file1 <- data
file1 <- file1 %>%dplyr::select(X, Value, Protrusiontype, name)

## define theme for plotting------

# Define the custom theme with the "sans" font, black text color, and smaller legend
custom_theme <- theme_classic()+
  theme(
  text = element_text(family = "sans", size = 7, color = "black"),
  legend.key.size = unit(1, "lines"),  # Adjust the size of the legend key
  legend.key.height = unit(0.5, "lines"),  # Adjust the height of the legend key
  # Customize the appearance of geom_text
  axis.text = element_text(size = 6),  # Customize the axis text size
  axis.title = element_text(size = 6),  # Customize the axis title text size
)


## make sense of column names & organize dataframe-----
rm(file)
rm(file2)
file <- separate(file1, (ncol(file1)), into=c("Identifier", "CellNR", "Orientation", 
                                              "LocationNR", "Imaging type", "RoiNR", "Fiber_type", "Linearity_index",
                                              "Channel", "Signal_type", "Branching_type", "Position_relative_to_fiber"), sep="_", remove = FALSE)
file$Value <- as.numeric(file$Value)
rm(file1)
rm(dirfiles)
rm(dirfiles2)

namelist <- unique(file$name)
uniquenames <- namelist[grep("col", namelist)]
uniques <- gsub("_col_parallel", "", uniquenames)


for(g in 1 :length(uniques)){
    sub <- file[grepl(uniques[g], file$name),]
    conditionname <- uniques[g]
    for(u in 1:length(unique(sub$Channel))){
      c <- unique(sub$Channel)
      c <- c[u]
      pp <- subset(sub, Channel==c)
      X <- pp$X
      pp_t <- pp$Value
      pp_t <- as.data.frame(pp_t)
      pp_t$Channel <- rep(c, nrow(pp_t))
      pp_t$RoiNR <- rep(unique(sub$RoiNR), nrow(pp_t))
      pp_t$X <- X
      pp_t$Value <- pp$Value
      pp_t$Type <- pp$Type
      colnames(pp_t)[1] <- "Value_z"
      df <- pp_t
      if(u==1){
        df1 <- df
      }else{
        df2 <- df
        df1 <- rbind(df1, df2)
      }
    }
    df1$name <- conditionname
    if(g==1){
      plotprofile <- df1
    }else{
      plotprofile2 <- df1
      plotprofile <- rbind(plotprofile, plotprofile2)
    }
  }



if(exists("B1int_GCX_ratio")==TRUE){
  analyzedrois <- length(unique(B1int_GCX_ratio$RoiNR))
}else{
  analyzedrois <- 0
}

file$Linearity_index <- as.numeric(file$Linearity_index)

file <- file %>%
  mutate(number = group_indices(., Orientation, CellNR, LocationNR, RoiNR))

file <- file %>%
  mutate(grouped_number = group_indices(., CellNR, LocationNR, RoiNR))


file <- file %>% 
  group_by(number, Channel, RoiNR, CellNR, Orientation) %>% 
  mutate(Value_bgcor = ifelse(Channel!="col", 
                              Value- mean(Value[Signal_type=="background"]), 
                                          Value))

bgcor = FALSE

if(bgcor==TRUE){
  file$Value <- file$Value_bgcor
}

file$Position_relative_to_fiber[is.na(file$Position_relative_to_fiber)] <- "aftercellstart"

## manually add analysis points and in the end generate dataframe with ratio B1int/gcx -----

 # B1int_GCX_ratio <- existingfile
  
file <- file %>%  group_by(Orientation, CellNR, RoiNR, Signal_type) %>%  mutate(X= ifelse(Position_relative_to_fiber=="beforecellstart",
                                                                                X-max(X[Position_relative_to_fiber=="beforecellstart"]),
                                                                                X))
file$X_micron <- file$X*pxsize
## plot a profile to determine that before/after cell-matrix works properly------

file <- file %>%  group_by(number, grouped_number, Channel, RoiNR, CellNR) %>%  
  mutate(Value_normalized = ifelse(Signal_type=="foreground",
                                   1/Value ,
                                   Value))

plot <- ggplot(subset(file, number==18), aes(x=X_micron, y=Value))+
  geom_line(aes(color=Channel, linetype=Signal_type))+
  facet_wrap(.~Channel)+theme_classic()

print(plot)


## manually select regions of focalized and non-focalized integrin along the plot profile------
  

  
  skip= "yes"
  
  if(skip!="yes"){
    
    
    
  for(r in (analyzedprofiles+1):length(unique(file$number))){
    print(r)
    
    plotprofile_sub <- subset(file, number==r)
    
    ppname <- unique(plotprofile_sub$name)
    ppname <- ppname[2]
    
    
    mx <- subset(plotprofile_sub,Channel=="B1int")
    mx <- max(mx$Value)
    mn <- subset(plotprofile_sub,Channel=="gcx")
    mn <- max(mn$Value)
    scale <- mx/mn
    
    p <- ggplot(plotprofile_sub, aes(x=X, y=Value_bgcor, color=Channel))+
      geom_line(data= plotprofile_sub %>% filter(Channel=="B1int" & Signal_type!="background"),aes(y=Value_bgcor), linewidth=1)+
    #  geom_line(data= plotprofile_sub %>% filter(Channel=="gcx"   & Signal_type!="background"),aes(y=Value_bgcor*scale),linewidth=1)+
      scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
      scale_x_continuous(breaks = seq(0, max(plotprofile_sub$X), by = 50),
                         minor_breaks = seq(0, max(plotprofile_sub$X), by = 10))+    
      theme(panel.grid.major = element_line(color = "white",
                                            linewidth = 1,
                                            linetype = 1))+
      
      ggtitle(paste(ppname, "\n", unique(plotprofile_sub$Orientation)))+
      
      scale_y_continuous(
        # Features of the first axis
        name = "B1int signal",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*scale, name="Gcx signal")
      )
    
    
    
    print(p)
    
    stop_loop <- "no"
    
    if(r==1){
      stop_loop <- dlg_message("Do you want to stop looping through all plot?", type = "yesno")$res
      if (stop_loop == "yes") {
        break  # Exit the loop
      }
    }
    x=1
    
    
    repeat{
      set <- subset(plotprofile_sub, Channel!="col")
      
      if(unique(set$Protrusiontype=="Bleb" | set$Protrusiontype=="Front")){
        print("Bleb/ retraction fiber detected -> length-based segregation analysis \n
              Leading edge protrusions -> focalization analysis")
        quant <- set
        quant$Type <- rep(NA, nrow(quant))
      }else{
        print("Leading edge protrusion detected -> focal adhesion-based segregation analysis")
        Xthr <- dlg_input("X value?")$res
        Typeq <- dlg_list(c("Focalized B1 integrin", "Control", "Cell start"))$res
        quant <- subset(set, X==Xthr)
        quant$Type <- rep(Typeq, nrow(quant))
      }
      
      
      
      Qualq <- dlg_message("Good?",type="yesno")$res
      if(Qualq=="yes"){
        quant$Quality <- rep("Good", nrow(quant))
      }else{
        quant$Quality <- rep("Bad", nrow(quant))
      }
      
      quant$number <- rep(unique(plotprofile_sub$number), nrow(quant))
      
      
      if(x==1){
        quant1 <- quant
      }else{
        quant2 <- quant
        quant1 <- rbind(quant1,quant2)
      }
      
      
      x=x+1
      
      again <- dlg_message("Analyze another point in graph?",type="yesno")$res
      if(again=="no"){ 
        break}
     
    }
    
    if(r==analyzedprofiles+1){
      stats <- quant1  
    }else{
      stats2 <- quant1
      stats <- rbind(stats,stats2)
    }
    
  }
    # Subset 'file' to get the non-common data
    non_common_data <- file %>%
      filter((Orientation %in% non_duplicates$Orientation & 
                CellNR %in% non_duplicates$CellNR & 
                RoiNR %in% non_duplicates$RoiNR &
                Signal_type %in% non_duplicates$Signal_type))
  stats <- subset(stats, Quality=="Good")
  
  }else{
    stats <- file
    stats$Type <- rep(NA, nrow(stats))
  }
  
  
## calculate ratio B1int/gcx for selected regions------


B1int_GCX_ratio <- stats %>% group_by(X, number, grouped_number, Protrusiontype, Identifier, CellNR, Orientation,
                                      LocationNR, "Imaging type", RoiNR, Fiber_type, Linearity_index, Signal_type,
                                      Branching_type, Position_relative_to_fiber, Type) %>%  
    summarise(ratio = Value[Channel=="B1int"]/Value[Channel=="gcx"])

ratio.bgcor <- stats %>% group_by(X, number, Protrusiontype, Identifier, CellNR, Orientation,
                                  LocationNR, "Imaging type", RoiNR, Fiber_type, Linearity_index, Signal_type,
                                  Branching_type, Position_relative_to_fiber, Type) %>%  
    summarise(ratio_bgcor = Value_bgcor[Channel=="B1int"]/Value_bgcor[Channel=="gcx"])

B1int_GCX_ratio$ratio_bgcor <- ratio.bgcor$ratio_bgcor
  
  # retrieve B1int data
  filtered_file <- stats[stats$Channel == "B1int", ]
  filtered_file <- filtered_file[c("number", "X", "Signal_type", "Value", "Position_relative_to_fiber")]
  colnames(filtered_file)[colnames(filtered_file) == "Value"] <- "B1int"
  
  
  # retrieve gcx data
  B1int_GCX_ratio <- merge(B1int_GCX_ratio, filtered_file, by = c("X", "number", "Signal_type", "Position_relative_to_fiber"), all.x = TRUE)
  
  filtered_file <- stats[stats$Channel == "gcx", ]
  filtered_file <- filtered_file[c("number", "X", "Signal_type", "Value", "Position_relative_to_fiber")]
  colnames(filtered_file)[colnames(filtered_file) == "Value"] <- "Gcx"
  
  # retrieve col data
  B1int_GCX_ratio <- merge(B1int_GCX_ratio, filtered_file, by = c("X", "number", "Signal_type", "Position_relative_to_fiber"), all.x = TRUE)
  
  
  filtered_file <- stats[stats$Channel == "col", ]
  filtered_file <- filtered_file[c("number", "X", "Signal_type", "Value", "Position_relative_to_fiber")]
  colnames(filtered_file)[colnames(filtered_file) == "Value"] <- "Col"
  
  ## merge B1int and Gcx data
  B1int_GCX_ratio <- merge(B1int_GCX_ratio, filtered_file, by = c("X", "number", "Signal_type", "Position_relative_to_fiber"), all.x = TRUE)
  
  

B1int_GCX_ratio$Type <- gsub("Control", "Non-focalized B1 integrin", B1int_GCX_ratio$Type)

B1int_GCX_ratio$Type[B1int_GCX_ratio$Type=="Focalized B1 integrin"] <- "Focalized"
B1int_GCX_ratio$Type[B1int_GCX_ratio$Type=="Non-focalized B1 integrin"] <- "Non-focalized"


B1int_GCX_ratio$X_fromcellstart[B1int_GCX_ratio$Orientation=="Bleb"] <- B1int_GCX_ratio$X[B1int_GCX_ratio$Orientation=="Bleb"]

B1int_GCX_ratio$X_fromcellstart_micron <- B1int_GCX_ratio$X_fromcellstart*pxsize
B1int_GCX_ratio$X_micron <- (B1int_GCX_ratio$X*pxsize)-pxsize

B1int_GCX_ratio$ratio[is.infinite(B1int_GCX_ratio$ratio)] <- 0

B1int_GCX_ratio$Orientation <- gsub("Front", "Leading edge protrusion", B1int_GCX_ratio$Orientation)
B1int_GCX_ratio$Orientation <- gsub("Rear", "Retraction fiber", B1int_GCX_ratio$Orientation)
B1int_GCX_ratio$Fiber_type[B1int_GCX_ratio$Linearity_index>0.95 & B1int_GCX_ratio$Orientation=="Leading edge protrusion"] <- "Tense"
B1int_GCX_ratio$Fiber_type[B1int_GCX_ratio$Linearity_index>0 & B1int_GCX_ratio$Linearity_index<=0.95& B1int_GCX_ratio$Orientation=="Leading edge protrusion"] <- "Relaxed"
B1int_GCX_ratio$Fiber_type[B1int_GCX_ratio$Linearity_index==0] <- "No fiber contact"
B1int_GCX_ratio$Protrusiontype[B1int_GCX_ratio$Fiber_type == "No fiber contact"] <- "No fiber contact"
B1int_GCX_ratio$Protrusiontype[B1int_GCX_ratio$Fiber_type != "No fiber contact"] <- "Fiber-contacting"
B1int_GCX_ratio$Linearity_index[is.na(B1int_GCX_ratio$Linearity_index)] <- 0
B1int_GCX_ratio$Type[is.na(B1int_GCX_ratio$Type)] <- "Unspecified"


B1int_GCX_ratio$ratio <- as.numeric(B1int_GCX_ratio$ratio)
B1int_GCX_ratio$B1int <- as.numeric(B1int_GCX_ratio$B1int)
B1int_GCX_ratio$X <- as.numeric(B1int_GCX_ratio$X)
B1int_GCX_ratio$Gcx <- as.numeric(B1int_GCX_ratio$Gcx)
B1int_GCX_ratio$X_fromcellstart <- as.numeric(B1int_GCX_ratio$X_fromcellstart)
B1int_GCX_ratio$Linearity_index <- as.numeric(B1int_GCX_ratio$Linearity_index)

B1int_GCX_ratio$X_fromcellstart_micron <- as.numeric(B1int_GCX_ratio$X_fromcellstart_micron)
B1int_GCX_ratio$X_micron <- as.numeric(B1int_GCX_ratio$X_micron)
B1int_GCX_ratio$number <- as.character(B1int_GCX_ratio$number)





B1int_GCX_ratio <- B1int_GCX_ratio %>% 
  mutate(number = group_indices(., Orientation, LocationNR, CellNR, RoiNR))



B1int_GCX_ratio$ratio <- as.numeric(B1int_GCX_ratio$ratio)
B1int_GCX_ratio$B1int <- as.numeric(B1int_GCX_ratio$B1int)
B1int_GCX_ratio$X <- as.numeric(B1int_GCX_ratio$X)
B1int_GCX_ratio$Gcx <- as.numeric(B1int_GCX_ratio$Gcx)
B1int_GCX_ratio$X_fromcellstart <- as.numeric(B1int_GCX_ratio$X_fromcellstart)
B1int_GCX_ratio$Linearity_index <- as.numeric(B1int_GCX_ratio$Linearity_index)

B1int_GCX_ratio$X_fromcellstart_micron <- as.numeric(B1int_GCX_ratio$X_fromcellstart_micron)
B1int_GCX_ratio$X_micron <- as.numeric(B1int_GCX_ratio$X_micron)

B1int_GCX_ratio$number <- as.character(B1int_GCX_ratio$number)

## data normalization------
  
  B1int_GCX_ratio <- protrusion_length_segregation(B1int_GCX_ratio, protrusion_length_segregation_threshold)

  ## for blebs, segregate Type as "Bleb"<1/3 length and "Cell body"<2/3 length
  
  #B1int_GCX_ratio <- bleb_length_segregation(B1int_GCX_ratio, bleb_length_segregation_threshold)

  ## for retraction fibers, segregate Type as "Retraction fiber"<1/5 length and "Cell body"<4/5 length
  
  #B1int_GCX_ratio <- retraction_fiber_length_segregation(B1int_GCX_ratio, retraction_fiber_length_segregation_threshold)
  
  ## convert necessary columns to numeric
  B1int_GCX_ratio$ratio <- as.numeric(B1int_GCX_ratio$ratio)
  B1int_GCX_ratio$B1int <- as.numeric(B1int_GCX_ratio$B1int)
  B1int_GCX_ratio$X <- as.numeric(B1int_GCX_ratio$X)
  B1int_GCX_ratio$Gcx <- as.numeric(B1int_GCX_ratio$Gcx)
  B1int_GCX_ratio$X_fromcellstart <- as.numeric(B1int_GCX_ratio$X_fromcellstart)
  B1int_GCX_ratio$Linearity_index <- as.numeric(B1int_GCX_ratio$Linearity_index)
  B1int_GCX_ratio$X_fromcellstart_micron <- as.numeric(B1int_GCX_ratio$X_fromcellstart_micron)
  B1int_GCX_ratio$X_micron <- as.numeric(B1int_GCX_ratio$X_micron)
  B1int_GCX_ratio$number <- as.character(B1int_GCX_ratio$number)
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%
    arrange(number, X) %>%
    group_by(number) %>%
    mutate(X_relative = X/max(X))
  
  
  ## ratio normalization for blebs 
  
  B1int_GCX_ratio <- ratio_normalization_blebs(B1int_GCX_ratio)
  
  ## calculate B1int positivity for blebs
  
    
  B1int_GCX_ratio <- B1int_GCX_ratio %>%
      group_by(number, Orientation, RoiNR, Signal_type, CellNR) %>% 
      mutate(
        b1int_bleb_mean = mean(B1int[Type == "Bleb"], na.rm=TRUE),
        b1int_cellbody_mean = mean(B1int[Type == "Cell body"], na.rm=TRUE),
        b1int_cellbody_sd = sd(B1int[Type == "Cell body"], na.rm=TRUE),
        thr_b1int_negative = b1int_cellbody_mean - b1int_cellbody_sd,
        thr_b1int_positive = b1int_cellbody_mean + b1int_cellbody_sd,
        B1int_status = ifelse(
          b1int_bleb_mean > thr_b1int_positive, "Positive",
          ifelse(b1int_bleb_mean < thr_b1int_negative, "Negative", "Unchanged")
        )
      )
  
  bleb_B1int_stats <- B1int_GCX_ratio %>% 
    distinct(number, CellNR, RoiNR, B1int_status, b1int_bleb_mean, b1int_cellbody_mean, b1int_cellbody_sd,
             thr_b1int_negative, thr_b1int_positive) %>%  na.omit() %>% filter(Orientation=="Bleb" & Signal_type=="foreground")
  
  bleb_B1int_count <- bleb_B1int_stats %>%
    group_by(B1int_status) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(percentage = (count / sum(count)) * 100)
  
  ## calculate GCX positivity for blebs
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>% 
    group_by(number, Orientation, RoiNR, Signal_type) %>%
    mutate(
      gcx_bleb_mean = mean(Gcx[Type == "Bleb"]),
      gcx_cellbody_mean = mean(Gcx[Type == "Cell body"]),
      gcx_cellbody_sd = sd(Gcx[Type == "Cell body"]),
      thr_gcx_negative = gcx_cellbody_mean - gcx_cellbody_sd,
      thr_gcx_positive = gcx_cellbody_mean + gcx_cellbody_sd,
      gcx_status = ifelse(
        gcx_bleb_mean > thr_gcx_positive, "Positive",
        ifelse(gcx_bleb_mean < thr_gcx_negative, "Negative", "Unchanged")
      )
    )
  
  summary_gcx_status <- B1int_GCX_ratio %>% 
    distinct(gcx_bleb_mean, gcx_cellbody_mean, gcx_cellbody_sd, 
                                          thr_gcx_negative, thr_gcx_positive, gcx_status) %>%  na.omit()%>% filter(Orientation=="Bleb" & Signal_type=="foreground")
    
    
  bleb_gcx_b1int_stats <- B1int_GCX_ratio %>%
    distinct(number, gcx_status, B1int_status, CellNR, Orientation) %>%  na.omit() %>% filter(Orientation=="Bleb" & Signal_type=="foreground")
  
  
  bleb_gcx_stats <- B1int_GCX_ratio %>% 
    distinct(number, gcx_status, CellNR) %>%  na.omit() %>% filter(Orientation=="Bleb" & Signal_type=="foreground") 
  
  bleb_gcx_b1int_count <- bleb_gcx_b1int_stats %>% 
    group_by(B1int_status) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(percentage = (count / sum(count)) * 100)
  
  bleb_gcx_b1int_count_no_unchanged <- bleb_gcx_b1int_stats %>% filter(B1int_status!="Unchanged") %>% 
    group_by(B1int_status, gcx_status) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(percentage = (count / sum(count)) * 100)
  
  
  write.csv(x=bleb_gcx_b1int_count, "Blebs_B1int_gcx_status.csv")
  write.csv(x=bleb_gcx_b1int_stats, "Blebs_B1int_status.csv")
  
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%  group_by(number, Signal_type) %>%  mutate(B1int_normalized = B1int / mean(B1int[Type=="Cell body"]))
  B1int_GCX_ratio <- B1int_GCX_ratio %>%  group_by(number, Signal_type) %>%  mutate(Gcx_normalized = Gcx / mean(Gcx[Type=="Cell body"]))
  
  
  X_thr_leadingedge = 20
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%group_by(number, Signal_type) %>% mutate(Type= ifelse(Orientation=="Leading edge protrusion" &
                                                                                                  X_micron>X_thr_leadingedge,
                                                                                   "Unspecified",
                                                                                   Type))
  
  ## calculate mean_ratio_leadingedgeprotrusion to get the mean gcx/B1int ratio on 1:
  
 # B1int_GCX_ratio <- define_focalized_b1int_cellbody_v2(B1int_GCX_ratio, Zdiff=0, smooth=FALSE)
  
  #B1int_GCX_ratio <- define_focalized_b1int_all_v2(B1int_GCX_ratio, Zdiff=0, smooth=FALSE)
  
  B1int_GCX_ratio <- define_focalized_b1int_leadingedge_v2(B1int_GCX_ratio, smooth=TRUE, smooth_span=0.30)
  
  #B1int_GCX_ratio <- define_focalized_b1int_retractionfiber(B1int_GCX_ratio, Xdiff=7)
  
  B1int_GCX_ratio <- ratio_normalization_leadingedgeprotrusion(B1int_GCX_ratio, Xthr=Xthr)
  
  B1int_GCX_ratio <- ratio_normalization_retractionfiber(B1int_GCX_ratio)
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%  group_by(number, Signal_type) %>% 
                      mutate(ratio_normalized_retractionfiber = ifelse(Orientation=="Retraction fiber",
                                                                       ratio/ mean(ratio[Fiber_type=="Cell body"]),
                                                                       NA))
  
  
  
  
  B1int_GCX_ratio <- B1int_GCX_ratio[!duplicated(B1int_GCX_ratio), ]
  
  B1int_GCX_ratio$B1int_status <- sub("Negative", "Dim", B1int_GCX_ratio$B1int_status)
  B1int_GCX_ratio$B1int_status <- sub("Positive", "Bright", B1int_GCX_ratio$B1int_status)
  
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%  
    group_by(number, Position_relative_to_fiber, Signal_type) %>%  mutate(Gcx_normalized_leadingedge = ifelse(Orientation=="Leading edge protrusion",
                                                                                                              Gcx/mean(Gcx[Length_label=="Contact-free"], na.rm=TRUE),
                                                                                                              NA),
                                                                          Gcx_mean_leadingedge = ifelse(Orientation=="Leading edge protrusion",
                                                                                                     Gcx/mean(Gcx, na.rm=TRUE),
                                                                                                     NA)
                                                                          )
  B1int_GCX_ratio <- B1int_GCX_ratio %>%  
    group_by(number, Position_relative_to_fiber, Signal_type) %>%  mutate(B1int_normalized_leadingedge = ifelse(Orientation=="Leading edge protrusion",
                                                                                                                B1int/mean(B1int[Length_label=="Contact-free"], na.rm=TRUE),
                                                                                                              NA),
                                                                          B1int_mean_leadingedge = ifelse(Orientation=="Leading edge protrusion",
                                                                                                          B1int/mean(B1int, na.rm=TRUE),
                                                                                                     NA)
                                                                          )
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%
    group_by(number, Position_relative_to_fiber, Signal_type) %>% 
    mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                           ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                           NA),
           ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                     ratio/mean(ratio, na.rm=TRUE),
                                                     NA)
    )
  
  
  
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>%  
    group_by(number, Position_relative_to_fiber, Signal_type) %>% 
    mutate(Xthr_pos = case_when(
      X_micron < Xthr ~ "Distal",
      X_micron >= Xthr & X_micron <= 20 ~ "Proximal",
      X_micron > 20 ~ "Exceeded"  # This line keeps the original value if none of the conditions above are met
    ))
  
  
  
  B1int_GCX_ratio <- B1int_GCX_ratio %>% 
    group_by(Orientation, number, Signal_type) %>% 
    mutate(Gcx_intracellular = Gcx/ mean(Gcx, na.rm=TRUE),
           Gcx_classification_intracellular = ifelse(Gcx < (mean(Gcx, na.rm=TRUE) - (0.5 * sd(Gcx, na.rm=TRUE))), "Low", NA),
           Gcx_classification_intracellular = ifelse(Gcx > (mean(Gcx, na.rm=TRUE) + (0.5 * sd(Gcx, na.rm=TRUE))), "High", Gcx_classification_intracellular)) %>% 
    ungroup() %>% group_by(Signal_type,CellNR) %>% 
           mutate(Gcx_intracellular_cellbody = Gcx/ mean (Gcx[Type=="Cell body"]))
  
  
  
  Individual_ROIs <- B1int_GCX_ratio %>% 
    distinct(Orientation, CellNR, RoiNR) %>%  filter(Signal_type=="foreground")
  
  print(Individual_ROIs)
  
  
  Individual_Cells <- B1int_GCX_ratio %>% 
    distinct(Orientation, CellNR) %>%  filter(Signal_type=="foreground")
  
  print(Individual_Cells)
  
  B1int_GCX_ratio_summarized_perROI <- B1int_GCX_ratio %>% 
    group_by(number, RoiNR, Orientation, CellNR, Type, Fiber_type, Signal_type, Position_relative_to_fiber, Length_label) %>%    
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>% 
    filter(!is.nan(ratio_normalized_leadingedgeprotrusion))
  
  write.csv2(B1int_GCX_ratio_summarized_perROI, "Data_summary_per_ROI.csv")
  
  B1int_GCX_ratio_summarized_percell <- B1int_GCX_ratio %>% 
    group_by(Orientation, CellNR, Type, Fiber_type, Signal_type, Position_relative_to_fiber, Length_label) %>%    
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
  
  write.csv2(B1int_GCX_ratio_summarized_percell, "Data_summary_per_Cell.csv")
  
      
  
  
## make corrections for  protrusions------

  
  protrusion_data <- B1int_GCX_ratio



if(any(protrusion_data$Orientation=="Cellbody")){
  cellbody_summary_percell <- protrusion_data  %>% filter(Orientation=="Cellbody" & Signal_type=="foreground") %>%  
    group_by(Orientation, CellNR, Signal_type) %>% 
    dplyr::summarize(B1int_cellbody_sd = sd(B1int, na.rm=TRUE),
                     B1int_cellbody_max = max(B1int, na.rm=TRUE),
                     across(where(is.numeric), mean, na.rm=TRUE)
    ) %>% ungroup() %>% 
    dplyr::select(c(CellNR, B1int, B1int_cellbody_sd, Gcx, B1int_cellbody_max)) %>% 
    mutate(ratio_cellbody_average = B1int/ Gcx,
    ) %>% 
    rename(Gcx_cellbody_average = Gcx,
           B1int_cellbody_average = B1int) %>% 
    mutate(peak_b1int_threshold = B1int_cellbody_max+1)
  
}
protrusion_data <- protrusion_data %>%
  group_by(Signal_type, CellNR) %>% 
  left_join(cellbody_summary_percell)



## copy from here to other Orientations if needed
protrusion_data <- protrusion_data %>%  
  mutate_at(vars(Fiber_type, Type, Orientation, Xthr_pos), factor)


#protrusion_data <- calc_local_segregation_leadingedge(protrusion_data, 3)
protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)



protrusion_data <- protrusion_data %>%
  group_by(Orientation, number, Position_relative_to_fiber, Signal_type) %>% 
  mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                         ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                         NA),
         ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                   ratio/mean(ratio, na.rm=TRUE),
                                                   NA)
  )


protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Position_relative_to_fiber, Signal_type) %>% 
  mutate( X_relative = X/ max(row_number()),
          Length_label = ifelse(Orientation=="Retraction fiber" ,"Retraction fiber tip", Length_label),
          Length_label = ifelse(Orientation=="Retraction fiber" & X_relative<retraction_fiber_length_segregation_threshold & X_relative>(0.5*retraction_fiber_length_segregation_threshold),
                                   "Retraction fiber intermediate zone", Length_label),
          Length_label = ifelse(Orientation=="Retraction fiber" & X_relative>=retraction_fiber_length_segregation_threshold, "Perifiber zone", Length_label),
          X = ifelse(Orientation=="Retraction fiber", -X, X),
          X = ifelse(Orientation=="Retraction fiber", X - min(X), X),
          X_micron = ifelse(Orientation=="Retraction fiber", -X_micron, X_micron),
          X_micron = ifelse(Orientation=="Retraction fiber", X_micron - min(X_micron), X_micron),
          X_relative = ifelse(Orientation == "Retraction fiber", -X_relative, X_relative),
          X_relative = ifelse(Orientation == "Retraction fiber", X_relative - min(X_relative), X_relative)
)



protrusion_data <- protrusion_data %>% 
  mutate(Length_label = ifelse(Orientation=="Cellbody", "Cell body", Length_label ),
         Length_label = ifelse(Orientation=="Perpendicular control", "Membrane border", Length_label),
         Length_label = ifelse(Orientation=="Outwardclusters", "Protrusion Tip", Length_label),
         Length_label = ifelse(Orientation=="Innerclusters", "Membrane border", Length_label)
         )

protrusion_data$Type[protrusion_data$Type=="Focalized" & is.na(protrusion_data$PeakID_perprofile)] <- "Unspecified"


protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
  filter(Signal_type=="foreground")


library(zoo)



protrusion_data <- protrusion_data %>% 
  group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
  arrange(X) %>% 
  mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
         WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
         Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
         WithinPeak = TRUE
  )

protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
  arrange(X) %>% 
  mutate(Clustersize = max(cumsum(WithinPeak==TRUE)*pxsize),
         Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
         Clustersize = ifelse(Clustersize < Clustersize_threshold, NA, Clustersize),
         new_Peak_ID = ifelse(is.na(Clustersize), NA, new_Peak_ID),
         Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
         Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
  )

filter_b1_peak_threshold = TRUE

if(filter_b1_peak_threshold==TRUE){
  protrusion_data <- protrusion_data %>% 
    group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
    arrange(X) %>% 
    mutate(peak_b1int_threshold = peak_b1int_threshold
    )
}
if(filter_b1_peak_threshold==TRUE){
protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
  arrange(X) %>% 
  mutate(
    new_Peak_ID = ifelse(!(Orientation %in% c("Perpendicular control", "Innerclusters", "Outwardclusters", "Cellbody")) & max(B1int, na.rm=TRUE) < peak_b1int_threshold, NA, new_Peak_ID),
    Clustersize = ifelse(!(Orientation %in% c("Perpendicular control", "Innerclusters", "Outwardclusters", "Cellbody")) & is.na(new_Peak_ID), NA, Clustersize),
    Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
  )
}


# determine local B1 integrin and Gcx background for the integrin clusters
px.clustersize.micron.local = 0.6
px.clustersize.local = floor(px.clustersize.micron.local/pxsize)


protrusion_data <- protrusion_data %>%
  group_by(Orientation, number, new_Peak_ID) %>%
  arrange(X) %>% 
  mutate(min_B1int_local = ifelse(X==min(X) | X==max(X), B1int, NA),
         min_Gcx_local =   ifelse(X==min(X) | X==max(X), Gcx, NA),
         interval.local = 1
  )%>%
  ungroup() %>% 
  group_by(Orientation, number, Signal_type, new_Peak_ID) %>%
  mutate(B1int_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID), 
                                         zoo::na.approx(min_B1int_local, na.rm=FALSE, rule=2), NA),
         Gcx_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID),
                                       zoo::na.approx(min_Gcx_local, na.rm=FALSE, rule=2), NA)
  )%>% 
  dplyr::select(-c(min_B1int_local, min_Gcx_local, interval.local))





## determine global B1 integrin and Gcx background for the integrin clusters

px.clustersize.micron = 3.5
px.clustersize = floor(px.clustersize.micron/pxsize)

protrusion_data <- protrusion_data %>%
  group_by(Orientation, number, Signal_type) %>%
  arrange(X) %>%
  mutate(
    interval.end = ifelse(!is.na(B1int_cluster_local_bg), cumsum(!is.na(B1int_cluster_local_bg) & lead(is.na(B1int_cluster_local_bg))), NA),
    interval = ifelse(!is.na(interval.end) & lead(is.na(interval.end)), NA, interval.end),
    interval = as.factor(interval)
  ) %>% 
  ungroup()


# Define the custom function
run_protrusion_code <- function(data, px_clustersize) {
  result <- data %>% 
    group_by(Orientation, number, Signal_type, interval)  %>% 
    mutate(
      composite.length = ifelse(!is.na(interval), max(row_number()), NA),
      dividecluster = ifelse(!is.na(interval) & composite.length > px_clustersize, TRUE, FALSE),
      clusterX = ifelse(!is.na(interval) & dividecluster, row_number(), NA),
      pastthrescoord = ifelse(!is.na(interval) & dividecluster & clusterX > px_clustersize, TRUE, FALSE),
      maxedge = ifelse(!is.na(interval) & pastthrescoord, max(Peak_position == "Edge"), NA),
      lastedge = ifelse(!is.na(interval) & cumsum(pastthrescoord & Peak_position == "Edge") == maxedge, TRUE, FALSE),
      lastedge = ifelse(lastedge==FALSE & lag(lastedge==TRUE), TRUE, lastedge),
      lastedge = ifelse(lastedge==TRUE & lead(lastedge==TRUE),FALSE, lastedge),
      dividecoord = ifelse(!is.na(interval) & length(unique(pastthrescoord)) > 1 & lastedge, TRUE, FALSE),
      dividecoord = ifelse(is.na(dividecoord), FALSE, dividecoord),
      subinterval = cumsum(dividecoord)
    ) %>% 
    ungroup() %>%
    mutate(interval = ifelse(!is.na(interval), group_indices(., interval, subinterval), NA))
  
  return(result)
}

repeated_samplingdistance = 10

# Apply the custom function to your data twice
for(p in 1:repeated_samplingdistance){
  print(p)
  protrusion_data <- run_protrusion_code(protrusion_data, px.clustersize)
}


protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type, interval) %>% 
  arrange(X) %>% 
  mutate(
    min_B1int = ifelse(X==min(X, na.rm=TRUE), B1int, NA),
    min_B1int = ifelse(X==max(X, na.rm=TRUE), B1int, min_B1int),
    # min_B1int = ifelse(dividecoord==TRUE, B1int, min_B1int),
    min_Gcx = ifelse(X==min(X, na.rm=TRUE) | X==max(X, na.rm=TRUE), Gcx, NA),
  )

protrusion_data <- protrusion_data %>% 
  ungroup() %>% 
  group_by(Orientation, number, Signal_type) %>%
  mutate(B1int_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_B1int, na.rm=FALSE, rule=2), NA),
         Gcx_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_Gcx, na.rm=FALSE, rule=2), NA),
         B1int_cluster_global_bg = ifelse(B1int_cluster_local_bg < B1int_cluster_global_bg, B1int_cluster_local_bg, B1int_cluster_global_bg)
  )



measuring_bandwidth_edge = 0
measuring_bandwidth_peak = 0


protrusion_data <- protrusion_data %>% 
  group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
  arrange(X) %>% 
  mutate(X.cluster =  ifelse(!is.na(new_Peak_ID), X - min(X)+1, NA),
         Peak_position = "Transition",
         Peak_position = ifelse(X.cluster == min(X.cluster, na.rm=TRUE) | 
                                  X.cluster==max(X.cluster, na.rm=TRUE), "Edge", "Transition"),
         Peak_position = ifelse(!is.na(Peaklocation_X), "Peak", Peak_position)
  )

if(measuring_bandwidth_peak>0){
  for(g in 1:(measuring_bandwidth_peak)){
    print(g)
    protrusion_data <- protrusion_data %>% 
      group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
      arrange(X) %>% 
      mutate(Peak_position = ifelse(lag(Peak_position=="Peak"), "Peak", Peak_position),
             Peak_position = ifelse(lead(Peak_position=="Peak"), "Peak", Peak_position)
      )
  }
}


protrusion_data <- protrusion_data %>% 
  ungroup() %>% 
  group_by(number, Signal_type, LocationNR, interval, new_Peak_ID) %>% 
  arrange(X) %>% 
  mutate(Edgenr = cumsum(X.cluster>1 & Peak_position=="Edge")
  )

protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)





interval_randompeaks = 14

random_edges_cellbody=FALSE
if(random_edges_cellbody==TRUE){
protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number) %>% 
  mutate(new_Peak_ID = as.numeric(new_Peak_ID),
         new_Peak_ID = ifelse(Orientation=="Cellbody", ceiling(row_number()/interval_randompeaks), new_Peak_ID),
         new_Peak_ID = as.factor(new_Peak_ID)) 
  
protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, new_Peak_ID) %>% 
  arrange(X) %>% 
  mutate(X_percluster = row_number(),
          edges = ifelse(X==min(X) | X==max(X), B1int, NA),
         Peaklocation_X = ifelse(Orientation=="Cellbody", NA, Peaklocation_X),
         Peaklocation_X = ifelse(Orientation=="Cellbody" & X_percluster==floor(max(row_number())/2), X, Peaklocation_X),
         B1int_cluster_global_bg = ifelse(!is.na(interval) & Orientation=="Cellbody", zoo::na.approx(edges, na.rm=FALSE, rule=2), NA),
         Peak_position = ifelse(Orientation=="Cellbody", NA, Peak_position),
         Peak_position = ifelse(Orientation=="Cellbody" & X==min(X) | Orientation=="Cellbody"& X==max(X), "Edge", Peak_position),
         Peak_position = ifelse(Orientation=="Cellbody" &  !is.na(Peaklocation_X), "Peak", Peak_position),
  ) %>% dplyr::select(-c(edges))
}


protrusion_data <- protrusion_data %>%  
  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>% 
  mutate(unique_Peak_ID = cur_group_id())

protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type, new_Peak_ID, unique_Peak_ID) %>% 
  mutate(B1int_global_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_global_bg, NA),
         B1int_local_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_local_bg, NA),
         Gcx_global_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_global_bg, NA),
         Gcx_local_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_local_bg, NA),
         Gcx_enrichment = ifelse(X==Peaklocation_X, mean(Gcx[Peak_position=="Peak"], na.rm=TRUE) / mean(Gcx[Peak_position=="Edge"], na.rm=TRUE), NA),
         B1int_enrichment = ifelse(X==Peaklocation_X, mean(B1int[Peak_position=="Peak"], na.rm=TRUE) / mean(B1int[Peak_position=="Edge"], na.rm=TRUE), NA),
        
  )

protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type) %>% 
  arrange(X) %>% 
  mutate(X_relative = X/ max(row_number(), na.rm=TRUE),
         Length_label= ifelse(Orientation=="Bleb" & X_relative>= (0/12*max(X_relative)) & X_relative<=(12/12*(max(X_relative))), "Peribleb region", Length_label),
         Length_label= ifelse(Orientation=="Bleb" & X_relative>= (3.5/12*max(X_relative)) & X_relative<=(8.5/12*(max(X_relative))), "Turning point", Length_label),
         Length_label= ifelse(Orientation=="Bleb" & X_relative>= (4.5/12*max(X_relative)) & X_relative<=(7.5/12*(max(X_relative))), "Bleb", Length_label),
         Length_label= ifelse(Orientation=="Bleb" & (X) == round(0.5*max(X)), "Apex", Length_label),
         Length_label= ifelse(Orientation=="Bleb" & lag(Length_label=="Apex") | lead(Length_label=="Apex"), "Apex", Length_label),
         Length_label= ifelse(Orientation=="Bleb" & is.na(Length_label), "Peribleb region", Length_label),
         Gcx_bleb_base_normalized = Gcx/ mean(Gcx[Length_label=="Peribleb region"]),
         B1int_bleb_base_normalized = B1int / mean(B1int[Length_label=="Peribleb region"]),
         Col_bleb_base_normalized = Col/ mean(Col[Length_label=="Peribleb region"]),
         Bleb_length = length(Gcx[Length_label=="Peribleb region"])*pxsize
  ) %>% arrange(X) %>% 
  mutate(Length_classification_nr = cumsum(X>1 & Length_label != lag(Length_label)),
         Length_classification_nr = ifelse(Orientation=="Perpendicular control"| Orientation=="Outwardclusters", 1, Length_classification_nr)
  )


protrusion_data <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type) %>% 
  mutate(B1int_cellbody_corrected = B1int/ B1int_cellbody_average,
         Gcx_cellbody_corrected = Gcx/ Gcx_cellbody_average)


protrusion_data_summarized_perROI <- protrusion_data %>% 
  group_by(Orientation, number, Signal_type, Length_label) %>% 
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))

protrusion_data_summarized_percell <- protrusion_data %>% 
  group_by(Orientation, CellNR, Type, Fiber_type, Signal_type) %>% 
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))

subset_protrusion <- protrusion_data %>%  filter(is.na(Peak_position)) 

summary_sds_protrusion <- subset_protrusion %>%  group_by(Orientation, number) %>%  summarize(B1int_sd= sd(B1int, na.rm=TRUE))




  ## copy to here for other Orientations
  
   ## all data -  save all lineplots to depict Gcx and B1 integrin intensity profiles------


  unique_numbers = unique(protrusion_data$number)

for(u in 1:length(unique_numbers)){
  print(u)
  unique_number <- unique_numbers[u]
  sub <- subset(protrusion_data, number==unique_number)
  
  legendsize=0.3
  
  
  Zdiff = 0.25
  
  sub_summary <- sub %>%
    filter(Signal_type=="foreground") %>%
    summarize(mean_B1int = mean(B1int_smooth),
              sd_B1int = sd(B1int_smooth),
              iv = mean_B1int+sd_B1int,
              thr = mean_B1int + sd_B1int*Zdiff)
  
  
  sub$Signal_type <- as.factor(sub$Signal_type)
  sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
  
  sub <- sub %>%  filter(Signal_type=="foreground")
  
  lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int, shape=NULL))+
    geom_line( color="magenta", alpha=1)+
    geom_line(data=sub, aes(y=B1int_cluster_global_bg),  color="magenta", alpha=1, linetype="dotted")+
    geom_line(data=sub, aes(y=B1int_cluster_local_bg), color="magenta", alpha=1, linetype="dashed")+
    theme_classic()+
    scale_y_continuous(limits=c(0,max(sub$B1int*1.1)))
  
  
  
  lengthlabels <- sub %>% ungroup() %>% 
    filter(Signal_type=="foreground") %>% filter(!is.na(Length_classification_nr)) %>% 
    group_by(Length_label, Length_classification_nr) %>% 
    summarize(
      xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf)
  
  
  lengthlabel.rectangles <- geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Length_label), 
                                      color= alpha("black", 0), inherit.aes = FALSE, alpha = 0.2, size=0.5)
  

  
  
  lineplot2 <- ggplot(data=sub, aes(x=X_micron, y=Gcx))+
    geom_line( color="gold3", alpha=1)+
    geom_line(data=sub, aes(y=Gcx_cluster_local_bg), color="gold3", alpha=1, linetype="dashed")+
    theme_classic()+
    scale_y_continuous(limits=c(0,max(sub$Gcx*1.2)))
  
  rect_data <- sub %>%
    filter(Signal_type=="foreground") %>%
    group_by(new_Peak_ID) %>% 
    summarize(Gcx_local_enrichment = mean(Gcx_local_enrichment, na.rm=TRUE),
              Gcx_global_enrichment = mean(Gcx_global_enrichment, na.rm=TRUE),
              B1int_global_enrichment = mean(B1int_global_enrichment, na.rm=TRUE),
              B1int_position = max(B1int, na.rm=TRUE) *1.05,
              Gcx_position = max(Gcx, na.rm=TRUE)*1.05,
              xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
              xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize-pxsize
    ) %>%  filter(!is.na(new_Peak_ID))
  
  clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("black", 0.3))
  
  cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", 
                                  fill=alpha("black",0), color=alpha("black", 0.3),
                                  inherit.aes = FALSE, alpha = 0, size=0.5)
  
  if(unique(sub$Orientation!="Outwardclusters")){
    bg <- geom_hline(yintercept=unique(sub$peak_b1int_threshold), color="magenta", alpha=0.5)
    lineplot <- lineplot + bg
  }else{
    lineplot <- lineplot + geom_line(data=sub, aes(y=Gcx), color="gold3")
  }
  
  annotations <-  geom_text(data=rect_data, aes(label = round(B1int_global_enrichment,digits=1), x=xintercept, y=B1int_position), vjust = -0.5, hjust = 0.5, inherit.aes=FALSE)
  
  lineplot <- lineplot + lengthlabel.rectangles + annotations
  
  annotations <-  geom_text(data=rect_data, aes(label = round(Gcx_local_enrichment,digits=1), x=xintercept, y=Gcx_position), vjust = -0.5, hjust = 0.5, inherit.aes=FALSE)
  
  
  lineplot2 <- lineplot2 +lengthlabel.rectangles + annotations
  
  
  
  
  lineplot <- lineplot + cluster.rectangles+ scale_y_continuous(expand=c(0,0.2))+
    labs(linetype="", x="Distance (\u03BCm)", y = " \u03B21int intensity")
  lineplot2 <- lineplot2 + cluster.rectangles+ scale_y_continuous(expand=c(0,0.2))+
    labs(linetype="", x="Distance (\u03BCm)", y="Gcx intensity")
  
  
  
  
  lineplot <- lineplot + clustermids+
    labs(linetype="", x="Distance (\u03BCm)")
  lineplot2 <- lineplot2 + clustermids+
    labs(linetype="", x="Distance (\u03BCm)")+
    theme(legend.position = "none")
  
  
  
  lineplot <- lineplot + 
    theme(legend.position="top")+
    theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
          legend.key.height = unit(legendsize, 'cm'), #change legend key height
          legend.key.width = unit(legendsize*1.5, 'cm'))
  
  
  if(unique(sub$Orientation)=="Innerclusters"){
    lineplot <- lineplot + geom_vline(xintercept = round(max(sub$X)/2)*pxsize, color="black")
    lineplot2 <- lineplot2 + geom_vline(xintercept = round(max(sub$X)/2)*pxsize, color="black")
    
  }
  
  
  lineplot <- plot_grid(lineplot, lineplot2, nrow=2, rel_heights = c(0.6,0.4))
  
  
  
  
  
  lineplot <-  plot_grid(lineplot, labels=paste0(unique(sub$Orientation), " - Number: ",unique_number), vjust=0)+
    theme(plot.margin = margin(t = 5, unit = "mm"))
  

  
  print(lineplot)
  

  pdf(paste0("Lineplot_", unique_number,"_", paste0(unique(sub$CellNR)), "_",paste0(unique(sub$RoiNR)),"_", paste0(unique(sub$Orientation)), ".pdf"), width=4, height=2.5)
  print(lineplot)
  dev.off()
  
  
  if(unique(sub$Orientation=="Bleb")){
    
    lengthlabels.bleb.example.table <- sub %>% ungroup() %>% 
      filter(Signal_type=="foreground") %>% filter(!is.na(Length_classification_nr)) %>% 
      group_by(Length_label, Length_classification_nr) %>% 
      summarize(
        xmin =0, xmax = 1, ymin = -Inf, ymax = Inf)
    
    
    lengthlabels.bleb.example <- geom_rect(data = lengthlabels.bleb.example.table, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Length_label), 
                                           color= alpha("black", 0), inherit.aes = FALSE, alpha = 0.2, size=0.5)
  }
  
  
  
  
  
}



  
  
  

   ## bleb data -  save plot of combined lines of Gcx and B1 integrin intensity profiles------

bleb_lines <- subset(protrusion_data, Orientation=="Bleb")

sub <- bleb_lines
legendsize=0.3


Zdiff = 0.25

sub_summary <- sub %>%
  filter(Signal_type=="foreground") %>%
  summarize(mean_B1int = mean(B1int_smooth),
            sd_B1int = sd(B1int_smooth),
            iv = mean_B1int+sd_B1int,
            thr = mean_B1int + sd_B1int*Zdiff)


sub$Signal_type <- as.factor(sub$Signal_type)
sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))

sub <- sub %>%  filter(Signal_type=="foreground")

lineplot <- ggplot(data=sub, aes(x=X_relative, shape=NULL))+
  geom_line(aes(y=B1int_cellbody_corrected, group=number), color="magenta", alpha=0.2)+
  geom_smooth(aes(y=B1int_cellbody_corrected), color="magenta", alpha=1)+
  theme_classic()


lengthlabels <- sub %>% ungroup() %>% 
  filter(Signal_type=="foreground") %>% filter(!is.na(Length_classification_nr)) %>% 
  group_by(Length_label, Length_classification_nr) %>% 
  summarize(
    xmin = min(X_relative), xmax = max(X_relative), ymin = -Inf, ymax = Inf)


lengthlabel.rectangles <- geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Length_label), 
                                    color= alpha("black", 0.05), inherit.aes = FALSE, alpha = 0.2, size=0.5)

lineplot <- lineplot + lengthlabel.rectangles

lineplot2 <- ggplot(data=sub, aes(x=X_relative, shape=NULL))+
  geom_smooth(aes(y=Gcx_bleb_base_normalized), color="gold3", alpha=1)+
  geom_line(aes(y=Gcx_bleb_base_normalized, group=number), color="gold3", alpha=0.5)+
  theme_classic()

lineplot2 <- lineplot2 + lengthlabel.rectangles



lineplot <- lineplot  +
  labs(linetype="", x="Relative distance", y = "Normalized \u03B21int")+ theme(legend.position="top")


lineplot2 <- lineplot2  +
  labs(linetype="", x="Relative distance", y = "Normalized Gcx")+ theme(legend.position="top")+
  theme(legend.position = "none")



lineplot <- lineplot + 
  theme(legend.position="top")+
  theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
        legend.key.height = unit(legendsize, 'cm'), #change legend key height
        legend.key.width = unit(legendsize*1.5, 'cm'))+
  labs(fill="")



lineplot.bleb <-  plot_grid(lineplot, lineplot2, ncol=1, vjust=0, rel_heights = c(0.6,0.4))+
  theme(plot.margin = margin(t = 5, unit = "mm"))

print(lineplot.bleb)
  
  pdf(paste0("Lineplot_blebs_combined.pdf"), width=8, height=5)
  print(lineplot.bleb)
  dev.off()
  









   ## all data -  save all lineplots to depict Gcx and B1 integrin intensity profiles with Gcx and B1int combined in one plot------


unique_numbers = unique(protrusion_data$number)



for(u in 1:length(unique_numbers)){
  print(u)
  unique_number <- unique_numbers[u]
  sub <- subset(protrusion_data, number==unique_number)
  
  cellnr <- unique(sub$CellNR)
  celldata <- subset(protrusion_data, CellNR==cellnr & Signal_type=="foreground")
  mingcx_cell <- min(celldata$Gcx)
  maxgcx_cell <- max(celldata$Gcx)
  minb1int_cell <- min(celldata$B1int)
  maxb1int_cell <- max(celldata$B1int)
  
  limits_conversion = maxb1int_cell/ maxgcx_cell
  
  lim = ifelse(limits_conversion <1 , maxb1int_cell, maxgcx_cell)

  
  legendsize=0.3
  
  
  Zdiff = 0.25
  
  sub_summary <- sub %>%
    filter(Signal_type=="foreground") %>%
    summarize(mean_B1int = mean(B1int_smooth),
              sd_B1int = sd(B1int_smooth),
              iv = mean_B1int+sd_B1int,
              thr = mean_B1int + sd_B1int*Zdiff)
  
  
  sub$Signal_type <- as.factor(sub$Signal_type)
  sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
  
  sub <- sub %>%  filter(Signal_type=="foreground")
  
  
  

  lineplot <- ggplot(data=sub, aes(x=X_micron, shape=NULL))+
    geom_line(aes(y=B1int_cellbody_corrected), color="magenta", alpha=1)+
    geom_line(aes(y=Gcx_cellbody_corrected* max(B1int_cellbody_corrected)/ max(Gcx_cellbody_corrected)),color="gold3", alpha=1)+
    theme_classic()+
    scale_y_continuous(expand=c(0,0.2))+
    scale_y_continuous(
      name = "B1int",
      sec.axis = sec_axis(~./max(sub$B1int_cellbody_corrected) * max(sub$Gcx_cellbody_corrected), name="Gcx")
    )
  
  lineplot <- ggplot(data=sub, aes(x=X_micron, shape=NULL))+
    geom_line(aes(y=B1int/1000), color="magenta", alpha=1)+
    geom_line(aes(y=Gcx/1000 * ((maxb1int_cell/1000)/ (maxgcx_cell/1000))),color="gold3", alpha=1)+
    theme_classic()+
    scale_y_continuous(expand=c(0,0.2))+
    scale_y_continuous(
      name = "\u03B21int (x1000)",
      limits = c(0, maxb1int_cell/1000),
      sec.axis = sec_axis(~./(maxb1int_cell/maxgcx_cell), name="Gcx (x1000)")
    )+
    labs(linetype="", x="Distance (\u03BCm)")

  lengthlabels <- sub %>% ungroup() %>% 
    filter(Signal_type=="foreground") %>% filter(!is.na(Length_classification_nr)) %>% 
    group_by(Length_label, Length_classification_nr) %>% 
    summarize(
      xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf)
  
  
  lengthlabel.rectangles <- geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Length_label), 
                                      color= alpha("black", 0), inherit.aes = FALSE, alpha = 0.2, size=0.5)
  
  
  
  lineplot2 <- ggplot(data=sub, aes(x=X_micron, y=Gcx))+
    geom_line( color="gold3", alpha=1)+
    geom_line(data=sub, aes(y=Gcx_cluster_local_bg), color="gold3", alpha=1, linetype="dashed")+
    theme_classic()
  
  rect_data <- sub %>%
    filter(Signal_type=="foreground") %>%
    group_by(new_Peak_ID) %>% 
    summarize(Gcx_local_enrichment = mean(Gcx_local_enrichment, na.rm=TRUE),
              Gcx_global_enrichment = mean(Gcx_global_enrichment, na.rm=TRUE),
              B1int_global_enrichment = mean(B1int_global_enrichment, na.rm=TRUE),
              B1int_position = max(B1int, na.rm=TRUE) *1.05,
              Gcx_position = max(Gcx, na.rm=TRUE)*1.05,
              xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
              xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize-pxsize
    ) %>%  filter(!is.na(new_Peak_ID))
  
  clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("black", 0.3))
  
  cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", 
                                  fill=alpha("black",0), color=alpha("black", 0.3),
                                  inherit.aes = FALSE, alpha = 0, size=0.5)
  
  annotations <-  geom_text(data=rect_data, aes(label = round(B1int_global_enrichment,digits=1), x=xintercept, y=B1int_position), vjust = -0.5, hjust = 0.5, inherit.aes=FALSE)
  
  lineplot <- lineplot + lengthlabel.rectangles 
  
  annotations <-  geom_text(data=rect_data, aes(label = round(Gcx_local_enrichment,digits=1), x=xintercept, y=Gcx_position), vjust = -0.5, hjust = 0.5, inherit.aes=FALSE)
  
  
  lineplot2 <- lineplot2 +lengthlabel.rectangles
  
  
  
  
  lineplot <- lineplot + cluster.rectangles+ clustermids
    labs(linetype="", x="Distance (\u03BCm)", y = "Norm. Intensity")


  
  lineplot <- lineplot + 
    theme(legend.position="top")+
    theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
          legend.key.height = unit(legendsize, 'cm'), #change legend key height
          legend.key.width = unit(legendsize*1.5, 'cm'))
  
  

  
  
  lineplot <-  plot_grid(lineplot, labels=paste0("Cell ", unique(sub$CellNR), 
                                                 "- ROI ",paste0(unique(sub$RoiNR)),
                                                 " - ", paste0(unique(sub$Orientation)), 
                                                 " - ID ",unique_number), vjust=0, label_size=10)+
    theme(plot.margin = margin(t = 5, unit = "mm"))
  
  print(lineplot)
  
  
  pdf(paste0("Lineplot_", unique_number,"_", paste0(unique(sub$CellNR)), "_",paste0(unique(sub$RoiNR)),"_", paste0(unique(sub$Orientation)), ".pdf"), width=2, height=1.7)
  print(lineplot)
  dev.off()
  
}







   ## example data -  save all lineplots to depict Gcx and B1 integrin intensity profiles with Gcx and B1int combined in one plot------

example <- subset(protrusion_data, CellNR==7 & RoiNR==1 & Orientation=="Leading edge protrusion")

unique_numbers = unique(example$number)



for(u in 1:length(unique_numbers)){
  print(u)
  unique_number <- unique_numbers[u]
  sub <- subset(example, number==unique_number)
  
  legendsize=0.3
  
  
  Zdiff = 0.25
  
  sub_summary <- sub %>%
    filter(Signal_type=="foreground") %>%
    summarize(mean_B1int = mean(B1int_smooth),
              sd_B1int = sd(B1int_smooth),
              iv = mean_B1int+sd_B1int,
              thr = mean_B1int + sd_B1int*Zdiff)
  
  
  sub$Signal_type <- as.factor(sub$Signal_type)
  sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
  
  
  lengthlabels <- sub %>% ungroup() %>% 
    filter(Signal_type=="foreground") %>% filter(!is.na(Length_classification_nr)) %>% 
    group_by(Length_label, Length_classification_nr) %>% 
    summarize(
      xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf)
  
  
  lengthlabel.rectangles <- ggplot(data=lengthlabels) + geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Length_label), 
                                                                  color= alpha("black", 0), inherit.aes = FALSE, alpha = 0.2, size=0.5)
  
  
  sub <- sub %>%  filter(Signal_type=="foreground")
  
  lineplot <- lengthlabel.rectangles+
    geom_line(data=sub, aes(y=B1int/ max(B1int), x=X_micron), color="magenta", alpha=1)+
    geom_line(data=sub, aes(y=Gcx/ max(Gcx), x=X_micron),color="gold3", alpha=1)+
    theme_classic()+
    scale_y_continuous(expand=c(0,0.2))
  
  
  
  

  rect_data <- sub %>%
    filter(Signal_type=="foreground") %>%
    group_by(new_Peak_ID) %>% 
    summarize(Gcx_local_enrichment = mean(Gcx_local_enrichment, na.rm=TRUE),
              Gcx_global_enrichment = mean(Gcx_global_enrichment, na.rm=TRUE),
              B1int_global_enrichment = mean(B1int_global_enrichment, na.rm=TRUE),
              B1int_position = max(B1int, na.rm=TRUE) *1.05,
              Gcx_position = max(Gcx, na.rm=TRUE)*1.05,
              xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
              xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize-pxsize
    ) %>%  filter(!is.na(new_Peak_ID))
  
  clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("black", 0.2))
  
  cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", 
                                  fill=alpha("black",0), color=alpha("black", 0.2),
                                  inherit.aes = FALSE, alpha = 0, size=0.5)
  
  annotations <-  geom_text(data=rect_data, aes(label = round(B1int_global_enrichment,digits=1), x=xintercept, y=B1int_position), vjust = -0.5, hjust = 0.5, inherit.aes=FALSE)
  
  lineplot <- lineplot  + clustermids
  
  annotations <-  geom_text(data=rect_data, aes(label = round(Gcx_local_enrichment,digits=1), x=xintercept, y=Gcx_position), vjust = -0.5, hjust = 0.5, inherit.aes=FALSE)
  

  
  
  
  lineplot <- lineplot + cluster.rectangles+
    labs(linetype="", x="Distance (\u03BCm)", y = "Norm. Intensity")
  
  
  
  lineplot <- lineplot + 
    theme(legend.position="top")+
    theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
          legend.key.height = unit(legendsize, 'cm'), #change legend key height
          legend.key.width = unit(legendsize*1.5, 'cm'))
  
  
  
  
  
  
  lineplot <-  plot_grid(lineplot, labels=paste0(unique(sub$Orientation), " - Number: ",unique_number), vjust=0)+
    theme(plot.margin = margin(t = 5, unit = "mm"))+
    theme(text = element_text(size = 7))
  
  print(lineplot)
  
  
  pdf(paste0("Lineplot_", unique_number,"_", paste0(unique(sub$CellNR)), "_",paste0(unique(sub$RoiNR)),"_", paste0(unique(sub$Orientation)), ".pdf"), width=3, height=2, family="Helvetica")
  print(lineplot)
  dev.off()
  
}








   ## outward clusters -  save all lineplots to depict Gcx and B1 integrin intensity profiles------

outwardprotrusion_data <- subset(protrusion_data, Orientation=="Outwardclusters" | Orientation == "Innerclusters")

outwardprotrusion_data <- outwardprotrusion_data %>%  group_by(grouped_number) %>%  mutate(Gcx_max = max(Gcx, na.rm=TRUE), 
                                                                                           B1int_max = max(B1int, na.rm=TRUE),
                                                                                           Gcx_min = min(Gcx,na.rm=TRUE),
                                                                                           B1int_min= min(B1int),na.rm=TRUE)

unique_numbers = unique(outwardprotrusion_data$number)


for(u in 1:length(unique_numbers)){
  print(u)
  unique_number <- unique_numbers[u]
  sub <- subset(outwardprotrusion_data, number==unique_number)
  
  legendsize=0.3
  
  
  Zdiff = 0.25
  
  sub <- sub %>%  filter(Signal_type=="foreground")
  
  lineplot <- ggplot(data=sub, aes(x=X_micron))+
    geom_line( color="magenta", aes(y=B1int), alpha=1)+
    geom_line( color="gold3", aes(y=Gcx), alpha=1)+
    theme_classic()
  

  
  
  rect_data <- sub %>%
    group_by(number) %>%
    summarize(B1int_peak_location = (X[which.max(B1int)])*pxsize-pxsize,
              Gcx_peak_location = (X[which.max(Gcx)])*pxsize-pxsize,
              Segregation_distance = (Gcx_peak_location - B1int_peak_location)
    ) 

  
  
  clustermids_b1int <-  geom_vline(data=rect_data, aes(xintercept=B1int_peak_location), linetype="longdash", color="magenta")
  clustermids_gcx <- geom_vline(data=rect_data, aes(xintercept=Gcx_peak_location), linetype="longdash", color="gold3")
  
  
  lineplot <- lineplot    + clustermids_b1int + clustermids_gcx + scale_x_continuous(breaks=c(0,1,2))
  

  
  
  
  lineplot <- lineplot +
    labs(linetype="", x="Distance (\u03BCm)", y = "Fraction of maximum intensity")
  
  
  lineplot <- lineplot + 
    theme(legend.position="top")+
    theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
          legend.key.height = unit(legendsize, 'cm'), #change legend key height
          legend.key.width = unit(legendsize*1.5, 'cm'))
  
  
  lineplot <-  plot_grid(lineplot, labels=paste0("Cell ", unique(sub$CellNR), 
                                                 "- ROI ",paste0(unique(sub$RoiNR)),
                                                 " - ", paste0(unique(sub$Orientation)), 
                                                 " - ID ",unique_number), vjust=0, label_size=10)+
    theme(plot.margin = margin(t = 5, unit = "mm"))
  
  print(lineplot)
  
  
  pdf(paste0("Lineplot_", unique_number,"_", paste0(unique(sub$CellNR)), "_",paste0(unique(sub$RoiNR)),"_", paste0(unique(sub$Orientation)), ".pdf"), width=2, height=1.7)
  print(lineplot)
  dev.off()
  
}











   ## inner clusters -  save all lineplots to depict Gcx and B1 integrin intensity profiles------

outwardprotrusion_data <- subset(protrusion_data, Orientation=="Innerclusters")

unique_numbers = unique(outwardprotrusion_data$number)


for(u in 1:length(unique_numbers)){
  print(u)
  unique_number <- unique_numbers[u]
  sub <- subset(outwardprotrusion_data, number==unique_number)
  
  legendsize=0.3
  
  
  Zdiff = 0.25
  
  sub <- sub %>%  filter(Signal_type=="foreground")
  
  lineplot <- ggplot(data=sub, aes(x=X_micron))+
    geom_line( color="magenta", aes(y=B1int/max(B1int)), alpha=1)+
    geom_line( color="magenta", aes(y=B1int_cluster_local_bg/max(B1int)), linetype="dashed")+
    geom_line( color="gold3", aes(y=Gcx/ max(Gcx)), alpha=1)+
    theme_classic()
  
  
  rect_data <- sub %>%
    group_by(number) %>%
    summarize(B1int_peak_location = (X[which.max(B1int)])*pxsize-pxsize,
              Gcx_peak_location = (X[which.max(Gcx)])*pxsize-pxsize,
              Segregation_distance = (Gcx_peak_location - B1int_peak_location)
    ) 
  
  
  
  clustermids_b1int <-  geom_vline(data=rect_data, aes(xintercept=max(sub$X_micron)/2), linetype="longdash", color="magenta")
 # clustermids_gcx <- geom_vline(data=rect_data, aes(xintercept=Gcx_peak_location), linetype="longdash", color="gold3")
  
  
  lineplot <- lineplot    + clustermids_b1int + clustermids_gcx + scale_x_continuous(breaks=c(0,1,2)) 
  
  
  
  
  
  lineplot <- lineplot +scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
    labs(linetype="", x="Distance (\u03BCm)", y = "Fraction of maximum intensity")
  
  
  lineplot <- lineplot + 
    theme(legend.position="top")+
    theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
          legend.key.height = unit(legendsize, 'cm'), #change legend key height
          legend.key.width = unit(legendsize*1.5, 'cm'))
  
  
  lineplot <-  plot_grid(lineplot, labels=paste0("Cell ", unique(sub$CellNR), 
                                                 "- ROI ",paste0(unique(sub$RoiNR)),
                                                 " - ", paste0(unique(sub$Orientation)), 
                                                 " - ID ",unique_number), vjust=0, label_size=10)+
    theme(plot.margin = margin(t = 5, unit = "mm"))
  
  print(lineplot)
  
  
  pdf(paste0("Lineplot_", unique_number,"_", paste0(unique(sub$CellNR)), "_",paste0(unique(sub$RoiNR)),"_", paste0(unique(sub$Orientation)), ".pdf"), width=2, height=1.7)
  print(lineplot)
  dev.off()
  
}












## micron-scale segregation data ----
   ## protrusions - plot of Gcx tip versus base in clusters -----
  
  clusterproperties <- protrusion_data %>%  group_by(number, CellNR, Orientation, new_Peak_ID, Signal_type, Clustersize, unique_Peak_ID, Length_label) %>%  
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))%>%  
      filter(Orientation=="Leading edge protrusion" | Orientation=="Retraction fiber") %>% 
    filter(Orientation=="Leading edge protrusion" & !is.na(Peaklocation_X) | 
           Orientation=="Retraction fiber" & !is.na(Peaklocation_X))



noncluster_data <- protrusion_data %>%  group_by(number, CellNR, Orientation, Signal_type, Length_label) %>%  
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Orientation=="Cellbody" | Orientation=="Bleb")


bleb_cluster_correction = FALSE

if(bleb_cluster_correction==TRUE){
  noncluster_data <- noncluster_data %>%  
    filter(is.nan(Peaklocation_X))
}

clusterproperties <- rbind(clusterproperties, noncluster_data)

clusterproperties <- clusterproperties %>%  mutate(unique_Peak_ID = ifelse(is.na(new_Peak_ID), NA, unique_Peak_ID))


clusterproperties <- clusterproperties %>%  filter(Length_label!="Cell body")

  clusterproperties$new_Peak_ID <- as.factor(clusterproperties$new_Peak_ID)
  clusterproperties$number <- as.factor(clusterproperties$number)
  clusterproperties$Length_label <- as.factor(clusterproperties$Length_label)
  clusterproperties$Length_label <- factor(clusterproperties$Length_label, levels=c("Peribleb region", "Turning point", "Bleb", "Apex", "Protrusion tip", 
                                                                                    "Intermediate zone", "Protrusion base", "Retraction fiber base", "Perifiber zone", 
                                                                                    "Retraction fiber intermediate zone", "Retraction fiber tip"))


  clusterproperties <- clusterproperties %>% 
    group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, unique_Peak_ID, Length_label) %>% 
    mutate(Gcx_cellbody_corrected = Gcx/ Gcx_cellbody_average,
           B1int_cellbody_corrected = B1int/ B1int_cellbody_average,
           ratio = B1int / Gcx,
           ratio_cellbody_corrected = ratio / ratio_cellbody_average)
  
  compare.data <- function(clusteproperties, subset){
  clusterproperties <- subset(clusterproperties, Orientation==subset)
  clusterproperties$Length_label <- droplevels(clusterproperties$Length_label)
  print(levels(clusterproperties$Length_label))
  print(unique(clusterproperties$Length_label))
  
  lengthlabel.vs.gcx <- custom_plot_with_statistics(data=clusterproperties, x_var=("Length_label"), group_vars="Orientation", alpha=0.5, size=2,
                               y_var="Gcx_cellbody_corrected", jitter=FALSE, annotate_signif = TRUE, fill_var=NULL,
                               annotate_only_signif = TRUE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
  
  lengthlabel.vs.gcx.plot <- lengthlabel.vs.gcx$plot
  lengthlabel.vs.gcx.plot <- lengthlabel.vs.gcx.plot + theme_classic()+
    labs(y="Normalized GCX", x="")
  print(lengthlabel.vs.gcx.plot)
  
  
  print(ggdraw(lengthlabel.vs.gcx.plot))
  
  return(lengthlabel.vs.gcx.plot)
  
  }
  
  
  
  clusterproperties <- clusterproperties %>%  mutate(CellNR = as.numeric(CellNR),
                                                     expnr = ifelse(CellNR<100, 1, 2),
                                                     expnr = ifelse(CellNR>=200, 3, expnr)
  )
  
  leadingedge.gcx <- compare.data(clusterproperties, "Leading edge protrusion") +   scale_x_discrete(labels = c("Tip", "Shaft", "Base")) 
  
  
  leadingedge.stats <- subset(clusterproperties, Orientation=="Retraction fiber")
  
  
  pairwise.results <- pairwise.wilcox.test(x=leadingedge.stats$Gcx_cellbody_corrected, g=leadingedge.stats$Length_label,p.adjust.method = "bonferroni")
  print(pairwise.results)
  
  library(rstatix)
    leadingedge.stats %>%  ungroup() %>% rstatix::kruskal_test(Gcx_cellbody_corrected ~ Length_label) 
    leadingedge.stats %>%  ungroup() %>% rstatix::kruskal_effsize(Gcx_cellbody_corrected ~ Length_label)
  
  rstatix::kruskal_effsize(leadingedge.stats, Gcx_cellbody_corrected ~ Length_label, ci = TRUE)
  
  leadingedge.stats <- leadingedge.stats %>%
    ungroup() %>% 
    group_by(expnr, Length_label) %>%
    arrange(Orientation) %>% 
    
    dplyr::summarise(
      n = n(),
      median = median(Gcx_cellbody_corrected, na.rm = TRUE),
      Q1 = quantile(Gcx_cellbody_corrected, 0.25, na.rm = TRUE),
      Q3 = quantile(Gcx_cellbody_corrected, 0.75, na.rm = TRUE),
      IQR = IQR(Gcx_cellbody_corrected, na.rm = TRUE),
      lower_whisker = max(min(Gcx_cellbody_corrected, na.rm = TRUE), Q1 - 1.5 * IQR),
      upper_whisker = min(max(Gcx_cellbody_corrected, na.rm = TRUE), Q3 + 1.5 * IQR)
    ) %>%
    ungroup() 
  
  print(leadingedge.stats)
  
  shapiro.test(leadingedge.stats$Gcx_cellbody_corrected[leadingedge.stats$Length_label=="Protrusion base"])
  

  pairwise.results <- pairwise.wilcox.test(x=leadingedge.stats$median, g=leadingedge.stats$Length_label,p.adjust.method = "bonferroni")
  print(pairwise.results)
  leadingedge.stats %>%  aov(median ~ Length_label, data = .) %>%  summary()
  leadingedge.stats %>%   tukey_hsd(median ~ Length_label) 
  dunn.test::dunn.test(x=leadingedge.stats$median, g=leadingedge.stats$Length_label,
                       method = "bonferroni")
  

  
  bleb.gcx <- compare.data(clusterproperties, "Bleb")+   scale_x_discrete(labels = c("Peribleb zone", "Neck", "Bleb", "Apex")) 
  retraction.gcx <- compare.data(clusterproperties, "Retraction fiber")+   scale_x_discrete(labels = c("Base", "Shaft", "Tip")) 
  data.together.gcx <- plot_grid(leadingedge.gcx, bleb.gcx, retraction.gcx, nrow=1, align="h")
  print(data.together.gcx)
  
  pdf("Gcx.violins.together.pdf", width=8, height=3)
  print(data.together.gcx)
  dev.off()
  
  
  individual_datapoints <- clusterproperties %>%  subset(Orientation=="Retraction fiber") %>% ungroup() %>% 
    dplyr::select(Length_label, Gcx_cellbody_corrected) %>% 
    filter(!is.nan(Gcx_cellbody_corrected))
  print(individual_datapoints, n=1000)

  
  
  individual_datapoints <- clusterproperties %>%  subset(Orientation=="Bleb") %>% ungroup() %>% 
    dplyr::select(Length_label, Gcx_cellbody_corrected) %>% 
    filter(!is.nan(Gcx_cellbody_corrected))
  print(individual_datapoints, n=1000)
  
  
  
  individual_datapoints <- clusterproperties %>%  subset(Orientation=="Leading edge protrusion") %>% ungroup() %>% 
    dplyr::select(Length_label, Gcx_cellbody_corrected) %>% 
    filter(!is.nan(Gcx_cellbody_corrected))
  print(individual_datapoints, n=1000)
  
  boxplot_stats <- clusterproperties %>%
    ungroup() %>% 
    group_by(Orientation, Length_label) %>%
    arrange(Orientation) %>% 
    
    dplyr::summarise(
      n = n(),
      median = median(Gcx_cellbody_corrected, na.rm = TRUE),
      Q1 = quantile(Gcx_cellbody_corrected, 0.25, na.rm = TRUE),
      Q3 = quantile(Gcx_cellbody_corrected, 0.75, na.rm = TRUE),
      IQR = IQR(Gcx_cellbody_corrected, na.rm = TRUE),
      lower_whisker = max(min(Gcx_cellbody_corrected, na.rm = TRUE), Q1 - 1.5 * IQR),
      upper_whisker = min(max(Gcx_cellbody_corrected, na.rm = TRUE), Q3 + 1.5 * IQR)
    ) %>%
    ungroup() 
  
  print(boxplot_stats)
  
  
  model <- aov(B1int_cellbody_corrected ~ Length_label, data=subset(clusterproperties, Orientation=="Retraction fiber"))
  residuals <- model$residuals
  norm <- shapiro.test(residuals)
  pv <- norm$p.value
  print(pv)
  
  setforplot <- subset(clusterproperties, Orientation=="Bleb")
  dunn.test::dunn.test(x=setforplot$B1int_cellbody_corrected, g= setforplot$Length_label, method= "bonferroni")
  
  
  
  Gcx.cellbody.protrusion.base <- median(clusterproperties$Gcx_cellbody_corrected[clusterproperties$Length_label=="Base"])
  Gcx.cellbody.protrusion.tip <- median(clusterproperties$Gcx_cellbody_corrected[clusterproperties$Length_label=="Tip"])
  B1int.cellbody.protrusion.base <- median(clusterproperties$B1int_cellbody_corrected[clusterproperties$Length_label=="Base"])
  B1int.cellbody.protrusion.tip <- median(clusterproperties$B1int_cellbody_corrected[clusterproperties$Length_label=="Tip"])
  ratio.cellbody.protrusion.base <- median(clusterproperties$ratio_cellbody_corrected[clusterproperties$Length_label=="Base"])
  ratio.cellbody.protrusion.tip <- median(clusterproperties$ratio_cellbody_corrected[clusterproperties$Length_label=="Tip"])
  
  
  
  
  
   ## protrusions - plot of B1int tip versus base in clusters -----
  
  
  
  compare.data <- function(clusteproperties, subset){
    clusterproperties <- subset(clusterproperties, Orientation==subset)
    clusterproperties$Length_label <- droplevels(clusterproperties$Length_label)
    
    lengthlabel.vs.gcx <- custom_plot_with_statistics(data=clusterproperties, x_var="Length_label", group_vars="Orientation", alpha=0.5, size=2,
                                                      y_var="B1int_cellbody_corrected", jitter=FALSE, annotate_signif = TRUE, fill_var=NULL,
                                                      annotate_only_signif = TRUE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
    
    lengthlabel.vs.gcx.plot <- lengthlabel.vs.gcx$plot
    lengthlabel.vs.gcx.plot <- lengthlabel.vs.gcx.plot + theme_classic()+
      labs(y="Normalized \u03B21int", x="")
    print(lengthlabel.vs.gcx.plot)
    
    print(ggdraw(lengthlabel.vs.gcx.plot))
    
    return(lengthlabel.vs.gcx.plot)
    
  }
  
  leadingedge.b1int <- compare.data(clusterproperties, "Leading edge protrusion")+   scale_x_discrete(labels = c("Tip", "Shaft", "Base")) +theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  clusterproperties <- clusterproperties %>%  mutate(CellNR = as.numeric(CellNR),
                                                     expnr = ifelse(CellNR<100, 1, 2),
                                                     expnr = ifelse(CellNR>=200, 3, expnr)
  )
  
  leadingedge.gcx <- compare.data(clusterproperties, "Leading edge protrusion") +   scale_x_discrete(labels = c("Tip", "Shaft", "Base")) 
  
  
  leadingedge.stats <- subset(clusterproperties, Orientation=="Leading edge protrusion")
  
  
  shapiro.test(leadingedge.stats$B1int_cellbody_corrected[leadingedge.stats$Length_label=="Protrusion base"])
  
  leadingedge.stats <- leadingedge.stats %>%
    ungroup() %>% 
    group_by(expnr, Length_label) %>%
    arrange(Orientation) %>% 
    
    dplyr::summarise(
      n = n(),
      median = median(B1int_cellbody_corrected, na.rm = TRUE),
      Q1 = quantile(B1int_cellbody_corrected, 0.25, na.rm = TRUE),
      Q3 = quantile(B1int_cellbody_corrected, 0.75, na.rm = TRUE),
      IQR = IQR(B1int_cellbody_corrected, na.rm = TRUE),
      lower_whisker = max(min(B1int_cellbody_corrected, na.rm = TRUE), Q1 - 1.5 * IQR),
      upper_whisker = min(max(B1int_cellbody_corrected, na.rm = TRUE), Q3 + 1.5 * IQR)
    ) %>%
    ungroup() 
  
  print(leadingedge.stats)
  
  
  
  pairwise.results <- pairwise.wilcox.test(x=leadingedge.stats$Gcx_cellbody_corrected, g=leadingedge.stats$Length_label,p.adjust.method = "bonferroni")
  print(pairwise.results)
  leadingedge.stats %>%  aov(median ~ Length_label, data = .) %>%  summary()
  leadingedge.stats %>%   tukey_hsd(median ~ Length_label) 
  dunn.test::dunn.test(x=leadingedge.stats$median, g=leadingedge.stats$Length_label,
                       method = "bonferroni")
  
  
  
  
  
  
   bleb.b1int <- compare.data(clusterproperties, "Bleb")+   scale_x_discrete(labels = c("Peribleb zone", "Neck", "Bleb", "Apex")) +theme(axis.text.x = element_text(angle = 45, hjust=1))
  retraction.b1int <- compare.data(clusterproperties, "Retraction fiber")+   scale_x_discrete(labels = c("Base", "Shaft", "Tip")) +theme(axis.text.x = element_text(angle = 45, hjust=1))
  

  
  clusterproperties <- protrusion_data %>%  group_by(number, CellNR, Orientation, new_Peak_ID, Signal_type, Clustersize, unique_Peak_ID, Length_label) %>%  
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
  
  clusterproperties_bleb <- clusterproperties %>%  filter(Orientation=="Bleb" | Orientation=="Cellbody") %>% ungroup() %>%  
    group_by(CellNR) %>%  mutate(B1int_cellbody_average = B1int/ mean(B1int[Orientation=="Cellbody"])) %>%  filter(Orientation=="Bleb") %>%  mutate(Length_label = as.factor(Length_label))
  
  clusterproperties_bleb$Length_label <- factor(clusterproperties_bleb$Length_label, levels = c("Peribleb region", "Turning point", "Bleb", "Apex"))
  clusterproperties_bleb$Length_label <- droplevels(clusterproperties_bleb$Length_label)

  
  model <- aov(B1int_cellbody_average ~ Length_label, data=clusterproperties_bleb)
  residuals <- model$residuals
  norm <- shapiro.test(residuals)
  pv <- norm$p.value
  print(pv)
  
  
  dunn_results <- dunn.test::dunn.test(clusterproperties_bleb$B1int_cellbody_average, 
                                       g = clusterproperties_bleb$Length_label, 
                                       method = "bonferroni", altp = TRUE)
  
  # Extract p-values
  p_values <- dunn_results$altP.adjusted
  comparisons <- dunn_results$comparisons
  
  # Create a data frame with comparisons and p-values
  p_values_df <- data.frame(comparisons, p_values)
  
  
  # Custom groups and p-value
  custom_comparison <- list(c("Peribleb region", "Apex"))
  custom_p_value <- round(p_values_df[2,2],3)  # Replace with your p-value
  
 
  
  
  
  bleb.b1int <- ggplot(subset(clusterproperties_bleb ), aes(x=Length_label, y=B1int_cellbody_average))+
    geom_violin(scale="width")+
    geom_boxplot(width=0.5, outlier.shape=NA)+
    theme_classic()+
    theme(legend.position="none")+
    labs(y="Normalized  \u03B21int", x="")+
    scale_y_continuous(breaks=c(0,1,2,3,4,5), expand=c(0,0.4))+
    facet_wrap(~Orientation)
  
  
  
  
  individual_datapoints <- clusterproperties_bleb %>%  subset(Orientation=="Bleb") %>% ungroup() %>% 
    dplyr::select(Length_label, B1int_cellbody_average) %>% 
    filter(!is.nan(B1int_cellbody_average))
  print(individual_datapoints, n=1000)
  
  
  
  boxplot_stats <- clusterproperties_bleb %>%
    ungroup() %>% 
    group_by(Orientation, Length_label) %>%
    arrange(Orientation) %>% 
    
    dplyr::summarise(
      n = n(),
      median = median(B1int_cellbody_average, na.rm = TRUE),
      Q1 = quantile(B1int_cellbody_average, 0.25, na.rm = TRUE),
      Q3 = quantile(B1int_cellbody_average, 0.75, na.rm = TRUE),
      IQR = IQR(B1int_cellbody_average, na.rm = TRUE),
      lower_whisker = max(min(B1int_cellbody_average, na.rm = TRUE), Q1 - 1.5 * IQR),
      upper_whisker = min(max(B1int_cellbody_average, na.rm = TRUE), Q3 + 1.5 * IQR)
    ) %>%
    ungroup() 
  
  print(boxplot_stats)
  
  
  # Add significance annotation
  bleb.b1int <- bleb.b1int + geom_signif(
    comparisons = custom_comparison,
    annotations = custom_p_value,
    map_signif_level = FALSE,  # Use the exact p-value provided
  )
  
  print(bleb.b1int) 
  
  
  
  
  pdf("Violins.b1int.bleb.pdf", width=3, height=2)
  print(bleb.b1int)
  dev.off()
  
  
  
  clusterproperties <- protrusion_data %>%  group_by(number, CellNR, Orientation, new_Peak_ID, Signal_type, Clustersize, unique_Peak_ID, Length_label) %>%  
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))%>%  
    filter(Orientation=="Leading edge protrusion" | Orientation=="Retraction fiber") %>% 
    filter(Orientation=="Leading edge protrusion" & !is.na(Peaklocation_X) | 
             Orientation=="Retraction fiber" & !is.na(Peaklocation_X))
  
  noncluster_data <- protrusion_data %>%  group_by(number, CellNR, Orientation, Signal_type, Length_label) %>%  
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Orientation=="Cellbody" | Orientation=="Bleb")
  
  
  bleb_cluster_correction = FALSE
  
  if(bleb_cluster_correction==TRUE){
    noncluster_data <- noncluster_data %>%  
      filter(is.nan(Peaklocation_X))
  }
  
  clusterproperties <- rbind(clusterproperties, noncluster_data)
  
  clusterproperties <- clusterproperties %>%  mutate(unique_Peak_ID = ifelse(is.na(new_Peak_ID), NA, unique_Peak_ID))
  
  
  clusterproperties <- clusterproperties %>%  filter(Length_label!="Cell body")
  
  
  
   ## cluster distance versus B1int -----
  
  
  
  compare.data.lineplot <- function(clusteproperties, subset, label){
    clusterproperties <- subset(clusterproperties, Orientation==subset)
    
    
    lengthlabels <- clusterproperties %>% ungroup() %>% 
      filter(Signal_type=="foreground")  %>% 
      group_by(Length_label) %>% 
      summarize()
    if(unique(clusterproperties$Orientation)=="Leading edge protrusion"){
      lengthlabels <- lengthlabels %>%  ungroup() %>%  mutate(xmin=c(0,0.25,0.75),
                                                              xmax=c(0.25, 0.75,1),
                                                              ymin= rep(-Inf,3),
                                                              ymax= rep(Inf,3))
      
      lengthlabel.rectangles <- geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=rev(Length_label)), 
                                          color= alpha("black", 0.05), inherit.aes = FALSE, alpha = 0.2, size=0.5)
    }else{
      lengthlabels <- lengthlabels %>%  ungroup() %>%  mutate(xmin=c(0,0.2,0.6),
                                                              xmax=c(0.2, 0.6,1),
                                                              ymin= rep(-Inf,3),
                                                              ymax= rep(Inf,3))
      lengthlabel.rectangles <- geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=(Length_label)), 
                                          color= alpha("black", 0.05), inherit.aes = FALSE, alpha = 0.2, size=0.5)
    }
  
  distvsb1int <- custom_lineplot_with_rsquared(clusterproperties, regression_type =  "linear", x_var=X_relative, group_vars = "Orientation",
                                             y_var=B1int_cellbody_corrected, size=1, alpha=0.5, add.extra.plot=lengthlabel.rectangles)
  
  distvsb1int.plot <-distvsb1int$plot + geom_line(data=clusterproperties, aes(x=X_relative, y=B1int_cellbody_corrected, group=number, color=number)) + guides(color="none")+ scale_x_continuous(breaks=c(0,0.5,1))
  distvsb1int.plot <- distvsb1int.plot + theme_classic() + theme(legend.position = "none")
    labs(x="Relative distance", y="Normalized \u03B21int")+ylim(c(0, max(clusterproperties$B1int_cellbody_corrected)))
  print(distvsb1int.plot)
  
  
  distvsb1int.plot <-  plot_grid(distvsb1int.plot, labels=label, vjust=0)+
    theme(plot.margin = margin(t = 5, unit = "mm"))
  
  print(ggdraw(distvsb1int.plot))
  
  return(distvsb1int.plot)
  }
  
  lineplot.b1int.leadingedge <- compare.data.lineplot(clusterproperties, "Leading edge protrusion", "a")
  lineplot.b1int.retraction <- compare.data.lineplot(clusterproperties, "Retraction fiber", "b")
  
  lineplots.b1int <- plot_grid(lineplot.b1int.leadingedge, lineplot.b1int.retraction)
  
  print(lineplots.b1int)
  
  pdf("Lineplots.b1int.pdf", width=5, height=2)
  print(lineplots.b1int)
  dev.off()
  
  individual.datapoints <- clusterproperties %>%  filter(Orientation=="Leading edge protrusion") %>% 
    ungroup() %>%  dplyr::select(X_relative, B1int_cellbody_corrected, Length_label) %>%  arrange(X_relative) 
  print(individual.datapoints, n=200)
  
  
  individual.datapoints <- clusterproperties %>%  filter(Orientation=="Retraction fiber") %>% 
    ungroup() %>%  dplyr::select(X_relative, B1int_cellbody_corrected, Length_label) %>%  arrange(X_relative) 
  print(individual.datapoints, n=200)
  
  
   ## cluster distance versus Gcx -----
  
  
  
  compare.data.lineplot <- function(clusteproperties, subset, label){
    clusterproperties <- subset(clusterproperties, Orientation==subset)
    
    
    lengthlabels <- clusterproperties %>% ungroup() %>% 
      filter(Signal_type=="foreground")  %>% 
      group_by(Length_label) %>% 
      summarize()
    if(unique(clusterproperties$Orientation)=="Leading edge protrusion"){
    lengthlabels <- lengthlabels %>%  ungroup() %>%  mutate(xmin=c(0,0.25,0.75),
                                             xmax=c(0.25, 0.75,1),
                                             ymin= rep(-Inf,3),
                                             ymax= rep(Inf,3))
    
    lengthlabel.rectangles <- geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=rev(Length_label)), 
                                        color= alpha("black", 0.05), inherit.aes = FALSE, alpha = 0.2, size=0.5)
    }else{
      lengthlabels <- lengthlabels %>%  ungroup() %>%  mutate(xmin=c(0,0.2,0.6),
                                                              xmax=c(0.2, 0.6,1),
                                                              ymin= rep(-Inf,3),
                                                              ymax= rep(Inf,3))
      lengthlabel.rectangles <- geom_rect(data = lengthlabels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=(Length_label)), 
                                          color= alpha("black", 0.05), inherit.aes = FALSE, alpha = 0.2, size=0.5)
    }
    
    


    distvsb1int <- custom_lineplot_with_rsquared(clusterproperties, regression_type =  "linear", x_var=X_relative, group_vars = "Orientation",
                                                 y_var=Gcx_cellbody_corrected, size=1, alpha=0.5, add.extra.plot=lengthlabel.rectangles)
    
    distvsb1int.plot <-distvsb1int$plot + geom_line(aes(x=X_relative, y=Gcx_cellbody_corrected, group=number, color=number)) + guides(color="none") + scale_x_continuous(breaks=c(0,0.5,1))
    distvsb1int.plot <- distvsb1int.plot + theme_classic() + theme(legend.position="none")+
      labs(x="Relative distance", y="Normalized Gcx")+ylim(c(0, max(clusterproperties$Gcx_cellbody_corrected)))
    print(distvsb1int.plot)
    
    
    distvsb1int.plot <-  plot_grid(distvsb1int.plot, labels=label, vjust=0)+
      theme(plot.margin = margin(t = 5, unit = "mm"))
    
    print(ggdraw(distvsb1int.plot))
    
    return(distvsb1int.plot)
  }
  
  lineplot.gcx.leadingedge <- compare.data.lineplot(individual_datapoints, "Leading edge protrusion", "a")
  lineplot.gcx.retraction <- compare.data.lineplot(clusterproperties, "Retraction fiber", "b")
  
  individual_datapoints <- clusterproperties %>%  filter(Orientation=="Leading edge protrusion") %>% ungroup() %>% 
    dplyr::select(X_relative, Gcx_cellbody_corrected)
  print(individual_datapoints, n=1000)
  
  
  individual_datapoints <- clusterproperties %>%  filter(Orientation=="Retraction fiber") %>% ungroup() %>% 
    dplyr::select(X_relative, Gcx_cellbody_corrected)
  print(individual_datapoints, n=1000)
  
  lineplots.gcx <- plot_grid(lineplot.gcx.leadingedge, lineplot.gcx.retraction)
  
  print(lineplots.gcx)
  
  pdf("Lineplots.gcx.pdf", width=5, height=2)
  print(lineplots.gcx)
  dev.off()
  
   ## blebs - scatter plot B1int vs Gcx at apex ----
  
  bleb_data <- subset(clusterproperties, Orientation=="Bleb")
  
  bleb_data <- subset(clusterproperties, Length_label=="Apex" | Length_label=="Bleb")
  
  
  breaks_log10 <- function(x) {
    library(scales)
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    
    10^(seq.int(low, high))
  }
  
  paired.blebvsapex <- ggplot(bleb_data, aes(x=B1int, y=Gcx, color=Length_label, shape=Length_label))+
    geom_point(size=2)+
    geom_line(aes(group=number), color="black")+
    annotation_logticks()+
    scale_x_log10(breaks = breaks_log10,
                  limits = c(10,10000),
                  labels = trans_format(log10, math_format(10^.x))) +
    scale_y_log10(breaks = breaks_log10,
                  limits = c(100,10000),
                  labels = trans_format(log10, math_format(10^.x)))+

    theme_classic()+
    theme(legend.position=c(0.2,0.8))+
    theme(legend.title = element_blank())+
    labs(x="\u03B21int", y="Gcx")
  
  print(paired.blebvsapex)
  
  individual.datapoints <- bleb_data %>%  ungroup() %>% 
    dplyr::select(B1int, Gcx, Length_label, number) %>%  arrange(number) %>%  
  print(individual.datapoints, n=200)
  
  bleb_data<- bleb_data %>%  ungroup() %>% dplyr::select(Orientation, CellNR, number, Length_label, X, B1int, Gcx)
  write.csv2(bleb_data, "Bleb_data_paired.csv")
  
  pdf("Blebs.paired.apexbleb.b1intvsgcx.pdf", width=3, height=3)
  print(paired.blebvsapex)
  dev.off()
  
  
  bleb_data <- bleb_data %>% ungroup() %>% group_by(number) %>% 
    mutate( Gcx_bleb_corrected=Gcx/ Gcx[Length_label=="Bleb"],
           Gcx_bleb_enrichment = ifelse(Gcx_bleb_corrected >1, TRUE, FALSE),
           B1int_bleb_corrected =  B1int/ B1int[Length_label=="Bleb"],
           B1int_bleb_enrichment = ifelse(B1int_bleb_corrected>1, TRUE, FALSE)
    ) %>%  filter(Length_label=="Apex")
  
  write.csv2(bleb_data, "Bleb_data_enrichment.csv")
  
  
  bleb_scatterplot <- custom_lineplot_with_rsquared(bleb_data, x_var=B1int_bleb_corrected, y_var=Gcx_bleb_corrected, 
                                                    group_vars="Orientation")
  

  
  
  
  bleb_scatterplot <- bleb_scatterplot$plot
  print(bleb_scatterplot)
  
  print(paste0("Percentage blebs with co-depletion of B1int and GCx: ",length(bleb_data$B1int_bleb_enrichment[bleb_data$B1int_bleb_enrichment==FALSE])/length(bleb_data$B1int_bleb_enrichment)*100))
  
## submicron segregation data -----
   ## outward clusters vs. paired lateral clusters-----
  outwardcluster_data <- subset(protrusion_data, Orientation=="Outwardclusters" | Orientation=="Innerclusters" | Orientation=="Cellbody")
  outwardcluster_data <- outwardcluster_data %>%  ungroup() %>%  group_by(CellNR, RoiNR, number, grouped_number, Signal_type) %>% 
    mutate(B1int_peak = ifelse(Orientation=="Outwardclusters", max(B1int), B1int[round(which.max(X)/2)]),
           Gcx_peak = ifelse(Orientation=="Outwardclusters", Gcx[which.max(B1int)], Gcx[round(which.max(X)/2)]),
           B1int_peak_location = ifelse(Orientation=="Outwardclusters", X[which.max(B1int)], X[round(which.max(X)/2)]),
           Gcx_peak_location = X[which.max(Gcx)]
  ) %>% 
    ungroup() %>%  group_by(CellNR, RoiNR, grouped_number, Signal_type) %>% 
    mutate(
           B1int_edge = ifelse(Orientation=="Innerclusters", B1int_cluster_local_bg[round(which.max(X[Orientation=="Innerclusters"])/2)], NA),
           B1int_edge = mean(B1int_edge, na.rm=TRUE),
           Gcx_edge = ifelse(Orientation=="Innerclusters", Gcx_cluster_local_bg[round(which.max(X[Orientation=="Innerclusters"])/2)], NA),
           Gcx_edge = mean(B1int_edge, na.rm=TRUE),
           Gcx_inner_peak = ifelse(Orientation!="Cellbody", Gcx[round(which.max(X[Orientation=="Innerclusters"])/2)], Gcx[which.max(Peaklocation_X)]),
           B1int_enrichment = B1int_peak/ B1int_edge,
           B1int_enrichment = ifelse(Orientation=="Cellbody", B1int_local_enrichment, B1int_enrichment),
           Gcx_enrichment = Gcx_peak/ Gcx_edge,
           Gcx_enrichment = ifelse(Orientation=="Cellbody", Gcx_local_enrichment, Gcx_enrichment)
             ) %>% 
    ungroup() %>% group_by(CellNR, RoiNR, number, grouped_number, Signal_type) %>% 
    mutate(new_Peak_ID = as.numeric(new_Peak_ID),
           Location_within_cluster = ifelse(Orientation=="Outwardclusters" & B1int==B1int_peak, "Peak", NA),
           Location_within_cluster = ifelse(Orientation=="Innerclusters" & B1int==B1int_peak, "Peak", Location_within_cluster),
           Location_within_cluster = ifelse(Orientation=="Innerclusters" & Peak_position=="Edge" & new_Peak_ID == mean(new_Peak_ID[Location_within_cluster=="Peak"], na.rm=TRUE), 
                                            Peak_position, Location_within_cluster),
           Location_within_cluster = ifelse(Orientation=="Outwardclusters" & X==X[which.max(Gcx)], "Edge", Location_within_cluster),
           Location_within_cluster = ifelse(Orientation=="Cellbody", Peak_position, Location_within_cluster),
           Segregation_distance = ifelse(Orientation=="Outwardclusters", ((Gcx_peak_location - B1int_peak_location))*pxsize, NA),
           min_Edge_X = ifelse(Orientation=="Innerclusters" | Orientation=="Cellbody", min(X[Location_within_cluster=="Edge" & new_Peak_ID == mean(new_Peak_ID[Location_within_cluster=="Peak"],na.rm=TRUE)  ]), NA),
           mid_manual = X[round(which.max(X)/2)],
           Segregation_distance = ifelse(Orientation=="Innerclusters", mid_manual- min_Edge_X*pxsize, Segregation_distance),
           new_Peak_ID = as.character(new_Peak_ID)
    ) %>%  dplyr::select(Orientation, CellNR, RoiNR, number, grouped_number, Signal_type, B1int, X, X_micron, min_Edge_X, mid_manual,
                         Gcx, B1int_edge, B1int_peak, Gcx_peak, Gcx_edge, Gcx_inner_peak, Gcx_enrichment, B1int_enrichment, Location_within_cluster, new_Peak_ID, Segregation_distance)

  
  
  clusters_summarized_positions <- outwardcluster_data %>%  group_by(CellNR, RoiNR, grouped_number, Orientation, Signal_type, Location_within_cluster) %>% 
    dplyr::filter(!is.na(Location_within_cluster)) %>% 
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  ungroup() %>% 
    group_by(CellNR, RoiNR, grouped_number, Orientation, Signal_type) %>% 
    mutate(Gcx_enrichment = mean( Gcx[Location_within_cluster=="Peak"])/ mean(Gcx[Location_within_cluster=="Edge"]),
           B1int_enrichment = mean( B1int[Location_within_cluster=="Peak"])/ mean(B1int[Location_within_cluster=="Edge"])
           )%>% ungroup() %>%  
    mutate(
           B1int_enrichment = B1int_enrichment / mean(B1int_enrichment[Orientation=="Cellbody"])
           )
  
  
  
  
  breaks_log10 <- function(x) {
    library(scales)
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    
    10^(seq.int(low, high))
  }
  library(scales)
  
  paired <- ggplot(clusters_summarized_positions, aes(x=Gcx, y= B1int))+
    geom_point(aes(shape=Location_within_cluster, color=as.factor(grouped_number)))+
    geom_line(aes(group=grouped_number, color=as.factor(grouped_number)))+
    facet_wrap(.~Orientation)+
    annotation_logticks()+
    scale_x_log10(breaks = breaks_log10,
                  limits = c(100,10000),
                  labels = trans_format(log10, math_format(10^.x))) +
    scale_y_log10(breaks = breaks_log10,
                  limits = c(300,10000),
                  labels = trans_format(log10, math_format(10^.x)))+
    theme_classic()+
    theme(legend.position="top")+
    scale_color_viridis_d(option="turbo", begin=0.1, end=0.9)+
    guides(color=FALSE)
  
  print(paired)
    
  clusters_summarized <-clusters_summarized_positions %>% dplyr:: filter(Location_within_cluster=="Edge") %>%  group_by(grouped_number) %>% 
    filter(!(any(is.na(Gcx_enrichment))))
  
  
  mean(clusters_summarized$B1int_enrichment[clusters_summarized$Orientation=="Cellbody"], na.rm=TRUE)
  mean(clusters_summarized$B1int_enrichment[clusters_summarized$Orientation=="Innerclusters"], na.rm=TRUE)
  mean(clusters_summarized$B1int_enrichment[clusters_summarized$Orientation=="Outwardclusters"], na.rm=TRUE)
  
  

  
  ## compare Gcx enrichment in membrane subcompartments
  
  enrichment <- ggplot(subset(clusters_summarized), aes(x=Orientation, y=Gcx_enrichment))+
    geom_violin(scale="width", aes(fill=Orientation))+
    geom_boxplot(width=0.5, aes(fill=Orientation), outlier.shape=NA)+
    #geom_point(alpha=0.3)+
    #geom_line(data=subset(clusters_summarized,Orientation!="Cellbody"), aes(group=grouped_number), alpha=0.3)+
     theme_classic()+
    geom_hline(yintercept=1, linetype="dashed")+
    theme(legend.position = "none")+
    labs(x="", y="Gcx enrichment")+
    scale_fill_manual(values=c("gray", "cyan", "white"))
  print(enrichment)
  
  
  boxplot_stats <- clusters_summarized %>%
    ungroup() %>% 
    group_by(Orientation) %>%
    
    dplyr::summarise(
      n = n(),
      median = median(Gcx_enrichment, na.rm = TRUE),
      Q1 = quantile(Gcx_enrichment, 0.25, na.rm = TRUE),
      Q3 = quantile(Gcx_enrichment, 0.75, na.rm = TRUE),
      IQR = IQR(Gcx_enrichment, na.rm = TRUE),
      lower_whisker = max(min(Gcx_enrichment, na.rm = TRUE), Q1 - 1.5 * IQR),
      upper_whisker = min(max(Gcx_enrichment, na.rm = TRUE), Q3 + 1.5 * IQR)
    ) %>%
    ungroup()
  
  print(boxplot_stats)
  
  clusters_summarized$Orientation <- droplevels(clusters_summarized$Orientation)
  stattest <- clusters_summarized %>%
    kruskal.test(Gcx_enrichment ~ Orientation, data = .)
  
  
    
  pairwise.results <- pairwise.wilcox.test(x=clusters_summarized$Gcx_enrichment, g=clusters_summarized$Orientation,
      p.adjust.method = "bonferroni"
    )
  
  cbvso <- pairwise.results$p.value[2,1]
  ivso <- pairwise.results$p.value[2,2]
  
  annotations <- geom_signif(comparisons = list(c("Cellbody", "Outwardclusters")), 
                             annotations = paste("p =", format(cbvso, scientific = TRUE, digits = 2)), 
                             textsize = 4, color="black", margin_top=0.15)
  
  
  annotations2 <- geom_signif(comparisons = list(c("Outwardclusters", "Innerclusters")), 
                             annotations = paste("p =", format(ivso, scientific = TRUE, digits = 2)), 
                             textsize = 4, color="black", margin_top=0.05)
  
  enrichment <- enrichment + annotations + annotations2 + ylim(c(FALSE, 1.8))
  print(enrichment)
  

  
  library(rstatix)
  clusters_summarized %>%  ungroup() %>% rstatix::kruskal_test(Gcx_enrichment ~ Orientation) 
  clusters_summarized %>%  ungroup() %>% rstatix::kruskal_effsize(Gcx_enrichment ~ Orientation) 
  

  
  pdf("violins.gcx.enrichment.pdf", width=2, height=1.5)
  print(enrichment)
  dev.off()
  
  
  ## compare b1int enrichment in submicron compartments
  enrichment <- ggplot(subset(clusters_summarized), aes(x=Orientation, y=B1int_enrichment))+
    geom_violin(scale="width", aes(fill=Orientation))+
    geom_boxplot(width=0.5, aes(fill=Orientation), outlier.shape=NA)+
    #geom_point(alpha=0.3)+
    #geom_line(data=subset(clusters_summarized,Orientation!="Cellbody"), aes(group=grouped_number), alpha=0.3)+
    theme_classic()+
    geom_hline(yintercept=1, linetype="dashed")+
    theme(legend.position = "none")+
    labs(x="", y="B1int enrichment")+
    scale_fill_manual(values=c("gray", "cyan", "white"))
  print(enrichment)
  
  
  
  boxplot_stats <- clusters_summarized %>%
    ungroup() %>% 
    group_by(Orientation) %>%
    
    dplyr::summarise(
      n = n(),
      median = median(B1int_enrichment, na.rm = TRUE),
      Q1 = quantile(B1int_enrichment, 0.25, na.rm = TRUE),
      Q3 = quantile(B1int_enrichment, 0.75, na.rm = TRUE),
      IQR = IQR(B1int_enrichment, na.rm = TRUE),
      lower_whisker = max(min(B1int_enrichment, na.rm = TRUE), Q1 - 1.5 * IQR),
      upper_whisker = min(max(B1int_enrichment, na.rm = TRUE), Q3 + 1.5 * IQR)
    ) %>%
    ungroup()
  
  print(boxplot_stats)
  
  
  clusters_summarized$Orientation <- droplevels(clusters_summarized$Orientation)
  stattest <- clusters_summarized %>%
    kruskal.test(Gcx_enrichment ~ Orientation, data = .)
  
  pairwise.results <- pairwise.wilcox.test(x=clusters_summarized$B1int_enrichment, g=clusters_summarized$Orientation,
                                           p.adjust.method = "bonferroni"
  )
  
  library(rstatix)
  clusters_summarized %>%  ungroup() %>% rstatix::kruskal_test(B1int ~ Orientation) 
  clusters_summarized %>%  ungroup() %>% rstatix::kruskal_effsize(B1int ~ Orientation) 
  
  
  cbvso <- pairwise.results$p.value[2,1]
  ivso <- pairwise.results$p.value[2,2]
  cbvsi <- pairwise.results$p.value[1,1]
  
  annotations <- geom_signif(comparisons = list(c("Cellbody", "Outwardclusters")), 
                             annotations = paste("p =", format(cbvso, scientific = TRUE, digits = 2)), 
                             textsize = 4, color="black", margin_top=0.3)
  
  
  annotations2 <- geom_signif(comparisons = list(c("Outwardclusters", "Innerclusters")), 
                              annotations = paste("p =", format(ivso, scientific = TRUE, digits = 2)), 
                              textsize = 4, color="black", margin_top=0.05)
  
  annotations3 <- geom_signif(comparisons = list(c("Cellbody", "Innerclusters")), 
                              annotations = paste("p =", format(cbvsi, scientific = TRUE, digits = 2)), 
                              textsize = 4, color="black", margin_top=0.15)
  
  enrichment <- enrichment + annotations + annotations2 + annotations3 + ylim(c(FALSE, 9))
  print(enrichment)
  
  
  pdf("violins.b1int.enrichment.pdf", width=2, height=1.5)
  print(enrichment)
  dev.off()
  
  
  
  
  
  
  
  
  clusters_summarized <-clusters_summarized_positions %>% dplyr:: filter(Location_within_cluster=="Edge") %>%  group_by(grouped_number) %>% 
    filter(!(any(is.na(Gcx_enrichment)))) %>% 
    filter(any(Orientation == "Cellbody") & any(Orientation == "Innerclusters"))
  
  clusters_summarized <-clusters_summarized %>%  dplyr::filter(Orientation!="Outwardclusters")
  
  
  enrichment <- ggplot(subset(clusters_summarized), aes(x=Orientation, y=B1int_enrichment))+
    geom_violin(scale="width", aes(fill=Orientation))+
    geom_boxplot(width=0.5, aes(fill=Orientation), outlier.shape=NULL)+
    geom_point(alpha=0.3)+
    geom_line(aes(group=grouped_number), alpha=0.5)+
    theme_classic()+
    geom_hline(yintercept=1, linetype="dashed")+
    theme(legend.position = "none")+
    labs(x="", y="B1int enrichment")+
    scale_fill_manual(values=c("gray", "cyan"))
  print(enrichment)
  
  
  stattest <- stats::t.test(B1int_enrichment ~ Orientation, data = clusters_summarized, paired = TRUE)
  p_value <- stattest$p.value
  
  
  
  annotations <- geom_signif(comparisons = list(c("Cellbody", "Innerclusters")), 
                             annotations = paste("p =", format(p_value, scientific = FALSE, digits = 2)), 
                             textsize = 4, color="black", margin_top=0.1)
  
  enrichment <- enrichment + annotations
  print(enrichment)
  
  
  pdf("Cellbody.vs.Innerclusters.integrinenrichment.pdf", width=5, height=3)
  print(enrichment)
  dev.off()
  
  p <- subset(clusters_summarized_positions, Orientation=="Outwardclusters"| Orientation=="Innerclusters")
  p <- p %>%  group_by(grouped_number) %>%  mutate(Segregation_distance = ifelse(is.na(Segregation_distance), Segregation_distance[Orientation=="Outwardclusters"], Segregation_distance))
  
ggplot(subset(p), aes(x=Segregation_distance, y= B1int))+
  geom_point(aes(shape=Location_within_cluster))+
  geom_line(aes(group=grouped_number))+
  theme_classic()+
  facet_wrap(.~Orientation)
  
  
   ## outward clusters - calculate segregation distance -------
  

  
  outwardcluster_data <- subset(protrusion_data, Orientation=="Outwardclusters" | Orientation=="Perpendicular control")
  outwardcluster_data <- outwardcluster_data %>%  ungroup() %>%  group_by(number, Orientation, Signal_type) %>% 
    mutate(B1int_peak = max(B1int),
           Gcx_peak = max(Gcx),
           B1int_peak_location = X[which.max(B1int)],
           Gcx_peak_location = X[which.max(Gcx)],
           Segregation_distance = ((Gcx_peak_location - B1int_peak_location))*pxsize
    ) 
  
  
  
  unpaired.segregationdistance.data <- outwardcluster_data %>%  group_by(number, CellNR, Signal_type, Orientation) %>%  
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  
    dplyr::select(c(number, CellNR, Orientation, Signal_type, B1int, Gcx, Col, peak_b1int_threshold, B1int_peak, Gcx_peak, B1int_peak_location, Gcx_peak_location, Segregation_distance)) 
  
  unpaired.segregationdistance.data$Orientation <- droplevels(unpaired.segregationdistance.data$Orientation)
  unpaired.segregationdistance.data$Orientation <- factor(unpaired.segregationdistance.data$Orientation, levels = c("Perpendicular control", "Outwardclusters"))
  
  violins.segregation.distance <- ggplot(data=unpaired.segregationdistance.data, aes(x=Orientation, y=Segregation_distance))+
    geom_violin(scale="width", width=1)+theme_classic()+
    geom_boxplot(width=0.5, outlier.shape=NA)+
    #geom_jitter(width=0.25)+
    
    labs(x="Integrin type" ,y ="Distance (\u03BCm)")+
    scale_x_discrete(labels = c("Non-focal", "Focal"))+
    scale_y_continuous(expand=(c(0,0.2)))
  
  
  boxplot_stats <- unpaired.segregationdistance.data %>%
    ungroup() %>% 
    group_by(Orientation) %>%
    
    dplyr::summarise(
      n = n(),
      median = median(Segregation_distance, na.rm = TRUE),
      Q1 = quantile(Segregation_distance, 0.25, na.rm = TRUE),
      Q3 = quantile(Segregation_distance, 0.75, na.rm = TRUE),
      IQR = IQR(Segregation_distance, na.rm = TRUE),
      lower_whisker = max(min(Segregation_distance, na.rm = TRUE), Q1 - 1.5 * IQR),
      upper_whisker = min(max(Segregation_distance, na.rm = TRUE), Q3 + 1.5 * IQR)
    ) %>%
    ungroup()
  
  boxplot_stats
  

  stattest <- wilcox.test(unpaired.segregationdistance.data$Segregation_distance[unpaired.segregationdistance.data$Orientation=="Perpendicular control"], 
                          unpaired.segregationdistance.data$Segregation_distance[unpaired.segregationdistance.data$Orientation=="Outwardclusters"],
                          exact=FALSE
                          )
  p_value <- stattest$p.value
  
  annotations <- geom_signif(comparisons = list(c("Outwardclusters", "Perpendicular control")), 
                             annotations = paste("p =", format(p_value, scientific = TRUE, digits = 2)), 
                             textsize = 4, color="black", margin_top=0.1)
  
  violins.segregation.distance <- violins.segregation.distance + annotations
  
  print(violins.segregation.distance)
  
  pdf("violins.segregation.distance.pdf", width=2, height=1)
  print(violins.segregation.distance)
  dev.off()
  
  
  unpaired.segregationdistance.data %>%  ungroup() %>% rstatix::kruskal_effsize(Segregation_distance ~ Orientation) 
  
  
  
   ## stratify for Gcx intensity: plot  - Gcx at edge vs. integrin enrichment------

  remove_outliers <- function(df, column) {
    Q1 <- quantile(df[[column]], 0.25)
    Q3 <- quantile(df[[column]], 0.75)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    df %>%
      filter(df[[column]] >= lower_bound & df[[column]] <= upper_bound)
  }
  
  
  clusters_summarized <-clusters_summarized_positions %>% dplyr:: filter(Location_within_cluster=="Edge") %>%  group_by(grouped_number) %>% 
    filter(!(any(is.na(Gcx_enrichment)))) %>% filter(Orientation=="Outwardclusters")
  
  remove_data_outliers=FALSE
  if(remove_data_outliers==TRUE){
  clusters_summarized <- clusters_summarized %>% ungroup () %>%   remove_outliers("B1int_enrichment")
  }
  
  # Sample data
  data <- clusters_summarized$B1int_enrichment
  
  # Calculate Q1, Q3, and IQR
  Q1 <- quantile(data, 0.25)
  Q3 <- quantile(data, 0.75)
  IQR <- Q3 - Q1
  
  # Determine outliers
  outliers <- data[data < (Q1 - 1.5 * IQR) | data > (Q3 + 1.5 * IQR)]
  print(outliers)
  
  quantiles <- quantile(clusters_summarized$Gcx_inner_peak, probs = c(0, 1/3, 2/3, 1))
  
  clusters_summarized$Gcx_content_classification_peak <- cut(clusters_summarized$Gcx_inner_peak, breaks = quantiles, labels = c("L", "M", "H"), include.lowest = TRUE)
  
  
  Gcx_content_threshold_peak_low =  max(clusters_summarized$Gcx_inner_peak[clusters_summarized$Gcx_content_classification_peak=="L"])
  Gcx_content_threshold_peak_high=  max(clusters_summarized$Gcx_inner_peak[clusters_summarized$Gcx_content_classification_peak=="M"])
  
  
  clusters_summarized$Gcx_content_classification_peak <- factor(clusters_summarized$Gcx_content_classification_peak, levels=c("L", "M", "H"))
  
  # Create a data frame with rectangle information
  rect_data_edge <- data.frame(
    xmin = c(0, Gcx_content_threshold_peak_low, Gcx_content_threshold_peak_high),
    xmax = c(Gcx_content_threshold_peak_low, Gcx_content_threshold_peak_high, max(clusters_summarized$Gcx_inner_peak, na.rm=TRUE)),
    ymin = rep(-Inf, 3),
    ymax = rep(Inf, 3),
    fill = c( "#00BA38","#1661a8","#bd3333")
  )
  
  gg.rect <- geom_rect(data = rect_data_edge, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                       inherit.aes = FALSE, alpha=0.5, size=0.5)
  
  clusters_summarized$CellNR <- as.factor(clusters_summarized$CellNR)
  
  
  Gcx.vs.B1intenrichment <- custom_lineplot_with_rsquared(clusters_summarized, regression_type =  "log", x_var=Gcx_inner_peak, y_var=B1int_enrichment, shape_var=NULL,
                                                          size=2, alpha=0.5, add.extra.plot=NULL, annotate_signif=FALSE, group_vars="Orientation")
  Gcx.vs.B1intenrichment.graph <- Gcx.vs.B1intenrichment$plot+
    theme_classic()+
    scale_x_continuous(n.breaks=3)+
    theme(legend.position="none")+
    labs(x="Local Gcx (MFI)", y="B1int enrichment")+
    theme(axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black"))
  print(Gcx.vs.B1intenrichment.graph)
  
  
  boxplot_stats <- clusters_summarized %>% ungroup() %>% 
    dplyr::select(Gcx_inner_peak, B1int_enrichment) %>%
    ungroup()
  
  print(boxplot_stats, n=100)
  
  
  
  pdf("Outerclusters.B1intenrichmentvsgcx.pdf", width=2, height=2)
  Gcx.vs.B1intenrichment.graph
  dev.off()
  
  
  custom_lineplot_with_rsquared(clusters_summarized, regression_type =  "log", x_var=Gcx_inner_peak, y_var=B1int_enrichment, shape_var=NULL,
                                size=2, alpha=0.5, add.extra.plot=NULL, annotate_signif=FALSE, group_vars="Orientation")
  
  
  violin.B1int.enrichment <-  custom_plot_with_statistics(data=clusters_summarized, x_var="Gcx_content_classification_peak", alpha=0.5, size=1,
                                                          y_var="B1int_enrichment", group_vars="Orientation", jitter=TRUE, annotate_signif = TRUE, fill_var="Gcx_content_classification_peak",
                                                          annotate_only_signif = FALSE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
  violin.B1int.enrichment.stats <- violin.B1int.enrichment$statistical_test
  violin.B1int.enrichment <- violin.B1int.enrichment$plot+
    theme_classic()+
    labs(x="Inner peak Gcx (MFI)", y="B1int enrichment")+
    theme(legend.position="none")
  
  print(violin.B1int.enrichment)
  
  
  
  B1int.enrichment.together <- plot_grid(Gcx.vs.B1intenrichment.graph, violin.B1int.enrichment, nrow=1, align="v", axis="b")
  print(B1int.enrichment.together)
  
  pdf("Perpendicular.data.pdf", width=5, height=2)
  B1int.enrichment.together
  dev.off()
  
  violin.B1int.enrichment.stats <- violin.B1int.enrichment.stats %>%  dplyr::select(-groups)
  
  write.csv2(outwardcluster_data, "Perpendicular.analysis.csv")
  write.csv2(as.data.frame(violin.B1int.enrichment.stats), "Perpendicular.analyis.glycocalyx.vs.b1intenrichment.csv")
  

  
  
   ## stratify for Gcx intensity: plot  - Gcx at inner clusters vs. integrin enrichment------
  
  remove_outliers <- function(df, column) {
    Q1 <- quantile(df[[column]], 0.25)
    Q3 <- quantile(df[[column]], 0.75)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    df %>%
      filter(df[[column]] >= lower_bound & df[[column]] <= upper_bound)
  }
  
  
  clusters_summarized <-clusters_summarized_positions %>% dplyr:: filter(Location_within_cluster=="Edge") %>%  group_by(grouped_number) %>% 
    filter(!(any(is.na(Gcx_enrichment)))) %>% filter(Orientation=="Innerclusters")
  
  remove_data_outliers=FALSE
  if(remove_data_outliers==TRUE){
    clusters_summarized <- clusters_summarized %>% ungroup () %>%   remove_outliers("B1int_enrichment")
  }
  
  # Sample data
  data <- clusters_summarized$B1int_enrichment
  
  # Calculate Q1, Q3, and IQR
  Q1 <- quantile(data, 0.25)
  Q3 <- quantile(data, 0.75)
  IQR <- Q3 - Q1
  
  # Determine outliers
  outliers <- data[data < (Q1 - 1.5 * IQR) | data > (Q3 + 1.5 * IQR)]
  print(outliers)
  
  quantiles <- quantile(clusters_summarized$Gcx_inner_peak, probs = c(0, 1/3, 2/3, 1))
  
  clusters_summarized$Gcx_content_classification_peak <- cut(clusters_summarized$Gcx_inner_peak, breaks = quantiles, labels = c("L", "M", "H"), include.lowest = TRUE)
  
  
  Gcx_content_threshold_peak_low =  max(clusters_summarized$Gcx_inner_peak[clusters_summarized$Gcx_content_classification_peak=="L"])
  Gcx_content_threshold_peak_high=  max(clusters_summarized$Gcx_inner_peak[clusters_summarized$Gcx_content_classification_peak=="M"])
  
  
  clusters_summarized$Gcx_content_classification_peak <- factor(clusters_summarized$Gcx_content_classification_peak, levels=c("L", "M", "H"))
  
  # Create a data frame with rectangle information
  rect_data_edge <- data.frame(
    xmin = c(0, Gcx_content_threshold_peak_low, Gcx_content_threshold_peak_high),
    xmax = c(Gcx_content_threshold_peak_low, Gcx_content_threshold_peak_high, max(clusters_summarized$Gcx_inner_peak, na.rm=TRUE)),
    ymin = rep(-Inf, 3),
    ymax = rep(Inf, 3),
    fill = c( "#00BA38","#1661a8","#bd3333")
  )
  
  gg.rect <- geom_rect(data = rect_data_edge, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                       inherit.aes = FALSE, alpha=0.5, size=0.5)
  
  clusters_summarized$CellNR <- as.factor(clusters_summarized$CellNR)
  
  
  Gcx.vs.B1intenrichment <- custom_lineplot_with_rsquared(clusters_summarized, regression_type =  "log", x_var=Gcx_inner_peak, y_var=B1int_enrichment, shape_var=NULL,
                                                          size=2, alpha=0.5, add.extra.plot=NULL, annotate_signif=FALSE, group_vars="Orientation")
  Gcx.vs.B1intenrichment.graph <- Gcx.vs.B1intenrichment$plot+
    theme_classic()+
    scale_x_continuous(n.breaks=3)+
    theme(legend.position="none")+
    labs(x="Local Gcx (MFI)", y="B1int enrichment")+
    theme(axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black"))
  print(Gcx.vs.B1intenrichment.graph)
  
  
  boxplot_stats <- clusters_summarized %>% ungroup() %>% 
    dplyr::select(Gcx_inner_peak, B1int_enrichment) %>%
    ungroup()
  
  print(boxplot_stats, n=100)
  
  
  
  pdf("Innerclusters.B1intenrichmentvsgcx.pdf", width=2, height=2)
  Gcx.vs.B1intenrichment.graph
  dev.off()
  
  
  custom_lineplot_with_rsquared(clusters_summarized, regression_type =  "log", x_var=Gcx_inner_peak, y_var=B1int_enrichment, shape_var=NULL,
                                size=2, alpha=0.5, add.extra.plot=NULL, annotate_signif=FALSE, group_vars="Orientation")
  
  
  violin.B1int.enrichment <-  custom_plot_with_statistics(data=clusters_summarized, x_var="Gcx_content_classification_peak", alpha=0.5, size=1,
                                                          y_var="B1int_enrichment", group_vars="Orientation", jitter=TRUE, annotate_signif = TRUE, fill_var="Gcx_content_classification_peak",
                                                          annotate_only_signif = FALSE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
  violin.B1int.enrichment.stats <- violin.B1int.enrichment$statistical_test
  violin.B1int.enrichment <- violin.B1int.enrichment$plot+
    theme_classic()+
    labs(x="Inner peak Gcx (MFI)", y="B1int enrichment")+
    theme(legend.position="none")
  
  print(violin.B1int.enrichment)
  
  
  
  B1int.enrichment.together <- plot_grid(Gcx.vs.B1intenrichment.graph, violin.B1int.enrichment, nrow=1, align="v", axis="b")
  print(B1int.enrichment.together)
  
  pdf("Perpendicular.data.pdf", width=5, height=2)
  B1int.enrichment.together
  dev.off()
  
  violin.B1int.enrichment.stats <- violin.B1int.enrichment.stats %>%  dplyr::select(-groups)
  
  write.csv2(outwardcluster_data, "Perpendicular.analysis.csv")
  write.csv2(as.data.frame(violin.B1int.enrichment.stats), "Perpendicular.analyis.glycocalyx.vs.b1intenrichment.csv")
  
  
  
  
  
   ## protrusions & retraction fibers - determine whether local and global background differ from Gcx (non-stratified) ------


clusterproperties <- filter(clusterproperties, Orientation!="Bleb")
  
  clusterproperties <- clusterproperties %>%  filter(Orientation=="Leading edge protrusion")
  


signal <-    clusterproperties  %>%  ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Gcx) %>%  mutate(Gcx.signal = "Signal", identifier = row_number()) 
global.bg <- clusterproperties %>% ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Gcx_cluster_global_bg) %>% mutate(Gcx.signal = "Global_bg", identifier = row_number()) %>%  rename(Gcx = Gcx_cluster_global_bg)
local.bg <- clusterproperties %>% ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Gcx_cluster_local_bg) %>% mutate(Gcx.signal = "Local_bg", identifier = row_number()) %>%  rename(Gcx = Gcx_cluster_local_bg)


signal <- rbind(signal, global.bg, local.bg)


pool_pixels = FALSE
grouped_pixels = 5


signal_pooled <- signal %>%
  group_by(Orientation, number, Signal_type, Gcx.signal) %>% arrange(X) %>% 
  mutate(X_group = ceiling(row_number()/grouped_pixels))



signal_pooled <- signal_pooled %>% group_by(Orientation, number, Signal_type, X_group, Gcx.signal) %>% 
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}")
  ) %>%
  ungroup()


if(pool_pixels==TRUE){
  signal_unpooled <- signal
  signal <- signal_pooled
}



signal$random_numbering <- runif(nrow(signal), min = 1, max = 100)
signal$random_numbering <- as.factor(signal$random_numbering)

subset <- signal %>%  subset(Gcx.signal != "Local_bg")

connected.signalvsbg <- ggplot(data=subset, aes(x = Gcx.signal, y = Gcx), color="black") +
  geom_line(aes(group = identifier), alpha = 0.2) +
  geom_point(size=1, stroke=0, alpha=0.2)+
  theme_classic() +
  scale_y_log10()+
  theme(legend.position = "none") +
  labs(y = "Gcx MFI", x = "") +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)+
  facet_wrap(.~Orientation)


p_value <- subset %>%  group_by(Orientation) %>%  tukey_hsd(Gcx ~ Gcx.signal) %>%  add_xy_position(y.trans=log10) %>%  mutate(y.position=y.position*1.05)
#p_value <- wilcox.test(Gcx ~ Gcx.signal, data = subset, paired = TRUE)$p.value


# Add line for comparison
connected.signalvsbg <- connected.signalvsbg +
  stat_pvalue_manual(data=p_value, label= "p.adj")+scale_y_log10()+ scale_x_discrete(labels=c("Global Gcx bg", "Gcx signal"))


# Print the updated plot
print(connected.signalvsbg)


connected.signalvsbg.global <- connected.signalvsbg


subset <- signal %>%  subset(Gcx.signal != "Global_bg")

connected.signalvsbg <- ggplot(data=subset, aes(x = Gcx.signal, y = Gcx), color="black") +
  geom_line(aes(group = identifier), alpha = 0.2) +
  geom_point(size=1, stroke=0, alpha=0.2)+
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Gcx MFI", x = "") +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)+
  facet_wrap(.~Orientation)


p_value <- subset %>%  group_by(Orientation) %>%  tukey_hsd(Gcx ~ Gcx.signal) %>%  add_xy_position(y.trans=log10) %>%  mutate(y.position=y.position*1.05)
#p_value <- wilcox.test(Gcx ~ Gcx.signal, data = subset, paired = TRUE)$p.value


# Add line for comparison
connected.signalvsbg <- connected.signalvsbg +
  stat_pvalue_manual(data=p_value, label= "p.adj")+scale_y_log10() + scale_x_discrete(labels=c("Local Gcx bg", "Gcx signal"))


# Print the updated plot
print(connected.signalvsbg)

connected.signalvsbg.local <- connected.signalvsbg



signalvsbg.combined <- plot_grid(connected.signalvsbg.global, connected.signalvsbg.local, nrow=1)
print(signalvsbg.combined)


Gcx.signal <- median(signal$Gcx[signal$Gcx.signal=="Signal" & !is.na(signal$Gcx)])
Gcx.global.bg <- median(signal$Gcx[signal$Gcx.signal=="Global_bg" & !is.na(signal$Gcx)])
Gcx.local.bg <- median(signal$Gcx[signal$Gcx.signal=="Local_bg" & !is.na(signal$Gcx)])
Fc.signalglobal <- Gcx.signal/Gcx.global.bg
Fc.signallocal <- Gcx.signal/Gcx.local.bg





   ## protrusions & retraction fibers - local enrichment of B1int and Gcx in clusters -----
 
 
 library(viridis)
 
 
 
 clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, CellNR, LocationNR, RoiNR, Signal_type, Clustersize, Peak_position, unique_Peak_ID) %>%  
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 clusterproperties <- clusterproperties %>%  filter(Signal_type=="foreground", !is.na(new_Peak_ID) & !is.na(Peak_position) & !is.na(Peaklocation_X) & 
                                                      Peak_position!="Transition" & Orientation!="Bleb" & Orientation !="Outwardclusters" & Orientation!="Perpendicular control")
 
 clusterproperties <- clusterproperties %>%  filter(Orientation=="Leading edge protrusion")
 
 outlier.removal=FALSE
 if(outlier.removal==TRUE){
 clusterproperties <- clusterproperties %>%  ungroup() %>%  remove_outliers("B1int_enrichment")
 }
 
 clusterproperties$new_Peak_ID <- as.factor(clusterproperties$new_Peak_ID)
 clusterproperties$number <-as.factor(clusterproperties$number)
 
 
 
 clusterproperties$Peak_position <- as.factor(clusterproperties$Peak_position)
 clusterproperties$Peak_position <- factor(clusterproperties$Peak_position, levels=c("Edge", "Peak"), labels=c("Adjacent", "Peak"))
 

 clusterproperties <- clusterproperties %>%  filter(Peak_position=="Peak")
 
 
 # Define thresholds based on the number of datapoints (1/3 with lowest Gcx_global_bg_adjacent is categorized "Low", etc)
 
 clusterproperties <- clusterproperties %>%  filter(!is.nan(Gcx_cluster_local_bg))
 
 
 quantiles <- quantile(clusterproperties$Gcx_cluster_local_bg, probs = c(0, 1/3, 2/3, 1))
 
 clusterproperties$Gcx_content_classification_localbg <- cut(clusterproperties$Gcx_cluster_local_bg, breaks = quantiles, labels = c("L", "M", "H"), include.lowest = TRUE)
 
 
 Gcx_content_threshold_localbg_low =  max(clusterproperties$Gcx_cluster_local_bg[clusterproperties$Gcx_content_classification_localbg=="L"])
 Gcx_content_threshold_localbg_high=  max(clusterproperties$Gcx_cluster_local_bg[clusterproperties$Gcx_content_classification_localbg=="M"])

 
 clusterproperties$Gcx_content_classification_localbg <- factor(clusterproperties$Gcx_content_classification_localbg, levels=c("L", "M", "H"))
 
 # Create a data frame with rectangle information
 rect_data_localbg <- data.frame(
   xmin = c(0, Gcx_content_threshold_localbg_low, Gcx_content_threshold_localbg_high),
   xmax = c(Gcx_content_threshold_localbg_low, Gcx_content_threshold_localbg_high, max(clusterproperties$Gcx_cluster_local_bg, na.rm=TRUE)),
   ymin = rep(-Inf, 3),
   ymax = rep(Inf, 3),
   fill = c( "#00BA38","#1661a8","#bd3333")
 )
 
 gg.rect <- geom_rect(data = rect_data_localbg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                      inherit.aes = FALSE, alpha=0.5, size=0.5)
 
 clusterproperties$CellNR <- as.factor(clusterproperties$CellNR)
 
 
 Gcx.vs.B1intenrichment <- custom_lineplot_with_rsquared(clusterproperties, regression_type =  "log", x_var=Gcx_cluster_local_bg, y_var=B1int_local_enrichment, shape_var=NULL,
                                                         size=2, alpha=0.5, add.extra.plot=gg.rect, annotate_signif=FALSE, group_vars="Orientation")
 Gcx.vs.B1intenrichment.graph <- Gcx.vs.B1intenrichment$plot+
   theme_classic()+
   scale_x_continuous(n.breaks=3)+ 
   theme(legend.position="none")+
   scale_x_log10()+
   labs(x="Local Gcx background", y="\u03B21int enrichment")
 
 print(Gcx.vs.B1intenrichment.graph)
 
 individual.datapoints <- clusterproperties %>% ungroup() %>% 
   dplyr::select(Gcx_cluster_local_bg, B1int_local_enrichment, Gcx_content_classification_localbg) 
 print(individual.datapoints, n=200)
   
 
 
 model <- aov(B1int_enrichment ~ Gcx_content_classification_localbg, data=clusterproperties)
 residuals <- model$residuals
 norm <- shapiro.test(residuals)
 pv <- norm$p.value
 print(pv)
 
 setforplot <- subset(clusterproperties)
 dunn.test::dunn.test(x=setforplot$B1int_enrichment, g= setforplot$Gcx_content_classification_localbg, method= "bonferroni")
 
 
 
 violin.B1int.enrichment <-  custom_plot_with_statistics(data=clusterproperties, x_var="Gcx_content_classification_localbg", alpha=0.5, size=1,
                                                         y_var="B1int_enrichment", group_vars="Orientation", jitter=FALSE, annotate_signif = TRUE, fill_var="Gcx_content_classification_localbg",
                                                         annotate_only_signif = FALSE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 violin.B1int.enrichment.stats <- violin.B1int.enrichment$statistical_test
 violin.B1int.enrichment <- violin.B1int.enrichment$plot+
   theme_classic()+
   labs(x="Local Gcx background", y="\u03B21int enrichment")+
   theme(legend.position="none")
 
 print(violin.B1int.enrichment)
 
 
 
 leadingedge.stats <- subset(clusterproperties)
 
 
 pairwise.results <- pairwise.wilcox.test(x=clusterproperties$B1int_enrichment, 
                                          g=clusterproperties$Gcx_content_classification_localbg,p.adjust.method = "bonferroni")
 print(pairwise.results)
 
 library(rstatix)
 leadingedge.stats %>%  ungroup() %>% rstatix::kruskal_test(B1int_enrichment ~ Gcx_content_classification_localbg) 
 leadingedge.stats %>%  ungroup() %>% rstatix::kruskal_effsize(B1int_enrichment ~ Gcx_content_classification_localbg)
 
 rstatix::kruskal_effsize(leadingedge.stats, Gcx_cellbody_corrected ~ Length_label, ci = TRUE)
 
 leadingedge.stats <- leadingedge.stats %>%
   ungroup() %>% 
   group_by(expnr, Length_label) %>%
   arrange(Orientation) %>% 
   
   dplyr::summarise(
     n = n(),
     median = median(Gcx_cellbody_corrected, na.rm = TRUE),
     Q1 = quantile(Gcx_cellbody_corrected, 0.25, na.rm = TRUE),
     Q3 = quantile(Gcx_cellbody_corrected, 0.75, na.rm = TRUE),
     IQR = IQR(Gcx_cellbody_corrected, na.rm = TRUE),
     lower_whisker = max(min(Gcx_cellbody_corrected, na.rm = TRUE), Q1 - 1.5 * IQR),
     upper_whisker = min(max(Gcx_cellbody_corrected, na.rm = TRUE), Q3 + 1.5 * IQR)
   ) %>%
   ungroup() 
 
 print(leadingedge.stats)
 
 shapiro.test(leadingedge.stats$Gcx_cellbody_corrected[leadingedge.stats$Length_label=="Protrusion base"])
 
 
 pairwise.results <- pairwise.wilcox.test(x=leadingedge.stats$median, g=leadingedge.stats$Length_label,p.adjust.method = "bonferroni")
 print(pairwise.results)
 leadingedge.stats %>%  aov(median ~ Length_label, data = .) %>%  summary()
 leadingedge.stats %>%   tukey_hsd(median ~ Length_label) 
 dunn.test::dunn.test(x=leadingedge.stats$median, g=leadingedge.stats$Length_label,
                      method = "bonferroni")
 
 
 
 
 
 
 
 boxplot_stats <- clusterproperties %>%
   ungroup() %>% 
   group_by(Gcx_content_classification_localbg) %>%
   arrange(Gcx_content_classification_localbg) %>% 
   
   dplyr::summarise(
     n = n(),
     median = median(B1int_enrichment, na.rm = TRUE),
     Q1 = quantile(B1int_enrichment, 0.25, na.rm = TRUE),
     Q3 = quantile(B1int_enrichment, 0.75, na.rm = TRUE),
     IQR = IQR(B1int_enrichment, na.rm = TRUE),
     lower_whisker = max(min(B1int_enrichment, na.rm = TRUE), Q1 - 1.5 * IQR),
     upper_whisker = min(max(B1int_enrichment, na.rm = TRUE), Q3 + 1.5 * IQR)
   ) %>%
   ungroup()
 
 print(boxplot_stats)
 
 
 
 B1int.enrichment.together <- plot_grid(Gcx.vs.B1intenrichment.graph, print(violin.B1int.enrichment), ncol=1, axis="b")+
   theme(plot.margin = margin(r = 5, unit = "mm"))
 print(B1int.enrichment.together)
 
 pdf("B1int.enrichment.lateral.pdf", width=6, height=6)
 print(B1int.enrichment.together)
 dev.off()
 
 Gcx.vs.Gcxenrichment <- custom_lineplot_with_rsquared(clusterproperties, regression_type =  "log", x_var=Gcx_cluster_local_bg, y_var=Gcx_local_enrichment, shape_var=NULL,
                                                         size=2, alpha=0.5, add.extra.plot=gg.rect, annotate_signif=FALSE, group_vars="Orientation")
 Gcx.vs.Gcxenrichment.graph <- Gcx.vs.Gcxenrichment$plot+
   theme_classic()+
   scale_x_continuous(n.breaks=3)+ 
   theme(legend.position="none")+
   labs(x="Local Gcx background", y="Gcx enrichment")+
   scale_x_log10()
 
 print(Gcx.vs.Gcxenrichment.graph)
 
 individual.datapoints <- clusterproperties %>%  ungroup() %>% 
   dplyr::select(Gcx_cluster_local_bg, Gcx_local_enrichment, Gcx_content_classification_localbg)
 print(individual.datapoints, n=200)
 
 
 boxplot_stats <- clusterproperties %>%
   ungroup() %>% 
   group_by(Gcx_content_classification_localbg) %>%
   arrange(Gcx_content_classification_localbg) %>% 
   
   dplyr::summarise(
     n = n(),
     median = median(Gcx_local_enrichment, na.rm = TRUE),
     Q1 = quantile(Gcx_local_enrichment, 0.25, na.rm = TRUE),
     Q3 = quantile(Gcx_local_enrichment, 0.75, na.rm = TRUE),
     IQR = IQR(Gcx_local_enrichment, na.rm = TRUE),
     lower_whisker = max(min(Gcx_local_enrichment, na.rm = TRUE), Q1 - 1.5 * IQR),
     upper_whisker = min(max(Gcx_local_enrichment, na.rm = TRUE), Q3 + 1.5 * IQR)
   ) %>%
   ungroup()
 
 print(boxplot_stats)
 
 
 violin.Gcx.enrichment <-  custom_plot_with_statistics(data=clusterproperties, x_var="Gcx_content_classification_localbg", alpha=0.5, size=1,
                                                         y_var="Gcx_local_enrichment", group_vars="Orientation", jitter=FALSE, annotate_signif = TRUE, fill_var="Gcx_content_classification_localbg",
                                                         annotate_only_signif = FALSE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 violin.Gcx.enrichment.stats <- violin.Gcx.enrichment$statistical_test
 violin.Gcx.enrichment <- violin.Gcx.enrichment$plot+
   theme_classic()+
   labs(x="Local Gcx background", y="Gcx enrichment")+
   theme(legend.position="none")
 
 print(violin.Gcx.enrichment)
 
 
 Gcx.enrichment.together <- plot_grid(Gcx.vs.Gcxenrichment.graph, print(violin.Gcx.enrichment), ncol=1, axis="b")+
   theme(plot.margin = margin(r = 5, unit = "mm"))
 print(Gcx.enrichment.together)
 
 pdf("Gcx.enrichment.lateral.pdf", width=6, height=6)
 print(Gcx.enrichment.together)
 dev.off()
 
 
 
 lateral.enrichment.together <- plot_grid(connected.signalvsbg.local, B1int.enrichment.together, Gcx.enrichment.together, axis="b", nrow=1)+
   theme(plot.margin = margin(r = 5, unit = "mm"))
 print(lateral.enrichment.together)
 
 pdf("Lateral.enrichment.graphs.pdf", width=8, height=5)
 print(lateral.enrichment.together)
 dev.off()
 
 mean(clusterproperties$B1int_global_enrichment[clusterproperties$Gcx_content_classification_localbg=="L" & clusterproperties$Orientation=="Leading edge protrusion"], na.rm=TRUE)
 mean(clusterproperties$B1int_global_enrichment[clusterproperties$Gcx_content_classification_localbg=="L" & clusterproperties$Orientation=="Cellbody"], na.rm=TRUE)
 
 mean(clusterproperties$B1int_global_enrichment[clusterproperties$Gcx_content_classification_localbg=="H" & clusterproperties$Orientation=="Leading edge protrusion"], na.rm=TRUE)

 mean(clusterproperties$B1int_global_enrichment[clusterproperties$Gcx_content_classification_localbg=="H" & clusterproperties$Orientation=="Cellbody"], na.rm=TRUE)
 
 mean(clusterproperties$B1int_global_enrichment[clusterproperties$Orientation=="Cellbody"], na.rm=TRUE)
 mean(clusterproperties$B1int_local_enrichment[clusterproperties$Orientation=="Cellbody"], na.rm=TRUE)

   ## protrusions & retraction fibers - cluster size ----
 
 
 
 clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, CellNR, LocationNR, RoiNR, Signal_type, Clustersize, Peak_position, unique_Peak_ID) %>%  
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 clusterproperties <- clusterproperties %>%  filter(Signal_type=="foreground", !is.na(new_Peak_ID) & !is.na(Peak_position) & !is.na(Peaklocation_X) & 
                                                      Peak_position!="Transition" & Orientation!="Bleb" & Orientation!="Perpendicular control")
 
 
 clustersize.violin <- ggplot(subset(clusterproperties, Orientation=="Leading edge protrusion"), aes(x=Orientation, y=Clustersize))+
   geom_violin(scale="width")+
   
   geom_jitter(size=1, width=0.5)+
   geom_boxplot(outlier.shape=NA, width=0.5)+
   theme_classic()+
   labs(x="", y="Cluster size (micron)")
 
 print(clustersize.violin)
 
 max(clusterproperties$Clustersize)
 min(clusterproperties$Clustersize)
 
 pdf("Clustersize.violin.pdf", width=1, height=2)
 print(clustersize.violin)
 dev.off()
 
 
 boxplot_stats <- clusterproperties %>%  subset(Orientation=="Leading edge protrusion") %>% 
   ungroup() %>% 
   group_by(Orientation) %>%
   
   dplyr::summarise(
     n = n(),
     median = median(Clustersize, na.rm = TRUE),
     Q1 = quantile(Clustersize, 0.25, na.rm = TRUE),
     Q3 = quantile(Clustersize, 0.75, na.rm = TRUE),
     IQR = IQR(Clustersize, na.rm = TRUE),
     lower_whisker = max(min(Clustersize, na.rm = TRUE), Q1 - 1.5 * IQR),
     upper_whisker = min(max(Clustersize, na.rm = TRUE), Q3 + 1.5 * IQR)
   ) %>%
   ungroup()
 
 boxplot_stats
 
 ## -----
 ## outward clusters - calculate segregation -----
 
 protrusion_data <- subset(protrusion_data, Orientation=="Outwardclusters" | Orientation=="Perpendicular control")
 protrusion_data <- protrusion_data %>%  ungroup() %>%  group_by(number, Orientation, Signal_type) %>% 
   mutate(firstcluster=ifelse(new_Peak_ID == min(new_Peak_ID, na.rm=TRUE), "Yes", "No"),
          firstcluster = as.factor(firstcluster),
   )
 
 outwardclusters <- protrusion_data %>%  group_by(number, CellNR, Signal_type, Orientation, firstcluster, Edgenr, Peak_position,new_Peak_ID) %>%  
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 outwardclusters <- outwardclusters %>%  filter(Signal_type=="foreground", !is.na(firstcluster) & 
                                                  !is.na(new_Peak_ID) &
                                                  firstcluster=="Yes" &
                                                  Peak_position!="Transition" &
                                                  !(Peak_position=="Edge" & Edgenr==0))
 
 outwardclusters <- outwardclusters %>%
   mutate(ratio_cellbody_corrected = ratio/ ratio_cellbody_average,
          Gcx_cellbody_corrected = Gcx/ Gcx_cellbody_average) %>% 
   ungroup() %>%  group_by(unique_Peak_ID, new_Peak_ID) %>%
   dplyr::filter(max(row_number())>1) %>% 
   mutate(segregation_length = ifelse(Peak_position=="Peak", X_micron[Peak_position=="Edge" & Edgenr!=0] - X_micron[Peak_position=="Peak"], NA),
          Gcx_edge = ifelse(Peak_position=="Peak", Gcx[Peak_position=="Edge"], Gcx),
          Gcx_enrichment_edge = Gcx/ Gcx[Peak_position=="Edge"],
          B1int_enrichment_edge = B1int / B1int[Peak_position=="Edge"]
   ) %>%  dplyr::select(number, CellNR, Signal_type, Orientation, firstcluster, Edgenr, Peak_position,new_Peak_ID, segregation_length,
                        Gcx_edge, Gcx_cellbody_corrected, ratio_cellbody_corrected, Gcx_enrichment_edge, B1int_enrichment_edge, Gcx, B1int, ratio)
 
 
 paired.test.graph.gcx <- function(data,subset,label){
   outwardclusters <- subset(outwardclusters, Orientation==subset)
   
   outwardclusters_connected <- ggplot(data=outwardclusters,aes(x=Peak_position, y=Gcx))+
     geom_point()+
     geom_line(aes(group=unique_Peak_ID))+
     theme_classic()+
     scale_y_log10(expand=c(0,0.2))+
     labs(x="", y="Gcx")+facet_wrap(.~Orientation)
   wilcoxtestresults <- wilcox.test(Gcx ~ Peak_position, data = outwardclusters, paired = TRUE)
   p_value <- wilcoxtestresults$p.value
   
   annotations <- geom_signif(comparisons = list(c("Edge", "Peak")), 
                              annotations = paste("p =", format(p_value, scientific = FALSE, digits = 2)), 
                              textsize = 4, color="black", margin_top=0.1)
   
   outwardclusters_connected <- outwardclusters_connected + annotations
   
   
   outwardclusters_connected.gcx<- plot_grid(outwardclusters_connected,labels=c(label), vjust=0)+ 
     theme(plot.margin = margin(t = 5, unit = "mm"))
   
   print(outwardclusters_connected.gcx)
   return(list(outwardclusters_connected.gcx, wilcoxtestresults))
 }
 outwardclusters_connected.gcx <- paired.test.graph.gcx(data=outwardclusters, subset="Outwardclusters", label="e")[[1]]
 
 paired.test.graph.b1int <- function(data,subset,label){
   outwardclusters <- subset(outwardclusters, Orientation==subset)
   
   outwardclusters_connected <- ggplot(data=outwardclusters,aes(x=Peak_position, y=B1int))+
     geom_point()+
     geom_line(aes(group=unique_Peak_ID))+
     theme_classic()+
     scale_y_log10(expand=c(0,0.2))+
     labs(x="", y="Gcx")+facet_wrap(.~Orientation)
   wilcoxtestresults <- wilcox.test(B1int ~ Peak_position, data = outwardclusters, paired = TRUE)
   p_value <- wilcoxtestresults$p.value
   
   annotations <- geom_signif(comparisons = list(c("Edge", "Peak")), 
                              annotations = paste("p =", format(p_value, scientific = FALSE, digits = 2)), 
                              textsize = 4, color="black", margin_top=0.1)
   
   outwardclusters_connected <- outwardclusters_connected + annotations
   
   
   outwardclusters_connected.gcx<- plot_grid(outwardclusters_connected,labels=c(label), vjust=0)+ 
     theme(plot.margin = margin(t = 5, unit = "mm"))
   
   print(outwardclusters_connected.gcx)
   return(list(outwardclusters_connected.gcx, wilcoxtestresults))
 }
 outwardclusters_connected.b1int <- paired.test.graph.b1int(data=outwardclusters, subset="Outwardclusters", label="f")[[1]]
 
 outwardclusters_connected <- plot_grid(outwardclusters_connected.b1int, outwardclusters_connected.gcx, nrow=1)
 print(outwardclusters_connected)
 
 perpendicularcontrol_connected.gcx <- paired.test.graph.gcx(data=outwardclusters, subset="Perpendicular control", label="e")[[1]]
 perpendicularcontrol_connected.b1int <- paired.test.graph.b1int(data=outwardclusters, subset="Perpendicular control", label="f")[[1]]
 
 
 perpendicularcontrol_connected <- plot_grid(perpendicularcontrol_connected.b1int, perpendicularcontrol_connected.gcx, nrow=1)
 print(perpendicularcontrol_connected)
 
 ## stratify for Gcx intensity: plot  - Gcx at edge vs. integrin enrichment------
 outwardclusters <- subset(outwardclusters, Orientation=="Outwardclusters")
 outwardclusters <- subset(outwardclusters, Peak_position=="Peak")
 
 
 quantiles <- quantile(outwardclusters$Gcx_edge, probs = c(0, 1/3, 2/3, 1))
 
 outwardclusters$Gcx_content_classification_edge <- cut(outwardclusters$Gcx_edge, breaks = quantiles, labels = c("L", "M", "H"), include.lowest = TRUE)
 
 
 Gcx_content_threshold_edge_low =  max(outwardclusters$Gcx_edge[outwardclusters$Gcx_content_classification_edge=="L"])
 Gcx_content_threshold_edge_high=  max(outwardclusters$Gcx_edge[outwardclusters$Gcx_content_classification_edge=="M"])
 
 
 outwardclusters$Gcx_content_classification_edge <- factor(outwardclusters$Gcx_content_classification_edge, levels=c("L", "M", "H"))
 
 # Create a data frame with rectangle information
 rect_data_edge <- data.frame(
   xmin = c(0, Gcx_content_threshold_edge_low, Gcx_content_threshold_edge_high),
   xmax = c(Gcx_content_threshold_edge_low, Gcx_content_threshold_edge_high, max(outwardclusters$Gcx_edge, na.rm=TRUE)),
   ymin = rep(-Inf, 3),
   ymax = rep(Inf, 3),
   fill = c( "#00BA38","#1661a8","#bd3333")
 )
 
 gg.rect <- geom_rect(data = rect_data_edge, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                      inherit.aes = FALSE, alpha=0.5, size=0.5)
 
 outwardclusters$CellNR <- as.factor(outwardclusters$CellNR)
 
 
 Gcx.vs.B1intenrichment <- custom_lineplot_with_rsquared(outwardclusters, regression_type =  "linear", x_var=Gcx_edge, y_var=B1int_enrichment_edge, shape_var=NULL,
                                                         size=2, alpha=0.5, add.extra.plot=gg.rect, annotate_signif=FALSE, group_vars="Orientation")
 Gcx.vs.B1intenrichment.graph <- Gcx.vs.B1intenrichment$plot+
   theme_classic()+
   scale_x_continuous(n.breaks=3)+
   theme(legend.position="none")+
   labs(x="Edge Gcx", y="\u03B21int enrichment")
 
 print(Gcx.vs.B1intenrichment.graph)
 
 Gcx.vs.B1intenrichment.graph <- plot_grid(Gcx.vs.B1intenrichment.graph,labels=c("g"), vjust=0)+
   theme(plot.margin = margin(t = 5, unit = "mm"))
 print(Gcx.vs.B1intenrichment.graph)
 
 
 
 violin.B1int.enrichment <-  custom_plot_with_statistics(data=outwardclusters, x_var="Gcx_content_classification_edge", alpha=0.5, size=1,
                                                         y_var="B1int_enrichment_edge", group_vars="Orientation", jitter=FALSE, annotate_signif = TRUE, fill_var="Gcx_content_classification_edge",
                                                         annotate_only_signif = FALSE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 violin.B1int.enrichment.stats <- violin.B1int.enrichment$statistical_test
 violin.B1int.enrichment <- violin.B1int.enrichment$plot+
   theme_classic()+
   labs(x="Edge Gcx", y=" \u03B21int enrichment")+
   theme(legend.position="none")
 
 print(violin.B1int.enrichment)
 
 
 violin.B1int.enrichment <- plot_grid(violin.B1int.enrichment,labels=c("h"), vjust=0)+
   theme(plot.margin = margin(t = 5, unit = "mm"))
 print(violin.B1int.enrichment)
 
 
 B1int.enrichment.together <- plot_grid(Gcx.vs.B1intenrichment.graph, violin.B1int.enrichment, ncol=2, axis="b")
 print(B1int.enrichment.together)
 
 
 
 ## protrusions - plot of Gcx MFI vs. adjacent / peak region (connected)-----
 
 clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, Peak_position, unique_Peak_ID) %>%  
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 clusterproperties <- clusterproperties %>%  filter(Signal_type=="foreground", !is.na(new_Peak_ID)  & !is.na(Peak_position) &
                                                      Orientation!="Bleb" & Orientation!="Outwardclusters")
 
 clusterproperties$new_Peak_ID <- as.factor(clusterproperties$new_Peak_ID)
 clusterproperties$number <-as.factor(clusterproperties$number)
 
 clusterproperties <- clusterproperties %>%  filter(Peak_position!="Transition" & !is.na(Peak_position))
 
 clusterproperties$Peak_position <- as.factor(clusterproperties$Peak_position)
 clusterproperties$Peak_position <- factor(clusterproperties$Peak_position, levels=c("Edge", "Peak"), labels=c("Adjacent", "Peak"))
 
 
 clusterproperties$random_numbering <- runif(nrow(clusterproperties), min = 1, max = 100)
 clusterproperties$random_numbering <- as.factor(clusterproperties$random_numbering)
 
 library(ggplot2)
 library(ggpubr)
 library(ggsignif)
 
 # Your original ggplot code
 peakposition.gcx.mfi <- ggplot(clusterproperties, aes(x = Peak_position, y = Gcx, color = random_numbering)) +
   geom_point(size = 1, stroke = 0) +
   geom_line(aes(group = unique_Peak_ID), alpha = 0.5) +
   geom_smooth(method="lm")+
   scale_y_log10(expand=c(0,0.2)) +
   theme_classic() +
   theme(legend.position = "none") +
   labs(y = "Gcx MFI", x = "Region") +
   scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)
 
 
 p_value <- wilcox.test(Gcx ~ Peak_position, data = clusterproperties, paired = TRUE)$p.value
 
 
 # Add line for comparison
 peakposition.gcx.mfi <- peakposition.gcx.mfi +
   
   geom_signif(comparisons = list(c("Adjacent", "Peak")), 
               annotations = paste("p =", format(p_value, scientific = FALSE, digits = 2)), 
               textsize = 4, color="black")
 
 # Print the updated plot
 print(peakposition.gcx.mfi)
 
 peakposition.gcx.mfi <- plot_grid(peakposition.gcx.mfi, labels="e", vjust=0)+
   theme(plot.margin = margin(t = 5, unit = "mm"))
 
 print(ggdraw(peakposition.gcx.mfi))
 
 Gcx.median.adj <-  median(clusterproperties$Gcx[clusterproperties$Peak_position=="Adjacent"])
 Gcx.median.peak <-  median(clusterproperties$Gcx[clusterproperties$Peak_position=="Peak"])
 
 
 
   ## protrusions - Absolute enrichment of B1int and Gcx in clusters ------
 
 library(viridis)
 
 
 
 clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, CellNR, LocationNR, RoiNR, Signal_type, Clustersize, Peak_position, unique_Peak_ID) %>%  
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 clusterproperties <- clusterproperties %>%  filter(Signal_type=="foreground", !is.na(new_Peak_ID) & !is.na(Peak_position) & Peak_position!="Transition")
 
 clusterproperties$new_Peak_ID <- as.factor(clusterproperties$new_Peak_ID)
 clusterproperties$number <-as.factor(clusterproperties$number)
 
 
 
 clusterproperties$Peak_position <- as.factor(clusterproperties$Peak_position)
 clusterproperties$Peak_position <- factor(clusterproperties$Peak_position, levels=c("Edge", "Peak"), labels=c("Adjacent", "Peak"))
 
 clusterproperties <- clusterproperties %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, unique_Peak_ID) %>%  
   mutate(B1int_adjacent = ifelse(Peak_position=="Peak", B1int[Peak_position=="Adjacent"], B1int),
          Gcx_adjacent = ifelse(Peak_position=="Peak", Gcx[Peak_position=="Adjacent"], Gcx),
          Gcx_local_bg_adjacent = ifelse(Peak_position=="Peak", Gcx_cluster_local_bg[Peak_position=="Adjacent"], Gcx_cluster_local_bg),
          Gcx_global_bg_adjacent = ifelse(Peak_position=="Peak", Gcx_cluster_global_bg[Peak_position=="Adjacent"], Gcx_cluster_global_bg)
   )
 
 clusterproperties <- clusterproperties %>%  filter(Peak_position=="Peak")
 
 
 # Define thresholds based on the number of datapoints (1/3 with lowest Gcx_global_bg_adjacent is categorized "Low", etc)
 
 
 quantiles <- quantile(clusterproperties$Gcx_adjacent, probs = c(0, 1/3, 2/3, 1))
 
 clusterproperties$category <- cut(clusterproperties$Gcx_adjacent, breaks = quantiles, labels = c("Low", "Medium", "High"), include.lowest = TRUE)
 
 
 Gcx_content_threshold_low =  max(clusterproperties$Gcx_adjacent[clusterproperties$category=="Low"])
 Gcx_content_threshold_high=  max(clusterproperties$Gcx_adjacent[clusterproperties$category=="Medium"])
 
 clusterproperties <- clusterproperties %>%   group_by(number, Orientation, new_Peak_ID, unique_Peak_ID, Signal_type, Clustersize) %>%
   mutate(Gcx_content_classification= ifelse(Gcx_adjacent < Gcx_content_threshold_low, "L", "M"),
          Gcx_content_classification = ifelse(Gcx_adjacent > Gcx_content_threshold_high, "H", Gcx_content_classification),
          Gcx_content_classification= as.factor(Gcx_content_classification))
 
 
 clusterproperties$Gcx_content_classification <- factor(clusterproperties$Gcx_content_classification, levels=c("L", "M", "H"))
 
 # Create a data frame with rectangle information
 rect_data <- data.frame(
   xmin = c(0, Gcx_content_threshold_low, Gcx_content_threshold_high),
   xmax = c(Gcx_content_threshold_low, Gcx_content_threshold_high, max(clusterproperties$Gcx_adjacent, na.rm=TRUE)),
   ymin = rep(-Inf, 3),
   ymax = rep(Inf, 3),
   fill = c( "#00BA38","#1661a8","#bd3333"),
   alpha = 0.5,
   size = 0.5
 )
 
 gg.rect <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                      inherit.aes = FALSE, alpha = rect_data$alpha, size = rect_data$size)
 
 
 Gcx.vs.Gcxenrichment <- custom_lineplot_with_rsquared(clusterproperties, regression_type =  "log", x_var=Gcx_adjacent, 
                                                       y_var=Gcx_enrichment, size=1, alpha=0.5, add.extra.plot=gg.rect, annotate_signif=FALSE)
 Gcx.vs.Gcxenrichment.graph <-  Gcx.vs.Gcxenrichment$plot+ 
   theme_classic()+
   scale_x_continuous(n.breaks=3)+
   theme(legend.position="none")+
   labs(x="Adjacent Gcx", y="Gcx enrichment")+
   geom_hline(yintercept=1, linetype="dashed")
 
 print(Gcx.vs.Gcxenrichment.graph)
 
 
 library(ggbreak)
 
 violin.Gcx.enrichment <- custom_plot_with_statistics(data=subset(clusterproperties,!is.na(Gcx_content_classification)), x_var="Gcx_content_classification", alpha=0.5, size=1,
                                                      y_var="Gcx_enrichment", group_vars=NULL, jitter=FALSE, annotate_signif = TRUE, fill_var="Gcx_content_classification",
                                                      annotate_only_signif = FALSE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 violin.Gcx.enrichment <- violin.Gcx.enrichment$plot+
   theme_classic()+
   labs(x="Adjacent Gcx", y="Gcx enrichment")+
   theme(legend.position = "none")+
   geom_hline(yintercept=1, linetype="dashed")
 
 
 print(violin.Gcx.enrichment)
 
 Gcx.enrichment.together <- plot_grid(Gcx.vs.Gcxenrichment.graph, print(violin.Gcx.enrichment), ncol=2, axis="b")+
   theme(plot.margin = margin(r = 5, unit = "mm"))
 print(Gcx.enrichment.together)
 
 
 
 Gcx.vs.B1intenrichment <- custom_lineplot_with_rsquared(clusterproperties, regression_type =  "log", x_var=Gcx_adjacent, y_var=B1int_enrichment,
                                                         size=1, alpha=0.5, add.extra.plot=gg.rect, annotate_signif=FALSE)
 Gcx.vs.B1intenrichment.graph <- Gcx.vs.B1intenrichment$plot+
   theme_classic()+
   scale_x_continuous(n.breaks=3)+
   theme(legend.position="none")+
   labs(x="Adjacent Gcx", y="\u03B21int enrichment")
 
 print(Gcx.vs.B1intenrichment.graph)
 
 
 
 violin.B1int.enrichment <-  custom_plot_with_statistics(data=subset(clusterproperties,!is.na(Gcx_content_classification)), x_var="Gcx_content_classification", alpha=0.5, size=1,
                                                         y_var="B1int_enrichment", group_vars=NULL, jitter=FALSE, annotate_signif = TRUE, fill_var="Gcx_content_classification",
                                                         annotate_only_signif = FALSE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 violin.B1int.enrichment.stats <- violin.B1int.enrichment$statistical_test
 violin.B1int.enrichment <- violin.B1int.enrichment$plot+
   theme_classic()+
   labs(x="Adjacent Gcx", y="\u03B21int enrichment")+
   theme(legend.position="none")
 
 print(violin.B1int.enrichment)
 
 
 B1int.enrichment.together <- plot_grid(Gcx.vs.B1intenrichment.graph, print(violin.B1int.enrichment), ncol=2, axis="b")+
   theme(plot.margin = margin(r = 5, unit = "mm"))
 print(B1int.enrichment.together)
 
 
 
 Enrichment_graphs_global <- plot_grid( B1int.enrichment.together, Gcx.enrichment.together, nrow=2,labels=c("f", "g"), vjust=0)+ 
   theme(plot.margin = margin(t = 5, unit = "mm"))
 
 
 clusterproperties <- clusterproperties %>%  ungroup() %>% arrange(Gcx_enrichment)
 overview <- clusterproperties %>%  select(number, Orientation, CellNR, LocationNR, RoiNR, unique_Peak_ID, Gcx_enrichment)
 
 
 print(ggdraw(Enrichment_graphs_global))
 
 
 clusterproperties <- clusterproperties %>%  arrange(Gcx_enrichment)
 
 
 
   ## protrusions - linearity index versus gcx enrichment -----
      linearity <- custom_lineplot_with_rsquared(clusterproperties, x_var=Linearity_index, y_var=Gcx_local_enrichment, size=1, regression_type="linear")
      linearity.plot <- linearity$plot
      print(linearity.plot)
   ## protrusions - scatter plot B1int enrichment factor vs. Gcx enrichment factor to indicate co-segregation------
 legendsize=0.3
 
 
 cosegregation <- custom_lineplot_with_rsquared(clusterproperties, regression_type =  "linear", x_var=B1int_global_enrichment, y_var=Gcx_local_enrichment,
                                                         color_var="Gcx_cluster_local_bg", group_var=NULL, size=1, alpha=1)
 cosegregation.graph <- cosegregation$plot+
   theme_classic()+
   scale_y_log10()+
   scale_x_log10()+
   theme(legend.position="top", legend.direction="horizontal")+
   labs(x="Cumulative \u03B21 int enrichment", y="Local Gcx enrichment", color= "Adjacent cumulative Gcx bg")+
   theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
         legend.key.height = unit(legendsize, 'cm'), #change legend key height
         legend.key.width = unit(legendsize*1.5, 'cm'))+
   guides(color = guide_colourbar(title.position="top", title.hjust = 0.5))+
   scale_color_viridis_c(option = "inferno", begin = 0, end = 0.8, direction=-1)

 
 box_area <- clusterproperties %>% ungroup() %>% 
   filter(Signal_type=="foreground") %>%
   summarize(ymin = min(Gcx_local_enrichment, na.rm=TRUE), ymax = 1, xmin = min(clusterproperties$B1int_enrichment, na.rm=TRUE), 
             xmax = max(clusterproperties$B1int_enrichment, na.rm=TRUE))
 
 
 box <- geom_rect(data = box_area, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill="blue", 
           inherit.aes = FALSE, alpha = 0, color=alpha("black",1), linetype="dotted")
 

 cosegregation.graph <- cosegregation.graph + box
 
 print(cosegregation.graph)
 
 
 cosegregation.graph <-  plot_grid(cosegregation.graph, labels="h", vjust=0)+
   theme(plot.margin = margin(t = 5, unit = "mm"))
 
 
 print(ggdraw(cosegregation.graph))
 
 
   ## protrusions - 3d scatter plot of Gcx MFI vs. Gcx enrichment factor and B1 integrin enrichment factor ----

 library(plotly)
 
 # Assuming clusterproperties is a data frame or tibble
 plot_ly(
   z = clusterproperties$Gcx_adjacent,
   y = log10(clusterproperties$Gcx_global_enrichment),
   x = log10(clusterproperties$B1int_global_enrichment),
   type = "scatter3d",
   mode = "markers",
   marker = list(
     color = clusterproperties$Gcx_adjacent,  # Color based on the x-axis values
     colorscale = 'Viridis',  # You can choose a different color scale
     opacity = 0.5,
     size = 5,
     line = list(color = "black", width = 2)  # Border color and width
   )
 ) %>%
   layout(
     scene = list(
       zaxis = list(title = 'Cluster-adjacent Gcx'),
       yaxis = list(title = 'Gcx global enrichment', autorange="reversed"),
       xaxis = list(title = 'B1int global enrichment')
     )
   )
 
 
   ## protrusions - clusters - combined graphs ----
 

 
 
  combined_graphs_clusters_efg <- plot_grid(peakposition.gcx.mfi, Enrichment_graphs_global, axis="b", rel_widths = c(0.34,0.66))
 print(ggdraw(combined_graphs_clusters_efg))
 
 
 combined_graphs_clusters_h <- cosegregation.graph
 
 combined_graphs_clusters_all <- plot_grid(combined_graphs_clusters_efg, combined_graphs_clusters_h, align="h", axis="b", rel_widths=c(0.65,0.35))
 print(ggdraw(combined_graphs_clusters_all))
 
 
 
   ## protrusions - global B1 int enrichment vs local B1 int enrichment
   ## protrusions - clusters Gcx vs. Gcx cluster_global_bg-----
 
      signal <-    data.frame(Gcx = clusterproperties$Gcx, Type = rep("Signal", length(clusterproperties$Gcx)), clusternr = 1:length(clusterproperties$Gcx))
      global.bg <- data.frame(Gcx = clusterproperties$Gcx_cluster_global_bg, Type = rep("Global_bg", length(clusterproperties$Gcx_cluster_global_bg)), clusternr = 1:length(clusterproperties$Gcx_cluster_global_bg))
      local.bg <- data.frame(Gcx = clusterproperties$Gcx_cluster_local_bg, Type = rep("Local_bg", length(clusterproperties$Gcx_cluster_local_bg)), clusternr = 1:length(clusterproperties$Gcx_cluster_local_bg))
      cluster.signalvsbg <- rbind(signal, global.bg, local.bg)
      cluster.signalvsbg$Type <- as.factor(cluster.signalvsbg$Type)
      cluster.signalvsbg$random_numbering <- runif(nrow(cluster.signalvsbg), min = 1, max = 100)
      cluster.signalvsbg$random_numbering <- as.factor(cluster.signalvsbg$random_numbering)
      
      subset <- subset(cluster.signalvsbg, Type!="Local_bg")
      connected.signalvsbg <- ggplot(data=subset, aes(x = Type, y = Gcx, color=random_numbering)) +
        #geom_point(size = 1, stroke = 0) +
        geom_line(aes(group = clusternr), alpha = 0.5) +
        scale_y_log10(expand=c(0,0.2)) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(y = "Gcx MFI", x = "") +
        scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)+
        scale_x_discrete(labels=c("Global_bg" = "Multi-peak \n background", "Signal" = "Peak  \nsignal"))
      
      
      p_value <- wilcox.test(Gcx ~ Type, data = subset, paired = TRUE)$p.value
      
      
      # Add line for comparison
      connected.signalvsbg <- connected.signalvsbg +
        geom_signif(comparisons = list(c("Global_bg", "Signal")), 
                    annotations = paste("p =", format(p_value, scientific = TRUE, digits = 2)), 
                    textsize = 4, color="black")
      
      # Print the updated plot
      print(connected.signalvsbg)
      
      
      
      connected.signalvsbg.global <- connected.signalvsbg
      
      
      subset <- subset(cluster.signalvsbg, Type!="Global_bg")
      
      connected.signalvsbg <- ggplot(data=subset, aes(x = Type, y = Gcx, color=random_numbering)) +
        #geom_point(size = 1, stroke = 0) +
        geom_line(aes(group = clusternr), alpha = 0.5) +
        scale_y_log10(expand=c(0,0.2)) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(y = "Gcx MFI", x = "") +
        scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)+
        scale_x_discrete(labels=c("Local_bg" = "Single-peak \n background", "Signal" = "Peak \nsignal"))
      
      
      p_value <- wilcox.test(Gcx ~ Type, data = subset, paired = TRUE)$p.value
      
      
      # Add line for comparison
      connected.signalvsbg <- connected.signalvsbg +
        
        geom_signif(comparisons = list(c("Local_bg", "Signal")), 
                    annotations = paste("p =", format(p_value, scientific = TRUE, digits = 2)), 
                    textsize = 4, color="black")
      
      # Print the updated plot
      print(connected.signalvsbg)
      
      
      connected.signalvsbg.local <- connected.signalvsbg
      
      
      
      connected.signalvsbg <- plot_grid(connected.signalvsbg.global, connected.signalvsbg.local, labels="e", vjust=0)+
        theme(plot.margin = margin(t = 5, unit = "mm"))
      
      print(ggdraw(connected.signalvsbg))
      
   ## protrusions - paired Gcx of global bg vs foreground for adjacent and peak regions in clusters ------
      
      
      clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, Peak_position, unique_Peak_ID) %>%  
        dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
      clusterproperties <- clusterproperties %>%  filter(Signal_type=="foreground", !is.na(new_Peak_ID)  & !is.na(Peak_position))
      clusterproperties$new_Peak_ID <- as.factor(clusterproperties$new_Peak_ID)
      clusterproperties$number <-as.factor(clusterproperties$number)
      clusterproperties$Peak_position <- as.factor(clusterproperties$Peak_position)
      clusterproperties$Peak_position <- factor(clusterproperties$Peak_position, levels=c("Edge", "Inside"), labels=c("Adjacent", "Peak"))
      
      
      signal <-    clusterproperties  %>%  ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Peak_position, Gcx) %>%  mutate(Gcx.signal = "Signal", identifier = row_number()) 
      global.bg <- clusterproperties %>% ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Peak_position, Gcx_cluster_global_bg) %>% mutate(Gcx.signal = "Global_bg", identifier = row_number()) %>%  rename(Gcx = Gcx_cluster_global_bg)
      local.bg <- clusterproperties %>% ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Peak_position, Gcx_cluster_local_bg) %>% mutate(Gcx.signal = "Local_bg", identifier = row_number()) %>%  rename(Gcx = Gcx_cluster_local_bg)
      
      
      signal <- rbind(signal, global.bg, local.bg)
      
      
      pool_pixels = TRUE
      grouped_pixels = 1
      
      
      signal_pooled <- signal %>%
        group_by(Orientation, number, Signal_type, Peak_position, Gcx.signal) %>% arrange(X) %>% 
        mutate(X_group = ceiling(row_number()/grouped_pixels))
      
      
      
      signal_pooled <- signal_pooled %>% group_by(Orientation, number, Signal_type, Peak_position, X_group, Gcx.signal) %>% 
        summarise(
          across(where(is.numeric), ~ mean(.), .names = "{.col}")
        ) %>%
        ungroup()
      
      
      if(pool_pixels==TRUE){
        signal_unpooled <- signal
        signal <- signal_pooled
      }
      
      
      
      signal$random_numbering <- runif(nrow(signal), min = 1, max = 100)
      signal$random_numbering <- as.factor(signal$random_numbering)
      
      subset <- signal %>%  subset(Gcx.signal != "Local_bg")
      
      connected.signalvsbg <- ggplot(data=subset, aes(x = Gcx.signal, y = Gcx), color="black") +
        geom_line(aes(group = identifier), alpha = 0.2) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(y = "Gcx MFI", x = "") +
        scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)+
        facet_wrap(.~Peak_position)
      
      
      tukey <-  subset %>%  ungroup() %>%  group_by(Orientation, Peak_position) %>% tukey_hsd(Gcx ~ Gcx.signal) %>% 
        add_significance() %>% 
        add_xy_position() %>% 
        mutate(y.position.log10= log10(y.position))
      
      
      # Add line for comparison
      connected.signalvsbg <- connected.signalvsbg +
        ggsignif::geom_signif(data = tukey, aes(xmin= group1, xmax= group2, y_position=y.position.log10, annotations=p.adj.signif), manual=TRUE)+
        scale_x_discrete(labels=c("Global_bg" = "Multi-peak \n background", "Signal" = "Peak  \nsignal"))+
        scale_y_log10()
      
      # Print the updated plot
      print(connected.signalvsbg)
      
     
      
 ## make corrections for  blebs------
 
 
 protrusion_data <- subset(B1int_GCX_ratio, Orientation =="Bleb")
      
      
      if(any(protrusion_data$Orientation=="Cellbody")){
        cellbody_summary_percell <- protrusion_data  %>% filter(Orientation=="Cellbody" & Signal_type=="foreground") %>%  
          group_by(Orientation, CellNR, Signal_type) %>% 
          dplyr::summarize(B1int_cellbody_sd = sd(B1int, na.rm=TRUE),
                           across(where(is.numeric), mean, na.rm=TRUE)
          ) %>% ungroup() %>% 
          dplyr::select(c(CellNR, B1int, B1int_cellbody_sd, Gcx)) %>% 
          mutate(ratio_cellbody_average = B1int/ Gcx,
          ) %>% 
          rename(Gcx_cellbody_average = Gcx,
                 B1int_cellbody_average = B1int) %>% 
          mutate(peak_b1int_threshold = B1int_cellbody_average + 4 * B1int_cellbody_sd)
        
      }
      protrusion_data <- protrusion_data %>%
        group_by(Signal_type, CellNR) %>% 
        left_join(cellbody_summary_percell)
      
      
      
      ## copy from here to other Orientations if needed
      protrusion_data <- protrusion_data %>%  
        mutate_at(vars(Fiber_type, Type, Orientation, Xthr_pos), factor)
      
      
      #protrusion_data <- calc_local_segregation_leadingedge(protrusion_data, 3)
      
      protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
      
      
      
      protrusion_data <- protrusion_data %>%
        group_by(number, Position_relative_to_fiber, Signal_type) %>% 
        mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                               ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                               NA),
               ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                         ratio/mean(ratio, na.rm=TRUE),
                                                         NA)
        )
      
      
      
      protrusion_data$Type[protrusion_data$Type=="Focalized" & is.na(protrusion_data$PeakID_perprofile)] <- "Unspecified"
      
      
      protrusion_data <- protrusion_data %>% 
        group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
        filter(Signal_type=="foreground")
      
      
      library(zoo)
      
      
      
      protrusion_data <- protrusion_data %>% 
        group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
        arrange(X) %>% 
        mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
               WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
               Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
               WithinPeak = TRUE
        )
      
      protrusion_data <- protrusion_data %>% 
        group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
        arrange(X) %>% 
        mutate(Clustersize = max(cumsum(WithinPeak==TRUE)*pxsize),
               Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
               Clustersize = ifelse(Clustersize < Clustersize_threshold, NA, Clustersize),
               new_Peak_ID = ifelse(is.na(Clustersize), NA, new_Peak_ID),
               new_Peak_ID = ifelse(max(B1int, na.rm=TRUE) < peak_b1int_threshold, NA, new_Peak_ID),
               Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
               Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
        )
      
      
      
      # determine local B1 integrin and Gcx background for the integrin clusters
      px.clustersize.micron.local = 0.6
      px.clustersize.local = floor(px.clustersize.micron.local/pxsize)
      
      
      protrusion_data <- protrusion_data %>%
        group_by(Orientation, number, new_Peak_ID) %>%
        arrange(X) %>% 
        mutate(min_B1int_local = ifelse(X==min(X) | X==max(X), B1int, NA),
               min_Gcx_local =   ifelse(X==min(X) | X==max(X), Gcx, NA),
               interval.local = 1
        )%>%
        ungroup() %>% 
        group_by(Orientation, number, Signal_type, new_Peak_ID) %>%
        mutate(B1int_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID), 
                                               zoo::na.approx(min_B1int_local, na.rm=FALSE, rule=2), NA),
               Gcx_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID),
                                             zoo::na.approx(min_Gcx_local, na.rm=FALSE, rule=2), NA)
        )%>% 
        dplyr::select(-c(min_B1int_local, min_Gcx_local, interval.local))
      
      
      
      
      
      ## determine global B1 integrin and Gcx background for the integrin clusters
      
      # WORK HERE WORK HERE WORK HERE
      px.clustersize.micron = 4
      px.clustersize = floor(px.clustersize.micron/pxsize)
      
      protrusion_data <- protrusion_data %>%
        group_by(Orientation, number, Signal_type) %>%
        arrange(X) %>%
        mutate(
          interval.end = ifelse(!is.na(B1int_cluster_local_bg), cumsum(!is.na(B1int_cluster_local_bg) & lead(is.na(B1int_cluster_local_bg))), NA),
          interval = ifelse(!is.na(interval.end) & lead(is.na(interval.end)), NA, interval.end),
          interval = as.factor(interval)
        ) %>% 
        ungroup()
      
      
      # Define the custom function
      run_protrusion_code <- function(data, px_clustersize) {
        result <- data %>% 
          group_by(Orientation, number, Signal_type, interval)  %>% 
          mutate(
            composite.length = ifelse(!is.na(interval), max(row_number()), NA),
            dividecluster = ifelse(!is.na(interval) & composite.length > px_clustersize, TRUE, FALSE),
            clusterX = ifelse(!is.na(interval) & dividecluster, row_number(), NA),
            pastthrescoord = ifelse(!is.na(interval) & dividecluster & clusterX > px_clustersize, TRUE, FALSE),
            maxedge = ifelse(!is.na(interval) & pastthrescoord, max(Peak_position == "Edge"), NA),
            lastedge = ifelse(!is.na(interval) & cumsum(pastthrescoord & Peak_position == "Edge") == maxedge, TRUE, FALSE),
            lastedge = ifelse(lastedge==FALSE & lag(lastedge==TRUE), TRUE, lastedge),
            lastedge = ifelse(lastedge==TRUE & lead(lastedge==TRUE),FALSE, lastedge),
            dividecoord = ifelse(!is.na(interval) & length(unique(pastthrescoord)) > 1 & lastedge, TRUE, FALSE),
            dividecoord = ifelse(is.na(dividecoord), FALSE, dividecoord),
            subinterval = cumsum(dividecoord)
          ) %>% 
          ungroup() %>%
          mutate(interval = ifelse(!is.na(interval), group_indices(., interval, subinterval), NA))
        
        return(result)
      }
      
      repeated_samplingdistance = 10
      
      # Apply the custom function to your data twice
      for(p in 1:repeated_samplingdistance){
        print(p)
        protrusion_data <- run_protrusion_code(protrusion_data, px.clustersize)
      }
      
      
      protrusion_data <- protrusion_data %>% 
        group_by(Orientation, number, Signal_type, interval) %>% 
        arrange(X) %>% 
        mutate(
          min_B1int = ifelse(X==min(X, na.rm=TRUE), B1int, NA),
          min_B1int = ifelse(X==max(X, na.rm=TRUE), B1int, min_B1int),
          # min_B1int = ifelse(dividecoord==TRUE, B1int, min_B1int),
          min_Gcx = ifelse(X==min(X, na.rm=TRUE) | X==max(X, na.rm=TRUE), Gcx, NA),
        )
      
      protrusion_data <- protrusion_data %>% 
        ungroup() %>% 
        group_by(Orientation, number, Signal_type) %>%
        mutate(B1int_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_B1int, na.rm=FALSE, rule=2), NA),
               Gcx_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_Gcx, na.rm=FALSE, rule=2), NA),
               B1int_cluster_global_bg = ifelse(B1int_cluster_local_bg < B1int_cluster_global_bg, B1int_cluster_local_bg, B1int_cluster_global_bg)
        )
      
      
      
      measuring_bandwidth_edge = 4
      measuring_bandwidth_peak = 2
      
      
      protrusion_data <- protrusion_data %>% 
        group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
        arrange(X) %>% 
        mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
               WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
               Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
               Peak_position  = ifelse(is.na(new_Peak_ID), NA, Peak_position),
               X.cluster =  ifelse(!is.na(new_Peak_ID), X - min(X), NA),
               Peak_position = ifelse(X.cluster <=(measuring_bandwidth_edge-1) | X.cluster >= (max(X.cluster)-(measuring_bandwidth_edge-1)), "Edge", "Transition"),
               Peak_position = ifelse(!is.na(Peaklocation_X), "Peak", Peak_position)
        )
      
      for(g in 1:(measuring_bandwidth_peak-1)){
        print(g)
        protrusion_data <- protrusion_data %>% 
          group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
          arrange(X) %>% 
          mutate(Peak_position = ifelse(lag(Peak_position=="Peak"), "Peak", Peak_position),
                 Peak_position = ifelse(lead(Peak_position=="Peak"), "Peak", Peak_position)
          )
      }
      
      protrusion_data <- protrusion_data %>% 
        group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
        arrange(X) %>% 
        mutate(Peak_position = ifelse(is.na(Peak_position), "Edge", Peak_position),
               Peak_position_number = ifelse(Peak_position=="Edge" & X.cluster < mean(X.cluster[!is.na(Peaklocation_X)]), 1, 2)
        )
      
      
      protrusion_data <- protrusion_data %>% 
        ungroup() %>% 
        group_by(number, Signal_type, LocationNR, interval) %>% 
        arrange(X) %>% 
        mutate(X.composite = ifelse(!is.na(new_Peak_ID), X - min(X), NA),
               Edgeclusternr = ifelse(!is.na(new_Peak_ID), cumsum(Peak_position=="Edge" & lead(Peak_position!="Edge")),NA),
               Edgeclusternr = ifelse(!is.na(Edgeclusternr) & lag(Edgeclusternr)!=Edgeclusternr, lag(Edgeclusternr), Edgeclusternr),
               Peak_position_composite = ifelse(!is.na(new_Peak_ID) & Peak_position=="Edge" & Edgeclusternr==0,"Edge", "Transition"),
               Peak_position_composite = ifelse(!is.na(new_Peak_ID) & Peak_position=="Edge" & Edgeclusternr==max(Edgeclusternr, na.rm=TRUE),"Edge", Peak_position_composite),
               Peak_position_composite = ifelse(Peak_position=="Peak", "Peak", Peak_position_composite),
               Peak_position_composite = ifelse(is.na(new_Peak_ID), NA, Peak_position_composite),
               Peak_position_composite_number = ifelse(Peak_position_composite=="Edge" & X.composite < (px.clustersize/2), 1, 2)
        )                                                 
      
      
      protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
      
      
      protrusion_data <- protrusion_data %>%  
        group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>% 
        mutate(unique_Peak_ID = cur_group_id())
      
      protrusion_data <- protrusion_data %>% ungroup() %>% 
        group_by(Orientation, number, Signal_type, interval) %>% 
        mutate(  Gcx_adjacent_composite = mean(Gcx[Peak_position_composite=="Edge"], na.rm=TRUE),
                 B1int_adjacent_composite = mean(B1int[Peak_position_composite=="Edge"], na.rm=TRUE)
        )
      
      protrusion_data <- protrusion_data %>% 
        group_by(Orientation, number, Signal_type, new_Peak_ID, unique_Peak_ID) %>% 
        mutate(B1int_global_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_global_bg, NA),
               B1int_local_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_local_bg, NA),
               Gcx_global_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_global_bg, NA),
               Gcx_local_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_local_bg, NA),
               Gcx_enrichment = ifelse(X==Peaklocation_X, mean(Gcx[Peak_position=="Peak"], na.rm=TRUE) / mean(Gcx[Peak_position=="Edge"], na.rm=TRUE), NA),
               B1int_enrichment = ifelse(X==Peaklocation_X, mean(B1int[Peak_position=="Peak"], na.rm=TRUE) / mean(B1int[Peak_position=="Edge"], na.rm=TRUE), NA),
               B1int_enrichment_composite = ifelse(X==Peaklocation_X, mean(B1int[Peak_position=="Peak"], na.rm=TRUE) / mean(B1int_adjacent_composite, na.rm=TRUE)),
               Gcx_enrichment_composite = ifelse(X==Peaklocation_X, mean(Gcx[Peak_position=="Peak"], na.rm=TRUE) / mean(Gcx_adjacent_composite, na.rm=TRUE))
        )
      
      
      protrusion_data <- protrusion_data %>% 
        group_by(Orientation, number, Signal_type) %>% 
        mutate(X_relative = X/ max(row_number(), na.rm=TRUE),
               Length_label= ifelse(Orientation=="Bleb" & X_relative> (3/12*max(X_relative)) & X_relative<(9/12*(max(X_relative))), "Neck", "Base"),
               Length_label= ifelse(Orientation=="Bleb" & X_relative> (5/12*max(X_relative)) & X_relative<(7/12*(max(X_relative))), "Bleb", Length_label),
               Gcx_bleb_base_normalized = Gcx/ mean(Gcx[Length_label=="Base"]),
               B1int_bleb_base_normalized = B1int / mean(B1int[Length_label=="Base"]),
               Col_bleb_base_normalized = Col/ mean(Col[Length_label=="Base"]),
               Bleb_length = length(Gcx[Length_label=="Bleb"])*pxsize
        ) %>% arrange(X) %>% 
        mutate(Length_classification_nr = cumsum(X>1 &Length_label != lag(Length_label)))
      
      
      protrusion_data_summarized_percluster <- protrusion_data %>% 
        group_by(Orientation, number, Signal_type, new_Peak_ID, 
                 Clustersize) %>% 
        dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Signal_type=="foreground") 
      
      
      protrusion_data_summarized_perROI <- protrusion_data %>% 
        group_by(Orientation, number, Signal_type, Length_label) %>% 
        dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
      
      protrusion_data_summarized_percell <- protrusion_data %>% 
        group_by(Orientation, CellNR, Type, Fiber_type, Signal_type, Xthr_pos) %>% 
        dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
      
      subset_protrusion <- protrusion_data %>%  filter(is.na(Peak_position)) 
      
      summary_sds_protrusion <- subset_protrusion %>%  group_by(Orientation, number) %>%  summarize(B1int_sd= sd(B1int, na.rm=TRUE))
      
      
   ## blebs - make combined lineplot-----
      
      
      sub <- protrusion_data
      
      #sub <- subset(protrusion_data, number==37)
      
      legendsize=0.3
      
      
      Zdiff = 0.25
      
      sub_summary <- sub %>%
        filter(Signal_type=="foreground") %>%
        summarize(mean_B1int = mean(B1int_smooth),
                  sd_B1int = sd(B1int_smooth),
                  iv = mean_B1int+sd_B1int,
                  thr = mean_B1int + sd_B1int*Zdiff)
      
      
      sub$Signal_type <- as.factor(sub$Signal_type)
      sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
      
      sub <- sub %>%  filter(Signal_type=="foreground")
      
      lineplot <- ggplot(data=sub, aes(x=X_relative, shape=NULL))+
        geom_line(aes(y=B1int_bleb_base_normalized, group=number), color="magenta", alpha=0.5)+
        geom_smooth(aes(y=B1int_bleb_base_normalized), color="magenta", alpha=1)+
          theme_classic()
      
      lineplot2 <- ggplot(data=sub, aes(x=X_relative, shape=NULL))+
        geom_smooth(aes(y=Gcx_bleb_base_normalized), color="gold3", alpha=1)+
        geom_line(aes(y=Gcx_bleb_base_normalized, group=number), color="gold3", alpha=0.5)+
        theme_classic()
      
      rect_data <- sub %>% ungroup() %>% 
        filter(Signal_type=="foreground") %>%
        group_by(Length_label, Length_classification_nr) %>% 
        summarize(
          xmin = min(X_relative), xmax = max(X_relative), ymin = -Inf, ymax = Inf)
      

      cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color=Length_label), fill=alpha("black",0), color=alpha("black", 0.3),
                                      inherit.aes = FALSE, alpha = 0, size=0.5)
  
      
      lineplot <- lineplot + cluster.rectangles +
        labs(linetype="", x="Relative distance (fraction)", y = "Normalized B1int")+ theme(legend.position="top")
      
      
      lineplot2 <- lineplot2 + cluster.rectangles +
        labs(linetype="", x="Relative distance (fraction)", y = "Normalized Gcx")+ theme(legend.position="top")+
        theme(legend.position = "none")
      
      
      
      lineplot <- lineplot + 
        theme(legend.position="top")+
        theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
              legend.key.height = unit(legendsize, 'cm'), #change legend key height
              legend.key.width = unit(legendsize*1.5, 'cm'))
      
      
      
      
      lineplot.bleb <-  plot_grid(lineplot, lineplot2, ncol=1, vjust=0, rel_heights = c(0.6,0.4))+
        theme(plot.margin = margin(t = 5, unit = "mm"))
      
      print(lineplot.bleb)
      
      
   ## blebs - save all lineplots to depict Gcx and B1 integrin intensity profiles------
 
 unique_numbers = unique(protrusion_data$number)
 
 for(u in 1:length(unique_numbers)){
   print(u)
   unique_number <- unique_numbers[u]
   sub <- subset(protrusion_data, number==unique_number)
   
   #sub <- subset(protrusion_data, number==37)
   
   legendsize=0.3
   
   
   Zdiff = 0.25
   
   sub_summary <- sub %>%
     filter(Signal_type=="foreground") %>%
     summarize(mean_B1int = mean(B1int_smooth),
               sd_B1int = sd(B1int_smooth),
               iv = mean_B1int+sd_B1int,
               thr = mean_B1int + sd_B1int*Zdiff)
   
   
   sub$Signal_type <- as.factor(sub$Signal_type)
   sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
   
   sub <- sub %>%  filter(Signal_type=="foreground")
  
   
   
   lineplot <- ggplot(data=sub, aes(x=X_micron, shape=NULL))+
     geom_line(aes(y=B1int_bleb_base_normalized), color="magenta", alpha=1)+
     geom_line(aes(y=Gcx_bleb_base_normalized), color="gold3", alpha=1)+
     geom_line(aes(y=Col_bleb_base_normalized), color="cyan4", alpha=1)+
     geom_point(data=subset(sub, !is.na(new_Peak_ID)& !is.na(Peaklocation_X)), aes(y=B1int_bleb_base_normalized))+
     theme_classic()
   
   #B1int_thr <- geom_hline(yintercept=unique(sub$peak_b1int_threshold), color="magenta", linetype="dashed")
   
   rect_data <- sub %>% ungroup() %>% 
     filter(Signal_type=="foreground") %>% filter(!is.na(Length_classification_nr)) %>% 
     group_by(Length_label, Length_classification_nr) %>% 
     summarize(
               xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf)
   

   cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, linetype=Length_label), fill=alpha("black",0), color=alpha("black", 0.3),
                                   inherit.aes = FALSE, alpha = 0, size=0.5)
   

   
   lineplot <- lineplot + cluster.rectangles 
     labs(linetype="", x="Distance (\u03BCm)", y = " Intensity")+ theme(legend.position="top")

   
   
   
   lineplot <- lineplot + 
     theme(legend.position="top")+
     theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
           legend.key.height = unit(legendsize, 'cm'), #change legend key height
           legend.key.width = unit(legendsize*1.5, 'cm'))
   
   
   
   lineplot.bleb <-  plot_grid(lineplot, labels=paste0(unique(sub$Orientation), " - Number: ",unique_number), vjust=0)+
     theme(plot.margin = margin(t = 5, unit = "mm"))
   
   print(lineplot.bleb)
   
   pdf(paste0("Lineplot_", unique_number, ".pdf"), width=8, height=5)
   print(lineplot.bleb)
   dev.off()
   
 }
 
 
 
 
 
 
 
 
 
 
   ## blebs - plot of Gcx and B1int bleb vs tip-----
 
 bleb_data <- protrusion_data %>%  group_by(number, Orientation, Signal_type, Length_label) %>%  
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
      
      
 bleb_cluster_correction = FALSE
 
 if(bleb_cluster_correction==TRUE){
   bleb_data <- bleb_data %>%  group_by(number, Orientation, Signal_type, Length_label) %>% 
     mutate(filterout = ifelse(is.na(Peaklocation_X) & Length_label=="Bleb", TRUE, FALSE)) %>% 
     group_by(number, Orientation, Signal_type) %>% 
     filter((any(filterout==TRUE)))
 }
 
bleb_data <- bleb_data %>%  filter(Length_label!="Base")
 
 bleb.lengthlabel.vs.gcx <- custom_plot_with_statistics(data=bleb_data, x_var="Length_label", alpha=0.5, size=2,
                                                   y_var="Gcx_bleb_base_normalized", group_vars=NULL, jitter=TRUE, annotate_signif = TRUE, fill_var=NULL,
                                                   annotate_only_signif = TRUE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 
 bleb.lengthlabel.vs.gcx.plot <- bleb.lengthlabel.vs.gcx$plot
 bleb.lengthlabel.vs.gcx.plot <- bleb.lengthlabel.vs.gcx.plot + theme_classic()+
   labs(y="Normalized Gcx", x="")
 print(bleb.lengthlabel.vs.gcx.plot)
 
 
 bleb.lengthlabel.vs.gcx.plot <-  plot_grid(bleb.lengthlabel.vs.gcx.plot, labels="d", vjust=0)+
   theme(plot.margin = margin(t = 5, unit = "mm"))
 
 print(ggdraw(bleb.lengthlabel.vs.gcx.plot))
 
 
 bleb.lengthlabel.vs.b1int <- custom_plot_with_statistics(data=bleb_data, x_var="Length_label", alpha=0.5, size=2,
                                                        y_var="B1int_bleb_base_normalized", group_vars=NULL, jitter=TRUE, annotate_signif = TRUE, fill_var=NULL,
                                                        annotate_only_signif = TRUE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 
 bleb.lengthlabel.vs.b1int.plot <- bleb.lengthlabel.vs.b1int$plot
 bleb.lengthlabel.vs.b1int.plot <- bleb.lengthlabel.vs.b1int.plot + theme_classic()+
   labs(y="Normalized B1int", x="")
 print(bleb.lengthlabel.vs.b1int.plot)
 
 
 bleb.lengthlabel.vs.b1int.plot <-  plot_grid(bleb.lengthlabel.vs.b1int.plot, labels="c", vjust=0)+
   theme(plot.margin = margin(t = 5, unit = "mm"))
 
 print(ggdraw(bleb.lengthlabel.vs.b1int.plot))
 
 bleb.graphs <- plot_grid(bleb.lengthlabel.vs.b1int.plot, bleb.lengthlabel.vs.gcx.plot, nrow=1)
 
 print(bleb.graphs)
 

   ## blebs - col vs b1int ----
      scatterplot <- custom_lineplot_with_rsquared(subset(bleb_data, Length_label=="Neck"), x_var=Col, 
                                                   y_var=B1int, annotate_signif = TRUE, regression_type = "linear", group_vars="Length_label")
 
 print(scatterplot)
 
 
 ###-----
 ### ----
### -----
 ## old plots etcs------
 
 
 
 
 
 ## cell body calculations-----
 cellbody_data <- subset(B1int_GCX_ratio, Orientation=="Cellbody" | Orientation=="Leading edge protrusion")
 
 
 protrusion_data <- cellbody_data
 
 
 if(any(protrusion_data$Orientation=="Cellbody")){
   cellbody_summary_percell <- protrusion_data  %>% filter(Orientation=="Cellbody" & Signal_type=="foreground") %>%  
     group_by(Orientation, CellNR, Signal_type) %>% 
     dplyr::summarize(B1int_cellbody_sd = sd(B1int, na.rm=TRUE),
                      across(where(is.numeric), mean, na.rm=TRUE)
     ) %>% ungroup() %>% 
     dplyr::select(c(CellNR, B1int, B1int_cellbody_sd, Gcx)) %>% 
     mutate(ratio_cellbody_average = B1int/ Gcx,
     ) %>% 
     rename(Gcx_cellbody_average = Gcx,
            B1int_cellbody_average = B1int) %>% 
     mutate(peak_b1int_threshold = B1int_cellbody_average + 4 * B1int_cellbody_sd)
   
 }
 protrusion_data <- protrusion_data %>%
   group_by(Signal_type, CellNR) %>% 
   left_join(cellbody_summary_percell)
 
 
 
 ## copy from here to other Orientations if needed
 protrusion_data <- protrusion_data %>%  
   mutate_at(vars(Fiber_type, Type, Orientation, Xthr_pos), factor)
 
 
 #protrusion_data <- calc_local_segregation_leadingedge(protrusion_data, 3)
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 
 protrusion_data <- protrusion_data %>%
   group_by(number, Position_relative_to_fiber, Signal_type) %>% 
   mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                          ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                          NA),
          ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                    ratio/mean(ratio, na.rm=TRUE),
                                                    NA)
   )
 
 
 
 protrusion_data$Type[protrusion_data$Type=="Focalized" & is.na(protrusion_data$PeakID_perprofile)] <- "Unspecified"
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   filter(Signal_type=="foreground")
 
 
 library(zoo)
 
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          WithinPeak = TRUE
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   arrange(X) %>% 
   mutate(Clustersize = max(cumsum(WithinPeak==TRUE)*pxsize),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Clustersize = ifelse(Clustersize < Clustersize_threshold, NA, Clustersize),
          new_Peak_ID = ifelse(is.na(Clustersize), NA, new_Peak_ID),
          new_Peak_ID = ifelse(max(B1int, na.rm=TRUE) < peak_b1int_threshold, NA, new_Peak_ID),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
   )
 
 
 
 # determine local B1 integrin and Gcx background for the integrin clusters
 px.clustersize.micron.local = 0.6
 px.clustersize.local = floor(px.clustersize.micron.local/pxsize)
 
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(min_B1int_local = ifelse(X==min(X) | X==max(X), B1int, NA),
          min_Gcx_local =   ifelse(X==min(X) | X==max(X), Gcx, NA),
          interval.local = 1
   )%>%
   ungroup() %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>%
   mutate(B1int_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID), 
                                          zoo::na.approx(min_B1int_local, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID),
                                        zoo::na.approx(min_Gcx_local, na.rm=FALSE, rule=2), NA)
   )%>% 
   dplyr::select(-c(min_B1int_local, min_Gcx_local, interval.local))
 
 
 
 
 
 ## determine global B1 integrin and Gcx background for the integrin clusters
 
 # WORK HERE WORK HERE WORK HERE
 px.clustersize.micron = 4
 px.clustersize = floor(px.clustersize.micron/pxsize)
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, Signal_type) %>%
   arrange(X) %>%
   mutate(
     interval.end = ifelse(!is.na(B1int_cluster_local_bg), cumsum(!is.na(B1int_cluster_local_bg) & lead(is.na(B1int_cluster_local_bg))), NA),
     interval = ifelse(!is.na(interval.end) & lead(is.na(interval.end)), NA, interval.end),
     interval = as.factor(interval)
   ) %>% 
   ungroup()
 
 
 # Define the custom function
 run_protrusion_code <- function(data, px_clustersize) {
   result <- data %>% 
     group_by(Orientation, number, Signal_type, interval)  %>% 
     mutate(
       composite.length = ifelse(!is.na(interval), max(row_number()), NA),
       dividecluster = ifelse(!is.na(interval) & composite.length > px_clustersize, TRUE, FALSE),
       clusterX = ifelse(!is.na(interval) & dividecluster, row_number(), NA),
       pastthrescoord = ifelse(!is.na(interval) & dividecluster & clusterX > px_clustersize, TRUE, FALSE),
       maxedge = ifelse(!is.na(interval) & pastthrescoord, max(Peak_position == "Edge"), NA),
       lastedge = ifelse(!is.na(interval) & cumsum(pastthrescoord & Peak_position == "Edge") == maxedge, TRUE, FALSE),
       lastedge = ifelse(lastedge==FALSE & lag(lastedge==TRUE), TRUE, lastedge),
       lastedge = ifelse(lastedge==TRUE & lead(lastedge==TRUE),FALSE, lastedge),
       dividecoord = ifelse(!is.na(interval) & length(unique(pastthrescoord)) > 1 & lastedge, TRUE, FALSE),
       dividecoord = ifelse(is.na(dividecoord), FALSE, dividecoord),
       subinterval = cumsum(dividecoord)
     ) %>% 
     ungroup() %>%
     mutate(interval = ifelse(!is.na(interval), group_indices(., interval, subinterval), NA))
   
   return(result)
 }
 
 repeated_samplingdistance = 10
 
 # Apply the custom function to your data twice
 for(p in 1:repeated_samplingdistance){
   print(p)
   protrusion_data <- run_protrusion_code(protrusion_data, px.clustersize)
 }
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, interval) %>% 
   arrange(X) %>% 
   mutate(
     min_B1int = ifelse(X==min(X, na.rm=TRUE), B1int, NA),
     min_B1int = ifelse(X==max(X, na.rm=TRUE), B1int, min_B1int),
     # min_B1int = ifelse(dividecoord==TRUE, B1int, min_B1int),
     min_Gcx = ifelse(X==min(X, na.rm=TRUE) | X==max(X, na.rm=TRUE), Gcx, NA),
   )
 
 protrusion_data <- protrusion_data %>% 
   ungroup() %>% 
   group_by(Orientation, number, Signal_type) %>%
   mutate(B1int_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_B1int, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_Gcx, na.rm=FALSE, rule=2), NA),
          B1int_cluster_global_bg = ifelse(B1int_cluster_local_bg < B1int_cluster_global_bg, B1int_cluster_local_bg, B1int_cluster_global_bg)
   )
 
 
 
 measuring_bandwidth_edge = 4
 measuring_bandwidth_peak = 2
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          Peak_position  = ifelse(is.na(new_Peak_ID), NA, Peak_position),
          X.cluster =  ifelse(!is.na(new_Peak_ID), X - min(X), NA),
          Peak_position = ifelse(X.cluster <=(measuring_bandwidth_edge-1) | X.cluster >= (max(X.cluster)-(measuring_bandwidth_edge-1)), "Edge", "Transition"),
          Peak_position = ifelse(!is.na(Peaklocation_X), "Peak", Peak_position)
   )
 
 for(g in 1:(measuring_bandwidth_peak-1)){
   print(g)
   protrusion_data <- protrusion_data %>% 
     group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
     arrange(X) %>% 
     mutate(Peak_position = ifelse(lag(Peak_position=="Peak"), "Peak", Peak_position),
            Peak_position = ifelse(lead(Peak_position=="Peak"), "Peak", Peak_position)
     )
 }
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(Peak_position = ifelse(is.na(Peak_position), "Edge", Peak_position),
          Peak_position_number = ifelse(Peak_position=="Edge" & X.cluster < mean(X.cluster[!is.na(Peaklocation_X)]), 1, 2)
   )
 
 
 protrusion_data <- protrusion_data %>% 
   ungroup() %>% 
   group_by(number, Signal_type, LocationNR, interval) %>% 
   arrange(X) %>% 
   mutate(X.composite = ifelse(!is.na(new_Peak_ID), X - min(X), NA),
          Edgeclusternr = ifelse(!is.na(new_Peak_ID), cumsum(Peak_position=="Edge" & lead(Peak_position!="Edge")),NA),
          Edgeclusternr = ifelse(!is.na(Edgeclusternr) & lag(Edgeclusternr)!=Edgeclusternr, lag(Edgeclusternr), Edgeclusternr),
          Peak_position_composite = ifelse(!is.na(new_Peak_ID) & Peak_position=="Edge" & Edgeclusternr==0,"Edge", "Transition"),
          Peak_position_composite = ifelse(!is.na(new_Peak_ID) & Peak_position=="Edge" & Edgeclusternr==max(Edgeclusternr, na.rm=TRUE),"Edge", Peak_position_composite),
          Peak_position_composite = ifelse(Peak_position=="Peak", "Peak", Peak_position_composite),
          Peak_position_composite = ifelse(is.na(new_Peak_ID), NA, Peak_position_composite),
          Peak_position_composite_number = ifelse(Peak_position_composite=="Edge" & X.composite < (px.clustersize/2), 1, 2)
   )                                                 
 
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 protrusion_data <- protrusion_data %>%  
   group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>% 
   mutate(unique_Peak_ID = cur_group_id())
 
 protrusion_data <- protrusion_data %>% ungroup() %>% 
   group_by(Orientation, number, Signal_type, interval) %>% 
   mutate(  Gcx_adjacent_composite = mean(Gcx[Peak_position_composite=="Edge"], na.rm=TRUE),
            B1int_adjacent_composite = mean(B1int[Peak_position_composite=="Edge"], na.rm=TRUE)
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, unique_Peak_ID) %>% 
   mutate(B1int_global_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_global_bg, NA),
          B1int_local_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_local_bg, NA),
          Gcx_global_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_global_bg, NA),
          Gcx_local_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_local_bg, NA),
          Gcx_enrichment = ifelse(X==Peaklocation_X, mean(Gcx[Peak_position=="Peak"], na.rm=TRUE) / mean(Gcx[Peak_position=="Edge"], na.rm=TRUE), NA),
          B1int_enrichment = ifelse(X==Peaklocation_X, mean(B1int[Peak_position=="Peak"], na.rm=TRUE) / mean(B1int[Peak_position=="Edge"], na.rm=TRUE), NA),
          B1int_enrichment_composite = ifelse(X==Peaklocation_X, mean(B1int[Peak_position=="Peak"], na.rm=TRUE) / mean(B1int_adjacent_composite, na.rm=TRUE)),
          Gcx_enrichment_composite = ifelse(X==Peaklocation_X, mean(Gcx[Peak_position=="Peak"], na.rm=TRUE) / mean(Gcx_adjacent_composite, na.rm=TRUE))
   )
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type) %>% 
   mutate(X_relative = X/ max(row_number(), na.rm=TRUE),
          Bleb_length_classification= ifelse(X_relative> (3/12*max(X_relative)) & X_relative<(9/12*(max(X_relative))), "Neck", "Base"),
          Bleb_length_classification= ifelse(X_relative> (5/12*max(X_relative)) & X_relative<(7/12*(max(X_relative))), "Bleb", Bleb_length_classification),
          Gcx_bleb_base_normalized = Gcx/ mean(Gcx[Bleb_length_classification=="Base"]),
          B1int_bleb_base_normalized = B1int / mean(B1int[Bleb_length_classification=="Base"])
   ) %>% arrange(X) %>% 
   mutate(Bleb_classification_nr = cumsum(X>1 &Bleb_length_classification != lag(Bleb_length_classification)))
 
 
 
 protrusion_data_summarized_percluster <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, 
            Clustersize) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Signal_type=="foreground") 
 
 
 protrusion_data_summarized_perROI <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, Length_label) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 protrusion_data_summarized_percell <- protrusion_data %>% 
   group_by(Orientation, CellNR, Type, Fiber_type, Signal_type, Xthr_pos) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 subset_protrusion <- protrusion_data %>%  filter(is.na(Peak_position)) 
 
 summary_sds_protrusion <- subset_protrusion %>%  group_by(Orientation, number) %>%  summarize(B1int_sd= sd(B1int, na.rm=TRUE))
 
 ## copy to here for other Orientations
 
 
 cellbody_data <- protrusion_data
 
 ## cell body - example plot to depict Gcx and B1 integrin intensity profiles-----
 
 sub <- subset(protrusion_data, CellNR==7 & LocationNR==1 & RoiNR==1 & Orientation=="Cellbody")
 
 #sub <- subset(protrusion_data, number==37)
 
 legendsize=0.3
 
 
 Zdiff = 0.25
 
 sub_summary <- sub %>%
   filter(Signal_type=="foreground") %>%
   summarize(mean_B1int = mean(B1int_smooth),
             sd_B1int = sd(B1int_smooth),
             iv = mean_B1int+sd_B1int,
             thr = mean_B1int + sd_B1int*Zdiff)
 
 
 sub$Signal_type <- as.factor(sub$Signal_type)
 sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
 
 sub <- sub %>%  filter(Signal_type=="foreground")
 
 lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int))+
   geom_line( color="magenta", alpha=1)+
   theme_classic()
 
 limits_conversion <- max(sub$B1int)/ max(sub$Gcx)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx*limits_conversion), color="gold3", alpha=1)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_global_bg), linetype="dotted", color="magenta", alpha=1)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_local_bg), linetype="dashed", color="magenta", alpha=1)
 
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_global_bg*limits_conversion),linetype="dotted", color="gold3", alpha=1)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_local_bg*limits_conversion),linetype="dashed", color="gold3", alpha=1)
 
 
 
 
 print(lineplot)
 
 
 rect_data <- sub %>%
   filter(Signal_type=="foreground") %>%
   group_by(new_Peak_ID) %>%
   filter(!is.na(Clustersize)) %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
             xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize
   )
 
 lineplot <- lineplot + geom_point(data=subset(sub, !is.na(Peaklocation_X)), size=1.5, stroke=0)
 
 lineplot <- lineplot+
   geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", fill="black", color=alpha("blue", 0),
             inherit.aes = FALSE, alpha = 0, size=0.5)
 
 lineplot <- lineplot+
   # geom_vline(data=rect_data,  linetype="dashed", aes(xintercept=xintercept), alpha=0.5)+
   scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
   scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
   scale_y_continuous(
     # Features of the first axis
     name = "B1 integrin (x1000)",
     labels= scales::number_format(scale=1e-3),
     # Add a second axis and specify its features
     sec.axis = sec_axis(~./limits_conversion, name="Gcx (x1000)",
                         labels= scales::number_format(scale=1e-3))
   )
 
 rect_data2 <- sub %>% ungroup() %>% 
   filter(Signal_type=="foreground") %>%
   group_by(Length_label) %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf) %>%  filter(Length_label!="Intermediate zone")
 
 lineplot <- lineplot + 
   geom_rect(data=rect_data2, aes(xmin= xmin, xmax=xmax, ymin=ymin, ymax=ymax, linetype=Length_label), fill=NA, color="black",
             size=0.5, inherit.aes = FALSE)+
   scale_linetype_manual(values=c("solid", "dashed"))+
   labs(linetype="", x="Distance (\u03BCm)")+
   theme(legend.position="top")+
   theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
         legend.key.height = unit(legendsize, 'cm'), #change legend key height
         legend.key.width = unit(legendsize*1.5, 'cm'))
 
 
 rect_data3 <- sub %>%  ungroup() %>% 
   filter(Signal_type=="foreground") %>% 
   group_by(CellNR) %>% 
   summarize(peak_b1int_threshold = mean(peak_b1int_threshold))
 
 lineplot <- lineplot + geom_hline(data=rect_data3, aes(yintercept=peak_b1int_threshold), linetype="dotted", inherit.aes=FALSE)
 
 lineplot <- lineplot
 print(lineplot)
 
 lineplot_cellbody <- lineplot
 
 
 ## cell body -  save all lineplots to depict Gcx and B1 integrin intensity profiles------
 
 unique_numbers = unique(protrusion_data$number)
 
 for(u in 1:length(unique_numbers)){
   print(u)
   unique_number <- unique_numbers[u]
   sub <- subset(protrusion_data, number==unique_number)
   
   #sub <- subset(protrusion_data, number==37)
   
   legendsize=0.3
   
   
   Zdiff = 0.25
   
   sub_summary <- sub %>%
     filter(Signal_type=="foreground") %>%
     summarize(mean_B1int = mean(B1int_smooth),
               sd_B1int = sd(B1int_smooth),
               iv = mean_B1int+sd_B1int,
               thr = mean_B1int + sd_B1int*Zdiff)
   
   
   sub$Signal_type <- as.factor(sub$Signal_type)
   sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
   
   sub <- sub %>%  filter(Signal_type=="foreground")
   
   
   
   lineplot <- ggplot(data=sub, aes(x=X_micron, shape=NULL))+
     geom_line(aes(y=B1int), color="magenta", alpha=1)+
     geom_line(aes(y=Gcx), color="gold3", alpha=1)+
     geom_point(data=subset(sub, !is.na(new_Peak_ID)& !is.na(Peaklocation_X)), aes(y=B1int))+
     theme_classic()
   
   B1int_thr <- geom_hline(yintercept=unique(sub$peak_b1int_threshold), color="magenta", linetype="dashed")
   
   rect_data <- sub %>% ungroup() %>% 
     filter(Signal_type=="foreground") %>% filter(!is.na(Bleb_classification_nr)) %>% 
     group_by(Bleb_length_classification, Bleb_classification_nr) %>% 
     summarize(
       xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf)
   
   if(as.character(unique(sub$Orientation))=="Bleb"){
     cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, linetype=Bleb_length_classification), fill=alpha("black",0), color=alpha("black", 0.3),
                                     inherit.aes = FALSE, alpha = 0, size=0.5)
     lineplot <- lineplot + cluster.rectangles
   }
   
   
   lineplot <- lineplot  + B1int_thr
   labs(linetype="", x="Distance (\u03BCm)", y = " Intensity")+ theme(legend.position="top")
   
   
   
   
   lineplot <- lineplot + 
     theme(legend.position="top")+
     theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
           legend.key.height = unit(legendsize, 'cm'), #change legend key height
           legend.key.width = unit(legendsize*1.5, 'cm'))
   
   
   
   lineplot.bleb <-  plot_grid(lineplot, labels=paste0(unique(sub$Orientation), " - Number: ",unique_number), vjust=0)+
     theme(plot.margin = margin(t = 5, unit = "mm"))
   
   print(lineplot.bleb)
   
   pdf(paste0("Lineplot_", unique_number, ".pdf"), width=8, height=5)
   print(lineplot.bleb)
   dev.off()
   
 }
 
 
 
 
 
 
 
 
 
 
 ## cell body - plot Gcx. vs Gcx global and local background-----
 
 cellbody_data <- cellbody_data %>%  mutate(Gcx_global_enrichment_allcoordinates = Gcx / Gcx_cluster_global_bg,
                                            Gcx_local_enrichment_allcoordinates = Gcx/ Gcx_cluster_local_bg)
 
 cellbody_data <- cellbody_data %>% mutate(Length_label = ifelse(Orientation=="Cellbody","Contact-free", Length_label))
 
 cellbody_data <- cellbody_data %>%  filter(!is.na(Length_label) & !is.na(Orientation))
 cellbody_data$Length_label <- as.factor(cellbody_data$Length_label)
 cellbody_data$Orientation <- as.factor(cellbody_data$Orientation)
 
 Gcx.vs.Gcx.globalbg.plot <- custom_plot_with_statistics(data=cellbody_data, x_var="Orientation", alpha=0.5, size=2,
                                                         y_var="Gcx_global_enrichment_allcoordinates", group_vars=NULL, jitter=TRUE, annotate_signif = TRUE,
                                                         annotate_only_signif = TRUE, log=FALSE, y_adjust_min=NULL, y_adjust_max=NULL, plot_mean=FALSE, violin=TRUE)
 
 
 Gcx.vs.Gcx.globalbg.plot <-Gcx.vs.Gcx.globalbg.plot$plot + theme_classic() + labs(y= "Multi cluster \n Gcx enrichment", x="Region")
 print(Gcx.vs.Gcx.globalbg.plot)
 
 cellbody.gcx.enrichment.global <- median(cellbody_data$Gcx_global_enrichment_allcoordinates[cellbody_data$Orientation=="Cellbody"])
 protrusion.gcx.enrichment.global <- median(cellbody_data$Gcx_global_enrichment_allcoordinates[cellbody_data$Orientation=="Leading edge protrusion"])
 
 cellbody.gcx.enrichment.local <- median(cellbody_data$Gcx_local_enrichment_allcoordinates[cellbody_data$Orientation=="Cellbody"])
 protrusion.gcx.enrichment.local <- median(cellbody_data$Gcx_local_enrichment_allcoordinates[cellbody_data$Orientation=="Leading edge protrusion"])
 
 
 ## cell body - plot paired Gcx global_bg vs signal, for cell body and protrusion ------
 
 
 
 
 
 signal <-    cellbody_data  %>%  ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Gcx) %>%  mutate(Gcx.signal = "Signal", identifier = row_number()) 
 global.bg <- cellbody_data %>% ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Gcx_cluster_global_bg) %>% mutate(Gcx.signal = "Global_bg", identifier = row_number()) %>%  rename(Gcx = Gcx_cluster_global_bg)
 local.bg <- cellbody_data %>% ungroup () %>%  dplyr::select(Orientation, number, X, Signal_type, Gcx_cluster_local_bg) %>% mutate(Gcx.signal = "Local_bg", identifier = row_number()) %>%  rename(Gcx = Gcx_cluster_local_bg)
 
 
 signal <- rbind(signal, global.bg, local.bg)
 
 
 pool_pixels = TRUE
 grouped_pixels = 5
 
 
 signal_pooled <- signal %>%
   group_by(Orientation, number, Signal_type, Gcx.signal) %>% arrange(X) %>% 
   mutate(X_group = ceiling(row_number()/grouped_pixels))
 
 
 
 signal_pooled <- signal_pooled %>% group_by(Orientation, number, Signal_type, X_group, Gcx.signal) %>% 
   summarise(
     across(where(is.numeric), ~ mean(.), .names = "{.col}")
   ) %>%
   ungroup()
 
 
 if(pool_pixels==TRUE){
   signal_unpooled <- signal
   signal <- signal_pooled
 }
 
 
 
 signal$random_numbering <- runif(nrow(signal), min = 1, max = 100)
 signal$random_numbering <- as.factor(signal$random_numbering)
 
 subset <- signal %>%  subset(Gcx.signal != "Local_bg")
 
 connected.signalvsbg <- ggplot(data=subset, aes(x = Gcx.signal, y = Gcx), color="black") +
   geom_line(aes(group = identifier), alpha = 0.2) +
   theme_classic() +
   theme(legend.position = "none") +
   labs(y = "Gcx MFI", x = "") +
   scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)+
   facet_wrap(.~Orientation)
 
 
 tukey <-  subset %>%  ungroup() %>%  group_by(Orientation) %>% tukey_hsd(Gcx ~ Gcx.signal) %>% 
   add_significance() %>% 
   add_xy_position() %>% 
   mutate(y.position.log10= log10(y.position))
 
 
 # Add line for comparison
 connected.signalvsbg <- connected.signalvsbg +
   ggsignif::geom_signif(data = tukey, aes(xmin= group1, xmax= group2, y_position=y.position, annotations=p.adj.signif), manual=TRUE)+
   scale_x_discrete(labels=c("Global_bg" = "Multi-peak \n background", "Signal" = "Peak  \nsignal"))
 
 # Print the updated plot
 print(connected.signalvsbg)
 
 
 
 
 
 
 ## OLD APPROACH ---make corrections for retraction fibers------
 
 
 protrusion_data <- subset(B1int_GCX_ratio, Orientation =="Retraction fiber")
 
 
 if(any(protrusion_data$Orientation=="Cellbody")){
   cellbody_summary_percell <- protrusion_data  %>% filter(Orientation=="Cellbody" & Signal_type=="foreground") %>%  
     group_by(Orientation, CellNR, Signal_type) %>% 
     dplyr::summarize(B1int_cellbody_sd = sd(B1int, na.rm=TRUE),
                      across(where(is.numeric), mean, na.rm=TRUE)
     ) %>% ungroup() %>% 
     dplyr::select(c(CellNR, B1int, B1int_cellbody_sd, Gcx)) %>% 
     mutate(ratio_cellbody_average = B1int/ Gcx,
     ) %>% 
     rename(Gcx_cellbody_average = Gcx,
            B1int_cellbody_average = B1int) %>% 
     mutate(peak_b1int_threshold = B1int_cellbody_average + 4 * B1int_cellbody_sd)
   
 }
 protrusion_data <- protrusion_data %>%
   group_by(Signal_type, CellNR) %>% 
   left_join(cellbody_summary_percell)
 
 
 
 ## copy from here to other Orientations if needed
 protrusion_data <- protrusion_data %>%  
   mutate_at(vars(Fiber_type, Type, Orientation, Xthr_pos), factor)
 
 
 #protrusion_data <- calc_local_segregation_leadingedge(protrusion_data, 3)
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 
 protrusion_data <- protrusion_data %>%
   group_by(number, Position_relative_to_fiber, Signal_type) %>% 
   mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                          ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                          NA),
          ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                    ratio/mean(ratio, na.rm=TRUE),
                                                    NA)
   )
 
 
 
 protrusion_data$Type[protrusion_data$Type=="Focalized" & is.na(protrusion_data$PeakID_perprofile)] <- "Unspecified"
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   filter(Signal_type=="foreground")
 
 
 library(zoo)
 
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          WithinPeak = TRUE
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   arrange(X) %>% 
   mutate(Clustersize = max(cumsum(WithinPeak==TRUE)*pxsize),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Clustersize = ifelse(Clustersize < Clustersize_threshold, NA, Clustersize),
          new_Peak_ID = ifelse(is.na(Clustersize), NA, new_Peak_ID),
          new_Peak_ID = ifelse(max(B1int, na.rm=TRUE) < peak_b1int_threshold, NA, new_Peak_ID),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
   )
 
 
 
 # determine local B1 integrin and Gcx background for the integrin clusters
 px.clustersize.micron.local = 0.6
 px.clustersize.local = floor(px.clustersize.micron.local/pxsize)
 
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(min_B1int_local = ifelse(X==min(X) | X==max(X), B1int, NA),
          min_Gcx_local =   ifelse(X==min(X) | X==max(X), Gcx, NA),
          interval.local = 1
   )%>%
   ungroup() %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>%
   mutate(B1int_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID), 
                                          zoo::na.approx(min_B1int_local, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID),
                                        zoo::na.approx(min_Gcx_local, na.rm=FALSE, rule=2), NA)
   )%>% 
   dplyr::select(-c(min_B1int_local, min_Gcx_local, interval.local))
 
 
 
 
 
 ## determine global B1 integrin and Gcx background for the integrin clusters
 
 # WORK HERE WORK HERE WORK HERE
 px.clustersize.micron = 4
 px.clustersize = floor(px.clustersize.micron/pxsize)
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, Signal_type) %>%
   arrange(X) %>%
   mutate(
     interval.end = ifelse(!is.na(B1int_cluster_local_bg), cumsum(!is.na(B1int_cluster_local_bg) & lead(is.na(B1int_cluster_local_bg))), NA),
     interval = ifelse(!is.na(interval.end) & lead(is.na(interval.end)), NA, interval.end),
     interval = as.factor(interval)
   ) %>% 
   ungroup()
 
 
 # Define the custom function
 run_protrusion_code <- function(data, px_clustersize) {
   result <- data %>% 
     group_by(Orientation, number, Signal_type, interval)  %>% 
     mutate(
       composite.length = ifelse(!is.na(interval), max(row_number()), NA),
       dividecluster = ifelse(!is.na(interval) & composite.length > px_clustersize, TRUE, FALSE),
       clusterX = ifelse(!is.na(interval) & dividecluster, row_number(), NA),
       pastthrescoord = ifelse(!is.na(interval) & dividecluster & clusterX > px_clustersize, TRUE, FALSE),
       maxedge = ifelse(!is.na(interval) & pastthrescoord, max(Peak_position == "Edge"), NA),
       lastedge = ifelse(!is.na(interval) & cumsum(pastthrescoord & Peak_position == "Edge") == maxedge, TRUE, FALSE),
       dividecoord = ifelse(!is.na(interval) & length(unique(pastthrescoord)) > 1 & lastedge, TRUE, FALSE),
       dividecoord = ifelse(is.na(dividecoord), FALSE, dividecoord),
       subinterval = cumsum(dividecoord)
     ) %>% 
     ungroup() %>%
     mutate(interval = ifelse(!is.na(interval), group_indices(., interval, subinterval), NA))
   
   return(result)
 }
 
 repeated_samplingdistance = 10
 
 # Apply the custom function to your data twice
 for(p in 1:repeated_samplingdistance){
   print(p)
   protrusion_data <- run_protrusion_code(protrusion_data, px.clustersize)
 }
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, interval) %>% 
   arrange(X) %>% 
   mutate(
     min_B1int = ifelse(X==min(X, na.rm=TRUE), B1int, NA),
     min_B1int = ifelse(X==max(X, na.rm=TRUE), B1int, min_B1int),
     # min_B1int = ifelse(dividecoord==TRUE, B1int, min_B1int),
     min_Gcx = ifelse(X==min(X, na.rm=TRUE) | X==max(X, na.rm=TRUE), Gcx, NA),
   )
 
 protrusion_data <- protrusion_data %>% 
   ungroup() %>% 
   group_by(Orientation, number, Signal_type) %>%
   mutate(B1int_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_B1int, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_Gcx, na.rm=FALSE, rule=2), NA),
          #B1int_cluster_local_bg = ifelse(B1int_cluster_local_bg > B1int_cluster_global_bg, B1int_cluster_local_bg, B1int_cluster_global_bg)
   )
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge")
   )
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 protrusion_data <- protrusion_data %>%  
   group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>% 
   mutate(unique_Peak_ID = cur_group_id())
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, unique_Peak_ID) %>% 
   mutate(B1int_global_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_global_bg, NA),
          B1int_local_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_local_bg, NA),
          Gcx_global_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_global_bg, NA),
          Gcx_local_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_local_bg, NA)
   ) 
 
 
 
 
 protrusion_data_summarized_percluster <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, 
            Clustersize) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Signal_type=="foreground") 
 
 
 protrusion_data_summarized_perROI <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, Length_label) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 protrusion_data_summarized_percell <- protrusion_data %>% 
   group_by(Orientation, CellNR, Type, Fiber_type, Signal_type, Xthr_pos) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 subset_protrusion <- protrusion_data %>%  filter(is.na(Peak_position)) 
 
 summary_sds_protrusion <- subset_protrusion %>%  group_by(Orientation, number) %>%  summarize(B1int_sd= sd(B1int, na.rm=TRUE))
 
 
 
 ## copy to here for other Orientations
 
 ## retraction fibers - example plot to depict Gcx and B1 integrin intensity profiles------
 
 sub <- subset(protrusion_data, CellNR==19 & LocationNR==1 & RoiNR==4)
 
 #sub <- subset(protrusion_data, number==37)
 
 legendsize=0.3
 
 
 Zdiff = 0.25
 
 sub_summary <- sub %>%
   filter(Signal_type=="foreground") %>%
   summarize(mean_B1int = mean(B1int_smooth),
             sd_B1int = sd(B1int_smooth),
             iv = mean_B1int+sd_B1int,
             thr = mean_B1int + sd_B1int*Zdiff)
 
 
 sub$Signal_type <- as.factor(sub$Signal_type)
 sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
 
 sub <- sub %>%  filter(Signal_type=="foreground")
 
 lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int, linetype="Unmodified signal"))+
   geom_line( color="magenta", alpha=1)+
   theme_classic()
 
 limits_conversion <- max(sub$B1int)/ max(sub$Gcx)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_global_bg, linetype="Global background"), color="magenta", alpha=1)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_local_bg ,linetype="Local background"), color="magenta", alpha=1)+
   labs(linetype="")+
   theme(legend.position="top")+
   scale_linetype_manual(values = c("dotted", "dashed", "solid"))
 
 
 lineplot2 <- ggplot(data=sub, aes(x=X_micron, y=Gcx, linetype="Unmodified signal"))+
   geom_line( color="gold3", alpha=1)+
   theme_classic()
 
 lineplot2 <- lineplot2 + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_global_bg, linetype="Global background"),color="gold3", alpha=1)
 lineplot2 <- lineplot2 + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_local_bg, linetype="Local background"),color="gold3", alpha=1)+
   labs(linetype="")+
   theme(legend.position="none")+
   scale_linetype_manual(values = c("dotted", "dashed", "solid"))
 
 rect_data <- sub %>%
   filter(Signal_type=="foreground") %>%
   group_by(new_Peak_ID) %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
             xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize
   ) %>%  filter(!is.na(new_Peak_ID))
 
 clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("blue", 0.5))
 
 cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", fill="black", color=alpha("blue", 0.3),
                                 inherit.aes = FALSE, alpha = 0.1, size=0.5)
 
 example_subset <-  geom_point(data=subset(sub, !is.na(Peaklocation_X)), size=1.5, stroke=0)
 example_subset_globalbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X)), aes(y=B1int_cluster_global_bg), size=1.5, stroke=0)
 example_subset_localbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X)), aes(y=B1int_cluster_local_bg), size=1.5, stroke=0)
 
 
 lineplot <- lineplot +clustermids +example_subset+ example_subset_globalbg+ example_subset_localbg
 
 
 example_subset <-  geom_point(data=subset(sub, !is.na(Peaklocation_X)), size=1.5, stroke=0)
 example_subset_globalbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X)), aes(y=Gcx_cluster_global_bg), size=1.5, stroke=0)
 example_subset_localbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X)), aes(y=Gcx_cluster_local_bg), size=1.5, stroke=0)
 example_subset_gcx_adj <-  geom_point(data=subset(sub, Peak_position=="Edge"), aes(y=Gcx), size=1.5, stroke=0)
 
 
 lineplot2 <- lineplot2 + clustermids + example_subset+ example_subset_globalbg+ example_subset_localbg
 
 
 
 
 lineplot <- lineplot + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 lineplot2 <- lineplot2 + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 
 lineplot <- lineplot + 
   theme(legend.position="top")+
   theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
         legend.key.height = unit(legendsize, 'cm'), #change legend key height
         legend.key.width = unit(legendsize*1.5, 'cm'))
 
 
 lineplot <- plot_grid(lineplot, lineplot2, nrow=2, rel_heights = c(0.6,0.4))
 
 print(lineplot)
 
 lineplot_retraction <- lineplot
 
 
 ## make corrections for  protrusions------
 
 
 protrusion_data <- subset(B1int_GCX_ratio, Orientation =="Leading edge protrusion")
 
 
 
 
 if(any(protrusion_data$Orientation=="Cellbody")){
   cellbody_summary_percell <- protrusion_data  %>% filter(Orientation=="Cellbody" & Signal_type=="foreground") %>%  
     group_by(Orientation, CellNR, Signal_type) %>% 
     dplyr::summarize(B1int_cellbody_sd = sd(B1int, na.rm=TRUE),
                      across(where(is.numeric), mean, na.rm=TRUE)
     ) %>% ungroup() %>% 
     dplyr::select(c(CellNR, B1int, B1int_cellbody_sd, Gcx)) %>% 
     mutate(ratio_cellbody_average = B1int/ Gcx,
     ) %>% 
     rename(Gcx_cellbody_average = Gcx,
            B1int_cellbody_average = B1int) %>% 
     mutate(peak_b1int_threshold = B1int_cellbody_average + 4 * B1int_cellbody_sd)
   
 }
 protrusion_data <- protrusion_data %>%
   group_by(Signal_type, CellNR) %>% 
   left_join(cellbody_summary_percell)
 
 
 
 ## copy from here to other Orientations if needed
 protrusion_data <- protrusion_data %>%  
   mutate_at(vars(Fiber_type, Type, Orientation, Xthr_pos), factor)
 
 
 #protrusion_data <- calc_local_segregation_leadingedge(protrusion_data, 3)
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 
 protrusion_data <- protrusion_data %>%
   group_by(number, Position_relative_to_fiber, Signal_type) %>% 
   mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                          ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                          NA),
          ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                    ratio/mean(ratio, na.rm=TRUE),
                                                    NA)
   )
 
 
 
 protrusion_data$Type[protrusion_data$Type=="Focalized" & is.na(protrusion_data$PeakID_perprofile)] <- "Unspecified"
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   filter(Signal_type=="foreground")
 
 
 library(zoo)
 
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          WithinPeak = TRUE
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   arrange(X) %>% 
   mutate(Clustersize = max(cumsum(WithinPeak==TRUE)*pxsize),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Clustersize = ifelse(Clustersize < Clustersize_threshold, NA, Clustersize),
          new_Peak_ID = ifelse(is.na(Clustersize), NA, new_Peak_ID),
          new_Peak_ID = ifelse(max(B1int, na.rm=TRUE) < peak_b1int_threshold, NA, new_Peak_ID),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
   )
 
 
 
 # determine local B1 integrin and Gcx background for the integrin clusters
 px.clustersize.micron.local = 0.6
 px.clustersize.local = floor(px.clustersize.micron.local/pxsize)
 
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(min_B1int_local = ifelse(X==min(X) | X==max(X), B1int, NA),
          min_Gcx_local =   ifelse(X==min(X) | X==max(X), Gcx, NA),
          interval.local = 1
   )%>%
   ungroup() %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>%
   mutate(B1int_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID), 
                                          zoo::na.approx(min_B1int_local, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID),
                                        zoo::na.approx(min_Gcx_local, na.rm=FALSE, rule=2), NA)
   )%>% 
   dplyr::select(-c(min_B1int_local, min_Gcx_local, interval.local))
 
 
 
 
 
 ## determine global B1 integrin and Gcx background for the integrin clusters
 
 # WORK HERE WORK HERE WORK HERE
 px.clustersize.micron = 4
 px.clustersize = floor(px.clustersize.micron/pxsize)
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, Signal_type) %>%
   arrange(X) %>%
   mutate(
     interval.end = ifelse(!is.na(B1int_cluster_local_bg), cumsum(!is.na(B1int_cluster_local_bg) & lead(is.na(B1int_cluster_local_bg))), NA),
     interval = ifelse(!is.na(interval.end) & lead(is.na(interval.end)), NA, interval.end),
     interval = as.factor(interval)
   ) %>% 
   ungroup()
 
 
 # Define the custom function
 run_protrusion_code <- function(data, px_clustersize) {
   result <- data %>% 
     group_by(Orientation, number, Signal_type, interval)  %>% 
     mutate(
       composite.length = ifelse(!is.na(interval), max(row_number()), NA),
       dividecluster = ifelse(!is.na(interval) & composite.length > px_clustersize, TRUE, FALSE),
       clusterX = ifelse(!is.na(interval) & dividecluster, row_number(), NA),
       pastthrescoord = ifelse(!is.na(interval) & dividecluster & clusterX > px_clustersize, TRUE, FALSE),
       maxedge = ifelse(!is.na(interval) & pastthrescoord, max(Peak_position == "Edge"), NA),
       lastedge = ifelse(!is.na(interval) & cumsum(pastthrescoord & Peak_position == "Edge") == maxedge, TRUE, FALSE),
       lastedge = ifelse(lastedge==FALSE & lag(lastedge==TRUE), TRUE, lastedge),
       lastedge = ifelse(lastedge==TRUE & lead(lastedge==TRUE),FALSE, lastedge),
       dividecoord = ifelse(!is.na(interval) & length(unique(pastthrescoord)) > 1 & lastedge, TRUE, FALSE),
       dividecoord = ifelse(is.na(dividecoord), FALSE, dividecoord),
       subinterval = cumsum(dividecoord),
     ) %>% 
     ungroup() %>%
     mutate(interval = ifelse(!is.na(interval), group_indices(., interval, subinterval), NA))
   
   return(result)
 }
 
 repeated_samplingdistance = 10
 
 # Apply the custom function to your data twice
 for(p in 1:repeated_samplingdistance){
   print(p)
   protrusion_data <- run_protrusion_code(protrusion_data, px.clustersize)
 }
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, interval) %>% 
   arrange(X) %>% 
   mutate(
     min_B1int = ifelse(X==min(X, na.rm=TRUE), B1int, NA),
     min_B1int = ifelse(X==max(X, na.rm=TRUE), B1int, min_B1int),
     # min_B1int = ifelse(dividecoord==TRUE, B1int, min_B1int),
     min_Gcx = ifelse(X==min(X, na.rm=TRUE) | X==max(X, na.rm=TRUE), Gcx, NA),
   )
 
 protrusion_data <- protrusion_data %>% 
   ungroup() %>% 
   group_by(Orientation, number, Signal_type) %>%
   mutate(B1int_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_B1int, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_Gcx, na.rm=FALSE, rule=2), NA),
          #B1int_cluster_local_bg = ifelse(B1int_cluster_local_bg > B1int_cluster_global_bg, B1int_cluster_local_bg, B1int_cluster_global_bg)
   )
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge")
   )
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 protrusion_data <- protrusion_data %>%  
   group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>% 
   mutate(unique_Peak_ID = cur_group_id())
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, unique_Peak_ID) %>% 
   mutate(B1int_global_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_global_bg, NA),
          B1int_local_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_local_bg, NA),
          Gcx_global_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_global_bg, NA),
          Gcx_local_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_local_bg, NA)
   ) 
 
 
 
 
 protrusion_data_summarized_percluster <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, 
            Clustersize) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Signal_type=="foreground") 
 
 
 protrusion_data_summarized_perROI <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, Length_label) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 protrusion_data_summarized_percell <- protrusion_data %>% 
   group_by(Orientation, CellNR, Type, Fiber_type, Signal_type, Xthr_pos) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 subset_protrusion <- protrusion_data %>%  filter(is.na(Peak_position)) 
 
 summary_sds_protrusion <- subset_protrusion %>%  group_by(Orientation, number) %>%  summarize(B1int_sd= sd(B1int, na.rm=TRUE))
 
 
 
 
 ## copy to here for other Orientations
 
 ## protrusions - example plot to depict Gcx and B1 integrin intensity profiles------
 
 sub <- subset(protrusion_data, CellNR==7 & LocationNR==1 & RoiNR==1)
 
 #sub <- subset(protrusion_data, number==37)
 
 legendsize=0.3
 
 
 Zdiff = 0.25
 
 sub_summary <- sub %>%
   filter(Signal_type=="foreground") %>%
   summarize(mean_B1int = mean(B1int_smooth),
             sd_B1int = sd(B1int_smooth),
             iv = mean_B1int+sd_B1int,
             thr = mean_B1int + sd_B1int*Zdiff)
 
 
 sub$Signal_type <- as.factor(sub$Signal_type)
 sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
 
 sub <- sub %>%  filter(Signal_type=="foreground")
 
 lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int))+
   geom_line( color="magenta", alpha=1)+
   theme_classic()
 
 limits_conversion <- max(sub$B1int)/ max(sub$Gcx)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx*limits_conversion), color="gold3", alpha=1)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_global_bg), linetype="dotted", color="magenta", alpha=1)
 #lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_local_bg), linetype="dashed", color="magenta", alpha=1)
 
 #lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_global_bg*limits_conversion),linetype="dotted", color="gold3", alpha=1)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_local_bg*limits_conversion),linetype="dashed", color="gold3", alpha=1)
 
 
 
 
 print(lineplot)
 
 
 rect_data <- sub %>%
   filter(Signal_type=="foreground") %>%
   group_by(new_Peak_ID) %>%
   filter(!is.na(Clustersize)) %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
             xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize
   )
 
 lineplot <- lineplot + geom_point(data=subset(sub, !is.na(Peaklocation_X)), size=1.5, stroke=0)
 
 lineplot <- lineplot+
   geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", fill="black", color=alpha("blue", 0),
             inherit.aes = FALSE, alpha = 0, size=0.5)
 
 lineplot <- lineplot+
   # geom_vline(data=rect_data,  linetype="dashed", aes(xintercept=xintercept), alpha=0.5)+
   scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
   scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
   scale_y_continuous(
     # Features of the first axis
     name = "B1 integrin (x1000)",
     labels= scales::number_format(scale=1e-3),
     # Add a second axis and specify its features
     sec.axis = sec_axis(~./limits_conversion, name="Gcx (x1000)",
                         labels= scales::number_format(scale=1e-3))
   )
 
 rect_data2 <- sub %>% ungroup() %>% 
   filter(Signal_type=="foreground") %>%
   group_by(Length_label) %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf) %>%  filter(Length_label!="Intermediate zone")
 
 lineplot <- lineplot + 
   geom_rect(data=rect_data2, aes(xmin= xmin, xmax=xmax, ymin=ymin, ymax=ymax, linetype=Length_label), fill=NA, color="black",
             size=0.5, inherit.aes = FALSE)+
   scale_linetype_manual(values=c("solid", "dashed"))+
   labs(linetype="", x="Distance (\u03BCm)")+
   theme(legend.position="top")+
   theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
         legend.key.height = unit(legendsize, 'cm'), #change legend key height
         legend.key.width = unit(legendsize*1.5, 'cm'))
 
 
 lineplot <- lineplot
 print(lineplot)
 
 lineplot_example1 <- lineplot
 
 plot_data <- ggplot_build(lineplot)
 
 ylimits <- plot_data$layout$panel_scales_y[[1]]$range$range
 
 lineplot_cellbody <- lineplot_cellbody + scale_y_continuous(limits=ylimits, labels=scales::number_format(scale=1e-3, suffix=""), 
                                                             sec.axis = sec_axis(~./limits_conversion, name="Gcx (x1000)",
                                                                                 labels= scales::number_format(scale=1e-3))) + ylab("B1int (x1000)")
 
 lineplots_together <- plot_grid(lineplot_cellbody, lineplot, ncol=2, axis="tb", align="h", rel_widths = c(0.3,0.7))
 print(lineplots_together)
 
 ## protrusions - second example plot to depict Gcx and B1 integrin intensity profiles------
 
 sub <- subset(protrusion_data, Orientation=="Leading edge protrusion"& CellNR==7 & LocationNR==1 & RoiNR==1)
 
 
 legendsize=0.3
 
 
 Zdiff = 0.25
 
 sub_summary <- sub %>%
   summarize(mean_B1int = mean(B1int_smooth),
             sd_B1int = sd(B1int_smooth),
             iv = mean_B1int+sd_B1int,
             thr = mean_B1int + sd_B1int*Zdiff)
 
 
 sub$Signal_type <- as.factor(sub$Signal_type)
 sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
 
 lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int, linetype="Unmodified signal"))+
   geom_line( color="magenta", alpha=1)+
   theme_classic()
 
 limits_conversion <- max(sub$B1int)/ max(sub$Gcx)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_global_bg, linetype="Global background"), color="magenta", alpha=1)
 lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_local_bg ,linetype="Local background"), color="magenta", alpha=1)+
   labs(linetype="")+
   theme(legend.position="top")+
   scale_linetype_manual(values = c("dotted", "dashed", "solid"))
 
 
 lineplot2 <- ggplot(data=sub, aes(x=X_micron, y=Gcx, linetype="Unmodified signal"))+
   geom_line( color="gold3", alpha=1)+
   theme_classic()
 
 lineplot2 <- lineplot2 + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_global_bg, linetype="Global background"),color="gold3", alpha=1)
 lineplot2 <- lineplot2 + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_local_bg, linetype="Local background"),color="gold3", alpha=1)+
   labs(linetype="")+
   theme(legend.position="none")+
   scale_linetype_manual(values = c("dotted", "dashed", "solid"))
 
 rect_data <- sub %>%
   filter(Signal_type=="foreground") %>%
   group_by(new_Peak_ID) %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
             xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize
   ) %>%  filter(!is.na(new_Peak_ID))
 
 clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("blue", 0.5))
 
 cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", fill="black", color=alpha("blue", 0.3),
                                 inherit.aes = FALSE, alpha = 0.1, size=0.5)
 
 
 
 lineplot <- lineplot +clustermids
 
 
 
 lineplot2 <- lineplot2 + clustermids 
 
 
 
 lineplot <- lineplot + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 lineplot2 <- lineplot2 + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 
 lineplot <- lineplot + 
   theme(legend.position="top")+
   theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
         legend.key.height = unit(legendsize, 'cm'), #change legend key height
         legend.key.width = unit(legendsize*1.5, 'cm'))
 
 
 lineplot <- plot_grid(lineplot, lineplot2, nrow=2, rel_heights = c(0.6,0.4))
 
 print(lineplot)
 
 lineplot_alternative <- lineplot
 
 
 
 
 
 
 
 
 ##  NEW APPROACH - make corrections for retraction fibers------
 
 
 protrusion_data <- subset(B1int_GCX_ratio, Orientation =="Retraction fiber")
 
 if(any(protrusion_data$Orientation=="Cellbody")){
   cellbody_summary_percell <- protrusion_data  %>% filter(Orientation=="Cellbody" & Signal_type=="foreground") %>%  
     group_by(Orientation, CellNR, Signal_type) %>% 
     dplyr::summarize(B1int_cellbody_sd = sd(B1int, na.rm=TRUE),
                      across(where(is.numeric), mean, na.rm=TRUE)
     ) %>% ungroup() %>% 
     dplyr::select(c(CellNR, B1int, B1int_cellbody_sd, Gcx)) %>% 
     mutate(ratio_cellbody_average = B1int/ Gcx,
     ) %>% 
     rename(Gcx_cellbody_average = Gcx,
            B1int_cellbody_average = B1int) %>% 
     mutate(peak_b1int_threshold = B1int_cellbody_average + 4 * B1int_cellbody_sd)
   
 }
 protrusion_data <- protrusion_data %>%
   group_by(Signal_type, CellNR) %>% 
   left_join(cellbody_summary_percell)
 
 
 
 ## copy from here to other Orientations if needed
 protrusion_data <- protrusion_data %>%  
   mutate_at(vars(Fiber_type, Type, Orientation, Xthr_pos), factor)
 
 
 #protrusion_data <- calc_local_segregation_leadingedge(protrusion_data, 3)
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 
 protrusion_data <- protrusion_data %>%
   group_by(number, Position_relative_to_fiber, Signal_type) %>% 
   mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                          ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                          NA),
          ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                    ratio/mean(ratio, na.rm=TRUE),
                                                    NA)
   )
 
 
 
 protrusion_data$Type[protrusion_data$Type=="Focalized" & is.na(protrusion_data$PeakID_perprofile)] <- "Unspecified"
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   filter(Signal_type=="foreground")
 
 
 library(zoo)
 
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          WithinPeak = TRUE
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   arrange(X) %>% 
   mutate(Clustersize = max(cumsum(WithinPeak==TRUE)*pxsize),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Clustersize = ifelse(Clustersize < Clustersize_threshold, NA, Clustersize),
          new_Peak_ID = ifelse(is.na(Clustersize), NA, new_Peak_ID),
          new_Peak_ID = ifelse(max(B1int, na.rm=TRUE) < peak_b1int_threshold, NA, new_Peak_ID),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
   )
 
 
 
 # determine local B1 integrin and Gcx background for the integrin clusters
 px.clustersize.micron.local = 0.6
 px.clustersize.local = floor(px.clustersize.micron.local/pxsize)
 
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(min_B1int_local = ifelse(X==min(X) | X==max(X), B1int, NA),
          min_Gcx_local =   ifelse(X==min(X) | X==max(X), Gcx, NA),
          interval.local = 1
   )%>%
   ungroup() %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>%
   mutate(B1int_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID), 
                                          zoo::na.approx(min_B1int_local, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID),
                                        zoo::na.approx(min_Gcx_local, na.rm=FALSE, rule=2), NA)
   )%>% 
   dplyr::select(-c(min_B1int_local, min_Gcx_local, interval.local))
 
 
 
 
 
 ## determine global B1 integrin and Gcx background for the integrin clusters
 
 # WORK HERE WORK HERE WORK HERE
 px.clustersize.micron = 4
 px.clustersize = floor(px.clustersize.micron/pxsize)
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, Signal_type) %>%
   arrange(X) %>%
   mutate(
     interval.end = ifelse(!is.na(B1int_cluster_local_bg), cumsum(!is.na(B1int_cluster_local_bg) & lead(is.na(B1int_cluster_local_bg))), NA),
     interval = ifelse(!is.na(interval.end) & lead(is.na(interval.end)), NA, interval.end),
     interval = as.factor(interval)
   ) %>% 
   ungroup()
 
 
 # Define the custom function
 run_protrusion_code <- function(data, px_clustersize) {
   result <- data %>% 
     group_by(Orientation, number, Signal_type, interval)  %>% 
     mutate(
       composite.length = ifelse(!is.na(interval), max(row_number()), NA),
       dividecluster = ifelse(!is.na(interval) & composite.length > px_clustersize, TRUE, FALSE),
       clusterX = ifelse(!is.na(interval) & dividecluster, row_number(), NA),
       pastthrescoord = ifelse(!is.na(interval) & dividecluster & clusterX > px_clustersize, TRUE, FALSE),
       maxedge = ifelse(!is.na(interval) & pastthrescoord, max(Peak_position == "Edge"), NA),
       lastedge = ifelse(!is.na(interval) & cumsum(pastthrescoord & Peak_position == "Edge") == maxedge, TRUE, FALSE),
       dividecoord = ifelse(!is.na(interval) & length(unique(pastthrescoord)) > 1 & lastedge, TRUE, FALSE),
       dividecoord = ifelse(is.na(dividecoord), FALSE, dividecoord),
       subinterval = cumsum(dividecoord)
     ) %>% 
     ungroup() %>%
     mutate(interval = ifelse(!is.na(interval), group_indices(., interval, subinterval), NA))
   
   return(result)
 }
 
 repeated_samplingdistance = 10
 
 # Apply the custom function to your data twice
 for(p in 1:repeated_samplingdistance){
   print(p)
   protrusion_data <- run_protrusion_code(protrusion_data, px.clustersize)
 }
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, interval) %>% 
   arrange(X) %>% 
   mutate(
     min_B1int = ifelse(X==min(X, na.rm=TRUE), B1int, NA),
     min_B1int = ifelse(X==max(X, na.rm=TRUE), B1int, min_B1int),
     # min_B1int = ifelse(dividecoord==TRUE, B1int, min_B1int),
     min_Gcx = ifelse(X==min(X, na.rm=TRUE) | X==max(X, na.rm=TRUE), Gcx, NA),
   )
 
 protrusion_data <- protrusion_data %>% 
   ungroup() %>% 
   group_by(Orientation, number, Signal_type) %>%
   mutate(B1int_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_B1int, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_Gcx, na.rm=FALSE, rule=2), NA),
          #B1int_cluster_local_bg = ifelse(B1int_cluster_local_bg > B1int_cluster_global_bg, B1int_cluster_local_bg, B1int_cluster_global_bg)
   )
 
 
 
 measuring_bandwidth = 7
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          Peak_position  = ifelse(is.na(new_Peak_ID), NA, Peak_position),
          X.cluster =  ifelse(!is.na(new_Peak_ID), X - min(X), NA),
          Peak_position = ifelse(X.cluster <=(measuring_bandwidth-1) | X.cluster >= (max(X.cluster)-(measuring_bandwidth-1)), "Edge", "Transition"),
          Peak_position = ifelse(!is.na(Peaklocation_X), "Peak", Peak_position)
   )
 
 for(g in 1:(measuring_bandwidth-1)){
   print(g)
   protrusion_data <- protrusion_data %>% 
     group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
     arrange(X) %>% 
     mutate(Peak_position = ifelse(lag(Peak_position=="Peak"), "Peak", Peak_position),
            Peak_position = ifelse(lead(Peak_position=="Peak"), "Peak", Peak_position)
     )
 }
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(Peak_position = ifelse(is.na(Peak_position), "Edge", Peak_position),
          Peak_position_number = ifelse(Peak_position=="Edge" & X.cluster < mean(X.cluster[!is.na(Peaklocation_X)]), 1, 2)
   )
 
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 protrusion_data <- protrusion_data %>%  
   group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>% 
   mutate(unique_Peak_ID = cur_group_id())
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, unique_Peak_ID) %>% 
   mutate(B1int_global_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_global_bg, NA),
          B1int_local_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_local_bg, NA),
          Gcx_global_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_global_bg, NA),
          Gcx_local_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_local_bg, NA)
   ) 
 
 
 
 
 protrusion_data_summarized_percluster <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, 
            Clustersize) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Signal_type=="foreground") 
 
 
 protrusion_data_summarized_perROI <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, Length_label) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 protrusion_data_summarized_percell <- protrusion_data %>% 
   group_by(Orientation, CellNR, Type, Fiber_type, Signal_type, Xthr_pos) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 subset_protrusion <- protrusion_data %>%  filter(is.na(Peak_position)) 
 
 summary_sds_protrusion <- subset_protrusion %>%  group_by(Orientation, number) %>%  summarize(B1int_sd= sd(B1int, na.rm=TRUE))
 
 
 
 ## copy to here for other Orientations
 
 ## retraction fibers - example plot to depict Gcx and B1 integrin intensity profiles------
 
 sub <- subset(protrusion_data, CellNR==19 & LocationNR==1 & RoiNR==4)
 
 #sub <- subset(protrusion_data, number==37)
 
 legendsize=0.3
 
 
 Zdiff = 0.25
 
 sub_summary <- sub %>%
   filter(Signal_type=="foreground") %>%
   summarize(mean_B1int = mean(B1int_smooth),
             sd_B1int = sd(B1int_smooth),
             iv = mean_B1int+sd_B1int,
             thr = mean_B1int + sd_B1int*Zdiff)
 
 
 sub$Signal_type <- as.factor(sub$Signal_type)
 sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
 
 sub <- sub %>%  filter(Signal_type=="foreground")
 
 lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int))+
   geom_line( color="magenta", alpha=1)+
   geom_line(data=sub, aes(y=B1int_smooth), linetype="dotted", color="magenta")+
   theme_classic()
 
 
 
 lineplot2 <- ggplot(data=sub, aes(x=X_micron, y=Gcx))+
   geom_line( color="gold3", alpha=1)+
   theme_classic()
 
 rect_data <- sub %>%
   filter(Signal_type=="foreground") %>%
   group_by(new_Peak_ID) %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
             xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize-pxsize
   ) %>%  filter(!is.na(new_Peak_ID))
 
 clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("black", 1))
 
 cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", fill=alpha("black",0), color=alpha("black", 1),
                                 inherit.aes = FALSE, alpha = 0, size=0.5)
 
 lineplot <- lineplot +clustermids 
 
 
 lineplot2 <- lineplot2 + clustermids
 
 
 
 lineplot <- lineplot + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 lineplot2 <- lineplot2 + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 
 
 
 
 rect_data <- sub %>%
   filter(Signal_type=="foreground") %>%
   group_by(new_Peak_ID, Peak_position, Peak_position_number) %>% filter(Peak_position!="Transition") %>% 
   summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
             xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize
   ) %>%  filter(!is.na(new_Peak_ID))
 
 
 
 cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Peak_position), 
                                 linetype="solid", color=alpha("blue", 0),
                                 inherit.aes=FALSE, alpha = 0.5, size=0.5)
 
 
 lineplot <- lineplot + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 lineplot2 <- lineplot2 + cluster.rectangles+
   labs(linetype="", x="Distance (\u03BCm)")
 
 
 
 lineplot <- lineplot + 
   theme(legend.position="top")+
   theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
         legend.key.height = unit(legendsize, 'cm'), #change legend key height
         legend.key.width = unit(legendsize*1.5, 'cm'))
 
 
 lineplot <- plot_grid(lineplot, lineplot2, nrow=2, rel_heights = c(0.6,0.4))
 
 print(lineplot)
 
 lineplot_retraction <- lineplot
 
 
 
 
 
 
 
 ## make corrections for  protrusions------
 
 
 protrusion_data <- B1int_GCX_ratio
 
 
 
 if(any(protrusion_data$Orientation=="Cellbody")){
   cellbody_summary_percell <- protrusion_data  %>% filter(Orientation=="Cellbody" & Signal_type=="foreground") %>%  
     group_by(Orientation, CellNR, Signal_type) %>% 
     dplyr::summarize(B1int_cellbody_sd = sd(B1int, na.rm=TRUE),
                      across(where(is.numeric), mean, na.rm=TRUE)
     ) %>% ungroup() %>% 
     dplyr::select(c(CellNR, B1int, B1int_cellbody_sd, Gcx)) %>% 
     mutate(ratio_cellbody_average = B1int/ Gcx,
     ) %>% 
     rename(Gcx_cellbody_average = Gcx,
            B1int_cellbody_average = B1int) %>% 
     mutate(peak_b1int_threshold = B1int_cellbody_average + 4 * B1int_cellbody_sd)
   
 }
 protrusion_data <- protrusion_data %>%
   group_by(Signal_type, CellNR) %>% 
   left_join(cellbody_summary_percell)
 
 
 
 ## copy from here to other Orientations if needed
 protrusion_data <- protrusion_data %>%  
   mutate_at(vars(Fiber_type, Type, Orientation, Xthr_pos), factor)
 
 
 #protrusion_data <- calc_local_segregation_leadingedge(protrusion_data, 3)
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, Position_relative_to_fiber, Signal_type) %>% 
   mutate(ratio_normalized_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                          ratio/mean(ratio[Length_label=="Contact-free"], na.rm=TRUE),
                                                          NA),
          ratio_mean_leadingedgeprotrusion = ifelse(Orientation=="Leading edge protrusion",
                                                    ratio/mean(ratio, na.rm=TRUE),
                                                    NA)
   )
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Position_relative_to_fiber, Signal_type) %>% 
   mutate(Length_label = ifelse(Orientation=="Retraction fiber" & X_relative<retraction_fiber_length_segregation_threshold, "Retraction fiber", Length_label),
          Length_label = ifelse(Orientation=="Retraction fiber" & X_relative>=retraction_fiber_length_segregation_threshold, "Retraction fiber base", Length_label),
          X = ifelse(Orientation=="Retraction fiber", -X + min(X), X),
          X_micron = ifelse(Orientation=="Retraction fiber", -X_micron + min(X_micron), X_micron),
          X_relative = ifelse(Orientation == "Retraction fiber", -X_relative + min(X_relative), X_relative)
   )
 
 
 protrusion_data$Type[protrusion_data$Type=="Focalized" & is.na(protrusion_data$PeakID_perprofile)] <- "Unspecified"
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   filter(Signal_type=="foreground")
 
 
 library(zoo)
 
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          WithinPeak = TRUE
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>% 
   arrange(X) %>% 
   mutate(Clustersize = max(cumsum(WithinPeak==TRUE)*pxsize),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Clustersize = ifelse(Clustersize < Clustersize_threshold, NA, Clustersize),
          new_Peak_ID = ifelse(is.na(Clustersize), NA, new_Peak_ID),
          new_Peak_ID = ifelse(max(B1int, na.rm=TRUE) < peak_b1int_threshold, NA, new_Peak_ID),
          Clustersize = ifelse(is.na(new_Peak_ID), NA, Clustersize),
          Peaklocation_X = ifelse(B1int == max(B1int,na.rm=TRUE) & !is.na(new_Peak_ID), X, NA)
   )
 
 
 
 # determine local B1 integrin and Gcx background for the integrin clusters
 px.clustersize.micron.local = 0.6
 px.clustersize.local = floor(px.clustersize.micron.local/pxsize)
 
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(min_B1int_local = ifelse(X==min(X) | X==max(X), B1int, NA),
          min_Gcx_local =   ifelse(X==min(X) | X==max(X), Gcx, NA),
          interval.local = 1
   )%>%
   ungroup() %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID) %>%
   mutate(B1int_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID), 
                                          zoo::na.approx(min_B1int_local, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_local_bg = ifelse(!is.na(interval.local) & !is.na(new_Peak_ID),
                                        zoo::na.approx(min_Gcx_local, na.rm=FALSE, rule=2), NA)
   )%>% 
   dplyr::select(-c(min_B1int_local, min_Gcx_local, interval.local))
 
 
 
 
 
 ## determine global B1 integrin and Gcx background for the integrin clusters
 
 # WORK HERE WORK HERE WORK HERE
 px.clustersize.micron = 4
 px.clustersize = floor(px.clustersize.micron/pxsize)
 
 protrusion_data <- protrusion_data %>%
   group_by(Orientation, number, Signal_type) %>%
   arrange(X) %>%
   mutate(
     interval.end = ifelse(!is.na(B1int_cluster_local_bg), cumsum(!is.na(B1int_cluster_local_bg) & lead(is.na(B1int_cluster_local_bg))), NA),
     interval = ifelse(!is.na(interval.end) & lead(is.na(interval.end)), NA, interval.end),
     interval = as.factor(interval)
   ) %>% 
   ungroup()
 
 
 # Define the custom function
 run_protrusion_code <- function(data, px_clustersize) {
   result <- data %>% 
     group_by(Orientation, number, Signal_type, interval)  %>% 
     mutate(
       composite.length = ifelse(!is.na(interval), max(row_number()), NA),
       dividecluster = ifelse(!is.na(interval) & composite.length > px_clustersize, TRUE, FALSE),
       clusterX = ifelse(!is.na(interval) & dividecluster, row_number(), NA),
       pastthrescoord = ifelse(!is.na(interval) & dividecluster & clusterX > px_clustersize, TRUE, FALSE),
       maxedge = ifelse(!is.na(interval) & pastthrescoord, max(Peak_position == "Edge"), NA),
       lastedge = ifelse(!is.na(interval) & cumsum(pastthrescoord & Peak_position == "Edge") == maxedge, TRUE, FALSE),
       lastedge = ifelse(lastedge==FALSE & lag(lastedge==TRUE), TRUE, lastedge),
       lastedge = ifelse(lastedge==TRUE & lead(lastedge==TRUE),FALSE, lastedge),
       dividecoord = ifelse(!is.na(interval) & length(unique(pastthrescoord)) > 1 & lastedge, TRUE, FALSE),
       dividecoord = ifelse(is.na(dividecoord), FALSE, dividecoord),
       subinterval = cumsum(dividecoord)
     ) %>% 
     ungroup() %>%
     mutate(interval = ifelse(!is.na(interval), group_indices(., interval, subinterval), NA))
   
   return(result)
 }
 
 repeated_samplingdistance = 10
 
 # Apply the custom function to your data twice
 for(p in 1:repeated_samplingdistance){
   print(p)
   protrusion_data <- run_protrusion_code(protrusion_data, px.clustersize)
 }
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, interval) %>% 
   arrange(X) %>% 
   mutate(
     min_B1int = ifelse(X==min(X, na.rm=TRUE), B1int, NA),
     min_B1int = ifelse(X==max(X, na.rm=TRUE), B1int, min_B1int),
     # min_B1int = ifelse(dividecoord==TRUE, B1int, min_B1int),
     min_Gcx = ifelse(X==min(X, na.rm=TRUE) | X==max(X, na.rm=TRUE), Gcx, NA),
   )
 
 protrusion_data <- protrusion_data %>% 
   ungroup() %>% 
   group_by(Orientation, number, Signal_type) %>%
   mutate(B1int_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_B1int, na.rm=FALSE, rule=2), NA),
          Gcx_cluster_global_bg = ifelse(!is.na(interval), zoo::na.approx(min_Gcx, na.rm=FALSE, rule=2), NA),
          B1int_cluster_global_bg = ifelse(B1int_cluster_local_bg < B1int_cluster_global_bg, B1int_cluster_local_bg, B1int_cluster_global_bg)
   )
 
 
 
 measuring_bandwidth_edge = 4
 measuring_bandwidth_peak = 2
 
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(WithinPeak = ifelse(X==min(X), FALSE, TRUE ),
          WithinPeak = ifelse(X==max(X), FALSE, WithinPeak),
          Peak_position = ifelse(WithinPeak==TRUE, "Inside", "Edge"),
          Peak_position  = ifelse(is.na(new_Peak_ID), NA, Peak_position),
          X.cluster =  ifelse(!is.na(new_Peak_ID), X - min(X), NA),
          Peak_position = ifelse(X.cluster <=(measuring_bandwidth_edge-1) | X.cluster >= (max(X.cluster)-(measuring_bandwidth_edge-1)), "Edge", "Transition"),
          Peak_position = ifelse(!is.na(Peaklocation_X), "Peak", Peak_position)
   )
 
 for(g in 1:(measuring_bandwidth_peak-1)){
   print(g)
   protrusion_data <- protrusion_data %>% 
     group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
     arrange(X) %>% 
     mutate(Peak_position = ifelse(lag(Peak_position=="Peak"), "Peak", Peak_position),
            Peak_position = ifelse(lead(Peak_position=="Peak"), "Peak", Peak_position)
     )
 }
 
 protrusion_data <- protrusion_data %>% 
   group_by(number, Signal_type, LocationNR, new_Peak_ID) %>%
   arrange(X) %>% 
   mutate(Peak_position = ifelse(is.na(Peak_position), "Edge", Peak_position),
          Peak_position_number = ifelse(Peak_position=="Edge" & X.cluster < mean(X.cluster[!is.na(Peaklocation_X)]), 1, 2)
   )
 
 
 protrusion_data <- protrusion_data %>% 
   ungroup() %>% 
   group_by(number, Signal_type, LocationNR, interval) %>% 
   arrange(X) %>% 
   mutate(X.composite = ifelse(!is.na(new_Peak_ID), X - min(X), NA),
          Edgeclusternr = ifelse(!is.na(new_Peak_ID), cumsum(Peak_position=="Edge" & lead(Peak_position!="Edge")),NA),
          Edgeclusternr = ifelse(!is.na(Edgeclusternr) & lag(Edgeclusternr)!=Edgeclusternr, lag(Edgeclusternr), Edgeclusternr),
          Peak_position_composite = ifelse(!is.na(new_Peak_ID) & Peak_position=="Edge" & Edgeclusternr==0,"Edge", "Transition"),
          Peak_position_composite = ifelse(!is.na(new_Peak_ID) & Peak_position=="Edge" & Edgeclusternr==max(Edgeclusternr, na.rm=TRUE),"Edge", Peak_position_composite),
          Peak_position_composite = ifelse(Peak_position=="Peak", "Peak", Peak_position_composite),
          Peak_position_composite = ifelse(is.na(new_Peak_ID), NA, Peak_position_composite),
          Peak_position_composite_number = ifelse(Peak_position_composite=="Edge" & X.composite < (px.clustersize/2), 1, 2)
   )                                                 
 
 
 protrusion_data$new_Peak_ID <- as.character(protrusion_data$new_Peak_ID)
 
 
 protrusion_data <- protrusion_data %>%  
   group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>% 
   mutate(unique_Peak_ID = cur_group_id())
 
 protrusion_data <- protrusion_data %>% ungroup() %>% 
   group_by(Orientation, number, Signal_type, interval) %>% 
   mutate(  Gcx_adjacent_composite = mean(Gcx[Peak_position_composite=="Edge"], na.rm=TRUE),
            B1int_adjacent_composite = mean(B1int[Peak_position_composite=="Edge"], na.rm=TRUE)
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, unique_Peak_ID) %>% 
   mutate(B1int_global_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_global_bg, NA),
          B1int_local_enrichment = ifelse(X==Peaklocation_X, B1int / B1int_cluster_local_bg, NA),
          Gcx_global_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_global_bg, NA),
          Gcx_local_enrichment = ifelse(X==Peaklocation_X, Gcx / Gcx_cluster_local_bg, NA),
          Gcx_enrichment = ifelse(X==Peaklocation_X, mean(Gcx[Peak_position=="Peak"], na.rm=TRUE) / mean(Gcx[Peak_position=="Edge"], na.rm=TRUE), NA),
          B1int_enrichment = ifelse(X==Peaklocation_X, mean(B1int[Peak_position=="Peak"], na.rm=TRUE) / mean(B1int[Peak_position=="Edge"], na.rm=TRUE), NA),
          B1int_enrichment_composite = ifelse(X==Peaklocation_X, mean(B1int[Peak_position=="Peak"], na.rm=TRUE) / mean(B1int_adjacent_composite, na.rm=TRUE)),
          Gcx_enrichment_composite = ifelse(X==Peaklocation_X, mean(Gcx[Peak_position=="Peak"], na.rm=TRUE) / mean(Gcx_adjacent_composite, na.rm=TRUE))
   )
 
 protrusion_data <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type) %>% 
   mutate(X_relative = X/ max(row_number(), na.rm=TRUE),
          Length_label= ifelse(Orientation=="Bleb" & X_relative>= (1/12*max(X_relative)) & X_relative<=(12/12*(max(X_relative))), "Bleb Base", Length_label),
          Length_label= ifelse(Orientation=="Bleb" & X_relative>= (3/12*max(X_relative)) & X_relative<=(9/12*(max(X_relative))), "Bleb Neck", Length_label),
          Length_label= ifelse(Orientation=="Bleb" & X_relative>= (5/12*max(X_relative)) & X_relative<=(7/12*(max(X_relative))), "Bleb", Length_label),
          Gcx_bleb_base_normalized = Gcx/ mean(Gcx[Length_label=="Base"]),
          B1int_bleb_base_normalized = B1int / mean(B1int[Length_label=="Base"]),
          Col_bleb_base_normalized = Col/ mean(Col[Length_label=="Base"]),
          Bleb_length = length(Gcx[Length_label=="Bleb"])*pxsize
   ) %>% arrange(X) %>% 
   mutate(Length_classification_nr = cumsum(X>1 &Length_label != lag(Length_label)))
 
 
 
 
 protrusion_data_summarized_percluster <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, new_Peak_ID, 
            Clustersize) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) %>%  filter(Signal_type=="foreground") 
 
 
 protrusion_data_summarized_perROI <- protrusion_data %>% 
   group_by(Orientation, number, Signal_type, Length_label) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 protrusion_data_summarized_percell <- protrusion_data %>% 
   group_by(Orientation, CellNR, Type, Fiber_type, Signal_type, Xthr_pos) %>% 
   dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
 
 subset_protrusion <- protrusion_data %>%  filter(is.na(Peak_position)) 
 
 summary_sds_protrusion <- subset_protrusion %>%  group_by(Orientation, number) %>%  summarize(B1int_sd= sd(B1int, na.rm=TRUE))
 
 
 
 
 ## copy to here for other Orientations
 
 ## protrusions - example plot to depict Gcx and B1 integrin intensity profiles------
 
 if(any(protrusion_data$Orientation=="Leading edge protrusion")){
   
   sub <- subset(protrusion_data, CellNR==7 & LocationNR==1 & RoiNR==1)
   
   #sub <- subset(protrusion_data, number==37)
   
   legendsize=0.3
   
   
   Zdiff = 0.25
   
   sub_summary <- sub %>%
     filter(Signal_type=="foreground") %>%
     summarize(mean_B1int = mean(B1int_smooth),
               sd_B1int = sd(B1int_smooth),
               iv = mean_B1int+sd_B1int,
               thr = mean_B1int + sd_B1int*Zdiff)
   
   
   sub$Signal_type <- as.factor(sub$Signal_type)
   sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
   
   sub <- sub %>%  filter(Signal_type=="foreground")
   
   lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int))+
     geom_line( color="magenta", alpha=1)+
     theme_classic()
   
   limits_conversion <- max(sub$B1int)/ max(sub$Gcx)
   lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx), color="gold3", alpha=1)
   lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_global_bg), linetype="dotted", color="magenta", alpha=1)
   #lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_local_bg), linetype="dashed", color="magenta", alpha=1)
   
   #lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_global_bg*limits_conversion),linetype="dotted", color="gold3", alpha=1)
   lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_local_bg),linetype="dashed", color="gold3", alpha=1)
   
   
   
   
   print(lineplot)
   
   
   rect_data <- sub %>%
     filter(Signal_type=="foreground") %>%
     group_by(new_Peak_ID) %>%
     filter(!is.na(Clustersize)) %>% 
     summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
               xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize
     )
   
   lineplot <- lineplot + geom_point(data=subset(sub, !is.na(Peaklocation_X)), size=1.5, stroke=0)
   
   lineplot <- lineplot+
     geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", fill="black", color=alpha("blue", 0),
               inherit.aes = FALSE, alpha = 0, size=0.5)
   
   lineplot <- lineplot+
     # geom_vline(data=rect_data,  linetype="dashed", aes(xintercept=xintercept), alpha=0.5)+
     scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
     scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
     scale_y_continuous(
       # Features of the first axis
       name = "B1 integrin (x1000)",
       labels= scales::number_format(scale=1e-3),
       # Add a second axis and specify its features
       sec.axis = sec_axis(~./limits_conversion, name="Gcx (x1000)",
                           labels= scales::number_format(scale=1e-3))
     )
   
   rect_data2 <- sub %>% ungroup() %>% 
     filter(Signal_type=="foreground") %>%
     group_by(Length_label) %>% 
     summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf) %>%  filter(Length_label!="Intermediate zone")
   
   lineplot <- lineplot + 
     geom_rect(data=rect_data2, aes(xmin= xmin, xmax=xmax, ymin=ymin, ymax=ymax, linetype=Length_label), fill=NA, color="black",
               size=0.5, inherit.aes = FALSE)+
     scale_linetype_manual(values=c("solid", "dashed"))+
     labs(linetype="", x="Distance (\u03BCm)")+
     theme(legend.position="top")+
     theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
           legend.key.height = unit(legendsize, 'cm'), #change legend key height
           legend.key.width = unit(legendsize*1.5, 'cm'))
   
   
   lineplot <- lineplot
   print(lineplot)
   
   lineplot_example1 <- lineplot
   
   plot_data <- ggplot_build(lineplot)
   
   ylimits <- plot_data$layout$panel_scales_y[[1]]$range$range
   
   lineplot_cellbody <- lineplot_cellbody + scale_y_continuous(limits=ylimits, labels=scales::number_format(scale=1e-3, suffix=""), 
                                                               sec.axis = sec_axis(~./limits_conversion, name="Gcx (x1000)",
                                                                                   labels= scales::number_format(scale=1e-3))) + ylab("B1int (x1000)")
   
   lineplots_together <- plot_grid(lineplot_cellbody, lineplot, ncol=2, axis="tb", align="h", rel_widths = c(0.3,0.7))
   print(lineplots_together)
 }
 
 ## protrusions - second example plot to depict Gcx and B1 integrin intensity profiles------
 
 if(any(protrusion_data$Orientation=="Leading edge protrusion")){
   sub <- subset(protrusion_data, Orientation=="Leading edge protrusion" & CellNR==7 & LocationNR==1 & RoiNR==1)
   
   #sub <- subset(protrusion_data, number==37)
   
   legendsize=0.3
   
   
   Zdiff = 0.25
   
   sub_summary <- sub %>%
     filter(Signal_type=="foreground") %>%
     summarize(mean_B1int = mean(B1int_smooth),
               sd_B1int = sd(B1int_smooth),
               iv = mean_B1int+sd_B1int,
               thr = mean_B1int + sd_B1int*Zdiff)
   
   
   sub$Signal_type <- as.factor(sub$Signal_type)
   sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
   
   sub <- sub %>%  filter(Signal_type=="foreground")
   
   lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int, linetype="Unmodified signal"))+
     geom_line( color="magenta", alpha=1)+
     theme_classic()
   
   limits_conversion <- max(sub$B1int)/ max(sub$Gcx)
   lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_global_bg, linetype="Global background"), color="magenta", alpha=1)
   lineplot <- lineplot + geom_line(data=sub, aes(x=X_micron, y=B1int_cluster_local_bg ,linetype="Local background"), color="magenta", alpha=1)+
     labs(linetype="")+
     theme(legend.position="top")+
     scale_linetype_manual(values = c("dotted", "dashed", "solid"))
   
   lineplot <- lineplot + geom_point(data=subset(sub, !is.na(min_B1int)))
   
   lineplot2 <- ggplot(data=sub, aes(x=X_micron, y=Gcx, linetype="Unmodified signal"))+
     geom_line( color="gold3", alpha=1)+
     theme_classic()
   
   lineplot2 <- lineplot2 + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_global_bg, linetype="Global background"),color="gold3", alpha=1)
   lineplot2 <- lineplot2 + geom_line(data=sub, aes(x=X_micron, y=Gcx_cluster_local_bg, linetype="Local background"),color="gold3", alpha=1)+
     labs(linetype="")+
     theme(legend.position="none")+
     scale_linetype_manual(values = c("dotted", "dashed", "solid"))
   
   rect_data <- sub %>%
     filter(Signal_type=="foreground") %>%
     group_by(new_Peak_ID) %>% 
     summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf,
               xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize
     ) %>%  filter(!is.na(new_Peak_ID))
   
   clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("blue", 0.5))
   
   cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), linetype="solid", fill="black", color=alpha("blue", 0.3),
                                   inherit.aes = FALSE, alpha = 0.1, size=0.5)
   
   example_subset <-  geom_point(data=subset(sub, !is.na(Peaklocation_X) & new_Peak_ID == 7), size=1.5, stroke=0)
   example_subset_globalbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X) & new_Peak_ID == 7), aes(y=B1int_cluster_global_bg), size=1.5, stroke=0)
   example_subset_localbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X) & new_Peak_ID == 7), aes(y=B1int_cluster_local_bg), size=1.5, stroke=0)
   
   
   lineplot <- lineplot +clustermids +example_subset+ example_subset_globalbg+ example_subset_localbg
   
   
   example_subset <-  geom_point(data=subset(sub, !is.na(Peaklocation_X) & new_Peak_ID == 7), size=1.5, stroke=0)
   example_subset_globalbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X) & new_Peak_ID == 7), aes(y=Gcx_cluster_global_bg), size=1.5, stroke=0)
   example_subset_localbg <-  geom_point(data=subset(sub, !is.na(Peaklocation_X) & new_Peak_ID == 7), aes(y=Gcx_cluster_local_bg), size=1.5, stroke=0)
   example_subset_gcx_adj <-  geom_point(data=subset(sub, Peak_position=="Edge" & new_Peak_ID == 7), aes(y=Gcx), size=1.5, stroke=0)
   
   
   lineplot2 <- lineplot2 + clustermids + example_subset+ example_subset_globalbg+ example_subset_localbg
   
   
   
   #lineplot <- lineplot + xlim(0,11)
   #lineplot2 <- lineplot2 + xlim(0,11)
   
   
   
   lineplot <- lineplot + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")
   lineplot2 <- lineplot2 + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")
   
   lineplot <- lineplot + 
     theme(legend.position="top")+
     theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
           legend.key.height = unit(legendsize, 'cm'), #change legend key height
           legend.key.width = unit(legendsize*1.5, 'cm'))
   
   
   lineplot <- plot_grid(lineplot, lineplot2, nrow=2, rel_heights = c(0.6,0.4))
   
   print(lineplot)
   
   lineplot_alternative <- lineplot
   
 }
 
 
 ## protrusions -  fourth example plot to depict Gcx and B1 integrin intensity profiles------
 if(any(protrusion_data$Orientation=="Leading edge protrusion")){
   
   sub <- subset(protrusion_data, CellNR==7 & LocationNR==1 & RoiNR==1)
   
   #sub <- subset(protrusion_data, number==37)
   
   legendsize=0.3
   
   
   Zdiff = 0.25
   
   sub_summary <- sub %>%
     filter(Signal_type=="foreground") %>%
     summarize(mean_B1int = mean(B1int_smooth),
               sd_B1int = sd(B1int_smooth),
               iv = mean_B1int+sd_B1int,
               thr = mean_B1int + sd_B1int*Zdiff)
   
   
   sub$Signal_type <- as.factor(sub$Signal_type)
   sub$Signal_type <- factor(sub$Signal_type, levels=c("foreground", "background"))
   
   sub <- sub %>%  filter(Signal_type=="foreground")
   
   lineplot <- ggplot(data=sub, aes(x=X_micron, y=B1int))+
     geom_line( color="magenta", alpha=1)+
     geom_line(data=sub, aes(y=B1int_smooth), linetype="dotted", color="magenta")+
     theme_classic()
   
   
   
   lineplot2 <- ggplot(data=sub, aes(x=X_micron, y=Gcx))+
     geom_line( color="gold3", alpha=1)+
     theme_classic()
   
   rect_data <- sub %>%
     filter(Signal_type=="foreground") %>%
     group_by(new_Peak_ID, Peak_position) %>% 
     summarize(xmin= min(X_micron), xmax=max(X_micron), ymin=-Inf, ymax=Inf,
               xintercept= mean(Peaklocation_X,na.rm=TRUE)*pxsize-pxsize
     ) %>%  filter(!is.na(new_Peak_ID) & Peak_position=="Peak")
   
   clustermids <-  geom_vline(data=rect_data, aes(xintercept=xintercept), linetype="longdash", color=alpha("black", 1))
   
   cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Peak_position), 
                                   linetype="solid",color=alpha("black", 0),
                                   inherit.aes = FALSE, alpha = 0.5, size=0.5)
   
   
   lineplot <- lineplot + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")+
     scale_fill_manual(values=c("red3", "turquoise4"))
   lineplot2 <- lineplot2 + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")+
     scale_fill_manual(values=c("red3", "turquoise4"))
   
   rect_data <- sub %>%
     filter(Signal_type=="foreground") %>%
     group_by(interval) %>% 
     summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf
     ) %>%  filter(!is.na(interval))
   
   cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                                   linetype="solid", fill=alpha("black",0), color=alpha("black", 1),
                                   inherit.aes = FALSE, alpha = 0, size=0.5)
   
   lineplot <- lineplot +clustermids 
   
   
   lineplot2 <- lineplot2 + clustermids
   
   
   
   lineplot <- lineplot + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")
   lineplot2 <- lineplot2 + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")
   
   
   
   
   rect_data <- sub %>%
     filter(Signal_type=="foreground" & !is.na(interval) & !is.na(Peak_position_composite)) %>%
     group_by(interval, Peak_position_composite, Peak_position_composite_number) %>% filter(Peak_position_composite=="Edge") %>% 
     summarize(xmin = min(X_micron), xmax = max(X_micron), ymin = -Inf, ymax = Inf)
   
   
   
   
   cluster.rectangles <- geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=Peak_position_composite), 
                                   linetype="solid", color=alpha("blue", 0),
                                   inherit.aes=FALSE, alpha = 0.3, size=0.5)
   
   
   lineplot <- lineplot + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")
   lineplot2 <- lineplot2 + cluster.rectangles+
     labs(linetype="", x="Distance (\u03BCm)")+
     theme(legend.position = "none")
   
   
   
   lineplot <- lineplot + 
     theme(legend.position="top")+
     theme(legend.key.size = unit(legendsize, 'cm'), #change legend key size
           legend.key.height = unit(legendsize, 'cm'), #change legend key height
           legend.key.width = unit(legendsize*1.5, 'cm'))
   
   
   lineplot <- plot_grid(lineplot, lineplot2, nrow=2, rel_heights = c(0.6,0.4))
   
   print(lineplot)
   
   lineplot_compcluster <- lineplot
   
 }
 
 
 
 
 
 
 
   ## protrusions - plot of Gcx classification vs. clustersize (datapoint = one cluster)-----
  gcxvsclustersize <- custom_lineplot_with_rsquared(protrusion_data_summarized_percluster, 
                                                    regression_type="linear", color_var="Fiber_type" ,x_var=Gcx_normalized_vscluster, y_var=Clustersize)
  
  gcxvsclustersize_plot <- gcxvsclustersize[[1]] +custom_theme+
    theme(legend.position="top", legend.justification = "center")
  gcxvsclustersize_plot <- gcxvsclustersize_plot + xlab("Gcx relative to cluster") + ylab("Cluster size (micron)")
  print(gcxvsclustersize_plot)
  
   ## protrusions - histograms to see if classification of Gcx_intracellular_classification makes sense (mean + * times SD)------
  
  foreground_data <- subset(protrusion_data, Signal_type=="foreground")
    histograms <- create_histogram(foreground_data, "Gcx_intracellular", "number", num_facets=4, num_bins=25, factor=0.5)
    histograms <- histograms + ylab("Count scaled to 1") + xlab("Normalized Gcx signal") + custom_theme+ theme(legend.position="none")
    print(histograms)

    
    
   ## protrusions - violin plots of cluster number versus ratio------
    
    

    onlypeaks <- protrusion_data_summarized_percluster %>% filter(Signal_type=="foreground")
    sub <- onlypeaks
    sub <- sub %>%  filter(!is.na(PeakID_perprofile_grouped))
    
    peaknrvsratio <- custom_violinplot_with_statistics(data=sub, 
                                      group_vars = NULL,
                                      x_var="PeakID_perprofile_grouped", 
                                      y_var="Gcx_mean_leadingedge",
                                      log=FALSE, jitter=TRUE, annotate_signif = TRUE)
    
    peaknrvsratio_plot <- peaknrvsratio$plot
    peaknrvsratio_stats <- peaknrvsratio$statistical_test
    peaknrvsratio_plot <- peaknrvsratio_plot
      peaknrvsratio_plot <- peaknrvsratio_plot+custom_theme+
      theme(legend.position="none")
    print(peaknrvsratio_plot)
    
   ## protrusions - scatter plot of cluster number versus ratio----
    
    subset <- filter(protrusion_data, !is.na(new_Peak_ID) & Type!="Unspecified")
    
    onlypeaks$PeakID_perprofile <- as.numeric(onlypeaks$PeakID_perprofile)
    peaknrvsratio.scatter <- custom_lineplot_with_rsquared(data=subset, x_var=PeakID_perprofile, y_var=ratio_normalized_leadingedgeprotrusion, 
                                  regression_type="log", color_var=NULL, alpha=0.5, size=0.5,
                                  Orientation=NULL, Type=NULL, Fiber_type=NULL, Branching_type=NULL)
    
    peaknrvsratio.scatter.plot <- peaknrvsratio.scatter$plot
    peaknrvsratio.scatter.plot <- peaknrvsratio.scatter.plot+ ylim(-1,10) +
      ylab("B1 integrin/ Gcx ratio") + xlab("Cluster number from ECM interface")+custom_theme
    print(peaknrvsratio.scatter.plot)
    
   ## protrusions - scatter plot of cluster number versus distance --------
    
    peaknrvsXrel <- custom_lineplot_with_rsquared(data=sub, x_var=PeakID_perprofile, y_var=X_relative, 
                                  regression_type="log", color_var=NULL, alpha=0.5, size=0.5,
                                  Orientation="Orientation", Type=NULL, Fiber_type=NULL, Branching_type=NULL)
    peaknrvsXrel_plot <- peaknrvsXrel$plot
    peaknrvsXrel_plot <- peaknrvsXrel_plot + 
      geom_rect(aes(ymin=0.66, ymax=1, xmin=-Inf, xmax=Inf), fill="gray", alpha=0.01)+custom_theme
    print(peaknrvsXrel_plot)
    
   ## protrusions - violin plot of Length_label versus # of clusters-----
    
    
    nrofclusters <- protrusion_data %>% filter(Signal_type=="foreground") %>% 
      complete(Length_label, number, fill = list(count = 0)) %>%
      filter(Type == "Focalized") %>%
      group_by(Length_label, number) %>%
      summarize(count = n(), .groupd="drop")
    
    bar <- ggplot(data=nrofclusters, aes(x=Length_label, y=count))+
      geom_bar(stat="identity", color="black", fill="white")
    print(bar)
    
    
    
   ## protrusions - violin plot of Length_label versus ratio  etc -----
    
  
    protrusion_data_summarized_perROI <- protrusion_data_summarized_perROI %>%  filter(Signal_type=="foreground" & Length_label!="Intermediate zone")
    protrusion_data_summarized_perROI$Length_label <- as.factor(protrusion_data_summarized_perROI$Length_label)
    
    
    violin.ratios.perlengthlabel <- custom_violinplot_with_statistics(data=protrusion_data_summarized_perROI, 
                                      group_vars = NULL,
                                      x_var="Length_label", 
                                      y_var="ratio_mean_leadingedgeprotrusion",
                                      log=FALSE, jitter=TRUE, annotate_signif=TRUE, annotate_only_signif = FALSE)
    violin.ratios.perlengthlabel_plot <- violin.ratios.perlengthlabel$plot + ylab("B1 int / GCX ratio") + xlab("")
    violin.ratios.perlengthlabel_plot <- violin.ratios.perlengthlabel_plot +custom_theme
    print(violin.ratios.perlengthlabel_plot)
    
    violin.gcx <- custom_violinplot_with_statistics(data=protrusion_data_summarized_perROI, 
                                                                      group_vars = NULL,
                                                                      x_var="Length_label", 
                                                                      y_var="Gcx_mean_leadingedge",
                                                                      log=FALSE, jitter=TRUE, annotate_signif=TRUE, annotate_only_signif = FALSE)
    violin.gcx.plot <- violin.gcx$plot + ylab(" Normalized GCX") + xlab("")
    violin.gcx.plot <- violin.gcx.plot +custom_theme
    print(violin.gcx.plot)
    
    
    violin.b1int <- custom_violinplot_with_statistics(data=protrusion_data_summarized_perROI, 
                                                    group_vars = NULL,
                                                    x_var="Length_label", 
                                                    y_var="B1int_mean_leadingedge",
                                                    log=FALSE, jitter=TRUE, annotate_signif=TRUE, annotate_only_signif = FALSE)
    violin.b1int.plot <- violin.b1int$plot + ylab("Normalized B1 int") + xlab("")
    violin.b1int.plot <- violin.b1int.plot +custom_theme
    print(violin.b1int.plot)
    
  
    violins.supp <- plot_grid(violin.b1int.plot, violin.gcx.plot, ncol=2)
    print(violins.supp)
    
    
   ## protrusions - unfocalized versus focalized------
    
    ratio_cellbody <- mean(B1int_GCX_ratio$B1int[B1int_GCX_ratio$Orientation=="Cellbody"]) / mean(B1int_GCX_ratio$Gcx[B1int_GCX_ratio$Orientation=="Cellbody"])
    gcx_cellbody <- mean(mean(B1int_GCX_ratio$Gcx[B1int_GCX_ratio$Orientation=="Cellbody"]))
    
    protrusion_data <- protrusion_data %>%  
      mutate(ratio_normalized_cellbody_nocontact = B1int/Gcx / ratio_cellbody,
             Gcx_normalized_cellbody_nocontact = Gcx/ gcx_cellbody)
    
    
    protrusion_data_summarized_perROI <- protrusion_data %>% 
      group_by(Orientation, number, Signal_type, Length_label, Type) %>% 
      dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
    
    protrusion_data_summarized_perROI <- protrusion_data_summarized_perROI %>%  filter(Type !="Non-focalized")
    protrusion_data_summarized_perROI$Type <- as.factor(protrusion_data_summarized_perROI$Type)
    #protrusion_data_summarized_perROI$Type <- factor(protrusion_data_summarized_perROI$Type, levels= c("Focalized", "Non-focalized"))
    
    protrusion_data_summarized_perROI <- protrusion_data_summarized_perROI %>%  filter(!is.infinite(ratio_normalized_cellbody_nocontact))
    
    type.vs.ratio <- custom_violinplot_with_statistics(data=subset(protrusion_data_summarized_perROI, Length_label!="Intermediate zone"), 
                                      group_vars = "Length_label",
                                      x_var="Type", 
                                      y_var="ratio_normalized_cellbody_nocontact",
                                      log=FALSE, jitter=TRUE, annotate_signif=TRUE, annotate_only_signif = FALSE)
    type.vs.ratio.plot <- type.vs.ratio$plot
    type.vs.ratio.plot <- type.vs.ratio.plot + custom_theme + ylab("B1int / Gcx") + xlab("")
    print(type.vs.ratio.plot)
  
    
    

    
    
   ## protrusions (clusters) - histograms of Gcx and presentation of threshold values for Gcx classification------
    multiplier = 0.5
    clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, Length_label) %>%  
      dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
    clusterproperties <- clusterproperties %>%  filter(!is.na(new_Peak_ID))
    
    
    Gcx_mean <- mean(B1int_GCX_ratio$Gcx[protrusion_data$Signal_type=="foreground"], na.rm=TRUE)
    Gcx_sd <-  sd(B1int_GCX_ratio$Gcx[protrusion_data$Signal_type=="foreground"], na.rm=TRUE)
    
    hist.gcx <- ggplot(data=subset(B1int_GCX_ratio, Signal_type=="foreground"), aes(x=Gcx))+
      geom_histogram(bins=100, size=0.5, alpha=0.5)+
      geom_vline(xintercept=Gcx_mean)+
      geom_vline(xintercept=Gcx_mean-Gcx_sd*multiplier, linetype="dashed")+
      geom_vline(xintercept=Gcx_mean+Gcx_sd*multiplier, linetype="dashed")+
      
      theme_classic()
    
    print(hist.gcx)
    
   ## plot of distance normalized for maximum cluster distance versus Gcx and B1int normalized------
    multiplier = 0.5
    
    clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, Length_label) %>%  
      dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
    
    
    clusterproperties <- clusterproperties %>%  filter(!is.na(new_Peak_ID))
    
    
    
    clusterproperties <- clusterproperties %>%  ungroup() %>% 
      mutate(Gcx_mean_lengthlabel = Gcx_mean,
             Gcx_sd_lengthlabel = Gcx_sd
      ) %>% 
      group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, Length_label) %>% 
      mutate(Gcx_bottom_threshold = Gcx_mean_lengthlabel - Gcx_sd_lengthlabel* multiplier,
             Gcx_upper_threshold =  Gcx_mean_lengthlabel + Gcx_sd_lengthlabel* multiplier,
             Gcx_classification = ifelse(Gcx < Gcx_bottom_threshold, "Low", "Mid"),
             Gcx_classification = ifelse(Gcx > Gcx_upper_threshold, "High", Gcx_classification)
      ) %>% 
      dplyr::select(-c(Gcx_mean_lengthlabel, Gcx_sd_lengthlabel, Gcx_bottom_threshold, Gcx_upper_threshold))
    
    clusterproperties$Gcx_classification <- as.factor(clusterproperties$Gcx_classification)
    clusterproperties$Gcx_classification <- factor(clusterproperties$Gcx_classification, levels=c("Low", "Mid", "High"))
    
    Gcx.class.percluster <- clusterproperties %>%  group_by(new_Peak_ID, Gcx_classification) %>%  summarize()
    protrusion_data <- merge(protrusion_data, Gcx.class.percluster, by = "new_Peak_ID", all.x = TRUE)
    colnames(protrusion_data)[colnames(protrusion_data) == "Gcx_classification.y"] <- "Gcx_classification"
    protrusion_data <- protrusion_data[, !(colnames(protrusion_data) %in% "Gcx_classification.x")]
    
    
    
    
    
    protrusion_data <- protrusion_data %>% 
      group_by(number, Signal_type, new_Peak_ID, Clustersize, Gcx_classification) %>%  
      mutate(X_rel_cluster = (X_micron - min(X_micron)) / (max(X_micron) - min(X_micron)))
    
 
   ## clusters bordering on borders of Length_label are taken out of the analysis 
    protrusion_data <- protrusion_data %>% 
      group_by(number, Signal_type, new_Peak_ID, Clustersize, Gcx_classification) %>%  
      mutate(new_Peak_ID = ifelse(any(X_rel_cluster %in% c(0, 1) & WithinPeak != FALSE), NA, new_Peak_ID))
    
    
    normclustersize <- normclustersize  %>%  filter(!is.na(new_Peak_ID) & Clustersize>0.2 & !is.na(Gcx_classification))
    
    
    
    clusterprofile <- ggplot(subset(normclustersize, Gcx_classification!="Mid"), aes(x=X_rel_cluster, y=B1int_normalized_clusteredge))+
      
      geom_line(aes(x=X_rel_cluster, y=B1int_normalized_clusteredge, group=factor(new_Peak_ID)), size=0.5, color="magenta", alpha=0.2)+
      #geom_point(aes(x=X_rel_cluster, y=B1int_normalized_clusteredge), stroke=0, size=1, alpha=0.1, color="magenta")+
      geom_smooth(aes(x=X_rel_cluster, y=B1int_normalized_clusteredge), color="black", method="loess", span=0.5)+
      #geom_point(data=subset(normclustersize, WithinPeak==FALSE), aes(x=X_rel_cluster, y=B1int_normalized_clusteredge))+
      
      
      geom_vline(aes(xintercept=1),linetype="dashed")+
      geom_vline(aes(xintercept=0),linetype="dashed")+
      geom_vline(aes(xintercept=0.5), linetype="solid")+
      
      ylab("Normalized B1 int.")+
      
      facet_wrap(.~Gcx_classification)
    
    clusterprofile_gcx <- ggplot(subset(normclustersize, Gcx_classification!="Mid"), aes(x=X_rel_cluster, y=Gcx_normalized_clusteredge))+
      
      geom_line(aes(y=Gcx_normalized_clusteredge, group=new_Peak_ID), size=0.5, color="yellow3", alpha=0.5)+
      # geom_point(aes(y=Gcx_normalized_clusteredge), stroke=0, size=1, alpha=0.1, color="yellow3")+
      #geom_point(data=subset(normclustersize, WithinPeak==FALSE & Gcx_classification!="Mid"), aes(x=X_rel_cluster, y=Gcx_normalized_clusteredge))+
      geom_smooth(aes(y=Gcx_normalized_clusteredge), color="black", method="loess", span=0.5)+
      
      geom_vline(aes(xintercept=1),linetype="dashed")+
      geom_vline(aes(xintercept=0),linetype="dashed")+
      geom_vline(aes(xintercept=0.5), linetype="solid")+
      
      ylab("Normalized Gcx")+
      
      facet_wrap(.~Gcx_classification)
    
    clusterprofile <- clusterprofile + xlab("Relative distance") + ylab("Normalized fluorescence") + custom_theme + scale_x_continuous(breaks=c(0,0.5,1))
    print(clusterprofile)
    clusterprofile_gcx <- clusterprofile_gcx + xlab("Relative distance") + ylab("Normalized fluorescence") + custom_theme + scale_x_continuous(breaks=c(0,0.5,1))
    print(clusterprofile)
    
    clusterprofiles <- plot_grid(clusterprofile, clusterprofile_gcx, nrow=2)
    print(clusterprofiles)
    
    
   ## protrusions - Inside cluster subset -  scatter plots B1int vs Gcx intracluster + cluste rsize vs Gcx intracluster ----
    clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, Peak_position) %>%  
      dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
    
    
    
    Gcx.class.percluster <- normclustersize %>%  group_by(new_Peak_ID, Gcx_classification) %>%  summarize()
    clusterproperties <- merge(clusterproperties, Gcx.class.percluster, by = "new_Peak_ID", all.x = TRUE)
    colnames(clusterproperties)[colnames(clusterproperties) == "Gcx_classification.y"] <- "Gcx_classification"
    clusterproperties <- clusterproperties[, !(colnames(clusterproperties) %in% "Gcx_classification.x")]
      
      clusterproperties <- clusterproperties %>%  filter(Peak_position=="Inside" & Clustersize>0.2 & Signal_type=="foreground")
      clusterproperties <- filter_bad_values(clusterproperties, "Gcx_normalized_clusteredge", "B1int_normalized_clusteredge")
    
      relevant.subset <- subset(clusterproperties, Gcx_classification!="Mid")
      
      
      clustersize.scatterplot <- custom_lineplot_with_rsquared(data=relevant.subset, x_var=Gcx_normalized_clusteredge, y_var=B1int_normalized_clusteredge, 
                                  group_var="Gcx_classification", regression_type="linear")
    
    clustersize.scatterplot.stastistics <- clustersize.scatterplot[[2]]
    clustersize.scatterplot <- clustersize.scatterplot[[1]]
    clustersize.scatterplot <- clustersize.scatterplot + xlab("Gcx normalized against integrin cluster edge") + ylab("B1int normalized against edge")+custom_theme + 
       geom_vline(xintercept=1, linetype="dashed") + xlim(c(0,3)) + ylim(c(0,6))
    B1intvsGcx.edgenormalized.scatterplot <- clustersize.scatterplot
    print(B1intvsGcx.edgenormalized.scatterplot)
    
    relevant.subset <- subset(clusterproperties, Gcx_classification!="Mid")
    
    clustersize.scatterplot <- custom_lineplot_with_rsquared(data=relevant.subset, x_var=Gcx_normalized_clusteredge, y_var=Clustersize, 
                                                             group_var="Gcx_classification", regression_type="linear")
    
    clustersize.scatterplot.stastistics <- clustersize.scatterplot[[2]]
    clustersize.scatterplot <- clustersize.scatterplot[[1]]
    clustersize.scatterplot <- clustersize.scatterplot + xlab("Gcx normalized against integrin cluster edge") + ylab("Cluster width (micron)")+custom_theme + 
       geom_vline(xintercept=1, linetype="dashed") + xlim(c(0,3))
    print(clustersize.scatterplot)
   
    
    scatterplots.combined <- plot_grid(B1intvsGcx.edgenormalized.scatterplot, clustersize.scatterplot, ncol=2)
    print(scatterplots.combined)
    
    
   ## protrusions - violin plots - peak vs adjacent region - Gcx and B1int and ratio  ----
    
    clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize, Peak_position, Length_label) %>%  
      dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
    clusterproperties <- clusterproperties %>%  filter(Signal_type=="foreground", !is.na(Peak_position))
    
    
    
    #clusterproperties <- clusterproperties %>%  filter(Length_label!="Intermediate zone" & Length_label=="Contact")
    
    clusterproperties$Peak_position <- as.factor(clusterproperties$Peak_position)
    clusterproperties$Peak_position <- factor(clusterproperties$Peak_position, levels=c("Edge", "Inside"), labels=c("Adjacent", "Peak"))
    
    
    clusterproperties$Length_label <- as.factor(clusterproperties$Length_label)
    
    clusterproperties <- filter(clusterproperties, !is.na(Gcx_normalized_peak))
    clusterproperties <- filter(clusterproperties, !is.na(Gcx_classification))
    
    
    peakposition.gcx.mfi <- ggplot(clusterproperties, aes(x=Peak_position, y=Gcx))+
      geom_point(size=1, stroke=0)+
      geom_line(aes(group=new_Peak_ID))+
      scale_y_log10()+
      custom_theme
    print(peakposition.gcx.mfi)
    
    violin.cluster.gcx <- custom_violinplot_with_rsquared(data=subset(clusterproperties, Gcx_classification!="Mid"), group_vars=c("Gcx_classification"), x_var="Peak_position", 
                                                            y_var="Gcx_normalized_peak", jitter=TRUE)
    violin.cluster.gcx.statistics <- violin.cluster.gcx[[2]]
    violin.cluster.gcx <- violin.cluster.gcx[[1]]
    violin.cluster.gcx <- violin.cluster.gcx + xlab("") + ylab("Normalized Gcx") + custom_theme
    print(violin.cluster.gcx)
    
    
    violin.cluster.b1int <- custom_violinplot_with_statistics(data=subset(clusterproperties ,Gcx_classification !="Mid"), group_vars=c("Gcx_classification"), x_var="Peak_position",
                                                              y_var="B1int_normalized_peak", jitter=TRUE, y_adjust_max=2)
    violin.cluster.b1int.statistics <- violin.cluster.b1int[[2]]
    violin.cluster.b1int <- violin.cluster.b1int[[1]]
    violin.cluster.b1int <- violin.cluster.b1int + xlab("") + ylab("Normalized B1int") + custom_theme
    print(violin.cluster.b1int)
    
    
    
    violins.vsclusters <- plot_grid(violin.cluster.b1int, violin.cluster.gcx)
    print(violins.vsclusters)
    
   ## protrusions - Gcx classification vs Gcx normalized vs edge ------

    
        
    violin.gcx.classification <- custom_violinplot_with_statistics(data=clusterproperties, group_vars=NULL, x_var="Gcx_classification", 
                                                            y_var="Gcx_normalized_clusteredge", jitter=TRUE, annotate_signif=FALSE,
                                                            annotate_only_signif = FALSE, y_adjust_max = 2)
    violin.gcx.classification <- violin.gcx.classification$plot + custom_theme
    
    print(violin.gcx.classification)
      
   ## protrusions - violin plots - sd_B1int------
    clusterproperties$Length_label <- as.factor(clusterproperties$Length_label)
    
    violinplot.sd.thr <- custom_violinplot_with_statistics(clusterproperties, x_var="Length_label", 
                                                           y_var="B1int_sd", group_vars=NULL, y_adjust_max=2, jitter=TRUE)
    
    violinplot.sd.thr <- violinplot.sd.thr[[1]]
    violinplot.sd.thr <- violinplot.sd.thr + geom_hline(yintercept=B1int_sd_thr)+ ylab("B1int Stdev") + xlab("")+custom_theme
    print(violinplot.sd.thr)
    
   ## protrusions - clusters - distance versus ratio----
    clusterproperties <- protrusion_data %>%  group_by(number, Orientation, new_Peak_ID, Signal_type, Clustersize) %>%  
      dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE))
    
    clusterproperties$X_micron <- clusterproperties$X_micron + 0.0001
    
    clusterproperties <-clusterproperties %>%  filter(!is.na(Gcx_mean_leadingedge))
    
    distvsratio <- custom_lineplot_with_rsquared(data=clusterproperties, x_var=X_micron, y_var=Gcx_normalized_cellbody_nocontact, 
                                  group_var=NULL, regression_type="linear")
    distvsratio.plot <- distvsratio$plot + custom_theme + ylab("B1int / Gcx") + xlab("Distance (micron)")
    print(distvsratio.plot)
    
    
    
   ## scatterplot with smoothened means of X_micron versus ratio_leadingedge normalized (datapoint = Xlocations outside/inside cluster)-------
  subset_protrusion <- filter_bad_values(subset_protrusion, "X_micron", "ratio_normalized_leadingedgeprotrusion")
  protrusion_data_ratio <- custom_lineplot_with_rsquared(data=subset_protrusion, x_var=X_micron, y_var=ratio_normalized_leadingedgeprotrusion, 
                                                         regression_type="log", color_var="Xthr_pos", alpha=0.5, size=0.5,
                                                         Orientation="Orientation", Type=NULL, Fiber_type=NULL, Branching_type=NULL)
  
  protrusion_data_ratio_plot <- protrusion_data_ratio$plot + xlim(0.001,19.9) + ylim(0,50)+custom_theme+
    theme(legend.position="none")+ xlab("Distance (micron)") + ylab("B1 integrin/ Gcx ratio") 
  protrusion_data_ratio_regression <- protrusion_data_ratio$rsquared_values
  print(protrusion_data_ratio_plot)
  
  ## scatterplot with smoothened means of X_micron versus Gcx_leadingedge normalized (datapoint = Xlocations outside/inside cluster)
  subset_protrusion <- filter_bad_values(subset_protrusion, "X_micron", "Gcx_normalized_leadingedge")
  protrusion_data_gcx <- custom_lineplot_with_rsquared(data=subset_protrusion, x_var=X_micron, y_var=Gcx_normalized_leadingedge, 
                                                       regression_type="log", color_var="Xthr_pos", alpha=0.5,size=0.5,
                                                       Orientation="Orientation", Type=NULL, Fiber_type=NULL, Branching_type=NULL)
  
  protrusion_lineplot_gcx <- protrusion_data_gcx$plot +xlim(0.001,19.9) + ylim(0,6)+custom_theme + theme(legend.position = "none") + ylab("Normalized Gcx") + xlab("Distance (micron)")
  protrustion_regression_gcx <- protrusion_data_gcx$rsquared_values
  print(protrusion_lineplot_gcx)
  
  ## scatterplot with smoothened means of X_micron versus B1int_leadingedge normalized (datapoint = Xlocations outside/inside cluster)
  subset_protrusion <- filter_bad_values(subset_protrusion, "X_micron", "B1int_normalized_leadingedge")
  protrusion_data_b1int <- custom_lineplot_with_rsquared(data=subset_protrusion, x_var=X_micron, y_var=B1int_normalized_leadingedge, 
                                                       regression_type="log", color_var="Xthr_pos", alpha=0.5,size=0.5,
                                                       Orientation="Orientation", Type=NULL, Fiber_type=NULL, Branching_type=NULL)
  
  protrusion_lineplot_b1int <- protrusion_data_b1int$plot +xlim(0.001,19.9) + ylim(0,6)+custom_theme + theme(legend.position = "none") + ylab("Normalized B1int") + xlab("Distance (micron)")
  protrustion_regression_b1int <- protrusion_data_gcx$rsquared_values
  print(protrusion_lineplot_b1int)
  
  
  supp.lineplots.protrusion <- plot_grid(protrusion_lineplot_b1int, protrusion_lineplot_gcx, ncol=2)
  print(supp.lineplots.protrusion)

  
  data <- subset(subset_protrusion,!is.na(Gcx_overall_classification))
  protrusion_stats_ratio <- custom_violinplot_with_statistics(data=protrusion_data_summarized_percell, group_vars=NULL, x_var="Xthr_pos", 
                                                              y_var="Gcx_intercellular_cellbody", jitter=TRUE)
  
  violin_protrusion_ratio <- protrusion_stats_ratio$plot +
   
    scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    theme(legend.position = "bottom", legend.direction = "horizontal")+
    scale_x_discrete(labels=NULL, name=NULL)+custom_theme+ylim(0,5)
  protrusion_stats_tukey <- protrusion_stats_ratio$tukey_values
  print(violin_protrusion_ratio)
  
  
  
  
  
  ## violin plots with datapoint per ROI
  protrusion_data_summarized_perROI <- protrusion_data_summarized_perROI %>% filter(Signal_type=="foreground")

  protrusion_data_summarized_perROI <- filter_bad_values(protrusion_data_summarized_perROI, "Gcx_classification_intracellular", "ratio_normalized_leadingedgeprotrusion")
  protrusion_data_summarized_perROI$Gcx_classification_intracellular <-  factor(protrusion_data_summarized_perROI$Gcx_classification_intracellular,
                                                                                levels=c("Low", "High"))
  
  
  
  protrusion_stats_ratio_summarized <- custom_violinplot_with_statistics(data=protrusion_data_summarized_perROI, 
                                                               group_vars = "Xthr_pos",
                                                               x_var="Gcx_classification_intracellular", 
                                                               y_var="ratio_normalized_leadingedgeprotrusion",
                                                               log=FALSE, jitter=TRUE)
  
  violin_protrusion_ratio_summarized <- protrusion_stats_ratio_summarized$plot +
    
    scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    theme(legend.position = "bottom", legend.direction = "horizontal")+custom_theme+ylab("B1 integrin/ Gcx ratio") + xlab("Gcx level")+theme(legend.position="none")
  protrusion_stats_tukey <- protrusion_stats_ratio_summarized$tukey_values
  print(violin_protrusion_ratio_summarized)
  


  
  protrusion_stats_B1int <- custom_violinplot_with_statistics(data=protrusion_data_summarized_perROI, 
                                                               group_vars = "Xthr_pos",
                                                               x_var="Gcx_classification_intracellular", 
                                                               y_var="B1int_normalized_leadingedge",
                                                               log=FALSE, jitter=TRUE)
  violin_protrusion_B1int <- protrusion_stats_B1int$plot + ylab("Normalized B1 integrin")+ custom_theme+theme(legend.position="none")
  print(violin_protrusion_B1int)
  
  protrusion_stats_gcx <- custom_violinplot_with_statistics(data=protrusion_data_summarized_perROI, 
                                                             group_vars = "Xthr_pos",
                                                             x_var="Gcx_classification_intracellular", 
                                                             y_var="Gcx_normalized_leadingedge",
                                                             log=FALSE, jitter=TRUE)
  
  violin_protrusion_gcx <- protrusion_stats_gcx$plot + ylab("Normalized Gcx") + custom_theme+theme(legend.position="none")
  print(violin_protrusion_gcx)
  
  violin_protrusion_combined <- plot_grid(violin_protrusion_B1int, violin_protrusion_gcx, ncol=2)
  print(violin_protrusion_combined)
  
  
   ## violin plots of Gcx_classification versus Gcx_normalized and B1int_normalized against cluster edge ------
  gcx.perpeak.vs.peakposition <- custom_violinplot_with_statistics(protrusion_data_summarized_percluster, group_vars=NULL, 
                                    "Gcx_classification", "Gcx_normalized_clusteredge", log=FALSE, jitter=TRUE)
  gcx.perpeak.vs.peakposition <- gcx.perpeak.vs.peakposition[[1]]+ xlab("") + 
    theme(legend.position="none")+ ylab("Normalized Gcx")
  print(gcx.perpeak.vs.peakposition)
  b1int.perpeak.vs.peakposition <- custom_violinplot_with_statistics(protrusion_data_summarized_percluster, group_vars=NULL, 
                                    "Gcx_classification", "B1int_normalized_clusteredge", log=FALSE, jitter=TRUE)
  b1int.perpeak.vs.peakposition <- b1int.perpeak.vs.peakposition[[1]] + xlab("") + 
    theme(legend.position="none")+ ylab("Normalized B1 integrin")
  print(b1int.perpeak.vs.peakposition)
  violin.rel.to.FA <- plot_grid(b1int.perpeak.vs.peakposition, gcx.perpeak.vs.peakposition, ncol=2)
  print(violin.rel.to.FA)
  
  # violin plots of Fiber_type (Tense/Relaxed) versus Gcx_normalized and B1int_normalized
  fibertype.B1int.clusteredge <- custom_violinplot_with_statistics(protrusion_data_summarized_percluster, group_vars=NULL, 
                                     "Fiber_type", "B1int_normalized_clusteredge", log=FALSE, jitter=TRUE)

  fibertype.B1int.clusteredge <- fibertype.B1int.clusteredge[[1]]
  
  fibertype.gcx.clusteredge <- custom_violinplot_with_statistics(protrusion_data_summarized_percluster, group_vars=NULL, 
                                                                    "Fiber_type", "Gcx_normalized_clusteredge", log=FALSE, jitter=TRUE)
  fibertype.gcx.clusteredge <- fibertype.gcx.clusteredge[[1]]  
  
  violin.fibertype.rel.FA <- plot_grid(fibertype.gcx.clusteredge, fibertype.B1int.clusteredge, ncol=2)
  print(violin.fibertype.rel.FA)

## plot distance from cell start vs. B1int/gcx ratio for blebs------

bleb_ratio <- subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Bleb" & Signal_type=="foreground")
bleb_ratio$B1int_status <- as.factor(bleb_ratio$B1int_status)
bleb_ratio$Fiber_type <- as.factor(bleb_ratio$Fiber_type)
bleb_ratio$Orientation <- as.factor(bleb_ratio$Orientation)
bleb_ratio$Type <- as.factor(bleb_ratio$Type)
bleb_ratio$B1int_status <- factor(bleb_ratio$B1int_status, levels = c("Dim", "Unchanged", "Bright"))

bleb_ratio <- bleb_ratio %>%  group_by(number, Signal_type) %>%  mutate(Gcx_normalized_bleb=Gcx/ mean(Gcx[Type=="Cell body"]),
                                                                        B1int_normalized_bleb=B1int/ mean(B1int[Type=="Cell body"]))




bleb_ratio_summary <- bleb_ratio %>%
  group_by(number,Type,B1int_status) %>%
  summarise(
    mean_B1int = mean(B1int, na.rm = TRUE), # calculate mean, removing NA values
    mean_Gcx = mean(Gcx, na.rm = TRUE),
    mean_ratio = mean(B1int)/mean(Gcx) # calculate mean, removing NA values
  ) %>%  ungroup () %>% 
    mutate(ratio_corrected = mean_ratio / mean(mean_ratio[Type=="Cell body"], na.rm=TRUE),
           B1int_normalized = mean_B1int / mean(mean_B1int[Type=="Cell body"], na.rm=TRUE),
           Gcx_normalized = mean_Gcx / mean(mean_Gcx[Type=="Cell body"]), na.rm=TRUE)
# violin plots



bleb_stats_ratio <- custom_violinplot_with_anova_tukey(data=bleb_ratio_, 
                                                 group_vars = c("B1int_status"),
                                                 x_var="Type", 
                                                 y_var="ratio_normalized_bleb",
                                                 log=FALSE)

violin_bleb_ratio <- bleb_stats_ratio$plot

bleb_stats_b1int <- custom_violinplot_with_anova_tukey(bleb_ratio, 
                                   group_vars = c("B1int_status"),
                                   x_var="Type", 
                                   y_var="B1int_normalized_bleb",
                                   log=FALSE,
                                   jitter=TRUE)

violin_bleb_b1int <- bleb_stats_b1int$plot
print(violin_bleb_b1int)

# violin plots with datapoints per bleb


bleb_ratio_sub <- subset(bleb_ratio_summary,B1int_status!="Unchanged")

bleb_stats_ratio_perbleb <- custom_violinplot_with_anova_tukey(bleb_ratio_sub,
                                                             group_vars = c("B1int_status"),
                                                             x_var="Type", 
                                                             y_var="ratio_corrected",
                                                             log=FALSE,
                                                             jitter=TRUE)
violin_bleb_ratio_perbleb <- bleb_stats_ratio_perbleb$plot + geom_jitter(width=0.2, size=1, stroke=0)
print(violin_bleb_ratio_perbleb)


bleb_stats_B1int_perbleb <- custom_violinplot_with_anova_tukey(bleb_ratio_sub,
                                                             group_vars = c("B1int_status"),
                                                             x_var="Type", 
                                                             y_var="B1int_normalized",
                                                             log=FALSE,
                                                             jitter=TRUE)
bleb_stats_B1int_perbleb_tukey <- bleb_stats_B1int_perbleb$tukey_values
violin_bleb_B1int_perbleb <- bleb_stats_B1int_perbleb$plot + geom_jitter(width=0.2, size=1, stroke=0)
print(violin_bleb_B1int_perbleb)


bleb_stats_gcx_perbleb <- custom_violinplot_with_anova_tukey(bleb_ratio_sub,
                                                       group_vars = c("B1int_status"),
                                                       x_var="Type", 
                                                       y_var="Gcx_normalized",
                                                       log=FALSE)
violin_bleb_gcx_perbleb <- bleb_stats_gcx_perbleb$plot + geom_jitter(width=0.2, size=1, stroke=0)
print(violin_bleb_gcx_perbleb)

violin_bleb_perbleb_combined <- plot_grid(violin_bleb_ratio_perbleb, violin_bleb_B1int_perbleb, violin_bleb_gcx_perbleb, ncol=3)
print(violin_bleb_perbleb_combined)




anova_result <- performTwoWayANOVA(bleb_ratio, responseVar="ratio_normalized_bleb", factor1="Type", factor2="B1int_status")
print(anova_result)
B1int.negative.bleb <- subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Bleb" & B1int_status=="Negative")
regression.data.bleb <- calculate_residuals_and_tests_limitcorrected(B1int.negative.bleb$X_relative, B1int.negative.bleb$B1int)
print(regression.data.bleb)

## plot distance from cell start vs. B1int/gcx ratio for retraction fibers------


retraction_data <- subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber" 
                          & Signal_type=="foreground" & Type!="Unspecified")

retraction_data <- retraction_data %>%  group_by(number, Signal_type) %>%  mutate(Gcx_normalized_retractionfiber=
                                                                                    Gcx/ mean(Gcx[Fiber_type=="Cell body"]),
                                                                                  B1int_normalized_retractionfiber=
                                                                                    B1int/ mean(B1int[Fiber_type=="Cell body"]))

retraction_data$Fiber_type <- as.character(retraction_data$Fiber_type)
retraction_data$Fiber_type <- gsub(" ", "\n", retraction_data$Fiber_type)
retraction_data$Fiber_type <- as.factor(retraction_data$Fiber_type)

## lineplots

retract.data.lineplot <- custom_lineplot_with_rsquared(data=retraction_data, 
                                                       x_var=X_relative, 
                                                       y_var=ratio_normalized_retractionfiber, 
                                                      Orientation="Orientation", 
                                                       Type="Type", 
                                                       Fiber_type=NULL, 
                                                       Branching_type = NULL)


retract.data.lineplot_rsq <- retract.data.lineplot[[1]]
retract.data.lineplot_rsq <- retract.data.lineplot_rsq + xlim(0.001,1) + ylim(0, 300) + 
  geom_point(aes(color=Fiber_type), size=1, stroke=0)+ custom_theme +theme(legend.position="none")
print(retract.data.lineplot_rsq)


retract.data.lineplot.gcx <- custom_lineplot_with_rsquared(data=retraction_data, 
                                                       x_var=X_relative, 
                                                       y_var=Gcx_normalized_retractionfiber, 
                                                       Orientation="Orientation", 
                                                       Type="Type", 
                                                       Fiber_type=NULL, 
                                                       Branching_type = NULL)

retract.data.lineplot_gcx <- retract.data.lineplot.gcx$plot
retract.data.lineplot_gcx <- retract.data.lineplot_gcx + xlim(0.0001,1) + ylim(0, 1.25)+
  geom_point(aes(color=Fiber_type), size=1, stroke=0)+ custom_theme+ theme(legend.position="none")
print(retract.data.lineplot_gcx)


retract.data.data.B1int <- custom_lineplot_with_rsquared(data=retraction_data, 
                                                           x_var=X_relative, 
                                                           y_var=B1int_normalized_retractionfiber, 
                                                           Orientation="Orientation", 
                                                           Type="Type", 
                                                           Fiber_type=NULL, 
                                                           Branching_type = NULL)

retract.data.lineplot.B1int <- retract.data.data.B1int$plot
retract.data.lineplot.B1int <- retract.data.lineplot.B1int+
  
  geom_point(aes(color=Fiber_type), size=1, stroke=0)+ custom_theme+ theme(legend.position="none")
print(retract.data.lineplot.B1int)


retract_data_line_supp <- plot_grid(retract.data.lineplot.B1int, retract.data.lineplot_gcx)
print(retract_data_line_supp)



#changePlotSize(retract.data.lineplot_rsq)


retraction_stats <- custom_violinplot_with_anova_tukey(data=retraction_data, 
                                                         group_vars = c("Orientation", "Type"),
                                                         x_var="Fiber_type", 
                                                         y_var="ratio_normalized_retractionfiber",
                                                         log=TRUE)

violins_retraction <- retraction_stats$plot + ylab("B1 integrin / Gcx ratio") + geom_jitter(alpha=0.5, stroke=0, size=1)
print(violins_retraction)
tukey_values <- retraction_stats$tukey_values
print(tukey_values)


retraction_stats_gcx <- custom_violinplot_with_anova_tukey(data=retraction_data, 
                                                       group_vars = c("Orientation", "Type"),
                                                       x_var="Fiber_type", 
                                                       y_var="Gcx_normalized_retractionfiber",
                                                       log=FALSE)

violins_retraction_gcx <- retraction_stats_gcx$plot + ylab("Gcx MFI") + geom_jitter(alpha=0.5, stroke=0, size=1)
print(violings_retraction_gcx)
tukey_values_gcx <- retraction_stats_gcx$tukey_values


retraction_stats_B1int <- custom_violinplot_with_anova_tukey(data=retraction_data, 
                                                           group_vars = c("Orientation", "Type"),
                                                           x_var="Fiber_type", 
                                                           y_var="B1int_normalized_retractionfiber",
                                                           log=FALSE)

violins_retraction_B1int <- retraction_stats_B1int$plot + ylab("B1 integrin MFI") + geom_jitter(alpha=0.5, stroke=0, size=1)
print(violins_retraction_B1int)
tukey_values_B1int <- retraction_stats_B1int$tukey_values

supp.fig.retraction <- plot_grid(
                                 violins_retraction_B1int,
                                 violins_retraction_gcx, ncol=2 ) 
print(supp.fig.retraction)

## make an overview figure which is distance versus ratio------


distvsratio <- plot_grid(distvsratio_protrusion, distvsratio_bleb_combined, distvsratio_retraction,ncol=3)

print(distvsratio)

ggsave("Distvsratio.pdf", distvsratio, width = 20, height = 10, units = "cm")
saveRDS(distvsratio, paste0(workdir,"/","Distvsratio.rdata"))

## create violin plots for leading edge protrusion ----

distvsratio_protrusion_violin_ratio <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), 
                                        aes(x=Type, y=ratio_normalized_leadingedgeprotrusion)) +
  geom_violin(aes(fill=Fiber_type),scale="width", position=position_dodge(.9))+
  geom_boxplot(aes(fill=Fiber_type),width=0.2,outlier.shape=NA, position=position_dodge(.9))+
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  theme(legend.position = c(0.8,0.8)) +
  #ylim(0,8)+
  #geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.85, 0.55),
        legend.justification = c(1, 0),
        legend.box.just = "right",
        legend.title=element_blank())

print(distvsratio_protrusion_violin_ratio)

distvsratio_protrusion_violin_gcx <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), 
                                        aes(x=Type, y=Gcx)) +
  geom_violin(aes(fill=Fiber_type),scale="width", position=position_dodge(.9))+
  geom_boxplot(aes(fill=Fiber_type),width=0.2,outlier.shape=NA, position=position_dodge(.9))+
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  theme(legend.position = c(0.8,0.8)) +
  #ylim(0,8)+
  #geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.85, 0.55),
        legend.justification = c(1, 0),
        legend.box.just = "right",
        legend.title=element_blank())+
  theme(legend.position="none")

print(distvsratio_protrusion_violin_gcx)

distvsratio_protrusion_violin_b1int <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), 
                                            aes(x=Type, y=B1int)) +
  geom_violin(aes(fill=Fiber_type),scale="width", position=position_dodge(.9))+
  geom_boxplot(aes(fill=Fiber_type),width=0.2,outlier.shape=NA, position=position_dodge(.9))+
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  theme(legend.position = c(0.8,0.8)) +
  #ylim(0,8)+
  #geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.85, 0.55),
        legend.justification = c(1, 0),
        legend.box.just = "right",
        legend.title=element_blank())+
  theme(legend.position="none")

print(distvsratio_protrusion_violin_b1int)


protrusion.violins <- plot_grid(distvsratio_protrusion_violin_ratio, distvsratio_protrusion_violin_gcx, distvsratio_protrusion_violin_b1int, ncol=3)

print(protrusion.violins)



## plot distance from cell start vs. Gcx------

distvsratio_protrusion <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), aes(x=X_fromcellstart_micron, y=Gcx)) +
  geom_point() +
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black") +
 # labs(y="B1 Integrin / Glycocalyx ratio") +
  labs(x="Distance from most distal point -> cell body (micron)") +
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  theme(legend.position = c(0.7,0.7)) +
  geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  
  scale_x_continuous(breaks = c(0, 5, 10), limits=c(0.1,10)) 



distvsratio_protrusion <- distvsratio_protrusion +
  facet_nested(. ~ Orientation  + Fiber_type + Type)


distvsratio_bleb <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Bleb"), aes(x=X_fromcellstart_micron, y=Gcx, color=Type))+
  
  geom_point()+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black")+
  
 # labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Distance from most distal point -> cell body (micron)")+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = c(0.8, 0.8),  # Adjust the position of the legend
        legend.justification = c(1, 0),  # Align the legend to the top right corner
        legend.box.just = "right",  # Justify the legend box to the right
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

distvsratio_bleb <- distvsratio_bleb +
  facet_nested(. ~ Orientation ) 

distvsratio_retraction <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber"), aes(x=X_fromcellstart_micron, y=Gcx))+
  
  geom_point(aes(color=Fiber_type))+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black")+
  
  #labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Distance from most distal point fiber -> cell body (micron)")+
  
  #ylim(0,75)+
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  
  geom_vline(xintercept=Xthreshold, linetype=2)+
  
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = c(0.9, 0.8),  # Adjust the position of the legend
        legend.justification = c(1, 0),  # Align the legend to the top right corner
        legend.box.just = "right",  # Justify the legend box to the right
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

distvsratio_retraction <- distvsratio_retraction+
  facet_nested(. ~ Orientation + Type)


distvsratio <- plot_grid(distvsratio_bleb, distvsratio_protrusion, distvsratio_retraction,ncol=3)

print(distvsratio)

ggsave("Distvsgcx.pdf", distvsratio, width = 30, height = 10, units = "cm")
saveRDS(distvsratio, paste0(workdir,"/","Distvsgcx.rdata"))


## plot distance from cell start vs. B1int------

distvsratio_protrusion <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), aes(x=X_fromcellstart_micron, y=B1int)) +
  geom_point() +
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black") +
  # labs(y="B1 Integrin / Glycocalyx ratio") +
  labs(x="Distance from most distal point -> cell body (micron)") +
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  theme(legend.position = c(0.7,0.7)) +
  geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  
  scale_x_continuous(breaks = c(0, 5, 10), limits=c(0.1,10)) 



distvsratio_protrusion <- distvsratio_protrusion +
  facet_nested(. ~ Orientation  + Fiber_type + Type)


distvsratio_bleb <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Bleb"), aes(x=X_fromcellstart_micron, y=B1int, color=Type))+
  
  geom_point()+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black")+
  
  # labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Distance from most distal point -> cell body (micron)")+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = c(0.8, 0.8),  # Adjust the position of the legend
        legend.justification = c(1, 0),  # Align the legend to the top right corner
        legend.box.just = "right",  # Justify the legend box to the right
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

distvsratio_bleb <- distvsratio_bleb +
  facet_nested(. ~ Orientation ) 

distvsratio_retraction <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber"), aes(x=X_fromcellstart_micron, y=B1int))+
  
  geom_point(aes(color=Fiber_type))+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black")+
  
  #labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Distance from most distal point fiber -> cell body (micron)")+
  
  #ylim(0,75)+
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  
  geom_vline(xintercept=Xthreshold, linetype=2)+
  
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = c(0.9, 0.8),  # Adjust the position of the legend
        legend.justification = c(1, 0),  # Align the legend to the top right corner
        legend.box.just = "right",  # Justify the legend box to the right
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

distvsratio_retraction <- distvsratio_retraction+
  facet_nested(. ~ Orientation + Type)


distvsratio <- plot_grid(distvsratio_bleb, distvsratio_protrusion, distvsratio_retraction,ncol=3)

print(distvsratio)

ggsave("DistvsB1int.pdf", distvsratio, width = 30, height = 10, units = "cm")
saveRDS(distvsratio, paste0(workdir,"/","DistvsB1int.rdata"))



## plot distance vs gcx_normalized and B1int_normalized in one plot----------------

distvsratio_protrusion <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), aes(x=X_fromcellstart_micron, y=ratio)) +
  geom_point() +
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black") +
  labs(y="B1 Integrin / Glycocalyx ratio") +
  labs(x="Distance from most distal point -> cell body (micron)") +
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  theme(legend.position = c(0.7,0.7)) +
  geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  
  scale_x_continuous(breaks = c(0, 5, 10), limits=c(0.1,10)) 



distvsratio_protrusion <- distvsratio_protrusion +
  facet_nested(. ~ Orientation+  Fiber_type+Type)

print(distvsratio_protrusion)


distvsratio_bleb <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Bleb"), aes(x=X_fromcellstart_micron, y=ratio, color=Type))+
  
  geom_point()+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black")+
  
  labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Distance from most distal point -> cell body (micron)")+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = c(0.8, 0.8),  # Adjust the position of the legend
        legend.justification = c(1, 0),  # Align the legend to the top right corner
        legend.box.just = "right",  # Justify the legend box to the right
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

distvsratio_bleb <- distvsratio_bleb +
  facet_nested(. ~ Orientation ) 

distvsratio_retraction <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber"), aes(x=X_fromcellstart_micron, y=ratio))+
  
  geom_point(aes(color=Fiber_type))+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black")+
  
  labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Distance from most distal point fiber -> cell body (micron)")+
  
  #ylim(0,75)+
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  
  geom_vline(xintercept=Xthreshold, linetype=2)+
  
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        legend.position = c(0.9, 0.8),  # Adjust the position of the legend
        legend.justification = c(1, 0),  # Align the legend to the top right corner
        legend.box.just = "right",  # Justify the legend box to the right
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

distvsratio_retraction <- distvsratio_retraction+
  facet_nested(. ~ Orientation + Type)


distvsratio <- plot_grid(distvsratio_bleb, distvsratio_protrusion, distvsratio_retraction,ncol=3)

print(distvsratio)

ggsave("Distvscombined.pdf", distvsratio, width = 30, height = 10, units = "cm")
saveRDS(distvsratio, paste0(workdir,"/","Distvscombined.rdata"))


## plot distance vs gcx_normalized and B1int_normalized in one plot_interesting data combined----------------

distvsgcx_normalized <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start"), aes(x=X_fromcellstart))+
  
  geom_point(aes(y=Gcx_normalized), color="yellow3", alpha=0.2)+
  geom_point(aes(y=B1int_normalized), color="magenta", alpha=0.2)+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="yellow3", aes(y=Gcx_normalized))+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="magenta", aes(y=B1int_normalized))+
  
  labs(y="Normalized Glycocalyx or B1int MFI")+
  labs(x="Distance from protrusion \n start (micron)")+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  #theme(legend.position = c(0.7,0.7)) +
  

  geom_vline(xintercept=Xthreshold, linetype=2)+
  
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        #legend.position = "none",
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

distvsgcx_normalized <- distvsgcx_normalized+
  facet_nested(. ~ Orientation + Protrusiontype)

distvsgcx_normalized <- distvsgcx_normalized +
  #xlim(0.1,xlim)+
  ylim(-1,5)
  
  print(distvsgcx_normalized)

ggsave("Distvsintensity_normalized_merged.pdf", distvsgcx, width = 20, height = 10, units = "cm")
saveRDS(distvsgcx, paste0(workdir,"/","Distvsintensity_normalized_merged.rdata"))




## plot linearity_index vs. B1int/gcx ratio in line plot------

linearityvsratio <- ggplot(subset(Linearity_data,Type!="Cell start" & Orientation=="Leading edge protrusion"), aes(x=Linearity_index, y=mean_ratio))+
  geom_point()+
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black")+

  labs(y="B1-integrin / GCX ratio")+
  scale_colour_viridis_c(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 8),
        #legend.position = "none",
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),  # Add facet title borders
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  
scale_y_log10()


linearityvsratio <- linearityvsratio +
  facet_nested(. ~ Orientation + Protrusiontype +Type)
  



print(linearityvsratio)

ggsave("Linearityvsratio.pdf", device="pdf", path=workdir, width = 20, height = 10, units = "cm")
saveRDS(linearityvsratio, paste0(workdir,"/","Linearityvsratio.pdf.rdata"))



## plot and save all the plot profiles with the manually selected regions-----

file$X_micron <- (file$X*pxsize)-pxsize

plotwidth = 10
plotheight = 15

for(q in 1:length(unique(B1int_GCX_ratio$name))){
allnames <- unique(B1int_GCX_ratio$name)
plotprofile_sub <- subset(file, grepl(allnames[q],file$name) &Channel!="col")

ppname <- unique(plotprofile_sub$name)
ppname <- ppname[2]

orient <- unique(plotprofile_sub$Orientation)


mx <- subset(plotprofile_sub,Channel=="B1int")
mx <- max(mx$Value)
mn <- subset(plotprofile_sub,Channel=="gcx")
mn <- max(mn$Value)
scale <- mx/mn

thr <- subset(B1int_GCX_ratio, name==allnames[q])
thr <- thr$X
txt <- subset(B1int_GCX_ratio, name==allnames[q])
txt <- txt$Type
thr <- cbind(thr,txt)
thr <- as.data.frame(thr)
colnames(thr) <- c("value", "Type")
thr$value <- as.numeric(thr$value)



p <- ggplot(plotprofile_sub, aes(x=X_micron, y=Value, color=Channel))+
  geom_line(data= plotprofile_sub %>% filter(Channel=="B1int"),aes(y=Value),linewidth=1)+
  geom_line(data= plotprofile_sub %>% filter(Channel=="gcx"),aes(y=Value*scale),linewidth=1)+
  
  scale_color_manual(
    values = c("B1int" = "magenta", "gcx" = "yellow2"))+
  theme(panel.grid.major = element_line(color = "white",
                                        size = 1,
                                        linetype = 1))+
  xlab(paste0("Distance"," (\u00B5","m)"))+
  
  theme(legend.position="bottom")+
  
  ggtitle(paste(unique(ppname), "\n", unique(plotprofile_sub$Orientation)))+
  
  scale_y_continuous(
    # Features of the first axis
    name = "B1int signal",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*scale, name="Gcx signal")
  )

print(p)

ggsave(paste0(allnames[q],"_plotprofile.pdf"), device="pdf", path=workdir, width = plotwidth, height = plotheight, units = "cm")


    if(orient!="Bleb"){
    h <- p+
      geom_vline(data=subset(thr, Type!="Cell start"), 
                 aes(xintercept=value*pxsize, linetype=Type), linewidth=0.5, color="black", alpha=0.5)+
      
      geom_vline(data=subset(thr, Type=="Cell start"), 
                 aes(xintercept=value*pxsize), linewidth=0.5, color="purple", alpha=0.5)
    
    
    print(h)
    
    ggsave(paste0(allnames[q],"_plotprofile_withthresholds.pdf"), device="pdf", path=workdir, width = plotwidth, height = plotheight, units = "cm")
    
   
    
    
    g <- ggplot(plotprofile_sub, aes(x=X_micron, y=Value, color=Channel))+
      geom_line(data= plotprofile_sub %>% filter(Channel=="B1int"),aes(y=Value),linewidth=1)+
      geom_line(data= plotprofile_sub %>% filter(Channel=="gcx"),aes(y=Value*scale),linewidth=1)+
      
      
      scale_color_manual(
        values = c("B1int" = "magenta", "gcx" = "yellow2"))+
      
      
      xlab(paste0("Distance"," (\u00B5","m)"))+
      
      theme(legend.position="bottom")+
      
      ggtitle(paste(unique(ppname), "\n", unique(plotprofile_sub$Orientation)))+
      
      scale_y_continuous(
        # Features of the first axis
        name = "B1int signal",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*scale, name="Gcx signal")
      )+
      geom_vline(xintercept=0, linewidth=0.5, color="darkgreen", alpha=0.5, linetype=2)+
      geom_vline(xintercept=Xthreshold, linewidth=0.5, color="blue", alpha=0.5, linetype=2)
      if(orient=="Retraction fiber"){
        g <- g +
          geom_vline(xintercept=max(plotprofile_sub$X)*retraction_fiber_length_segregation_threshold*pxsize, 
                     linewidth=0.5, color="orange", alpha=0.5, linetype=2)
      }
      
    
    print(g)
    
    ggsave(paste0(allnames[q],"_plotprofile_cellstart.pdf"), device="pdf", path=workdir, width = plotwidth, height = plotheight, units = "cm")
    
    }else{
      h <- p+
        geom_vline(data=subset(thr, Type!="Cell start"), 
                   aes(xintercept=max(plotprofile_sub$X)*bleb_length_segregation_threshold*pxsize), 
                   linewidth=0.5, color="black", alpha=0.5, linetype=2)
      print(h)
      ggsave(paste0(allnames[q],"_plotprofile_withthresholds.pdf"), device="pdf", path=workdir, width = plotwidth, height = plotheigth, units = "cm")
    }

}




## calculate summarizing stats to pool for blebs------

pooled_data_bleb <- B1int_GCX_ratio %>%  subset(Orientation =="Bleb") %>%  group_by(Orientation, Type, number) %>%  summarize(mean_ratio = mean(ratio), stdev_ratio = sd(ratio), mean_gcx = mean(Gcx), mean_B1int = mean(B1int))
filtered_data <- pooled_data_bleb %>%
  filter(Type %in% c('Bleb', 'Cell body'))
ratio_data <- filtered_data %>%
  group_by(number, Orientation) %>%
  summarize(gcx_ratio = mean(mean_gcx[Type == 'Bleb']) / mean(mean_gcx[Type == 'Cell body']),
            B1int_ratio = mean(mean_B1int[Type == 'Bleb']) / mean(mean_B1int[Type == 'Cell body']))
bleb_integrinpos_thr <- 0.25

ratio_data$B1int_positive[ratio_data$B1int_ratio>bleb_integrinpos_thr] <- "Positive"
ratio_data$B1int_positive[ratio_data$B1int_ratio<bleb_integrinpos_thr] <- "Negative"

## plot summarized data

bleb_plot <- ggplot(ratio_data, aes(x=B1int_ratio, y = gcx_ratio))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)+
  geom_vline(xintercept=bleb_integrinpos_thr, linetype=2)

bleb_viol <- ggplot(ratio_data, aes(B1int_positive, y = gcx_ratio))+
  geom_point()+
  geom_violin(alpha=0.5)+
  geom_boxplot(alpha=0.5, width=0.2)

 

bleb_plot <- plot_grid(bleb_plot, bleb_viol)

print(bleb_plot)

ggsave("Bleb_plot_summary.pdf", bleb_plot, width = 20, height = 10, units = "cm")
saveRDS(bleb_plot, paste0(workdir,"/","Distvsratio.rdata"))


## plot examples of rel dist vs B1 int MFI for blebs----


distvsratio_bleb <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Bleb" & number == 8 | number ==9), aes(x=X_relative, y=ratio))+
  geom_point(size=2, stroke=0, alpha=0.5,aes(color=Type))+
  geom_line(size=0.5, aes(group=number, color=Type))+
  
  #geom_smooth(method="lm", formula = y~(x), fullrange=TRUE, se=TRUE) +
  
  #labs(y="B1 integrin / Glycocalyx ratio") +
  
  labs(x="Relative distance (fraction)")+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.95, 0.8),
        legend.title=element_blank(),
        strip.background = element_rect(color = "black", fill = NA, size = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.background = element_rect(fill = "transparent"))


distvsratio_bleb <- distvsratio_bleb +
  facet_nested(. ~ Orientation+B1int_status) +
  geom_vline(xintercept=0.33, linetype=2, size=0.5)


print(distvsratio_bleb)

## calculate summarizing stats to pool for retraction fibers------

pooled_data_bleb <- B1int_GCX_ratio %>%  subset(Orientation =="Retraction fiber") %>%  group_by(Orientation, Type, number) %>%  summarize(mean_ratio = mean(ratio), stdev_ratio = sd(ratio), mean_gcx = mean(Gcx), mean_B1int = mean(B1int))
filtered_data <- pooled_data_bleb %>%
  filter(Type %in% c('Retraction fiber', 'Cell body'))
ratio_data <- filtered_data %>%
  group_by(number, Orientation) %>%
  summarize(gcx_ratio = mean(mean_gcx[Type == 'Bleb']) / mean(mean_gcx[Type == 'Cell body']),
            B1int_ratio = mean(mean_B1int[Type == 'Bleb']) / mean(mean_B1int[Type == 'Cell body']))
bleb_integrinpos_thr <- 0.25

ratio_data$B1int_positive[ratio_data$B1int_ratio>bleb_integrinpos_thr] <- "Positive"
ratio_data$B1int_positive[ratio_data$B1int_ratio<bleb_integrinpos_thr] <- "Negative"

## plot summarized data

bleb_plot <- ggplot(ratio_data, aes(x=B1int_ratio, y = gcx_ratio))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)+
  geom_vline(xintercept=bleb_integrinpos_thr, linetype=2)

bleb_viol <- ggplot(ratio_data, aes(B1int_positive, y = gcx_ratio))+
  geom_point()+
  geom_violin(alpha=0.5)+
  geom_boxplot(alpha=0.5, width=0.2)



bleb_plot <- plot_grid(bleb_plot, bleb_viol)

print(bleb_plot)

ggsave("Bleb_plot_summary.pdf", bleb_plot, width = 20, height = 10, units = "cm")
saveRDS(bleb_plot, paste0(workdir,"/","Distvsratio.rdata"))



## delete unnecessary files-------


rm(ratio)
rm(ratio1)
rm(ratio2)
rm(set)
rm(sub)
rm(pp)
rm(pp_t)
rm(plotprofile_sub)
rm(plotprofile2)
rm(quant)
rm(quant1)
rm(quant2)
rm(stats)
rm(stats_ratio)
rm(coord_val)
rm(df)
rm(B1int_GCX_ratio2)
rm(df1)
rm(df2)
rm(unique_numbers)
rm(thr)
rm(stats2)
rm(stats)
##rm(result)
rm(h)
rm(p)
##rm(non_duplicates)
rm(newdata)



## plot profile example leading edge protrusion----

file <- file %>%  mutate(X_micron = X * pxsize)


overview <- file %>%  distinct(name)

library(stringr)

sub_protrusion <- file %>% 
  filter(str_detect(name, "Cell_7_Front_1_Airyscan Processing_3"))


# Calculate the maximum values for each channel
max_values_protrusion <- sub_protrusion %>%
  group_by(Channel) %>%
  summarize(max_value = max(Value))

# Create a new data frame with separate columns for each channel
sub_protrusion_wide <- sub_protrusion %>%
  pivot_wider(names_from = Channel, values_from = Value) 


# Create the plot with two y-axes and move the legend above
p_leadingedge <- ggplot(sub_protrusion_wide, aes(x = X_micron)) +
  geom_line(data = subset(sub_protrusion_wide, !is.na(B1int)), aes(y = B1int, color = "B1int"), size = 0.5) +
  geom_line(data = subset(sub_protrusion_wide, !is.na(gcx)), aes(y = gcx * (max_values_protrusion$max_value[1] / max_values_protrusion$max_value[3]), color = "gcx"), size = 0.5) +
  scale_color_manual(values = c("B1int" = "magenta", "gcx" = "yellow3")) +
  labs(
    y = "B1 integrin (x10^3)",
    color = "Channel"
  ) +
  scale_y_continuous(
    limits=c(0,2500),
    breaks=c(0,1000,2000,3000),
    labels=c(0,1,2,3),
    sec.axis = sec_axis(~./(max_values_protrusion$max_value[1] / max_values_protrusion$max_value[3])/1000, name = "GCX (x10^3)")
  ) +
  scale_x_continuous(
    breaks=c(0,1,2,3,4,5,6,7,8,9), 
    limits=c(3,9))+
  theme(legend.position = "top",  # Set legend position to "top"
        legend.box = "horizontal")  # Place legend items horizontally


p_leadingedge <- p_leadingedge + custom_theme
# Print the plot
print(p_leadingedge)





## plot profile example bleb----

## save examples for positive and negative beta integrin bleb -------

stats$X_micron <- stats$X * pxsize

B1int_GCX_ratio <- B1int_GCX_ratio %>%  group_by(number, Signal_type) %>%  mutate(X_relative = X/max(X)) 

sub_bleb <- B1int_GCX_ratio %>% 
  filter(CellNR==7 & RoiNR==5)

sub_bleb <- sub_bleb %>%  mutate(Type = ifelse(X_micron < (1/3* max(X_micron)), "Bleb", "Cell body") )

sub_bleb <- sub_bleb %>%  subset(Signal_type=="foreground")

sub_bleb <- sub_bleb %>%   mutate(B1int_normalized = B1int/mean(B1int[Type=="Cell body"])) %>% mutate(Gcx_normalized = Gcx/mean(Gcx[Type=="Cell body"])) 

standard_dev_mean <- sub_bleb %>%
  group_by(Type) %>%
  summarize(sd = sd(B1int_normalized), mean= mean(B1int_normalized), Xmin=min(X_micron), Xmax=max(X_micron),
            sd_gcx = sd(Gcx_normalized), mean_gcx = mean(Gcx_normalized))




# Create the plot with two y-axes and move the legend above
p_bleb  <- ggplot(sub_bleb) +
  geom_line(aes(x=X_relative, y=B1int_normalized),color="magenta")+
  geom_line(aes(x=X_relative, y=Gcx_normalized),color="yellow3")+
  geom_vline(xintercept = max(sub_bleb$X_relative)*0.33, linetype=2, color="black")+
  geom_hline(yintercept=1)+
  geom_rect(data = standard_dev_mean, aes(xmin = 0, 
                                 xmax =  1, 
                                 ymin = 1 - sd[Type=="Cell body"], 
                                 ymax = 1 + sd[Type=="Cell body"]), 
            fill = "black", alpha = 0.1) +
 
  labs(
    y = "Normalized Fluorescence",
  )+
  ylim(0,3)
  

p_bleb_negative <- p_bleb  + custom_theme
# Print the plot
print(p_bleb_negative)


sub_bleb <- B1int_GCX_ratio %>% 
  filter(CellNR==100, RoiNR==4)

sub_bleb <- sub_bleb %>%  mutate(Type = ifelse(X_micron < (1/3* max(X_micron)), "Bleb", "Cell body") )

sub_bleb <- sub_bleb %>%  subset(Signal_type=="foreground")

sub_bleb <- sub_bleb %>%   mutate(B1int_normalized = B1int/mean(B1int[Type=="Cell body"])) %>% mutate(Gcx_normalized = Gcx/mean(Gcx[Type=="Cell body"])) 

standard_dev_mean <- sub_bleb %>%
  group_by(Type) %>%
  summarize(sd = sd(B1int_normalized), mean= mean(B1int_normalized), Xmin=min(X_micron), Xmax=max(X_micron),
            sd_gcx = sd(Gcx_normalized), mean_gcx = mean(Gcx_normalized))




# Create the plot with two y-axes and move the legend above
p_bleb  <- ggplot(sub_bleb) +
  geom_line(aes(x=X_relative, y=B1int_normalized),color="magenta")+
  geom_line(aes(x=X_relative, y=Gcx_normalized),color="yellow3")+
  geom_vline(xintercept = max(sub_bleb$X_relative)*0.33, linetype=2, color="black")+
  geom_hline(yintercept=1)+
  geom_rect(data = standard_dev_mean, aes(xmin = 0, 
                                          xmax =  1, 
                                          ymin = 1 - sd[Type=="Cell body"], 
                                          ymax = 1 + sd[Type=="Cell body"]), 
            fill = "black", alpha = 0.1) +
  
  labs(
    y = "Normalized Fluorescence",
  )+
  ylim(0,3)


p_bleb_positive <- p_bleb  + custom_theme
# Print the plot
print(p_bleb_positive)

combined_bleb_examples <- plot_grid(p_bleb_positive, p_bleb_negative)
print(combined_bleb_examples)

## plot profile example retraction fiber ----



sub_retraction <- file %>% 
  filter(str_detect(name, "Cell_7_Retraction fiber_1_Airyscan Processing_2"))

sub_retraction <- sub_retraction %>%  mutate(Type = ifelse(X_micron < (4/5* max(X_micron)), "Retraction fiber", "Cell body") )


# Calculate the maximum values for each channel
max_values_retraction <- sub_retraction %>%
  group_by(Channel) %>%
  summarize(max_value = max(Value))

# Create a new data frame with separate columns for each channel
sub_retraction_wide <- sub_retraction %>%
  pivot_wider(names_from = Channel, values_from = Value) 


# Create the plot with two y-axes and move the legend above
p_retraction <- ggplot(sub_retraction_wide, aes(x = X_micron)) +
  
  geom_vline(xintercept = max(sub_retraction_wide$X_micron)*(4/5), linetype=2, color="black")+
  geom_line(data = subset(sub_retraction_wide, !is.na(B1int)), aes(y = B1int, color = "B1int"), size = 0.5) +
  geom_line(data = subset(sub_retraction_wide, !is.na(gcx)), aes(y = gcx * (max_values_retraction$max_value[1] / max_values_retraction$max_value[3]), color = "gcx"), size = 0.5) +
  scale_color_manual(values = c("B1int" = "magenta", "gcx" = "yellow3")) +
  
  
  labs(
    y = "B1 integrin (x10^3)",
    color = "Channel"
  ) +
  
  scale_y_continuous(
    labels=c(0,2,4,6),
    sec.axis = sec_axis(~./(max_values_retraction$max_value[1] / max_values_retraction$max_value[3])/1000, name = "GCX (x10^3)")
  ) +
  theme(legend.position = "top",  # Set legend position to "top"
        legend.box = "horizontal")  # Place legend items horizontally


p_retraction <- p_retraction + custom_theme
# Print the plot
print(p_retraction)

## combine profiles----
combined.profiles <- plot_grid(p_leadingedge, p_bleb, p_retraction,nrow=1, rel_widths = c(0.25,0.25,0.5))
print(combined.profiles)


## plot distance from cell start vs. B1int/gcx ratio for protrusion------

distvsratio_protrusion <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), 
                                 aes(x=X_fromcellstart_micron, y=Gcx)) +
  
  geom_point(size=2, stroke=0, alpha=0.5) +
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=TRUE, color="black", size=0.5) +
  
  labs(x="Distance from most distal point -> cell body (micron)") +
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  
  ylim(c(0, 3000))+
  #geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.85, 0.55),
        legend.justification = c(1, 0),
        legend.box.just = "right",
        legend.title=element_blank())+
  
  scale_x_continuous(breaks = c(0, 6,12), limits=c(0.1,15)) 



distvsratio_protrusion <- distvsratio_protrusion +
  facet_nested(. ~ Orientation+  Fiber_type+Type)

print(distvsratio_protrusion)

distvsratio_protrusion_gcx <- distvsratio_protrusion+ custom_theme

distvsratio_protrusion <- ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion"), 
                                 aes(x=X_fromcellstart_micron, y=B1int)) +
  
  geom_point(size=2, stroke=0, alpha=0.5) +
  geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=TRUE, color="black", size=0.5) +
  
  labs(x="Distance from most distal point -> cell body (micron)") +
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  
  ylim(c(0, 3000))+
  #geom_vline(xintercept=Xthreshold, linetype=2) +
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.85, 0.55),
        legend.justification = c(1, 0),
        legend.box.just = "right",
        legend.title=element_blank())+
  
  scale_x_continuous(breaks = c(0, 6,12), limits=c(0.1,15)) 



distvsratio_protrusion <- distvsratio_protrusion +
  facet_nested(. ~ Orientation+  Fiber_type+Type)

print(distvsratio_protrusion)

distvsratio_protrusion_b1int <- distvsratio_protrusion+custom_theme

distvsratio_rawvalues <- plot_grid(distvsratio_protrusion_gcx, distvsratio_protrusion_b1int,nrow=2)
print(distvsratio_rawvalues)


# Assuming you have already created your ggplot object distvsratio_protrusion_gcx
p_build <- ggplot_build(distvsratio_protrusion_gcx)

# Extract individual formulas used in geom_smooth layers for each facet
for (layer in p_build$plot$layers) {
  if ("geom_smooth" %in% layer$geom_params$geom) {
    formula <- layer$data[[1]]$formula
    print(formula)
  }
}

## plot distance from cell start vs. B1int/gcx ratio for bleb ----

bleb_ratio <- subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Bleb")

bleb_ratio$B1int_status <- as.factor(bleb_ratio$B1int_status)



distvsratio_bleb_gcx <-  ggplot(bleb_ratio, aes(x=X_relative, y=Gcx))+
  #geom_point(size=2, stroke=0, alpha=0.5,aes(color=Type))+
  geom_line(size=0.5, aes(group=number, color=Type))+
  
  #geom_smooth(method="lm", formula = y~(x), fullrange=TRUE, se=TRUE) +
  
  #labs(y="B1 integrin / Glycocalyx ratio") +
  
  labs(x="Relative distance (fraction)")+
  scale_x_continuous(breaks=c(0,0.5,1))+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.95, 0.8),
        legend.title=element_blank())


distvsratio_bleb_gcx <- distvsratio_bleb_gcx +xlab(NULL)+
  facet_nested(. ~ Orientation+B1int_status) + 

  geom_vline(xintercept=0.33, linetype=2, size=0.5)+custom_theme


print(distvsratio_bleb_gcx)

distvsratio_bleb_b1int <-  ggplot(bleb_ratio, aes(x=X_relative, y=B1int))+
  #geom_point(size=2, stroke=0, alpha=0.5,aes(color=Type))+
  geom_line(size=0.5, aes(group=number, color=Type))+
  
  #geom_smooth(method="lm", formula = y~(x), fullrange=TRUE, se=TRUE) +
  
  #labs(y="B1 integrin / Glycocalyx ratio") +
  
  labs(x="Relative distance (fraction)")+
  scale_x_continuous(breaks=c(0,0.5,1))+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.95, 0.8),
        legend.title=element_blank())


distvsratio_bleb_b1int <- distvsratio_bleb_b1int +
  facet_nested(. ~ Orientation+B1int_status) +theme(strip.text = element_blank(), strip.background = element_blank())+
  geom_vline(xintercept=0.33, linetype=2, size=0.5)+custom_theme


print(distvsratio_bleb_b1int)

combined_plot_bleb <- plot_grid(distvsratio_bleb_gcx, distvsratio_bleb_b1int, nrow=2)
print(combined_plot_bleb)


## plot distance from cell start vs. B1int/gcx ratio for retraction fibers------


retraction_data <- subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber")

distvsratio_retraction_gcx <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber"), 
                                  aes(x=X_relative, y = Gcx, color=Fiber_type))+
  geom_point(size=1, stroke=0, alpha=0.5) +
  
  geom_smooth(data=subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber" & Type=="Focalized"),
              method="lm", formula=y~x, fullrange=TRUE, color="black") +
  
  geom_smooth(data=subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber" & Type=="Non-focalized"),
              method="lm", formula=y~x, fullrange=TRUE, color="black") +
  
  #labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Relative distance (fraction)")+
  
  #ylim(0,50)+
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  scale_x_continuous(breaks=c(0,0.5,1))+
  scale_color_manual(values=c("#1661a8","#bd3333"))+
  
  #geom_vline(xintercept=Xthreshold, linetype=2)+
  
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.5, 0.35),
        legend.title=element_blank())

distvsratio_retraction_gcx <- distvsratio_retraction_gcx+
  facet_nested(.~Orientation+Type)+custom_theme

print(distvsratio_retraction_gcx)


distvsratio_retraction_b1int <-  ggplot(subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber"), 
                                  aes(x=X_relative, y = B1int, color=Fiber_type))+
  geom_point(size=1, stroke=0, alpha=0.5) +
  
  geom_smooth(data=subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber" & Type=="Focalized"),
              method="lm", formula = y~(x), fullrange=TRUE, color="black") +
  
  geom_smooth(data=subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Retraction fiber" & Type=="Non-focalized"),
              method="lm", formula = y~(x), fullrange=TRUE, color="black") +
  
  #labs(y="B1 Integrin / Glycocalyx ratio")+
  labs(x="Relative distance (fraction)")+
  
  #ylim(0,50)+
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  
  scale_x_continuous(breaks=c(0,0.5,1))+
  scale_color_manual(values=c("#1661a8","#bd3333"))+
  
  #geom_vline(xintercept=Xthreshold, linetype=2)+
  
  
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.5, 0.35),
        legend.title=element_blank())

distvsratio_retraction_b1int <- distvsratio_retraction_b1int+
  facet_nested(.~Orientation+Type)+custom_theme

print(distvsratio_retraction_b1int)


combined_plot <- plot_grid(distvsratio_retraction_gcx, distvsratio_retraction_b1int, ncol=2)
print(combined_plot)

## supplementary figures (raw values) -----
## calculate increase in gcx over distance ------
paste(print(max(leadingedge_data$Gcx[leadingedge_data$X_fromcellstart_micron>12 & leadingedge_data$X_fromcellstart_micron<15])/
      min(leadingedge_data$Gcx[leadingedge_data$X_fromcellstart_micron<1]))," is the difference of gcx abundance between begin and end")

## blebs violin plots with raw values ------

distvsratio_bleb <-  ggplot(bleb_ratio, aes(x=B1int_status, y=B1int, fill=Type))+
  #geom_jitter(size=2, stroke=0, position=position_jitterdodge(jitter.width = 0.2))+
  
  geom_violin(scale="width", position=position_dodge(.9))+
  geom_boxplot(width=0.2,outlier.shape=NA, position=position_dodge(.9))+
  
  #labs(y="B1 integrin / Glycocalyx ratio") +
  
  labs(x=" ")+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.95, 0.8),
        legend.title=element_blank())


distvsratio_bleb_b1int <- distvsratio_bleb +
  facet_nested(. ~ Orientation) +
  geom_vline(xintercept=0.33, linetype=2)+custom_theme


print(distvsratio_bleb_b1int)


distvsratio_bleb <-  ggplot(bleb_ratio, aes(x=B1int_status, y=Gcx, fill=Type))+
 # geom_jitter(size=2, stroke=0, position=position_jitterdodge(jitter.width = 0.2))+
  geom_violin(widht=0.2, scale="width", position=position_dodge(.9))+
  geom_boxplot(width=0.2,outlier.shape=NA, position=position_dodge(.9))+
  
  #labs(y="B1 integrin / Glycocalyx ratio") +
  
  labs(x=" ")+
  
  
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  theme(axis.text.y = element_text(size = 8),
        text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
        axis.text = element_text(size = 8),
        legend.position = c(0.95, 0.8),
        legend.title=element_blank())


distvsratio_bleb_gcx <- distvsratio_bleb +
  facet_nested(. ~ Orientation) +
  geom_vline(xintercept=0.33, linetype=2)+custom_theme


print(distvsratio_bleb_gcx)

violinplots_bleb <- plot_grid(distvsratio_bleb_b1int, distvsratio_bleb_gcx, nrow=2)

print(violinplots_bleb)



## statistics for violin plots blebs -----



library(multcomp)
anova_result <- aov(B1int ~ Type * B1int_status, data = bleb_ratio)
library(emmeans)
pairwise <- emmeans(anova_result, pairwise ~ Type*B1int_status)
aov.tukey.results <- summary(pairwise)
aov.tukey.results <- aov.tukey.results[[2]]
aov.tukey.results$Val <- rep("B1int", nrow(aov.tukey.results))
aov.tukey.results.bleb.b1int <- aov.tukey.results
print(aov.tukey.results.bleb.b1int)
write.csv2(aov.tukey.results.bleb.b1int, "anova.tukey.results.bleb.b1int")

library(multcomp)
anova_result <- aov(Gcx ~ Type * B1int_status, data = bleb_ratio)
library(emmeans)
pairwise <- emmeans(anova_result, pairwise ~ Type*B1int_status)
aov.tukey.results <- summary(pairwise)
aov.tukey.results <- aov.tukey.results[[2]]
aov.tukey.results$Val <- rep("Gcx", nrow(aov.tukey.results))
aov.tukey.results.bleb.gcx <- aov.tukey.results
print(aov.tukey.results.bleb.gcx)
write.csv2(aov.tukey.results.bleb.gcx, "anova.tukey.results.bleb.gcx")

library(multcomp)
anova_result <- aov(ratio_normalized_bleb ~ Type * B1int_status, data = bleb_ratio)
library(emmeans)
pairwise <- emmeans(anova_result, pairwise ~ Type*B1int_status)
aov.tukey.results <- summary(pairwise)
aov.tukey.results <- aov.tukey.results[[2]]
aov.tukey.results$Val <- rep("B1int-gcx_ratio", nrow(aov.tukey.results))
aov.tukey.results.bleb.gcx <- aov.tukey.results
print(aov.tukey.results.bleb.gcx)
write.csv2(aov.tukey.results.bleb.gcx, "anova.tukey.results.bleb.ratio")

## statistics for leading edge protrusions -----



ratio_leadingedge <- subset(B1int_GCX_ratio, Type!="Cell start" & Orientation =="Leading edge protrusion")


library(multcomp)
anova_result <- aov(ratio_normalized_leadingedgeprotrusion ~ Fiber_type * Type, data = ratio_leadingedge)
library(emmeans)
pairwise <- emmeans(anova_result, pairwise ~ Fiber_type * Type)
aov.tukey.results <- summary(pairwise)
aov.tukey.results <- aov.tukey.results[[2]]
aov.tukey.results$Val <- rep("B1int-gcx_ratio", nrow(aov.tukey.results))
aov.tukey.results.leadingedge.ratio <- aov.tukey.results
print(aov.tukey.results.leadingedge.ratio)
write.csv2(aov.tukey.results.leadingedge.ratio, "anova.tukey.results.leadingedge.ratio")


library(multcomp)
anova_result <- aov(Gcx ~ Fiber_type * Type, data = ratio_leadingedge)
library(emmeans)
pairwise <- emmeans(anova_result, pairwise ~ Fiber_type * Type)
aov.tukey.results <- summary(pairwise)
aov.tukey.results <- aov.tukey.results[[2]]
aov.tukey.results$Val <- rep("Gcx", nrow(aov.tukey.results))
aov.tukey.results.leadingedge.gcx <- aov.tukey.results
print(aov.tukey.results.leadingedge.gcx)
write.csv2(aov.tukey.results.leadingedge.gcx, "anova.tukey.results.leadingedge.gcx")

## statistical (regression) analysis for leading edge b1int/gcx ratio:------

leadingedge_data <- subset(B1int_GCX_ratio, Orientation=="Leading edge protrusion"& Type!="Cell start" & Signal_type=="foreground")

# Create a data frame containing all unique combinations of "Type" and "Fiber_type"
combinations <- expand.grid(Type = unique(leadingedge_data$Type), Fiber_type=unique(leadingedge_data$Fiber_type))


# Initialize a list to store regression results
regression_results <- list()

# Initialize a data frame to store problematic rows
problematic_rows <- data.frame()

# Loop through each combination and calculate regression
for (i in 1:nrow(combinations)) {
  subset_data <- leadingedge_data %>%
    filter(Type == combinations$Type[i] & Fiber_type == combinations$Fiber_type[i])
  
  print(i)
  
  if (nrow(subset_data) > 1) {  # Ensure there's enough data for regression
    # Extract x and y data
    x <- subset_data$X_micron
    y <- subset_data$ratio_normalized_leadingedgeprotrusion
    df_clean <- data.frame(X=x,
                           Y=y)
    columns_to_check <- c("X", "Y")
    df_clean <- df_clean <- df_clean[complete.cases(df_clean[columns_to_check]), ]
    x  <- df_clean$X
    y <- df_clean$Y
    data_clean <- data.frame(x, y)
    data_clean <- data_clean[data_clean$x > 0, ]
    
    # Check for missing and infinite values in x and y
    if (any(is.na(x)) || any(is.infinite(x)) || any(is.na(y)) || any(is.infinite(y))) {
      # Append problematic rows to the data frame
      problematic_subset <- subset_data[which(is.na(x) | is.infinite(x) | is.na(y) | is.infinite(y)), ]
      problematic_rows <- rbind(problematic_rows, problematic_subset)
      
      warning("Skipping regression due to missing or infinite values in x or y.")
    } else {
      # Call the modified function to perform the regression analysis
      regression_result <- calculate_residuals_and_tests_limitcorrected_noexponential(x, y)
      
      # Store the results in the list
      regression_results[[paste0(combinations$Type[i],"-", paste0(combinations$Fiber_type[i]))]] <- regression_result
    }
  } else {
    warning("Skipping regression due to insufficient data points.")
  }
}

# Initialize an empty list to store the combined dataframes
combined_dfs_list <- list()

# Iterate through each combination
for (comb_name in names(regression_results)) {
  # Get the regression results for the current combination
  regression_result <- regression_results[[comb_name]]
  
  # Create a dataframe with combination information
  combination_info <- data.frame(Combination = comb_name)
  
  # Bind the combination information to the regression results
  combined_df <- cbind(combination_info, regression_result)
  
  # Add the combined dataframe to the list
  combined_dfs_list[[comb_name]] <- combined_df
}

# Combine the list of dataframes into a single dataframe
regression.analysis.combined.leadingedgeprotrusions <- do.call(rbind, combined_dfs_list)

# Reset row names (if needed)
rownames(regression.analysis.combined.leadingedgeprotrusions) <- NULL

regression.analysis.combined.leadingedgeprotrusions$Var <- rep("B1int-gcx_ratio", nrow(regression.analysis.combined.leadingedgeprotrusions))

print(regression.analysis.combined.leadingedgeprotrusions)
write.csv(x=regression.analysis.combined.leadingedgeprotrusions, "Regression.analysis.results.leadingedgeprotrusion.ratio.csv")

## statistical (regression) analysis for leading edge b1int:------

leadingedge_data <- subset(B1int_GCX_ratio, Orientation=="Leading edge protrusion"& Type!="Cell start")

# Create a data frame containing all unique combinations of "Type" and "Fiber_type"
combinations <- expand.grid(Type = unique(leadingedge_data$Type), Fiber_type=unique(leadingedge_data$Fiber_type))


# Initialize a list to store regression results
regression_results <- list()

# Initialize a data frame to store problematic rows
problematic_rows <- data.frame()

# Loop through each combination and calculate regression
for (i in 1:nrow(combinations)) {
  subset_data <- leadingedge_data %>%
    filter(Type == combinations$Type[i] & Fiber_type == combinations$Fiber_type[i])
  
  print(i)
  
  if (nrow(subset_data) > 1) {  # Ensure there's enough data for regression
    # Extract x and y data
    x <- subset_data$X_fromcellstart_micron
    y <- subset_data$B1int
    df_clean <- data.frame(X=x,
                           Y=y)
    columns_to_check <- c("X", "Y")
    df_clean <- df_clean <- df_clean[complete.cases(df_clean[columns_to_check]), ]
    x  <- df_clean$X
    y <- df_clean$Y
    data_clean <- data.frame(x, y)
    data_clean <- data_clean[data_clean$x > 0, ]
    
    # Check for missing and infinite values in x and y
    if (any(is.na(x)) || any(is.infinite(x)) || any(is.na(y)) || any(is.infinite(y))) {
      # Append problematic rows to the data frame
      problematic_subset <- subset_data[which(is.na(x) | is.infinite(x) | is.na(y) | is.infinite(y)), ]
      problematic_rows <- rbind(problematic_rows, problematic_subset)
      
      warning("Skipping regression due to missing or infinite values in x or y.")
    } else {
      # Call the modified function to perform the regression analysis
      regression_result <- calculate_residuals_and_tests_limitcorrected_noexponential(x, y)
      
      # Store the results in the list
      regression_results[[paste0(combinations$Type[i],"-", paste0(combinations$Fiber_type[i]))]] <- regression_result
    }
  } else {
    warning("Skipping regression due to insufficient data points.")
  }
}

# Initialize an empty list to store the combined dataframes
combined_dfs_list <- list()

# Iterate through each combination
for (comb_name in names(regression_results)) {
  # Get the regression results for the current combination
  regression_result <- regression_results[[comb_name]]
  
  # Create a dataframe with combination information
  combination_info <- data.frame(Combination = comb_name)
  
  # Bind the combination information to the regression results
  combined_df <- cbind(combination_info, regression_result)
  
  # Add the combined dataframe to the list
  combined_dfs_list[[comb_name]] <- combined_df
}

# Combine the list of dataframes into a single dataframe
regression.analysis.combined.leadingedgeprotrusions <- do.call(rbind, combined_dfs_list)

# Reset row names (if needed)
rownames(regression.analysis.combined.leadingedgeprotrusions) <- NULL

regression.analysis.combined.leadingedgeprotrusions$Var <- rep("B1int", nrow(regression.analysis.combined.leadingedgeprotrusions))

print(regression.analysis.combined.leadingedgeprotrusions)
write.csv(x=regression.analysis.combined.leadingedgeprotrusions, "Regression.analysis.results.leadingedgeprotrusion.b1int.csv")



## statistical (regression) analysis for leading edge gcx:------

leadingedge_data <- subset(B1int_GCX_ratio, Orientation=="Leading edge protrusion"& Type!="Cell start")

# Create a data frame containing all unique combinations of "Type" and "Fiber_type"
combinations <- expand.grid(Type = unique(leadingedge_data$Type), Fiber_type=unique(leadingedge_data$Fiber_type))


# Initialize a list to store regression results
regression_results <- list()

# Initialize a data frame to store problematic rows
problematic_rows <- data.frame()

# Loop through each combination and calculate regression
for (i in 1:nrow(combinations)) {
  subset_data <- leadingedge_data %>%
    filter(Type == combinations$Type[i] & Fiber_type == combinations$Fiber_type[i])
  
  print(i)
  
  if (nrow(subset_data) > 1) {  # Ensure there's enough data for regression
    # Extract x and y data
    x <- subset_data$X_fromcellstart_micron
    y <- subset_data$Gcx
    df_clean <- data.frame(X=x,
                           Y=y)
    columns_to_check <- c("X", "Y")
    df_clean <- df_clean <- df_clean[complete.cases(df_clean[columns_to_check]), ]
    x  <- df_clean$X
    y <- df_clean$Y
    data_clean <- data.frame(x, y)
    data_clean <- data_clean[data_clean$x > 0, ]
    
    # Check for missing and infinite values in x and y
    if (any(is.na(x)) || any(is.infinite(x)) || any(is.na(y)) || any(is.infinite(y))) {
      # Append problematic rows to the data frame
      problematic_subset <- subset_data[which(is.na(x) | is.infinite(x) | is.na(y) | is.infinite(y)), ]
      problematic_rows <- rbind(problematic_rows, problematic_subset)
      
      warning("Skipping regression due to missing or infinite values in x or y.")
    } else {
      # Call the modified function to perform the regression analysis
      regression_result <- calculate_residuals_and_tests_limitcorrected_noexponential(x, y)
      
      # Store the results in the list
      regression_results[[paste0(combinations$Type[i],"-", paste0(combinations$Fiber_type[i]))]] <- regression_result
    }
  } else {
    warning("Skipping regression due to insufficient data points.")
  }
}

# Initialize an empty list to store the combined dataframes
combined_dfs_list <- list()

# Iterate through each combination
for (comb_name in names(regression_results)) {
  # Get the regression results for the current combination
  regression_result <- regression_results[[comb_name]]
  
  # Create a dataframe with combination information
  combination_info <- data.frame(Combination = comb_name)
  
  # Bind the combination information to the regression results
  combined_df <- cbind(combination_info, regression_result)
  
  # Add the combined dataframe to the list
  combined_dfs_list[[comb_name]] <- combined_df
}

# Combine the list of dataframes into a single dataframe
regression.analysis.combined.leadingedgeprotrusions <- do.call(rbind, combined_dfs_list)

# Reset row names (if needed)
rownames(regression.analysis.combined.leadingedgeprotrusions) <- NULL

regression.analysis.combined.leadingedgeprotrusions$Var <- rep("Gcx", nrow(regression.analysis.combined.leadingedgeprotrusions))

print(regression.analysis.combined.leadingedgeprotrusions)
write.csv(x=regression.analysis.combined.leadingedgeprotrusions, "Regression.analysis.results.leadingedgeprotrusion.gcx.csv")




## statistical (regression) analysis for retraction fiber b1int/gcx ratio:------

retractionfiber_data <- subset(B1int_GCX_ratio, Orientation=="Retraction fiber"& Type!="Cell start")

# Create a data frame containing all unique combinations of "Type" and "Fiber_type"
combinations <- expand.grid(Type = unique(retractionfiber_data$Type))


# Initialize a list to store regression results
regression_results <- list()

# Initialize a data frame to store problematic rows
problematic_rows <- data.frame()

# Loop through each combination and calculate regression
for (i in 1:nrow(combinations)) {
  subset_data <- retractionfiber_data %>%
    filter(Type == combinations$Type[i])
  
  print(i)
  
  if (nrow(subset_data) > 1) {  # Ensure there's enough data for regression
    # Extract x and y data
    x <- subset_data$X_fromcellstart_micron
    y <- subset_data$ratio_normalized_retractionfiber
    df_clean <- data.frame(X=x,
                           Y=y)
    columns_to_check <- c("X", "Y")
    df_clean <- df_clean <- df_clean[complete.cases(df_clean[columns_to_check]), ]
    x  <- df_clean$X
    y <- df_clean$Y
    data_clean <- data.frame(x, y)
    data_clean <- data_clean[data_clean$x > 0, ]
    
    # Check for missing and infinite values in x and y
    if (any(is.na(x)) || any(is.infinite(x)) || any(is.na(y)) || any(is.infinite(y))) {
      # Append problematic rows to the data frame
      problematic_subset <- subset_data[which(is.na(x) | is.infinite(x) | is.na(y) | is.infinite(y)), ]
      problematic_rows <- rbind(problematic_rows, problematic_subset)
      
      warning("Skipping regression due to missing or infinite values in x or y.")
    } else {
      # Call the modified function to perform the regression analysis
      regression_result <- calculate_residuals_and_tests_limitcorrected_noexponential(x, y)
      
      # Store the results in the list
      regression_results[[paste0(combinations$Type[i])]] <- regression_result
    }
  } else {
    warning("Skipping regression due to insufficient data points.")
  }
}

# Initialize an empty list to store the combined dataframes
combined_dfs_list <- list()

# Iterate through each combination
for (comb_name in names(regression_results)) {
  # Get the regression results for the current combination
  regression_result <- regression_results[[comb_name]]
  
  # Create a dataframe with combination information
  combination_info <- data.frame(Combination = comb_name)
  
  # Bind the combination information to the regression results
  combined_df <- cbind(combination_info, regression_result)
  
  # Add the combined dataframe to the list
  combined_dfs_list[[comb_name]] <- combined_df
}

# Combine the list of dataframes into a single dataframe
regression.analysis.combined.retractionfibers <- do.call(rbind, combined_dfs_list)

# Reset row names (if needed)
rownames(regression.analysis.combined.retractionfibers) <- NULL

regression.analysis.combined.retractionfibers$Var <- rep("B1int-Gcx_ratio", nrow(regression.analysis.combined.retractionfibers))


write.csv(x=regression.analysis.combined, "Regression.analysis.results.retractionfiber.ratio.csv")


retractionfiber_data <- subset(B1int_GCX_ratio, Orientation=="Retraction fiber"& Type!="Cell start")

# Create a data frame containing all unique combinations of "Type" and "Fiber_type"
combinations <- expand.grid(Type = unique(retractionfiber_data$Type))


# Initialize a list to store regression results
regression_results <- list()

# Initialize a data frame to store problematic rows
problematic_rows <- data.frame()

# Loop through each combination and calculate regression
for (i in 1:nrow(combinations)) {
  subset_data <- retractionfiber_data %>%
    filter(Type == combinations$Type[i])
  
  print(i)
  
  if (nrow(subset_data) > 1) {  # Ensure there's enough data for regression
    # Extract x and y data
    x <- subset_data$X_fromcellstart_micron
    y <- subset_data$ratio_normalized_retractionfiber
    df_clean <- data.frame(X=x,
                           Y=y)
    columns_to_check <- c("X", "Y")
    df_clean <- df_clean <- df_clean[complete.cases(df_clean[columns_to_check]), ]
    x  <- df_clean$X
    y <- df_clean$Y
    data_clean <- data.frame(x, y)
    data_clean <- data_clean[data_clean$x > 0, ]
    
    # Check for missing and infinite values in x and y
    if (any(is.na(x)) || any(is.infinite(x)) || any(is.na(y)) || any(is.infinite(y))) {
      # Append problematic rows to the data frame
      problematic_subset <- subset_data[which(is.na(x) | is.infinite(x) | is.na(y) | is.infinite(y)), ]
      problematic_rows <- rbind(problematic_rows, problematic_subset)
      
      warning("Skipping regression due to missing or infinite values in x or y.")
    } else {
      # Call the modified function to perform the regression analysis
      regression_result <- calculate_residuals_and_tests_limitcorrected_noexponential(x, y)
      
      # Store the results in the list
      regression_results[[paste0(combinations$Type[i])]] <- regression_result
    }
  } else {
    warning("Skipping regression due to insufficient data points.")
  }
}

# Initialize an empty list to store the combined dataframes
combined_dfs_list <- list()

# Iterate through each combination
for (comb_name in names(regression_results)) {
  # Get the regression results for the current combination
  regression_result <- regression_results[[comb_name]]
  
  # Create a dataframe with combination information
  combination_info <- data.frame(Combination = comb_name)
  
  # Bind the combination information to the regression results
  combined_df <- cbind(combination_info, regression_result)
  
  # Add the combined dataframe to the list
  combined_dfs_list[[comb_name]] <- combined_df
}

# Combine the list of dataframes into a single dataframe
regression.analysis.combined.retractionfibers <- do.call(rbind, combined_dfs_list)

# Reset row names (if needed)
rownames(regression.analysis.combined.retractionfibers) <- NULL

write.csv(x=regression.analysis.combined.retractionfibers, "Regression.analysis.results.retractionfiber.csv")


## statistical (regression) analysis for retraction fiber gcx:------

retractionfiber_data <- subset(B1int_GCX_ratio, Orientation=="Retraction fiber"& Type!="Cell start")

# Create a data frame containing all unique combinations of "Type" and "Fiber_type"
combinations <- expand.grid(Type = unique(retractionfiber_data$Type))


# Initialize a list to store regression results
regression_results <- list()

# Initialize a data frame to store problematic rows
problematic_rows <- data.frame()

# Loop through each combination and calculate regression
for (i in 1:nrow(combinations)) {
  subset_data <- retractionfiber_data %>%
    filter(Type == combinations$Type[i])
  
  print(i)
  
  if (nrow(subset_data) > 1) {  # Ensure there's enough data for regression
    # Extract x and y data
    x <- subset_data$X_fromcellstart_micron
    y <- subset_data$Gcx
    df_clean <- data.frame(X=x,
                           Y=y)
    columns_to_check <- c("X", "Y")
    df_clean <- df_clean <- df_clean[complete.cases(df_clean[columns_to_check]), ]
    x  <- df_clean$X
    y <- df_clean$Y
    data_clean <- data.frame(x, y)
    data_clean <- data_clean[data_clean$x > 0, ]
    
    # Check for missing and infinite values in x and y
    if (any(is.na(x)) || any(is.infinite(x)) || any(is.na(y)) || any(is.infinite(y))) {
      # Append problematic rows to the data frame
      problematic_subset <- subset_data[which(is.na(x) | is.infinite(x) | is.na(y) | is.infinite(y)), ]
      problematic_rows <- rbind(problematic_rows, problematic_subset)
      
      warning("Skipping regression due to missing or infinite values in x or y.")
    } else {
      # Call the modified function to perform the regression analysis
      regression_result <- calculate_residuals_and_tests_limitcorrected_noexponential(x, y)
      
      # Store the results in the list
      regression_results[[paste0(combinations$Type[i])]] <- regression_result
    }
  } else {
    warning("Skipping regression due to insufficient data points.")
  }
}

# Initialize an empty list to store the combined dataframes
combined_dfs_list <- list()

# Iterate through each combination
for (comb_name in names(regression_results)) {
  # Get the regression results for the current combination
  regression_result <- regression_results[[comb_name]]
  
  # Create a dataframe with combination information
  combination_info <- data.frame(Combination = comb_name)
  
  # Bind the combination information to the regression results
  combined_df <- cbind(combination_info, regression_result)
  
  # Add the combined dataframe to the list
  combined_dfs_list[[comb_name]] <- combined_df
}

# Combine the list of dataframes into a single dataframe
regression.analysis.combined.retractionfibers <- do.call(rbind, combined_dfs_list)

# Reset row names (if needed)
rownames(regression.analysis.combined.retractionfibers) <- NULL

regression.analysis.combined.retractionfibers$Var <- rep("Gcx", nrow(regression.analysis.combined.retractionfibers))

write.csv(x=regression.analysis.combined.retractionfibers, "Regression.analysis.results.retractionfiber.gcx.csv")



## statistical (regression) analysis for retraction fiber B1int:------

retractionfiber_data <- subset(B1int_GCX_ratio, Orientation=="Retraction fiber"& Type!="Cell start")

# Create a data frame containing all unique combinations of "Type" and "Fiber_type"
combinations <- expand.grid(Type = unique(retractionfiber_data$Type))


# Initialize a list to store regression results
regression_results <- list()

# Initialize a data frame to store problematic rows
problematic_rows <- data.frame()

# Loop through each combination and calculate regression
for (i in 1:nrow(combinations)) {
  subset_data <- retractionfiber_data %>%
    filter(Type == combinations$Type[i])
  
  print(i)
  
  if (nrow(subset_data) > 1) {  # Ensure there's enough data for regression
    # Extract x and y data
    x <- subset_data$X_fromcellstart_micron
    y <- subset_data$B1int
    df_clean <- data.frame(X=x,
                           Y=y)
    columns_to_check <- c("X", "Y")
    df_clean <- df_clean <- df_clean[complete.cases(df_clean[columns_to_check]), ]
    x  <- df_clean$X
    y <- df_clean$Y
    data_clean <- data.frame(x, y)
    data_clean <- data_clean[data_clean$x > 0, ]
    
    # Check for missing and infinite values in x and y
    if (any(is.na(x)) || any(is.infinite(x)) || any(is.na(y)) || any(is.infinite(y))) {
      # Append problematic rows to the data frame
      problematic_subset <- subset_data[which(is.na(x) | is.infinite(x) | is.na(y) | is.infinite(y)), ]
      problematic_rows <- rbind(problematic_rows, problematic_subset)
      
      warning("Skipping regression due to missing or infinite values in x or y.")
    } else {
      # Call the modified function to perform the regression analysis
      regression_result <- calculate_residuals_and_tests_limitcorrected_noexponential(x, y)
      
      # Store the results in the list
      regression_results[[paste0(combinations$Type[i])]] <- regression_result
    }
  } else {
    warning("Skipping regression due to insufficient data points.")
  }
}

# Initialize an empty list to store the combined dataframes
combined_dfs_list <- list()

# Iterate through each combination
for (comb_name in names(regression_results)) {
  # Get the regression results for the current combination
  regression_result <- regression_results[[comb_name]]
  
  # Create a dataframe with combination information
  combination_info <- data.frame(Combination = comb_name)
  
  # Bind the combination information to the regression results
  combined_df <- cbind(combination_info, regression_result)
  
  # Add the combined dataframe to the list
  combined_dfs_list[[comb_name]] <- combined_df
}

# Combine the list of dataframes into a single dataframe
regression.analysis.combined.retractionfibers <- do.call(rbind, combined_dfs_list)

# Reset row names (if needed)
rownames(regression.analysis.combined.retractionfibers) <- NULL

regression.analysis.combined.retractionfibers$Var <- rep("B1int", nrow(regression.analysis.combined.retractionfibers))

write.csv(x=regression.analysis.combined.retractionfibers, "Regression.analysis.results.retractionfiber.b1int.csv")


