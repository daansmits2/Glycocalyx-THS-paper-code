## load packages-----
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)
library(reshape2)
library(zoo)
library(purrr)
library(cowplot)
library(fpc)
library(ggalt)
library(dbscan)
library(factoextra)
## set working directory-----
workdir <- choose.dir(getwd(), "Choose folder containing subdirectories with .csv files")
parentdirname <- basename(workdir)


## set parameters ----

Timessigma = 3
pool_bypixels = TRUE
grouped_pixels = 12
grouped_pixels = 4

pool_relative = FALSE
grouped_fraction = 0.1

select_linethickness=T

grouped_frames = 1

## Load functions-------
read_and_clean <- function(file, file_id) {
  df <- read.csv(file)
  df$FILE_ID <- file_id
  return(df)
}


# Define your read_and_clean function
read_and_clean <- function(file, file_id) {
  # Read the CSV file into a data frame
  df <- read.csv(file)
  df$FILE_ID <- file_id
  df$BASENAME <- basename(file)
  print(unique(df$BASENAME))
  
  return(df)
}


calc.threshold <- function(data, threshold_var, min, max, stepsize, log=FALSE){
  breaks_log10 <- function(x) {
    library(scales)
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    
    10^(seq.int(low, high))
  }
  
  threshold_var <- enquo(threshold_var)
  
  # Initialize a data frame to store the results
  results <- data.frame(ET_ratio = numeric(0), 
                        Thresh = numeric(0), 
                        Event_Count = numeric(0), 
                        Event_Percentage = numeric(0),  # Add this column for percentage
                        TRACK_ID_UNIQUE = character(0))
  
  # Get unique ET_ratio values
  et_ratios <- unique(data$ET)
  
  # Loop through each ET_ratio value
  for (et_ratio in et_ratios) {
    # Filter the data for the current ET_ratio
    data_filtered <- data %>% filter(ET == et_ratio)
    
    # Calculate the total number of events for this ET_ratio
    total_events <- nrow(data_filtered)
    
    # Initialize list to store TRACK_ID_UNIQUE for each threshold
    filtered_track_ids <- list()
    
    # Loop through the threshold values from 3 to 6 in increments of 0.01
    for (threshold in seq(min, max, by = stepsize)) {
      # Filter the data based on the current threshold
      filtered_data <- data_filtered %>%
        group_by(TRACK_ID_UNIQUE) %>%
        filter(!!threshold_var > threshold)
      
      # Count the number of events remaining
      event_count <- nrow(filtered_data)
      
      # Calculate the percentage of total events
      event_percentage <- (event_count / total_events) * 100
      
      # Collect TRACK_ID_UNIQUE values
      track_ids <- pull(filtered_data, TRACK_ID_UNIQUE)
      filtered_track_ids[[as.character(threshold)]] <- track_ids
      
      # Add the results to the results data frame
      results <- rbind(results, 
                       data.frame(ET_ratio = et_ratio, 
                                  Thresh = threshold, 
                                  Event_Count = event_count, 
                                  Event_Percentage = event_percentage,  # Add the percentage here
                                  TRACK_ID_UNIQUE = paste(track_ids, collapse = ", ")))
    }
  }
  
  results <- results %>%  ungroup() %>% group_by(Thresh) %>% 
    mutate(event_count_ratio = Event_Count / mean(Event_Count[ET_ratio=="0"], na.rm=TRUE),
           event_score = event_count_ratio * Event_Count)
  
  eventsperthr <- ggplot(results, aes(x=Thresh, y=Event_Count, color=ET_ratio)) +
    theme_classic() +
    geom_point() +
    geom_line() +
    ggtitle(threshold_var)+
    labs(x="Threshold", y= "Event #")+
    
    theme(legend.position = c(0.8,0.8))
  
  if(log==TRUE){
    eventsperthr <- eventsperthr +
      annotation_logticks(sides="l")+
      scale_y_log10(breaks = breaks_log10)
  }
  
  
  eventcountratio <- ggplot(results, aes(x=Thresh, y=event_count_ratio, color=ET_ratio)) +
    theme_classic() +
    geom_point() +
    geom_line() +
    ggtitle(threshold_var)+
    
    labs(x="Threshold", y= "Event Ratio")
  if(log==TRUE){
    eventcountratio <- eventcountratio +
      annotation_logticks(sides="l")+
      scale_y_log10(breaks = breaks_log10)
  }
  
  eventscore <- ggplot(results, aes(x=Thresh, y=event_score, color=ET_ratio)) +
    theme_classic() +
    geom_point() +
    geom_line() +
    ggtitle(threshold_var)+
    
    labs(x="Threshold", y= "Event Ratio * Event count")
  
  if(log==TRUE){
    eventscore <- eventscore +
      annotation_logticks(sides="l")+
      scale_y_log10(breaks = breaks_log10)
  }
  
  grid <- plot_grid(eventsperthr, eventcountratio, eventscore, nrow=1)
  print(grid)
  
  df <- results  %>% filter(ET_ratio!="0")  %>%  filter(!is.infinite(event_count_ratio)) %>%  filter(!is.nan(event_count_ratio))
  chosen_thr <- max(df$Thresh)
  print(chosen_thr)
  ratio <- max(df$event_count_ratio[df$Thresh==chosen_thr])
  print(ratio)
  
  
  return(list(results, eventsperthr, eventcountratio, eventscore))
}

# Define the function to fit an exponential model
fit_exponential <- function(df, x_col, y_col) {
  tryCatch({
    # Fit the model y = a * exp(b * x) + c using dynamically specified columns
    fit <- nls(as.formula(paste(y_col, "~ a * exp(b *", x_col, ") + c")), 
               data = df, 
               start = list(a = 1, b = 0.5, c = 1))
    
    # Extract coefficients
    coefs <- coef(fit)
    data.frame(a = coefs["a"], b = coefs["b"], c = coefs["c"])
  }, error = function(e) {
    # Return NA in case of error
    data.frame(a = NA, b = NA, c = NA)
  })
}



# Define a custom theme with black text
my_theme <- function() {
  theme_classic() +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black")
    )
}
# Set the custom theme as the default theme for ggplot
theme_set(my_theme())
## set working directory -----

setwd(workdir)
## Import data ----

# Define file name
output_file <- "Import_data.csv"
save=F

# Check if the file already exists in the working directory
if (file.exists(output_file)) {
  # If the file exists, load the data
  df <- read.csv(output_file)
  message("File 'Import_data.csv' already exists. Data loaded successfully.")
} else {
  
  roimanager_files <- list.files(workdir, pattern=paste0("_Profiles.csv"), full.names = TRUE)
  data_measurement <- bind_rows(lapply(seq_along(roimanager_files), function(i) read_and_clean(roimanager_files[i], i)))
  
  # If the file does not exist, process the data
  df <- data_measurement %>%
    separate(Segment, into = c("Bleached", "Roi_region", "Segment"), sep = "_") %>%
    dplyr::select(-X.1) %>%
    separate(BASENAME, into = c("Total_region", "A", "B", "C", "D", "E", "Channel"), sep = "_", remove = FALSE) %>%
    dplyr::select(-c(FILE_ID, B, C, D, E)) %>%
    mutate(Channel = ifelse(Channel == "Profiles.csv", "Gcx", Channel)) %>%
    mutate(
      ImageNR = gsub("Transitionzone_", "", Image),
      ImageNR = gsub("Granule_", "", ImageNR),
      ImageNR = gsub("Single-cell-segment_", "", ImageNR),
      ImageNR = gsub("Cell-cell-contact_", "", ImageNR),
      ImageNR = gsub(".czi", "", ImageNR),
      ImageNR = gsub("C1-", "", ImageNR),
      ImageNR = gsub("C2-", "", ImageNR),
      Total_region = ifelse(Total_region == "Single-cell-membrane", "Single-membrane-segment", Total_region),
      Total_region = ifelse(Total_region == "Cell", "Cell-cell-contact", Total_region),
      Total_region = as.factor(Total_region),
      Total_region = factor(Total_region, levels = c("Granule", "Single-membrane-segment", "Cell-cell-contact", "Transitionzone")),
      label = A
    ) %>%
    group_by(A) %>%
    mutate(
      ImageNR = as.numeric(ImageNR),
      ExpNR = ifelse(ImageNR >= 100, 2, 1),
      ExpNR = ifelse(ImageNR >= 200, 3, ExpNR),
      ExpNR = ifelse(ImageNR >= 300, 4, ExpNR),
      ImageNR = cur_group_id()
    )
  
  df <- df %>%  filter(ExpNR!=3)
  
  df <- df %>%  mutate(Total_region = Segment)
  
  
  
  
  df <- df %>%  group_by(ImageNR, CellNR, Total_region, Channel, Linethickness, Bleached, label, A, Image) %>%  mutate(RoiNR = cur_group_id(), RoiNR = as.factor(RoiNR))
  
  if(save==T){
  # Save the processed data to a CSV file
  write.csv(df, output_file, row.names = FALSE)
  message("File 'Import_data.csv' has been created and data saved successfully.")
  }
}

## Organize data------
data <- df

data <- data %>% ungroup() %>% group_by(ImageNR, CellNR, RoiNR, Frame, Image) %>% mutate(Bleached = ifelse(Frame<=10, "Unbleached", Bleached),
                                                                                  Bleached = as.factor(Bleached),
                                                                                  Bleached = factor(Bleached, levels=c("Unbleached", "Bleached"))
                                                                                  )



bleachingframe <- max(data$Frame[data$Bleached=="Unbleached"], na.rm=T)

data <- data %>%  mutate(X_micron = X* Pixelsize_micron, 
                         Frame = Frame - bleachingframe,
                         timeinterval = Time, 
                         time = Frame * timeinterval - timeinterval) %>% 
  dplyr::select(-Time)

data <- data %>%  arrange(time) %>%  group_by(RoiNR) %>% mutate(time.interval = ifelse(Bleached=="Bleached", time[Frame==1] - time[Frame==0], NA))

data <- data %>% group_by(RoiNR) %>%  mutate(Intensity_bgcor = Intensity - Background)



data <- data %>%  group_by(RoiNR, Total_region) %>%  mutate(I_minmax_cor = Intensity_bgcor - min(Intensity_bgcor),
                                              I_minmax_cor = I_minmax_cor / max(I_minmax_cor)
                                              )

data <- data %>% ungroup() %>%  group_by(RoiNR, Total_region) %>% mutate(radius = (max(X_micron[Bleached=="Bleached"]) - min (X_micron[Bleached=="Bleached"]))/2)


data <- data %>%  mutate(Symmetry = ifelse(Total_region !="Transitionzone", "Symmetrical", "Asymmetrical"))

data <- data %>% mutate(Segment = ifelse(Segment=="Cell-cell-contact", "Cell-cell contact", Segment))

## Calculate average signal and stdev per line thickness ------
thicknesses <- data %>%  group_by(ImageNR, CellNR, RoiNR, Bleached, Linethickness, Channel) %>% summarise(Mean = mean(Intensity),
                                                                                                 Stdev = sd(Intensity),
                                                                                                 ratio = Mean/Stdev
                                                                                                 ) %>% 
  filter(Bleached=="Unbleached") %>%  filter(Channel=="Gcx")

ggplot(thicknesses, aes(x=Linethickness, y=ratio))+
  geom_violin(scale="width", aes(group=as.factor(Linethickness)))+
  #geom_jitter(width=0.5)+
  geom_boxplot(outlier.shape=NA, width=0.5, aes(group=as.factor(Linethickness)))+
  
  #geom_smooth(color="black", span=1)+
  scale_color_viridis(option="turbo")+
  scale_x_continuous(breaks=1:10)


boxplot_stats <- thicknesses %>%
  ungroup() %>%
  group_by(Linethickness) %>%
  summarise(
    n = n(),
    mean = mean(ratio, na.rm = TRUE),
    sd = sd(ratio, na.rm = TRUE),
    median = median(ratio, na.rm = TRUE),
    Q1 = quantile(ratio, 0.25, na.rm = TRUE),
    Q3 = quantile(ratio, 0.75, na.rm = TRUE),
    IQR = IQR(ratio, na.rm = TRUE),
    lower_whisker = max(min(ratio, na.rm = TRUE), Q1 - 1.5 * IQR),
    upper_whisker = min(max(ratio, na.rm = TRUE), Q3 + 1.5 * IQR),
    ci_lower = mean - qt(0.975, df = n - 1) * sd / sqrt(n),
    ci_upper = mean + qt(0.975, df = n - 1) * sd / sqrt(n)
  ) %>%
  ungroup()

print(boxplot_stats)

ggplot(thicknesses, aes(x=Mean, y=Stdev, color=Linethickness, group=as.factor(Linethickness)))+
  geom_point()+
  geom_line()+
  scale_color_viridis(option="turbo")


smry_thickness <- thicknesses %>%  group_by(Linethickness) %>% summarise(ratio  = median(ratio)) 

Linethickness_pixel <- filter(smry_thickness, Linethickness == which.max(ratio))$Linethickness


## Filter on specific line thickness -------
if(select_linethickness==T){
data <- data %>%  dplyr::filter(Linethickness==Linethickness_pixel)
}
## Display graph of example ROI per frame -----

display_example  = FALSE

if(display_example==TRUE){
roi <- droplevels(subset(data, RoiNR==1))

lines <- ggplot(roi, aes(x=X_micron, y= Intensity_bgcor))+
  geom_line(aes(group=Frame, color=Segment), size=0.5)+
  facet_wrap(.~Frame)

print(lines)
}

## Determine Bleached region computationally ---- if asked for ------

Bleached_computational = F
smooth_span=0.003

data <- data %>% mutate(RoiNR, CellNR, ImageNR, Frame, Total_region) %>% 
  mutate(Intensity_bgcor_smooth = loess(Intensity_bgcor ~ seq_along(Intensity_bgcor), span=smooth_span, family="gaussian")$fitted)

if(Bleached_computational == TRUE){
  data <- data %>% mutate(Average_prebleach_region = mean(Intensity_bgcor_smooth[Bleached=="Unbleached"]),
                          Intensity_vsunbleached = Intensity_bgcor_smooth/ Average_prebleach_region
                          )
}



## Bin Pixels ----




segment_size = grouped_pixels




if(pool_bypixels==TRUE){
data <- data %>%
  group_by(RoiNR, CellNR, ImageNR, Bleached, Frame, Total_region, Channel, Image) %>% 
  mutate(X_group = floor(X/grouped_pixels)*grouped_pixels,
         X_group_range = paste0(X_group, "/", X_group + grouped_pixels - 1)
         ) %>% 
  group_by(RoiNR, CellNR, ImageNR, X_group, Bleached, Frame, Total_region, Channel, Image) %>% 
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}")
  ) %>%
  ungroup()


# Calculate the midpoint of the X values
mid_X <- (min(data$X) + max(data$X)) / 2

# Assign segments using dplyr
assign_sides = FALSE

  if(assign_sides==TRUE){
    data <- data %>%
      group_by(RoiNR, CellNR, ImageNR, Bleached, Frame, Total_region) %>%
      mutate(
        X_min = min(X),    # Find the minimum X
        X_max = max(X),    # Find the maximum X
        mid_X = (X_min + X_max) / 2,  # Calculate the midpoint (center of X)
        
        # Calculate the distance from both sides
        Distance_from_left = (X - X_min) %/% segment_size,    # Segment number from the left
        Distance_from_right = (X_max - X) %/% segment_size,   # Segment number from the right
        
        # Bin start and end positions for the current segment
        Bin_start = pmin(Distance_from_left, Distance_from_right) * segment_size + 1,
        Bin_end = Bin_start + segment_size - 1,
        
        # Compute the average X of the bin (the center of the bin)
        Bin_middle = (Bin_start + Bin_end) / 2,
        
        # Adjust X_group to mirror around the midpoint
        X_group = case_when(
          Bin_middle <= mid_X ~ Bin_middle,                    # For left side, use the bin middle directly
          Bin_middle > mid_X ~ mid_X - (Bin_middle - mid_X)    # For right side, mirror the bin middle around the center
        ),
        
        # Assign side labels (left or right)
        Side = case_when(
          Distance_from_left < Distance_from_right ~ "left",   # If closer to the left, assign "left"
          TRUE ~ "right"                                        # Otherwise, assign "right"
        ),
        Side = case_when(Distance_from_left==Distance_from_right ~ NA,
                         TRUE ~ Side
                         )
      ) %>%
      ungroup()
    
    data <- data %>% mutate(X_group = X_group - 0.5)
    
  }else{
    data <- data %>%
      group_by(RoiNR, CellNR, ImageNR, Bleached, Frame, Total_region) %>%
      mutate(X_group = floor(X / segment_size),
             Side = "no_side")  # Divide by bin size (10) and floor to get group
  }

}




if(pool_relative==TRUE){
  
  data <- data %>%
    group_by(RoiNR, CellNR, ImageNR, Segment, Bleached, Frame, Total_region, Channel) %>% 
    mutate(Bleached_segment_count = ifelse(Bleached=="Bleached", n(), NA),
           Relative_position = row_number()/ Bleached_segment_count,
           Subsegment = ifelse(Relative_position <=0.5, "Start-Mid", "Mid-End"),
           Subsegment = as.factor(Subsegment),
           Subsegment = factor(Subsegment, levels = c("Start-Mid", "Mid-End"))
    )
  
  
  data <- data %>%
    group_by(RoiNR, CellNR, ImageNR, Bleached, Frame, Total_region, Channel) %>%
    mutate(
      Relative_group = cut(Relative_position, breaks = seq(0, 1, by = grouped_fraction), include.lowest = TRUE, labels = FALSE),
    ) %>%
    group_by(RoiNR, CellNR, ImageNR, Relative_group, Subsegment, Bleached, Frame, Total_region) %>%
    summarise(
      across(where(is.numeric), ~ mean(.), .names = "{.col}")
    ) %>%
    ungroup()
  
}


## Bin Frames  -----


# Parameters
pool_frames = T

# Get the data for group 1 before summarizing
group_1_data <- data %>%
  ungroup() %>%
  group_by(RoiNR, CellNR, ImageNR, Bleached, Total_region, Channel, Side, X_micron, Image) %>%
  mutate(Frame_bin = ceiling(Frame / grouped_frames) * grouped_frames) %>%
  group_by(RoiNR, CellNR, ImageNR, Bleached, Total_region, Channel, Frame_bin, Side, X_micron) %>% 
  filter(RoiNR == 1) %>%  # Filter for group 1 (adjust condition as needed)
  mutate(rows = row_number()) %>% 
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}")
  ) %>%
  ungroup() %>%
  mutate(Frame = Frame_bin)

# Proceed with summarizing
file_averaged <- data %>%
  ungroup() %>%
  group_by(RoiNR, CellNR, ImageNR, Bleached, Total_region, Channel, Side, X_micron, Image) %>%
  mutate(Frame_bin = ceiling(Frame / grouped_frames)) %>%
  group_by(RoiNR, CellNR, ImageNR, Bleached, Total_region, Channel, Frame_bin, Side, X_micron, Image) %>%
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}")
  ) %>%
  ungroup() %>%
  mutate(Frame = Frame_bin)

if (pool_frames == TRUE) {
  file_unpooled <- data
  data <- file_averaged
}



## smoothen pixel trajectories if chosen ----

smoothen = F
smooth_span = 1
if(smoothen==T){
  data <- data %>% 
    group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, Frame, Channel) %>%
    mutate(Intensity_bgcor = loess(Intensity_bgcor ~ seq_along(Intensity_bgcor), span=smooth_span, family="gaussian")$fitted
    )
  
}



## normalize for max and min----
data <- data %>% ungroup() %>% 
  group_by(RoiNR, X_group, Side, Total_region, Channel, Image) %>%
  mutate(I_minmax_cor = Intensity_bgcor - min(Intensity_bgcor)) %>%
  mutate(I_minmax_cor = I_minmax_cor / max(I_minmax_cor)) %>% 
  ungroup() %>%
  group_by(RoiNR, Side, Total_region) %>%
  mutate(I_minmax_cor_allregions = Intensity_bgcor - min(Intensity_bgcor)) %>%
  mutate(I_minmax_cor_allregions = I_minmax_cor_allregions / max(I_minmax_cor_allregions))




## plot a profile per X_group ----


plot_xgroup= F

if(plot_xgroup==T){
s <- subset(data)

sum <- s %>%  group_by(X_group, Side) %>% summarise(
  across(where(is.numeric), ~ mean(.), .names = "{.col}")
) 

ggplot(s,aes(x=Frame, y=I_minmax_cor))+
  geom_line(aes(group=interaction(RoiNR, Side),color=Side), size=1)+
  geom_point()+
  facet_wrap(Bleached+Side~X_group)
}

## fit exponential curve ----
DEopt = F

if(DEopt==T){
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(broom)
# Fit the data



library(DEoptim)

# Define the cost function
cost_function <- function(params, data) {
  yf <- params[1]
  y0 <- params[2]
  log_alpha <- params[3]
  
  # Calculate predicted values
  predicted <- yf + (y0 - yf) * exp(log_alpha * data$time)
  
  # Sum of squared residuals
  sum((data$I_minmax_cor - predicted)^2)
}


# Use DEoptim for optimization
safe_DEoptim <- possibly(
  ~ {min_intensity <- min(data$I_minmax_cor)
  max_intensity <- max(data$I_minmax_cor)
  
  result <- DEoptim(
    cost_function, 
    lower = c(min_intensity, min_intensity, -20),  # Widen bounds
    upper = c(max_intensity, max_intensity, 20),
    data = .
  )
  
  list(params = result$optim$bestmem)
  },
  otherwise = NA
)

safe_DEoptim <- possibly(
  ~ {
    min_intensity <- min(data$I_minmax_cor)
    max_intensity <- max(data$I_minmax_cor)
    
    result <- DEoptim(
      cost_function,
      lower = c(0, 0, -20),
      upper = c(1, 1, 20),
      data = .,
      DEoptim.control(
        strategy = 2,
        itermax = 50,          # Increased iterations
        NP = 50,                # Larger population size
        F = 0.9,                 # Differential weight
        CR = 0.9                 # Crossover probability
      )
    )
    list(params = result$optim$bestmem)
  },
  otherwise = NA
)


  fitted <- data %>%
    group_by(RoiNR, X_group, Bleached, Total_region, Side, Channel) %>%
    nest() %>% 
    mutate(
      fit = map(data, safe_DEoptim),
      yf = map_dbl(fit, ~.x$params[1]),
      y0 = map_dbl(fit, ~.x$params[2]),
      log_alpha = map_dbl(fit, ~.x$params[3]),
      t_half = log(0.5) / log_alpha,
      t_half = round(t_half, 3),
      alpha = exp(log_alpha)
    )
  
  
  # Merge fitted parameters back into the original data
  data_with_fit <- data %>%
    left_join(
      fitted %>% select(RoiNR, X_group, yf, y0, log_alpha, t_half),
      by = c("RoiNR", "X_group", "Bleached", "Total_region", "Side", "Channel")
    ) %>%
    # Calculate predicted values based on the fitted parameters
    mutate(
      Intensity_bgcor_expfit = yf + (y0 - yf) * exp(log_alpha * Frame),
      Intensity_delta = yf - y0,
      Intensity_fit_minmaxcorrected = Intensity_bgcor_expfit - y0,
    ) 
  
  
  
  # calculate diffusion constant
  
  data_with_fit <- data_with_fit %>%  group_by(ImageNR, CellNR, X_group, Bleached, RoiNR, ExpNR, Total_region, Side, Channel, Image) %>% 
    mutate(Diffusion_coefficient = radius^2 / (4*t_half))

  
}



## nls------



# Define the model formula
exponential_model <- I_minmax_cor ~ yf + (y0 - yf) * exp(log_alpha * time)
# Safe fitting with nls
safe_nls <- possibly(
  ~ {
    # Provide reasonable starting values based on the current group's data
    start_vals <- list(
      yf = 1,
      y0 = 0,
      log_alpha = -0.1
    )
    
    
    # Calculate weights: increase emphasis on early points by using a higher exponent (k)
    k <- 0.1 # Adjust as needed for stronger emphasis
    weights <- 1 / (1 + .x$time^k)
    
    # Fit the exponential model using nls with weights
    fit <- nls(
      formula = I_minmax_cor ~ yf + (y0 - yf) * exp(log_alpha * time),
      data = .x,
      start = start_vals,
      weights = weights,  # Add weights
      control = nls.control(maxiter = 100, minFactor = 1e-8)
    )
    
    # Extract coefficients as a named list
    as.list(coef(fit))
  },
  otherwise = list(yf = NA_real_, y0 = NA_real_, log_alpha = NA_real_)
)


# Perform group-wise fitting and calculations
if (pool_bypixels == TRUE) {
  fitted <- data %>%
    group_by(RoiNR, X_group, Bleached, Total_region, Side, Channel) %>%
    nest() %>% 
    mutate(
      # Apply safe_nls to each group
      fit = map(data, safe_nls),
      yf = map_dbl(fit, "yf"),
      y0 = map_dbl(fit, "y0"),
      log_alpha = map_dbl(fit, "log_alpha"),
      t_half = ifelse(!is.na(log_alpha), round(log(0.5) / log_alpha, 3), NA_real_),
      t_half = ifelse(
        !is.na(log_alpha) & (y0 != yf),
        log((0.5 * y0 - yf) / (y0 - yf)) / log_alpha,
        NA_real_
      ),
      alpha = ifelse(!is.na(log_alpha), exp(log_alpha), NA_real_)
    )
  
  # Merge fitted parameters back into the original data
  data_with_fit <- data %>%
    left_join(
      fitted %>% select(RoiNR, X_group, yf, y0, log_alpha, t_half),
      by = c("RoiNR", "X_group", "Bleached", "Total_region", "Side", "Channel")
    ) %>%
    # Calculate predicted values based on the fitted parameters
    mutate(
      Intensity_bgcor_expfit = ifelse(
        !is.na(yf) & !is.na(y0) & !is.na(log_alpha),
        yf + (y0 - yf) * exp(log_alpha * Frame),
        NA_real_
      ),
      Intensity_delta = ifelse(!is.na(yf) & !is.na(y0), yf - y0, NA_real_),
      Intensity_fit_minmaxcorrected = ifelse(
        !is.na(Intensity_bgcor_expfit) & !is.na(y0),
        Intensity_bgcor_expfit - y0,
        NA_real_
      )
    )
  
}


data_with_fit <- data_with_fit %>%  mutate(yf_half = yf*0.5)

# Calculate t_half by finding the closest value to yf * 0.5
data_with_fit <- data_with_fit %>%
  group_by(RoiNR, X_group, Bleached, Total_region, Side, Channel) %>%
  arrange(Frame) %>% 
  mutate(
    I0 = first(Intensity_bgcor_expfit),
    yf_half = (yf +I0)*0.5,
    # Calculate the absolute difference between fit values and yf_half
    diff_to_yf_half = abs(Intensity_bgcor_expfit - yf_half),
    min_diff_yf_half = min(diff_to_yf_half, na.rm=T),
    t_half_coord = ifelse(diff_to_yf_half == min_diff_yf_half & Bleached=="Bleached", T,F),
    t_half = mean(time[t_half_coord==T], na.rm=T)
  )

  
## filter out negative exponential fits and too shallow fits------
data_with_fit <- data_with_fit %>%  mutate(X_group_micron = X_group * Pixelsize_micron)

filter=T

if(filter==TRUE){
data_with_fit <- data_with_fit %>% 
  group_by(ImageNR, CellNR, X_group, Bleached, RoiNR, ExpNR, Total_region, Side, Channel, Image) %>% 
  mutate(
    Intensity_bgcor_expfit = case_when(
      Intensity_bgcor_expfit >= first(Intensity_bgcor_expfit) ~ Intensity_bgcor_expfit,
      TRUE ~ NA_real_
    ),
    Intensity_bgcor_expfit = ifelse(yf < I0, NA, Intensity_bgcor_expfit),
    Intensity_bgcor_expfit = ifelse(yf > 1 , NA, Intensity_bgcor_expfit)
  )
}

calc.groups <- function(df){
groups <- df %>% 
  group_by(ImageNR, CellNR, X_group, Bleached, RoiNR, ExpNR, Total_region, Side) %>%  tally()

print(groups)
}

## redefine bleached area-----

recoverythr = 0.8
recoverybelowthr = 0.4

  data_with_fit <- data_with_fit %>%  ungroup() %>% 
    group_by(ImageNR, CellNR, X_group, RoiNR, ExpNR, Total_region, Side, Channel, Image) %>% 
    mutate(
      First_Bleached_Frame = min(Frame[Bleached == "Bleached"], na.rm = TRUE)+1,
      First_Bleached_Frame_fit = mean(Intensity_bgcor_expfit[Frame<=First_Bleached_Frame & Bleached=="Bleached"],na.rm=T),
      Bleaching_delta = mean(I_minmax_cor[Bleached=="Unbleached"] & Frame <0)-First_Bleached_Frame_fit,
      
      Unbleach_mean = mean(I_minmax_cor[Bleached == "Unbleached" & Frame < 0], na.rm = TRUE),
      Sd = sd(I_minmax_cor[Bleached == "Unbleached" & Frame < First_Bleached_Frame], na.rm = TRUE),
      Bleaching_threshold = Unbleach_mean - Timessigma*Sd,
      Include = ifelse(First_Bleached_Frame_fit < Bleaching_threshold, TRUE, FALSE),
      Include = ifelse(Bleached=="Unbleached", F, Include),
      Include = ifelse(yf >1, F, Include),
      Include = ifelse(yf < I0, F, Include),
      Include = ifelse(yf > (Unbleach_mean * recoverythr), F, Include),
      Include = ifelse(yf < (Unbleach_mean * recoverybelowthr), F, Include),
      Include = ifelse(is.na(Include), F, Include),
      t_half = ifelse(Include ==F, NA, t_half),
      y0 = ifelse(Include==F, NA, y0),
      yf = ifelse(Include==F, NA, yf),
      yf_half = ifelse(Include==F, NA, yf_half),
      I0 = ifelse(Include==F, NA, I0)
    
    ) 
  
  
  
  # Calculate diffusion constant
  data_with_fit <- data_with_fit %>% ungroup () %>% 
    group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, Channel, Image) %>% 
    mutate(Bleaching_midpoint = (min(X_micron[Include==T], na.rm=T) + max(X_micron[Include==T], na.rm=T))/2
           ) %>% 
  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, Channel, X_group, Image) %>% 
    mutate(
           X_micron_bleaching_midpoint = X_micron - Bleaching_midpoint,
           radius = abs(X_micron_bleaching_midpoint)+0.001,
           Diffusion_coefficient = ifelse(!is.na(t_half), (radius^2) / (4 * t_half), NA_real_)
    )
  
  data_with_fit %>%  select(Unbleach_mean, Sd, Bleaching_threshold, Include) %>% print()
  
  sub <- data_with_fit
  
  sd <- data_with_fit %>%  
    group_by(ImageNR, CellNR, X_group, RoiNR, ExpNR, Total_region, Side, Channel, Image) %>% 
    summarise(stdev = sd(I_minmax_cor),
    )
  
  printallplots = F
  
  if(printallplots==TRUE){
    uniques <- unique(data_with_fit$RoiNR)
    for(u in 1:length(uniques)){
      
      set <- subset(data_with_fit,RoiNR== uniques[u])
      image <- unique(set$ImageNR)
      
      plot <- ggplot(data = subset(set, Total_region=="Transitionzone"), aes(x=Frame, y= I_minmax_cor))+
        geom_line(data=set, aes(group=RoiNR))+
        geom_line(data=subset(set, Frame>10 & Include==T), aes(y=Intensity_bgcor_expfit), color="blue")+
        facet_wrap(X_group~Include~RoiNR+Side)+
        scale_color_viridis_c(option="turbo")+
        ggtitle(paste0("Roi ", u, "_ImageNR ", image))
      print(plot)
    }
  }
  
  
  printallplots = F
  
  if(printallplots==TRUE){
    uniques <- unique(data_with_fit$RoiNR[data_with_fit$Total_region=="Granule"])
    for(u in 1:length(uniques)){
      
      set <- subset(data_with_fit,RoiNR== uniques[u])
      image <- unique(set$ImageNR)
      
      plot <- ggplot(data = subset(set, Total_region=="Transitionzone"), aes(x=Frame, y= I_minmax_cor))+
        geom_line(data=set, aes(group=RoiNR))+
        geom_line(data=subset(set, Frame>10 & Include==T), aes(y=Intensity_bgcor_expfit), color="blue")+
        facet_wrap(X_group~Include~RoiNR+Side)+
        scale_color_viridis_c(option="turbo")+
        ggtitle(paste0("Roi ", u, "_ImageNR ", image))
      print(plot)
    }
  }
  
  
## plot frap curves -----
  
  datatoplot <- data_with_fit %>%  filter(Channel=="Gcx")
  
  trans <- unique(datatoplot$RoiNR[datatoplot$Total_region=="Transitionzone"])
  print(length(trans))
  contacts <- unique(datatoplot$RoiNR[datatoplot$Total_region=="Cell-cell-contact"])
  print(length(contacts))
  single <- unique(datatoplot$RoiNR[datatoplot$Total_region=="Single-cell-segment"])
  print(length(single))
  
  
  nr = 709
  #datatoplot <- subset(data_with_fit, RoiNR == nr)
  datatoplot <- subset(datatoplot, Image == "Transitionzone_6.czi" & CellNR==1)
  datatoplot$Bleached <- droplevels(datatoplot$Bleached)
  datatoplot$time <- as.numeric(as.character(datatoplot$time))
  datatoplot$I_minmax_cor <- as.numeric(as.character(datatoplot$I_minmax_cor))
  datatoplot$X_group <- as.factor(datatoplot$X_group)
  
  
  plot= T
  
  if(plot==T){
    sub <- subset(datatoplot,X_group==32)
    
  individual_frapcurves <- ggplot(sub, aes(x=time, y= I_minmax_cor))+
    #geom_point()+
    geom_line(aes(group=interaction(X_group)))+
    facet_wrap(.~X_group)+
    geom_line(data=sub, aes(y=Intensity_bgcor_expfit, group=interaction(X_group)), size=1)+
    geom_vline(data=sub, aes(xintercept=t_half))+
    geom_hline(aes(yintercept=yf_half))+
    geom_hline(aes(yintercept=I0))+
    geom_hline(aes(yintercept=yf))+
    ggtitle(unique(datatoplot$Image))+
    ylim(0,1)
  print(individual_frapcurves)
  }
  
  
  
  
## define transition zone (v2) ----

print(unique(data_with_fit$ImageNR))

print_groupswithouttransition <- function(file){
  # Filter groups where Section is not "Transition"
  wine <- file %>% filter(Total_region=="Transitionzone") %>% 
    group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, Channel) %>% 
    mutate(Transitionzone_present = ifelse(any(Section=="Transition"), T, F)
    )
  
  # Print the non-"Transition" groups
  
  print(unique(wine$Transitionzone_present))
  print(unique(wine$ImageNR[wine$Transitionzone_present==F & wine$Total_region=="Transitionzone"]))
  
  groups <- file %>% filter(Total_region=="Transitionzone") %>% ungroup()%>% 
    group_by(Total_region, ImageNR, Section) %>% 
    summarize(single_cell_mean= mean(single_cell_mean),
              cell_contact_mean = mean(cell_contact_mean)
    ) %>%  arrange(ImageNR) %>%  
    mutate(FCcontact = cell_contact_mean / single_cell_mean)
  
  print(groups)
  return(groups)
  
}


file <- data_with_fit


file <- file %>% ungroup() %>%   group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, Channel, Bleached ) %>% 
  mutate(Section = ifelse(X_micron > (max(X_micron)-5), "Cell-cell-contact", "Single-cell segment"),
         Value_relative = Intensity_bgcor/mean(Intensity_bgcor[Section=="Single-cell segment"]),
         Lineprofiletype = "Contact",
         X_micron_normalized = X_micron - min(X_micron[Section=="Cell-cell-contact"]),
         RoiNR = cur_group_id(),
         Value_bgcor = Intensity_bgcor
  ) %>% group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, Bleached, Frame) %>% 
  mutate(number = cur_group_id())

upperXthr = 0
lowerXthr = -4




## Define Segment=Transition based on the distance between the mean of the single-cell segment and cell-cell contact
file <- file %>% 
  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, number, Channel, Bleached) %>%
  mutate(
    single_cell_mean = mean(Value_relative[Section == "Single-cell segment" & Lineprofiletype=="Contact"&
                                             X_micron_normalized < lowerXthr]),
    single_cell_mfi_mean = mean(Value_bgcor[Section == "Single-cell segment" & Lineprofiletype=="Contact"&
                                              X_micron_normalized < lowerXthr]),
    cell_contact_mean = mean(Value_relative[Section == "Cell-cell-contact"& Lineprofiletype=="Contact"]),
    cell_contact_mfi_mean = mean(Value_bgcor[Section == "Cell-cell-contact"& Lineprofiletype=="Contact"]),
    single_cell_sd = sd(Value_relative[Section == "Single-cell segment"& Lineprofiletype=="Contact"&
                                         X_micron_normalized < lowerXthr]),
    cell_contact_sd = sd(Value_relative[Section== "Cell-cell-contact"& Lineprofiletype=="Contact"])
  ) %>% mutate(
    Section = case_when(
      # First condition for Transition
      X_micron_normalized > lowerXthr & 
        X_micron_normalized <= upperXthr & 
        Value_relative > (single_cell_mean + 0 * single_cell_sd) & 
        Value_relative <= (cell_contact_mean + 0 * cell_contact_sd) ~ "Transition",
      
      # Default case
      TRUE ~ Section
    )
  )




example <- print_groupswithouttransition(file)

file <- file %>% 
  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, number, Channel, Bleached) %>% 
  arrange(X_micron) %>% 
  mutate(
    delta_relative = Value_relative  - lag(Value_relative),
    NoTransition = ifelse(all(Section!="Transition"), T, F),
    maxdeltarelative= max(delta_relative, na.rm=T),
    Section = ifelse(NoTransition==T & delta_relative == maxdeltarelative, "Transition", Section),
    Section = ifelse(NoTransition==T & lead(Section=="Transition"), "Transition", Section),
    consec_tran = ifelse(Section == "Transition", 1, 0),
    consec_tran = ifelse(consec_tran ==1 & lead(consec_tran ==1), consec_tran  + lead(consec_tran), consec_tran),
    consec_tran = ifelse(lag(consec_tran)> consec_tran & consec_tran !=0, lag(consec_tran), consec_tran),
    
    # Identifying the start of a new "Transition" event
    group_start = ifelse(Section == "Transition" & lag(Section, default = "") != "Transition", 1, 0),
    # Apply cumulative sum to assign a cumulative group number to each event
    cumulative_nr = cumsum(group_start),
    group_id = ifelse(consec_tran>0 & cumulative_nr>0, cumulative_nr, NA),
    group_id = ifelse(is.na(group_id) & !is.na(lag(group_id)) , lag(group_id), group_id),
    group_id = ifelse(is.na(group_id) & !is.na(lead(group_id)), lead(group_id), group_id),
    
  ) %>%  
  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, number, Channel, Bleached, group_id) %>% 
  mutate(delta = Value_relative[which.max(X_micron_normalized)] /
           Value_relative[which.min(X_micron_normalized)]
         ) %>% 
  ungroup() %>% 
  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, number, Channel) %>% 
  mutate(maxdelta = ifelse(delta == max(delta, na.rm=T), T, F),
         midpoint = mean(X_micron_normalized[Section=="Transition"], na.rm=T),
         midX_transition = midpoint,
         Section = ifelse(Section=="Transition" & maxdelta == F, "Temp", Section),
         Section = ifelse(Section=="Temp" & all(Section!="Transition"), "Transition", Section)
         )



file <- file %>%  group_by(number) %>% 
  mutate(Section = ifelse(X_micron_normalized == max(X_micron_normalized, na.rm=T), "Cell-cell-contact", Section),
         Section = ifelse(X_micron_normalized == min(X_micron_normalized, na.rm=T), "Single-cell segment", Section)
  )

file <- file %>%  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Side, Channel) %>% 
  mutate(Transition_mid_original_X = mean(X[Section=="Transition" & Frame<=(10/grouped_frames)], na.rm=T),
         Transition_mid_original_X_micron = mean(X_micron[Section=="Transition" & Frame<=(10/grouped_frames)], na.rm=T),
         X_micron_transitionmid = X_micron - Transition_mid_original_X_micron,
         Section = ifelse(!is.na(Section) & Section == "Temp" & X_micron_transitionmid < 0, 
                          "Single-cell segment", Section),
         Section = ifelse(!is.na(Section) & Section== "Temp" & X_micron_transitionmid > 0, 
                          "Cell-cell-contact", Section),
         Section = ifelse(is.na(Section) & X_micron_transitionmid < 0, "Single-cell segment", Section),
         Section = ifelse(is.na(Section) & X_micron_transitionmid > 0, "Cell-cell-contact", Section)
  )

             

file <- file %>%  group_by(number) %>% 
  mutate(Section = ifelse(Section!="Transition" & Section!="Background" & maxdelta == F &
                            X_micron_normalized > midX_transition, "Cell-cell-contact", Section),
         Section = ifelse(Section!="Transition" & Section !="Background" & maxdelta == F &
                            X_micron_normalized < midX_transition, "Single-cell segment", Section)
  )


recalc_value_relative=TRUE
if(recalc_value_relative==TRUE){
  file <- file %>%  group_by(number) %>%  mutate(Value_relative = Value_bgcor / mean(Value_bgcor[Section=="Single-cell segment"], na.rm=T))
  file <- file %>%  group_by(number) %>%  mutate(Value_vsmean = Value_bgcor / mean(Value_bgcor))
}



# define transition zone length

file <- file %>%  group_by(number) %>%  mutate(Transitionzone_length = row_number(Section=="Transition") * grouped_pixels * Pixelsize_micron)

## plot example transition zones -----


spatialprofile <- ggplot(subset(file,Bleached=="Unbleached" & Channel=="Gcx" & RoiNR<15), aes(x=X_micron, y= Value_bgcor))+
        geom_path(aes(group=number, color=Section), alpha=0.3, size=2)+
  #geom_smooth()+
  geom_point(aes(color=Section))+
  facet_wrap(.~number)
print(spatialprofile)

framesmry <- file %>%  group_by(ImageNR, CellNR, Side, Bleached, RoiNR, ExpNR, Total_region, X_group, Include, Channel, Section, X_micron) %>%  
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}")
  ) %>% filter(Bleached=="Unbleached") %>%  filter(Total_region!="Granule")

profiles <- ggplot(framesmry, aes(x=X_micron_transitionmid , y= Value_relative))+
  geom_path(aes(group=number, color=Section), alpha=0.5, size=1)+
  geom_smooth(color="black", fill="black", aes(group=interaction(Channel, Total_region)))+
  facet_wrap(Total_region~Channel)+
  xlim(-7,7)
print(profiles)
  

## define further transition zone mid in unbleached sections ----









## if transition zone is NA, define based on fixed value ----

correction = F


file <- file %>%  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Channel) %>% 
  mutate(Automated_detection = ifelse(is.na(X_micron_transitionmid)& Total_region=="Transitionzone", F, T),
         Automated_detection = ifelse(Total_region!="Transitionzone", T, Automated_detection)
  )


if(correction==T){

file <- file %>%  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Channel) %>% 
  mutate(X_micron_transitionmid = ifelse(Automated_detection==F, X_micron_normalized + 2, X_micron_transitionmid)
  )

}

file <- file %>%  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Channel) %>% 
  mutate(Section = ifelse(Total_region=="Single membrane", "Single-cell segment", Section),
         Section = ifelse(Total_region=="Cell-cell contact", "Cell-cell contact", Section),
         Section = ifelse(Total_region=="Granule", "Granule", Section),
         Mid_transitionzone_coord = mean(X_micron[Section=="Transition"], na.rm=T),
          X_micron_transitionmid = ifelse(Total_region=="Transitionzone", X_micron - Mid_transitionzone_coord, X_micron_transitionmid)
         )




file <- file %>%  group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Channel) %>% 
  mutate(X_micron_transitionmid = ifelse(Total_region!="Transitionzone", X_micron - mean(X_micron[Bleached=="Bleached"], na.rm=T), X_micron_transitionmid)
  )

## split trajectories into two Sides ----

Transitionborder = 1.5

file <- file %>% ungroup() %>% group_by(ImageNR, CellNR, RoiNR, ExpNR, Total_region, Channel, Image) %>%  
  mutate(X_midpoint_abs = abs(X_micron_transitionmid),
         Side = ifelse(X_micron_transitionmid<=0 & Total_region=="Transitionzone", "Single membrane", "Cell-cell-contact"),
         Side = ifelse(Total_region=="Single-cell-segment","Single membrane", Side),
         Side = ifelse(Total_region=="Cell-cell-contact", "Cell-cell-contact", Side),
         Side = ifelse(Total_region=="Granule", "Granule", Side),
         Insidethreshold = ifelse(X_micron_transitionmid<(Transitionborder*-1) , F, T),
         Insidethreshold = ifelse(X_micron_transitionmid>Transitionborder, F, Insidethreshold),
         #Insidethreshold = ifelse(Total_region!="Transitionzone", F, Insidethreshold)
         ) 




 
## plot transition zone length vs value relative ----


summary <- file %>%  group_by(ImageNR, CellNR, Bleached, RoiNR, ExpNR, Channel,Section) %>%  
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}"),
    count = n()
  ) %>% 
  filter(Channel=="Gcx") %>% 
  filter(Section=="Transition") %>% 
  filter(!any(Value_relative<0))
  

transitionlength <- ggplot(summary, aes(x=Transitionzone_length, y= Value_relative))+
  geom_point()+
  geom_smooth(method="lm", formula = "y~x")
print(transitionlength)

## plots -----



HiDeltathr = 1.5
HiSpeedthr = 1

file <- file %>%  ungroup () %>% group_by(ImageNR, CellNR, Bleached, ExpNR, Total_region, X_group, Include, Image) %>% 
  mutate(Side = ifelse(Channel=="Membrane", Side[Channel=="Gcx"], Side),
         X_midpoint_abs = ifelse(Channel=="Membrane", X_midpoint_abs[Channel=="Gcx"], X_midpoint_abs),
         Side = as.factor(Side)
         )

# filter out granules 


## summarize without taking into account the transition zone

summary <- file %>%  group_by(ImageNR, Image, CellNR, Side, Bleached, RoiNR, ExpNR, Total_region, X_group, Include, Channel, Insidethreshold) %>%  
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}"),
    count = n()
  ) %>%   
  filter(Bleached=="Bleached") %>% 
  filter(Total_region!="Granule") %>% 
  filter(Include==T) %>% 
    ungroup() %>% 
  group_by(ImageNR, Image, CellNR, Bleached, RoiNR, ExpNR, Total_region, Include, Channel) %>% 
  mutate(segments_included = max(row_number()))
  

perc0.25 <- quantile(summary$t_half[summary$Total_region=="Transitionzone" & summary$Side!="Granule"], 0.025)
perc0.975 <- quantile(summary$t_half[summary$Total_region=="Transitionzone" & summary$Side!="Granule"], 0.975)

within.ci = F
if(within.ci==T){
summary <- summary %>% 
 mutate(
    within_conf_interval = ifelse(t_half > perc0.25  & t_half < perc0.975,  T,F)
  ) %>% 
   filter(within_conf_interval==T)
}

trim = T

if(trim==T){
summary <- summary %>% group_by(ImageNR, CellNR, Bleached, RoiNR, ExpNR, X_group, Total_region, number, Include, Channel) %>% 
  filter(X_micron_transitionmid > -4.5 & X_micron_transitionmid < 4.5 ) %>% 
  ungroup() %>% 
  group_by(ImageNR, CellNR, Total_region, Channel) %>% 
  mutate(count = max(row_number())) 
}


summary <- summary %>% 
  filter(!(Total_region=="Transitionzone" & Insidethreshold==F)) %>% 
  filter(!(Total_region=="Single-cell-segment" & Insidethreshold==F)) %>% 
  filter(!(Total_region=="Cell-cell-contact" & Insidethreshold==F))


#summary <- summary %>%  filter(Channel=="Gcx")
summary <- summary %>%  mutate(Total_region = factor(Total_region, levels=c("Single-cell-segment", "Transitionzone", "Cell-cell-contact")))


rolling_median <- function(formula, data, n_roll = 75, ...) {
  x <- data$x[order(data$x)]
  y <- data$y[order(data$x)]
  y <- zoo::rollmedian(y, n_roll, na.pad = TRUE)
  structure(list(x = x, y = y, f = approxfun(x, y)), class = "rollmed")
}

predict.rollmed <- function(mod, newdata, ...) {
  setNames(mod$f(newdata$x), newdata$x)
}

median_data <- summary %>%
  group_by(X_micron_transitionmid, Channel) %>%
  summarize(Median_t_half = median(t_half, na.rm = TRUE), .groups = 'drop')

trajectory_oneplot <- ggplot(subset(summary), aes(x=X_micron_transitionmid, y=t_half))+
  #geom_point()+
  annotate("rect",xmin = -Transitionborder, xmax= 0, ymin = -Inf, ymax=Inf, fill="green", alpha=0.2)+
  annotate("rect",xmin = 0, xmax= Transitionborder, ymin = -Inf, ymax=Inf, fill="red", alpha=0.2)+
  #geom_line(data=median_data, aes(x=X_micron_transitionmid, y=Median_t_half, group=X_micron_transitionmid))+
  geom_vline(xintercept=0, linetype="dashed")+
geom_line(aes(group=RoiNR, color=as.factor(ExpNR)))+
  #geom_hline(yintercept=HiDeltathr)+
  geom_smooth(formula = y ~ x, method = "rolling_median", se = FALSE, size=2, color="black") +
  #geom_smooth(span=0.5)+
  #geom_quantile(quantiles=c(0.5))+
  #scale_color_viridis(option="turbo")+
  scale_fill_viridis_d()+
  facet_wrap(.~Total_region)+
  theme(legend.position="none")
print(trajectory_oneplot)



# Create interaction term
summary$Group <- interaction(summary$Side, summary$Total_region, summary$Insidethreshold)

set <- summary



summary_all <- set %>% group_by(Total_region, Side, Channel) %>%   summarise(
  across(where(is.numeric), ~ mean(.), .names = "{.col}")
)


# Dunn's test for pairwise comparisons
library(FSA)  # For Dunn's test
kruskal.test(t_half ~ Total_region, data = summary_all)

dunn.results.perexp <- dunnTest(t_half ~ Group, data = summary_all, method = "sidak")

print(dunn.results.perexp$res)



perexp <- ggplot(subset(summary_all), aes(x=interaction(Total_region, Insidethreshold), y=t_half, fill=interaction(Side)))+
  #geom_violin(scale = "width", position= position_dodge(.9), trim=T)+
  
  geom_boxplot(width=0.4, outlier.shape=NA, position=position_dodge(.9))+
  geom_hline(yintercept=(median(summary$t_half[summary$Group=="Single membrane.Single-cell-segment.FALSE"])))+
  geom_point(position=position_dodge(.9))+
  theme(legend.position=c(0.8,0.7))+
  facet_wrap(.~Channel)

print(perexp)

#shapiro.test(summary$t_half[summary$Group=="Single membrane.Single-cell-segment.FALSE" & summary$ExpNR==4])
kruskal.test(t_half ~ Group, data = summary)
dunn.results <- dunnTest(t_half ~ Group, data = summary, method = "sidak")

print(dunn.results$res)
write.csv(dunn.results$res, "DunnTest.csv")


violins <- ggplot(subset(set, Total_region=="Transitionzone"), aes(x=interaction(Group), y=t_half, fill=Side))+
  geom_violin(scale = "width", position= position_dodge(.9), trim=T)+
  geom_boxplot(width=0.4, outlier.shape=NA, position=position_dodge(.9))+
  geom_hline(yintercept=(median(summary$t_half[summary$Group=="Single membrane.Single-cell-segment.FALSE"])))+
# geom_point(position=position_dodge(.9))+
  facet_wrap(.~Channel)
print(violins)


summary_percell <- set %>%  group_by(ImageNR, Image, CellNR, Side, Bleached, RoiNR, ExpNR, Total_region, Include, Channel, Insidethreshold, Group) %>%  
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}"),
    count = n()
  ) 


kruskal.test(t_half ~ Total_region, data = summary_percell)

dunn.results.perexp <- dunnTest(t_half ~ Total_region, data = summary_percell, method = "sidak")

print(dunn.results.perexp$res)



violins_percell <- ggplot(subset(summary_percell), aes(x=interaction(Total_region), y=t_half, fill=Total_region))+
  geom_violin(scale = "width", position= position_dodge(.9), trim=T)+
  geom_boxplot(width=0.4, outlier.shape=NA, position=position_dodge(.9))+
  #geom_point(position=position_dodge(.9))+
  facet_wrap(.~Channel)+
  theme(legend.position = "none")
print(violins_percell)

selection <- summary_percell %>% ungroup() %>% arrange(Total_region) %>% dplyr::select(Total_region, t_half)
print(selection,n=130)


boxplot_stats <- selection %>%
  ungroup() %>%
  group_by(Total_region) %>%
  summarise(
    n = n(),
    mean = mean(t_half, na.rm = TRUE),
    sd = sd(t_half, na.rm = TRUE),
    median = median(t_half, na.rm = TRUE),
    Q1 = quantile(t_half, 0.25, na.rm = TRUE),
    Q3 = quantile(t_half, 0.75, na.rm = TRUE),
    IQR = IQR(t_half, na.rm = TRUE),
    lower_whisker = max(min(t_half, na.rm = TRUE), Q1 - 1.5 * IQR),
    upper_whisker = min(max(t_half, na.rm = TRUE), Q3 + 1.5 * IQR),
    ci_lower = mean - qt(0.975, df = n - 1) * sd / sqrt(n),
    ci_upper = mean + qt(0.975, df = n - 1) * sd / sqrt(n)
  ) %>%
  ungroup()

print(boxplot_stats)



trajectory_oneplot <- ggplot(subset(summary), aes(x=X_micron_bleaching_midpoint, y=t_half))+
  geom_point()+
  #geom_line(aes(group=interaction(RoiNR)))+
  geom_vline(xintercept=0, linetype="dashed")+
  #geom_hline(yintercept=HiDeltathr)+
  geom_smooth(se=T, aes(group=interaction(Channel)), span=1, color="black", fill="black")+
  labs(x="Distance from bleaching midpoint(um)", y="Half-max recovery time (s)")+
  scale_color_viridis(option="turbo")
print(trajectory_oneplot)


superviolins <- summary %>%  ungroup() %>% dplyr::select(Group, ExpNR, Diffusion_coefficient)
write.csv(superviolins,"superviolins.csv" )

ggplot(subset(summary, Bleached=="Bleached"  & !is.na(Side) & !is.na(Intensity_bgcor_expfit) & Include==T ), aes(x=Channel, y=Diffusion_coefficient, color=Value_relative))+
  geom_violin(scale="width", bw=0.3, draw_quantiles = T)+
  geom_jitter(width=0.2)



trajectory_oneplot <- ggplot(subset(summary, Bleached=="Bleached"  & !is.na(Side) & !is.na(Intensity_bgcor_expfit) & Include==T ), aes(x=X_micron_transitionmid, y=Value_relative, color=(Diffusion_coefficient)))+
  geom_point()+
  
  geom_line(aes(group=interaction(ImageNR, CellNR)))+
  geom_vline(xintercept=0, linetype="dashed")+
  facet_wrap(Channel~Total_region)+
  geom_smooth(se=T, aes(group=interaction(Channel)), span=1)+
  labs(x="Distance (um)", y="Max diffusion speed (um/sec)")+
  scale_color_viridis(option="turbo")

print(trajectory_oneplot)





reg <- subset(summary, Bleached=="Bleached"  & !is.na(Side) & !is.na(Intensity_bgcor_expfit) & Include==T & Channel=="Gcx" & Total_region!="Granule")

model <- lm(Value_relative ~ Diffusion_coefficient, data = reg)
summary_model <- summary(model)
p_value <- summary_model$coefficients["Diffusion_coefficient", "Pr(>|t|)"]
print(p_value)


selection <- subset(summary, Linethickness == 1 & Section=="Cell-cell-contact")

trajectory_oneplot <- ggplot(subset(selection, Bleached=="Bleached"  & !is.na(Side) & !is.na(Intensity_bgcor_expfit) & Include==T & Channel=="Gcx"), aes(x=Value_relative, y=Diffusion_coefficient))+
  geom_point(aes(color=Side))+
  geom_line(aes(group=interaction(RoiNR, Side), color=Side))+
  facet_wrap(Channel~Total_region)+
  geom_smooth(se=T, aes(group=Channel), color="black", span=2)


print(trajectory_oneplot)






## plot of interest -----



nr = 19
datatoplot <- subset(data_with_fit, RoiNR == nr)
datatoplot$Bleached <- droplevels(datatoplot$Bleached)
datatoplot$time <- as.numeric(as.character(datatoplot$time))
datatoplot$I_minmax_cor <- as.numeric(as.character(datatoplot$I_minmax_cor))
datatoplot$X_group <- as.factor(datatoplot$X_group)

plot= T

if(plot==T){
  individual_frapcurves <- ggplot(datatoplot, aes(x=time, y= I_minmax_cor, color=Include))+
    #geom_point()+
    geom_line(aes(group=interaction(X_group)))+
    facet_wrap(.~X_group)+
    geom_line(data=datatoplot, aes(y=Intensity_bgcor_expfit, group=interaction(X_group)), linetype=2, size=1)+
    geom_vline(data=datatoplot, aes(xintercept=t_half))+
    geom_hline(aes(yintercept=yf_half))+
    geom_hline(aes(yintercept=I0))+
    geom_hline(aes(yintercept=yf))+
    ggtitle(unique(datatoplot$Image))+
    ylim(0,1)
  print(individual_frapcurves)
}

trajectory_oneplot <- ggplot(subset(summary, Total_region=="Transitionzone" & RoiNR== nr),
                             aes(x=X_micron_transitionmid, y=t_half))+
  geom_point()+
  annotate("rect",xmin = -Transitionborder, xmax= Transitionborder, ymin = -Inf, ymax=Inf, fill="blue", alpha=0.2)+
  annotate("rect",xmin = -Inf, xmax= -Transitionborder, ymin = -Inf, ymax=Inf, fill="green", alpha=0.2)+
  annotate("rect",xmin = Transitionborder, xmax= Inf, ymin = -Inf, ymax=Inf, fill="red", alpha=0.2)+
  #geom_line(data=median_data, aes(x=X_micron_transitionmid, y=Median_t_half, group=X_micron_transitionmid))+
  geom_vline(xintercept=0, linetype="dashed")+
  # geom_line(aes(group=RoiNR))+
  #geom_hline(yintercept=HiDeltathr)+
 # geom_smooth(formula = y ~ x, method = "rolling_median", se = FALSE, size=2) +
  #geom_smooth(span=0.5)+
  #geom_quantile(quantiles=c(0.5))+
  scale_color_viridis(option="turbo")+
  facet_wrap(.~Total_region)
print(trajectory_oneplot)




int <- subset(set, Image==nr)

violins <- ggplot(subset(int), aes(x=interaction(Total_region), y=t_half, fill=Side))+
  geom_violin(scale = "width", position= position_dodge(.9), trim=T)+
  geom_boxplot(width=0.4, outlier.shape=NA, position=position_dodge(.9))+
  geom_hline(yintercept=(median(summary$t_half[summary$Group=="Single membrane.Single-cell-segment.FALSE"])))+
  geom_point(position=position_dodge(.9))+
  theme(legend.position=c(0.8,0.7))+
  facet_wrap(.~Channel)
print(violins)



## statistics------

summary_perexp <- summary %>%  group_by(ImageNR, CellNR, Bleached, RoiNR, ExpNR, Total_region, Side, Channel) %>% 
 filter(!is.na(t_half)) %>% 
 filter(X_midpoint_abs >2) %>% 
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}"),
    count = n()
  ) 

b <-subset(summary, Total_region!="Granule" &ImageNR ==3)
ggplot(subset(summary, Total_region!="Granule"), 
       aes(x=Value_relative, y=Diffusion_coefficient)) +
  geom_point(aes(color=Total_region))+
  geom_line(aes(group=RoiNR))+
  geom_smooth(span=2)+
  facet_wrap(RoiNR~Channel)

ggplot(subset(summary, X_midpoint_abs>2), aes(x=Side, y=Diffusion_coefficient))+
  geom_violin(scale="width")+
  geom_boxplot(width=0.5, outlier.shape = NA)+
  
  geom_point()+
  facet_wrap(.~Total_region)

summary_all<- summary_perexp %>%  group_by(Total_region, Side, Bleached, X_group, ExpNR, Channel) %>% 
  summarise(
    across(where(is.numeric), ~ mean(.), .names = "{.col}"),
    count = n()
  )




# compare individual conditions

paired_test_results <- summary_perexp  %>% filter(!is.na(Intensity_bgcor_expfit)) %>%  filter(!is.na(Bleached)) %>%  
  filter(!is.na(Side)) %>% 
  ungroup() %>% 
  group_by(X_group, Total_region, Channel, ImageNR) %>% 
  mutate(count = n()) %>% 
  group_by(X_group, Total_region,RoiNR) %>% 
  filter(Channel=="Gcx") %>% 
  summarise(
    p_value = wilcox.test(t_half[Side == "Left Single membrane"], t_half[Side == "Right cell-cell contact"], paired = TRUE)$p.value
  )

plot_data <- subset(summary_perexp, Bleached == "Bleached" & Total_region != "Granule") %>%
  left_join(paired_test_results, by = c("X_group", "Total_region")) %>%  
  filter(Total_region=="Transitionzone")



violin <- ggplot(subset(plot_data, Bleached=="Bleached" & Total_region!="Granule"), aes(x=X_group, y=t_half))+
  geom_violin()+
  geom_boxplot(outlier.shape=NA, width=0.4)+
  geom_line(aes(group=RoiNR))+
  geom_point()+
  facet_wrap(.~Symmetry + Total_region)

print(violin)
library(ggsignif)


# Subset paired_test_results for 'Transitionzone' data
p_values_transitionzone <- droplevels(subset(paired_test_results, Total_region == "Transitionzone"))
max_y <- subset(plot_data, Bleached=="Bleached" & Total_region=="Transitionzone")
max_y <- droplevels(max_y)
                                    
                                      

# Create the violin plot with annotation
violin <- ggplot(subset(plot_data, Bleached == "Bleached" & Total_region == "Transitionzone"), 
                 aes(x = Side, y = t_half)) +
  geom_line(aes(group = RoiNR)) +
  geom_point() +
  facet_wrap(.~X_group, nrow=1) +
  # Annotate the p-values on the plot
  geom_text(data = p_values_transitionzone, aes(label = paste("p =", round(p_value, 3)), x = 1.5, y = max(max_y$t_half) * 1.05), 
            inherit.aes = FALSE, color = "black", size = 5)+
  
  labs(x="Distance from unbleached region (micron)", y="Half-recovery time (sec)")

# Display the plot
violin