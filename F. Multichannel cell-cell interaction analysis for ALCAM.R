## install necessary packages------
library(ggplot2)
library(gridExtra)
library(tidyr)
library(matrixStats)
library(plyr)
library(ggpubr)
library(Rmisc)
library(svDialogs)
library(cowplot)
library(ggh4x)
library(ggpmisc)
library(nortest)
library(moments)
library(stats)
library(grid)
library(Hmisc)
library(zoo)
library(dplyr)
library(multcomp)


## chosen parameters ----


selected_optimallinethickness=6
library(tcltk)

# Function to create a pop-up dialog
get_pxsize <- function() {
  dialog <- tktoplevel()
  tkwm.title(dialog, "Input Pixel Size")
  
  input_var <- tclVar("") # Variable to store user input
  
  # Create entry box
  tkgrid(tklabel(dialog, text = "Enter the pixel size:"))
  entry <- tkentry(dialog, textvariable = input_var)
  tkgrid(entry)
  
  # Button to confirm
  on_ok <- function() {
    tclvalue(input_var) <<- tclvalue(input_var) # Update the variable
    tkdestroy(dialog) # Close the dialog
  }
  tkgrid(tkbutton(dialog, text = "OK", command = on_ok))
  
  # Wait for the user to close the dialog
  tkwait.window(dialog)
  
  # Return the numeric input
  return(as.numeric(tclvalue(input_var)))
}

# Call the function
pxsize <- get_pxsize()

Xthreshold = 3



upperXthr = 1
lowerXthr = -3


xlim = 15

contact_length_segregation_threshold = 1/2

## define theme for plotting------


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

contact_length_segregation_v6 <- function(df) {
  condition_subset <- df
  
  condition_subset$number <- as.numeric(condition_subset$number)
  
  max_values <- vector("numeric", length(unique(condition_subset$number)))
  
  for (i in 1:length(unique(condition_subset$number))) {
    print(paste0("number: ", i))
    
    condition_subset_nr <- subset(condition_subset, number == i)
    
    for (b in 1:length(unique(condition_subset_nr$ROINR))) {
      cond_name <- unique(condition_subset_nr$ROINR)
      cond_name <- cond_name[b]
      print(paste0("RoiNR: ", cond_name))
      
      profilelength <- max(condition_subset_nr$X[condition_subset_nr$ROINR == cond_name])
      print(paste("Profile length:", profilelength, "Pixels"))
      
      ##profilethreshold <- profilelength * contact_length_segregation_threshold
      ##print(paste("Single cell/cell-cell contact threshold: ", profilethreshold, "Pixels"))
      
  
      
      # Calculate X_relative for each subset (per ROINR)
      condition_subset$X_relative[condition_subset$number == i &
                                    condition_subset$ROINR == cond_name] <- condition_subset$X[condition_subset$number == i &
                                                                                                 condition_subset$ROINR == cond_name] / profilelength
      
      
    }
    
    
    condition_subset$Value_bgcor[condition_subset$number==i] <- condition_subset$Value[condition_subset$number==i] - mean(condition_subset$Value[condition_subset$number==i & condition_subset$Section=="Background"])
    
    
    
    
    condition_subset$Value_relative[condition_subset$number==i] <-       condition_subset$Value_bgcor[condition_subset$number==i] / 
                                                                    mean(condition_subset$Value_bgcor[condition_subset$Section=="Single-cell segment"])
    
    
  
    
    col_names <- colnames(condition_subset)
    data_frame <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
    colnames(data_frame) <- col_names
    
    data_frame <- rbind(data_frame,condition_subset)
  }
  
  
  return(data_frame)
}

contact_length_segregation_v5 <- function(df) {
  condition_subset <- df
  
  condition_subset$number <- as.numeric(condition_subset$number)
  
  max_values <- vector("numeric", length(unique(condition_subset$number)))
  
  for (i in 1:length(unique(condition_subset$number))) {
    print(paste0("number: ", i))
    
    condition_subset_nr <- subset(condition_subset, number == i)
    
    for (b in 1:length(unique(condition_subset_nr$ROINR))) {
      cond_name <- unique(condition_subset_nr$ROINR)
      cond_name <- cond_name[b]
      print(paste0("RoiNR: ", cond_name))
      
      profilelength <- max(condition_subset_nr$X[condition_subset_nr$ROINR == cond_name])
      print(paste("Profile length:", profilelength, "Pixels"))
      
      ##profilethreshold <- profilelength * contact_length_segregation_threshold
      ##print(paste("Single cell/cell-cell contact threshold: ", profilethreshold, "Pixels"))
      
      condition_subset$Section[condition_subset$number == i &
                                 condition_subset$ROINR == cond_name &
                                 condition_subset$X < (profilelength * 3.5/16)] <- "Single-cell"
      
      condition_subset$Section[condition_subset$number == i &
                                 condition_subset$ROINR == cond_name &
                                 condition_subset$X >= (profilelength * 3.5/16) & 
                                 condition_subset$X < (profilelength * 4.5/16)] <- "Transition"
      
      condition_subset$Section[condition_subset$number == i &
                                 condition_subset$ROINR == cond_name &
                                 condition_subset$X >= (profilelength * 4.5/16) & 
                                 condition_subset$X < (profilelength * 8/16)] <- "Cell-cell-contact"
      
      condition_subset$Section[condition_subset$number == i &
                                 condition_subset$ROINR == cond_name &
                                 condition_subset$X >= (profilelength * 8/16) & 
                                 condition_subset$X < (profilelength * 10/16)] <- "Transition-to-bg"
      
      condition_subset$Section[condition_subset$number == i &
                                 condition_subset$ROINR == cond_name &
                                 condition_subset$X >= (profilelength * 10/16)] <- "Background"
      
      condition_subset$Section[condition_subset$number == i &
                                 condition_subset$Lineprofiletype == "Control"] <- "Single-cell"
      
      # Calculate X_relative for each subset (per ROINR)
      condition_subset$X_relative[condition_subset$number == i &
                                    condition_subset$ROINR == cond_name] <- condition_subset$X[condition_subset$number == i &
                                                                                                 condition_subset$ROINR == cond_name] / profilelength
      
      
    }
    
    
    Max_value <- mean(condition_subset_nr$Value)
    print(paste0("Max value for number: ", i, " --> ", Max_value))
    
    condition_subset$Value_bgcor[condition_subset$number==i] <- condition_subset$Value[condition_subset$number==i] - mean(condition_subset$Value[condition_subset$number==i & condition_subset$Section=="Background"])
    
    condition_subset$Value_relative[condition_subset$number==i] <- condition_subset$Value_bgcor[condition_subset$number==i] / mean(condition_subset$Value_bgcor[condition_subset$number==i & condition_subset$Section=="Single-cell"])
    
    print(colnames(condition_subset))
    
    col_names <- colnames(condition_subset)
    data_frame <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
    colnames(data_frame) <- col_names
    
    data_frame <- rbind(data_frame,condition_subset)
  }
  
  
  return(data_frame)
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

changePlotSize <- function(plot) {
  
  cat("Enter 'q' to quit without saving or provide a plot name:\n")
  plotname <- readline(prompt = "Enter plot name: ")
  
  # Check if the user wants to quit
  if (tolower(plotname) == "q") {
    cat("Operation aborted.\n")
    return(invisible())
  }
  
  # Open a new graphics device
  dev.new()
  
  # Set the size of the plot
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
  
  # Apply your custom theme to the plot (assuming custom_theme is defined)
  plot <- plot + custom_theme  # Apply your custom theme here
  
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
  saveRDS(plot, file = paste0(plotname, ".RData"))
  dev.off()
}

calculate_pv_t.test <- function(expected_value, actual_value, standard_deviation, n) {
  # Calculate the t-statistic
  t_statistic <- (actual_value - expected_value) / (standard_deviation / sqrt(n))
  
  # Calculate the degrees of freedom
  df <- n - 1
  
  # Calculate the two-tailed p-value
  p_value <- 2 * pt(abs(t_statistic), df, lower.tail = FALSE)
  
  return(p_value)
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
  
  if(is.null(group_vars)){
    # Check the number of levels in x_var
    num_levels <- length(unique(data[[x_var]]))
    
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
      geom_point(size = size, stroke = 0, alpha = alpha, aes_string(color = color_var))+
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
    dplyr::summarize(
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
  if (regression_type %in% c("linear", "log", "exp") && any(rsquared_p_values$significance_coefficient != "N.S")) {
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
              hjust = 1, vjust = 1, size = 3)
  
  # Return both the plot and the R-squared values dataframe
  return(list(plot = final_plot, rsquared_values = rsquared_values_df))
}






plot_number <- function(df, nr){
  distvssignal <- ggplot(subset(file,number==nr), aes(x=X_relative, y=Value, color=Section))+
    geom_line()+
    geom_point(size=2, stroke=0, alpha=0.5) +
    
    #geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black") +
    
    labs(y="Normalized Glycocalyx MFI") +
    xlim(0,1)+
    
    labs(x="Relative distance (fraction)")+
    
    scale_x_continuous(breaks=c(0,0.5,1))+
    
    facet_nested(.~number+Lineprofiletype+CellNR)+
    
    
    scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    theme(axis.text.y = element_text(size = 8),
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
          axis.text = element_text(size = 8),
          legend.position = c(0.95, 0.2),
          legend.title=element_blank(),
          strip.background = element_rect(color = "black", fill = NA, size = 0.5),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          legend.background = element_rect(fill = "transparent"))
  
  print(distvssignal)
  
  
  return(distvssignal)
}

plot_number_relative <- function(df, nr){
  distvssignal <- ggplot(subset(file,number==nr), aes(x=X_relative, y=Value_relative, color=Section))+
    geom_line()+
    geom_point(size=2, stroke=0, alpha=0.5) +
    
    #geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black") +
    
    labs(y="Normalized Glycocalyx MFI") +
    xlim(0,1)+
    
    labs(x="Relative distance (fraction)")+
    
    scale_x_continuous(breaks=c(0,0.5,1))+
    
    facet_nested(.~number+Lineprofiletype+CellNR)+
    
    
    scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    theme(axis.text.y = element_text(size = 8),
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(c(2, 1), "mm"),  # Specify the size of the legend key (width, height)
          axis.text = element_text(size = 8),
          legend.position = c(0.95, 0.2),
          legend.title=element_blank(),
          strip.background = element_rect(color = "black", fill = NA, size = 0.5),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          legend.background = element_rect(fill = "transparent"))
  
  print(distvssignal)
  
  
  return(distvssignal)
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
  
  if (any(x > log(.Machine$integer.max))) {
    warning("Exponential regression skipped due to high values.")
  } else {
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
  
  # Create a data frame to store the results
  results_df <- data.frame(
    Formula = c(linear_formula, exponential_formula, logarithmic_formula),
    JB_p_value = c(linear_jarque$p.value, ifelse(is.null(exponential_jarque), NA, exponential_jarque$p.value), logarithmic_jarque$p.value),
    SW_p_value = c(linear_sw$p.value, ifelse(is.null(exponential_sw), NA, exponential_sw$p.value), logarithmic_sw$p.value),
    R_squared = c(linear_r_squared, ifelse(is.null(exponential_r_squared), NA, exponential_r_squared), logarithmic_r_squared),
    Adj_R_squared = c(linear_adj_r_squared, ifelse(is.null(exponential_adj_r_squared), NA, exponential_adj_r_squared), logarithmic_adj_r_squared),
    AIC = c(linear_aic, ifelse(is.null(exponential_aic), NA, exponential_aic), logarithmic_aic),
    BIC = c(linear_bic, ifelse(is.null(exponential_bic), NA, exponential_bic), logarithmic_bic),
    Coefficient = c(linear_slope, ifelse(is.null(exponential_slope), NA, exponential_slope), logarithmic_slope),
    Intercept = c(linear_intercept, ifelse(is.null(exponential_intercept), NA, exponential_intercept), logarithmic_intercept),
    Coefficient_Significance = c(linear_slope_sig, ifelse(is.null(exponential_slope_sig), NA, exponential_slope_sig), logarithmic_slope_sig),
    Intercept_Significance = c(linear_intercept_sig, ifelse(is.null(exponential_intercept_sig), NA, exponential_intercept_sig), logarithmic_intercept_sig),
    Coefficient_p_value = c(linear_slope_p_value, ifelse(is.null(exponential_slope_p_value), NA, exponential_slope_p_value), logarithmic_slope_p_value),
    Intercept_p_value = c(linear_intercept_p_value, ifelse(is.null(exponential_intercept_p_value), NA, exponential_intercept_p_value), logarithmic_intercept_p_value),
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
  linear_std_error <- summary(linear_model)$coefficients[2, "Std. Error"]
  linear_stdev <- sd(linear_residuals)  # Standard deviation for linear relationship
  linear_n <- length(y)  # Sample size for linear relationship
  
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
  exponential_std_error <- NULL
  exponential_stdev <- NULL
  exponential_n <- NULL
  
  if (any(x > log(.Machine$integer.max))) {
    warning("Exponential regression skipped due to high values.")
  } else {
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
    exponential_std_error <- summary(exponential_model)$coefficients[2, "Std. Error"]
    exponential_stdev <- sd(exponential_residuals)  # Standard deviation for exponential relationship
    exponential_n <- length(y)  # Sample size for exponential relationship
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
  logarithmic_std_error <- summary(logarithmic_model)$coefficients[2, "Std. Error"]
  logarithmic_stdev <- sd(logarithmic_residuals)  # Standard deviation for logarithmic relationship
  logarithmic_n <- length(y)  # Sample size for logarithmic relationship
  
  # Create a data frame to store the results
  results_df <- data.frame(
    Formula = c(linear_formula, exponential_formula, logarithmic_formula),
    JB_p_value = c(linear_jarque$p.value, ifelse(is.null(exponential_jarque), NA, exponential_jarque$p.value), logarithmic_jarque$p.value),
    SW_p_value = c(linear_sw$p.value, ifelse(is.null(exponential_sw), NA, exponential_sw$p.value), logarithmic_sw$p.value),
    R_squared = c(linear_r_squared, ifelse(is.null(exponential_r_squared), NA, exponential_r_squared), logarithmic_r_squared),
    Adj_R_squared = c(linear_adj_r_squared, ifelse(is.null(exponential_adj_r_squared), NA, exponential_adj_r_squared), logarithmic_adj_r_squared),
    AIC = c(linear_aic, ifelse(is.null(exponential_aic), NA, exponential_aic), logarithmic_aic),
    BIC = c(linear_bic, ifelse(is.null(exponential_bic), NA, exponential_bic), logarithmic_bic),
    Coefficient = c(linear_slope, ifelse(is.null(exponential_slope), NA, exponential_slope), logarithmic_slope),
    Intercept = c(linear_intercept, ifelse(is.null(exponential_intercept), NA, exponential_intercept), logarithmic_intercept),
    Coefficient_Significance = c(linear_slope_sig, ifelse(is.null(exponential_slope_sig), NA, exponential_slope_sig), logarithmic_slope_sig),
    Intercept_Significance = c(linear_intercept_sig, ifelse(is.null(exponential_intercept_sig), NA, exponential_intercept_sig), logarithmic_intercept_sig),
    Coefficient_p_value = c(linear_slope_p_value, ifelse(is.null(exponential_slope_p_value), NA, exponential_slope_p_value), logarithmic_slope_p_value),
    Intercept_p_value = c(linear_intercept_p_value, ifelse(is.null(exponential_intercept_p_value), NA, exponential_intercept_p_value), logarithmic_intercept_p_value),
    Std_Error = c(linear_std_error, ifelse(is.null(exponential_std_error), NA, exponential_std_error), logarithmic_std_error),
    Stdev = c(linear_stdev, ifelse(is.null(exponential_stdev), NA, exponential_stdev), logarithmic_stdev),
    Sample_Size = c(linear_n, ifelse(is.null(exponential_n), NA, exponential_n), logarithmic_n),
    stringsAsFactors = FALSE
  )
  
  return(results_df)
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

# Function to calculate mean squared error of "loess" smoothing
loess_mse <- function(data, span) {
  smooth <- predict(loess(y ~ x, data = data, span = span))
  mse <- mean((data$y - smooth)^2)
  return(mse)
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
  
  
  existingfile_alldata <- read.csv2(paste0(workdir,"/",existingfilename_alldata), header=T, sep = ",")
  existingfile_alldata <- existingfile_alldata[2:ncol(existingfile_alldata)]
  
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
      file <- read.csv(filepath, header = TRUE, sep = ",", quote = "")
      print(filepath)
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
  colnames(data)[ncol(data)] <- "Subfolder"
  }else{
    
    print(protrusionbasename)
    data2 <- cbind(file1, rep(protrusionbasename, nrow(file1)))
    colnames(data2)[ncol(data2)] <- "Subfolder"
    
    data <- rbind(data,data2)
  }
}

file1 <- data
file1 <- dplyr::select(file1, "X", "Value", "Subfolder", "name")





## fix naming problems -----
library(stringr)

file1 <- file1 %>%
  mutate(num_underscores = str_count(name, "_"))

file1 <- file1 %>%
  mutate(name = if_else(
    str_count(name, "_") == 10,
    str_replace(name, "^((?:[^_]*_){2})", "\\1Unpermeabilized_"),
    name
  ))


file1 <- file1 %>%
  mutate(num_underscores = str_count(name, "_"))

## make sense of column names & organize dataframe-----
rm(file)
file <- separate(file1, name, into=c("Celltype", "Channel", "Permeabilization_status",  "ImageNR", "Acquisitiontype", "ROINR", "CellNR", 
                                              "Lineprofiletype", "Section", "Linethickness", "Reference", "Interactiontype"), sep="_", remove=FALSE)
file$Value <- as.numeric(file$Value)
#rm(file1)
rm(dirfiles)
rm(dirfiles2)



file$Linethickness <- as.numeric(file$Linethickness)

file$Acquisitiontype <- as.numeric(file$Acquisitiontype)

file <- file %>%  mutate(Exposure = "1")

file <- file %>%  mutate(expnr = match(Subfolder, unique(Subfolder)))

file <- file %>% mutate(number= group_indices(.,Exposure, Linethickness, Subfolder, expnr, ImageNR , ROINR, CellNR, Channel, Permeabilization_status))

file <- file %>% mutate(analysisnumber= group_indices(.,expnr, Lineprofiletype, CellNR, ImageNR, Subfolder, Section, ROINR, Exposure, Channel, Permeabilization_status))

file <- file %>% mutate(ID= group_indices(.,expnr, ImageNR, Subfolder , Exposure, Linethickness))


file$number <- as.character(file$number)


file$X_micron <- (file$X * pxsize)-pxsize


file <- contact_length_segregation_v6(file)


file <- file %>%  group_by(number) %>%  mutate(Value_bgcor = Value-mean(Value[Section=="Background"]))

file <- file %>%  group_by(number) %>%  mutate(Value_relative = Value_bgcor / mean(Value_bgcor[Section=="Single-cell segment"]))


file <- file %>%  group_by(Exposure, Linethickness, Section) %>%  mutate(Value_vsmean = Value_bgcor / mean(Value_bgcor))

file <- file %>%
  group_by(analysisnumber) %>%
  mutate(X_micron_normalized = ifelse(Section!="Background",
                                      X_micron - max(X_micron[Section != "Background"]),
                                      X_micron + max(X_micron[Section != "Background "]))) %>%
  ungroup()

file <- file %>%
  group_by(analysisnumber) %>% 
  mutate(X_micron_normalized = ifelse(Section=="Cell-cell-contact", 
                                      X_micron_normalized-min(X_micron_normalized),
                                      X_micron_normalized)) %>% 
  ungroup()




file$analysisnumber <- as.character(file$analysisnumber)



recalc_transition_zone=TRUE

file <- file %>%  mutate(Channel = ifelse(Channel=="ALCAM93", "ALCAM", Channel))

## pool the X coordinates by desired number of consecutive pixels------
# Calculate the means using the sliding window of three consecutive rows of "X"

pool_pixels = TRUE
grouped_pixels = ceiling(0.21 / pxsize)



file_averaged <- file %>%
  group_by(number, Section) %>% 
  mutate(X_group = ceiling(row_number()/grouped_pixels))

file_averaged <- file_averaged %>% 
  group_by(Exposure, Linethickness, Lineprofiletype, name, Subfolder, ImageNR, ROINR, Channel,
           CellNR, number, analysisnumber, Section, X_group, Permeabilization_status) %>% 
  dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.), .names = "{.col}")) %>%
  ungroup()



if(pool_pixels==TRUE){
  file_unpooled <- file
  file <- file_averaged
}



file$number <- as.character(file$number)





## Define Transition zone and Transition zone length ----

## Define Segment=Transition based on the distance between the mean of the single-cell segment and cell-cell contact
file <- file %>% 
  group_by(Exposure, Lineprofiletype, number, ImageNR, expnr) %>%
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
  ) %>%
  mutate(
    Section = case_when(
      X_micron_normalized > lowerXthr & 
        X_micron_normalized <= upperXthr &
        Value_relative > (single_cell_mean + 0 * single_cell_sd) & 
        Value_relative <= (cell_contact_mean + 0 * cell_contact_sd) ~ "Transition",
      TRUE ~ Section
    )
  ) 

## Define the length and an id for each individual segment where Section==Transition for further filtering
## since the transition zone should be a single, uninterrupted section
file <- file %>% group_by(number) %>% 
  arrange(X_micron_normalized) %>%
  mutate(consecutive_transitions = ifelse(Section == "Transition", ave(Section == "Transition", 
                                                                       cumsum(Section != lag(Section, default = first(Section))), FUN = cumsum), 0))
num_iterations <- 100

for (i in 1:num_iterations) {
  file <- file %>%
    group_by(number) %>%  arrange(X_micron_normalized) %>% 
    mutate(consecutive_transitions = ifelse(lead(consecutive_transitions + 1) & lead(consecutive_transitions != 0) & consecutive_transitions !=0, 
                                            lead(consecutive_transitions), consecutive_transitions),
           consecutive_transitions = ifelse(is.na(consecutive_transitions) & lag(!is.na(consecutive_transitions)), 
                                            lag(consecutive_transitions), consecutive_transitions)
    )
}


file <- file %>%
  group_by(number) %>%  arrange(X_micron_normalized) %>% 
  mutate(
    unique_section_id = cumsum(consecutive_transitions == 0),
    unique_section_id = ifelse(consecutive_transitions == 0 , NA, unique_section_id),
    unique_section_id = unique_section_id +1,
    unique_section_id = as.factor(unique_section_id)
  )



Xnormvalue <- file %>% group_by(number) %>% 
  filter(Section == "Transition") %>%
  dplyr::summarize(minX = min(X_micron_normalized[consecutive_transitions==max(consecutive_transitions)]),
            maxX = max(X_micron_normalized[consecutive_transitions==max(consecutive_transitions)]))
Xnormvalue <- Xnormvalue %>% 
  mutate(midX_transition = (maxX + minX)/2) %>% 
  dplyr::select(number, midX_transition, maxX, maxX)

file <- file %>%
  left_join(Xnormvalue, by = "number") 



## calculate the delta of each section and only include the section with the largest delta per line profile
minmax_persection <- file %>% group_by(number, unique_section_id) %>% 
  filter(Section == "Transition")
minmax_persection <- minmax_persection %>% 
  dplyr::summarize(minValue = min(Value_relative),
            maxValue = max(Value_relative)) %>% 
  mutate(Section_Delta_rel = maxValue - minValue) %>% ungroup () %>% group_by(number) %>% 
  mutate(maxDelta = ifelse(Section_Delta_rel== max(Section_Delta_rel), TRUE, FALSE)) %>% 
  dplyr::select(number, unique_section_id, maxDelta)

file <- file %>%
  left_join(minmax_persection, by = c("number", "unique_section_id")) 

## redefine Single-cell segment and Cell-cell contact based on their position respective to the newly defined transition zone
file <- file %>% group_by(number) %>% 
  mutate(consecutive_transitions = as.numeric(consecutive_transitions),
         Section = ifelse(Section=="Transition" & maxDelta==FALSE & X_micron_normalized > midX_transition , 
                          "Cell-cell-contact", Section),
         Section = ifelse(Section=="Transition" & maxDelta==FALSE & X_micron_normalized < midX_transition , 
                          "Single-cell segment", Section),
  ) 


file <- file %>%  group_by(number) %>% 
  mutate(Section = ifelse(consecutive_transitions>0 &
    Section!="Transition" & Section!="Background" &
                            X_micron_normalized > midX_transition, "Cell-cell-contact", Section),
         Section = ifelse(consecutive_transitions>0 &
           Section!="Transition" & Section !="Background" &
                            X_micron_normalized < midX_transition, "Single-cell segment", Section)
  )


## calculate the delta of each section and only include the section with the largest delta per line profile
minmax_persection <- file %>% group_by(number, unique_section_id) %>% 
  filter(Section == "Transition") %>%
  dplyr::summarize(minValue = min(Value_relative),
            maxValue = max(Value_relative) ,
            minValue = Value_relative[which.min(X)],
            maxValue = Value_relative[which.max(X)])%>% 
  mutate(Section_Delta_rel = maxValue - minValue) %>% ungroup () %>% group_by(number) %>% 
  mutate(maxDelta = ifelse(Section_Delta_rel== max(Section_Delta_rel), TRUE, FALSE)) %>% 
  dplyr::select(number, unique_section_id, maxDelta)


file <- file %>%
  left_join(minmax_persection, by = c("number", "unique_section_id")) 

file <- file %>%  group_by(number) %>% 
  mutate(Section = ifelse(consecutive_transitions>0 &
                            Section!="Transition" & Section!="Background" &
                            X_micron_normalized > midX_transition, "Cell-cell-contact", Section),
         Section = ifelse(consecutive_transitions>0 &
                            Section!="Transition" & Section !="Background" &
                            X_micron_normalized < midX_transition, "Single-cell segment", Section)
  )


Xnormvalue <- file %>% group_by(number) %>% 
  filter(Section == "Transition") %>%
  dplyr::summarize(minX = min(X_micron_normalized[consecutive_transitions==max(consecutive_transitions)]),
            maxX = max(X_micron_normalized[consecutive_transitions==max(consecutive_transitions)]))
Xnormvalue <- Xnormvalue %>% 
  mutate(midX_transition = (maxX + minX)/2) %>% 
  dplyr::select(number, midX_transition, minX, maxX)

file <- file %>%
  left_join(Xnormvalue, by = "number") 


file <- file %>%  dplyr::rename(midX_transition = midX_transition.x)

file <- file %>%  group_by(number) %>% 
  mutate(Section = ifelse(consecutive_transitions>0 &
                            Section!="Transition" & Section!="Background" &
                            X_micron_normalized > midX_transition, "Cell-cell-contact", Section),
         Section = ifelse(consecutive_transitions>0 &
                            Section!="Transition" & Section !="Background" &
                            X_micron_normalized < midX_transition, "Single-cell segment", Section)
  )


calcmid <- file %>% group_by(number) %>% 
  filter(Section == "Transition") %>%
  dplyr::summarize(minX = min(X_micron_normalized),
            maxX = max(X_micron_normalized),
            Transition_mid = mean(X_micron_normalized)) %>% 
  dplyr::select(number, Transition_mid)
file <- file %>%
  left_join(calcmid, by = "number") 

#file <- file %>%  rename(Transition_mid = Transition_mid.x)

file <- file %>%  group_by(number) %>% 
  mutate(Section = ifelse(consecutive_transitions>0 &
                            Section!="Transition" & Section!="Background" &
                            X_micron_normalized > midX_transition, "Cell-cell-contact", Section),
         Section = ifelse(consecutive_transitions>0 &
                            Section!="Transition" & Section !="Background" &
                            X_micron_normalized < midX_transition, "Single-cell segment", Section)
  )



file <- file %>%  group_by(number) %>% 
  mutate(Section= ifelse(Lineprofiletype=="Control", "Single-cell segment", Section))



# Now, add the "Transitionzone_length" column
file <- file %>%
  group_by(number) %>%
  mutate(
    Transitionzone_length = sum(Section == "Transition")*pxsize * grouped_pixels
  ) %>%
  ungroup()



# Calculate Transitionzone_length per number
transition_lengths <- file %>%
  dplyr::group_by(number) %>%
  dplyr::summarise(Transitionzone_length = sum(Section == "Transition") * pxsize* grouped_pixels)

# Join the calculated values back to the original dataframe
file <- dplyr::left_join(file, transition_lengths, by = "number")

file <- file %>%  dplyr::select(-Transitionzone_length.y)
file <- file %>%  dplyr::rename(Transitionzone_length = Transitionzone_length.x)

file <- file %>% 
  group_by(Exposure, Linethickness, number) %>%  
  mutate(
    Lineprofiletype = ifelse(
      Lineprofiletype == "Contact" && any(Section == "Transition"),
      "Contact",  
      "Control"  
    )
  )


file <- file %>% group_by(Exposure, Linethickness, number) %>%  
  mutate(X_fromtransitionzone = ifelse(Lineprofiletype=="Contact",
                                       X_micron_normalized - min(X_micron_normalized[Section=="Transition"], na.rm=TRUE), NA))





recalc_value_relative=TRUE
if(recalc_value_relative==TRUE){
  file <- file %>%  group_by(number) %>%  mutate(Value_relative = Value_bgcor / mean(Value_bgcor[Section=="Single-cell segment"]))
  file <- file %>%  group_by(Exposure, Linethickness, Section) %>%  mutate(Value_vsmean = Value_bgcor / mean(Value_bgcor))
  file <- file %>%  group_by(number) %>% 
    mutate(Contact_enrichment = ifelse(mean(Value_relative[Section=="Cell-cell-contact"] > 
                                              mean(Value_relative[Section=="Single-cell segment"])),
                                       TRUE,
                                       FALSE))
  file <- file %>%  group_by(number) %>%  mutate(Transitionzone_length=ifelse(Lineprofiletype=="Contact" & Contact_enrichment==TRUE,
                                                                              Transitionzone_length,
                                                                              NA))
}


file$Exposure[file$Exposure==0.13] <- 0.125




## recalculate Value_relative based on the newly-defined Single-cell and cell-cell-contact sections -----


file <- file %>%  group_by(number) %>%  mutate(Value_relative = Value_bgcor / mean(Value_bgcor[Section=="Single-cell segment"]))



## save the file ----

write.csv(x=file, "Cell-cell-contact-data.csv")



## dplyr::summarize the data----

file_summarized <- file %>%
  group_by(Subfolder, expnr , number, ID, ImageNR , CellNR, ROINR, Lineprofiletype, Exposure , Linethickness, Section, Transitionzone_length, Channel, Permeabilization_status) %>%
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) 


data_optimallinethickness <- subset(file_summarized, Linethickness==selected_optimallinethickness)


file_optimallinethickness <- subset(file, Linethickness==selected_optimallinethickness)


  all_summarized <- data_optimallinethickness %>%
    group_by(Section) %>%
    dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) 
  


transition_zone_coef=FALSE

if(transition_zone_coef==TRUE){
print("Regression analysis calculated for transition zone")

# Apply the function to subsets of the data
     subset_results <- file_optimallinethickness %>%
       filter(Lineprofiletype=="Contact") %>%
       group_by(number, Section) %>% 
       filter(Section=="Transition") %>% 
       
       
    summarise(results = calculate_residuals_and_tests_limitcorrected(subset_results$X_micron_normalized, subset_results$Value_relative)) %>% 
    ungroup()
#  
    subset_results <- subset_results %>%  filter(results$Formula== "y ~ x")
    subset_results <- subset_results %>%
     unnest(results) %>%
      filter(Section=="Transition") %>% 
     select(number, Section, Coefficient, Coefficient_Significance,
            Coefficient_p_value)
  
 
      #Print the results
     print(subset_results)
      
    data_optimallinethickness <- data_optimallinethickness %>%
      left_join(subset_results, by = c("number", "Section"))
  
}

transitionzone_uniques <- file_optimallinethickness %>%  group_by(Transitionzone_length,number) %>%  summarise(Transitionzone_length= mean(Transitionzone_length))

## calculate some statistics------



n_value_relative_percell <- length(data_optimallinethickness$Value_relative[data_optimallinethickness$Section=="Cell-cell-contact"])
mean_value_relative_percell <- mean(data_optimallinethickness$Value_relative[data_optimallinethickness$Section=="Cell-cell-contact"])
sd_value_relative_percell <- sd(data_optimallinethickness$Value_relative[data_optimallinethickness$Section=="Cell-cell-contact"])
expected_value_relative_percell <- 2

pv_t.test <- calculate_pv_t.test(expected_value_relative_percell,
                    mean_value_relative_percell,
                    sd_value_relative_percell,
                    n_value_relative_percell)

write.csv2(pv_t.test, "Cell-cell-contact.rel.vs.2.csv")


nr_contact_enriched <- file_optimallinethickness

## save summarizing files------
write.csv(x=file_summarized, "Cell-cell-contact-data_summarized.csv")




## plot the data (choose which to plot - they require some tweaking to run)-------
  


## 1. plot profiles versus mfi------

f <- ggplot(data=subset(file_optimallinethickness, Lineprofiletype=="Contact" & Section=="Single-cell segment"),aes(x=X_micron))+
  geom_line(aes(y=Value, color=Section, group=Section))+
  
  facet_nested(.~expnr + ImageNR  + CellNR, labeller=label_both)+
  ylab("MFI")

print(f)

g <- ggplot(subset(file_optimallinethickness, Lineprofiletype=="Contact" & Section!="Background"),aes(x=X_micron))+
  geom_line(aes(y=Value_relative, group=Section, color=Section))
  
  #facet_nested(.~expnr+ Subfolder+ ImageNR  + CellNR + number, labeller=label_both)


print(g)


plot.specific.profile <- function(data){
library(svDialogs)
  
print(unique(data$number[data$Lineprofiletype=="Contact" & data$Section!="Transition-to-bg" & data$Section!="Background"]))
  
  nr <- dlg_input("Number?")$res
  
  g <- ggplot(data=subset(file_optimallinethickness, Lineprofiletype=="Contact" & number==nr),aes(x=X_micron))+
    geom_line(aes(y=Value_bgcor, group=number, color=Section))+
    
    facet_nested(.~expnr+ ImageNR  + CellNR, labeller=label_both)+
   # scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
    ylab("Background corrected MFI")
  
  return(g)
}


## 2. plot linethickness vs MFI -----



a <- ggplot(data=subset(file_summarized, Section=="Single-cell segment"), aes(x=Linethickness, y = Value_bgcor))+
  geom_smooth(aes(color=Exposure, group=Exposure), method="loess", span=0.1,na.rm=TRUE)+
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  ylab("MFI")+
  xlab("Line thickness (px)")+
  facet_wrap(.~Section)


print(a)

a <- a+ custom_theme

print(a)

## 3a. Plot histograms without mean and sd of Value_relative to see variability  ------

set.to.plot <- subset(file_optimallinethickness, Section != "Background" & Lineprofiletype=="Contact")


histo <- ggplot(set.to.plot, 
                aes(x = Value_relative, fill = Section)) +
  geom_histogram(data=subset(file_optimallinethickness, Section =="Cell-cell-contact"& Lineprofiletype=="Contact"), bins=20, 
                 aes(y = ..ncount..), alpha=0.5)+
  geom_histogram(data=subset(file_optimallinethickness, Section =="Single-cell segment"& Lineprofiletype=="Contact"), bins=20, 
                 aes(y = ..ncount..), alpha=0.5)+
  geom_density(data=subset(set.to.plot, Section!="Transition"), aes(Value_relative, after_stat(scaled),color=Section), alpha=0)

histo <- histo + facet_wrap(.~number)

print(histo)

library(ggplot2)
library(dplyr)

## 3b. Plot histograms of Value_relative with mean and sd to see variabilty & calculate parameters Gcx_enrichment_contact and Gcx_density ------



file_optimallinethickness <- file_optimallinethickness %>% 
  group_by(number) %>%  mutate(ID = group_indices())

# Calculate mean SD and CV for each facet group
summary_data <- file_optimallinethickness %>%
  filter(Section %in% c("Cell-cell-contact", "Transition", "Single-cell segment")) %>%
  group_by(Lineprofiletype, number, Section, ID) %>%
  dplyr::summarize(mean_value = mean(Value_bgcor),
                   sd_value = sd(Value_bgcor), 
                   mean_value_relative = mean(Value_relative),
                   sd_value_relative = sd(Value_relative),
                   cv = sd_value / mean_value * 100,
                   Transitionzone_length= mean(Transitionzone_length)
  )

summary_data$Section <- as.factor(summary_data$Section)
summary_data$Section <- factor(summary_data$Section, levels=c("Single-cell segment", "Transition", "Cell-cell-contact"))

summary_data <- summary_data %>% ungroup() %>%  group_by(number, ID) %>%  
  mutate(Foldchange = ifelse(Section=="Cell-cell-contact", mean_value[Section=="Cell-cell-contact"] / mean_value[Section=="Single-cell segment"], NA),
         Delta = Transitionzone_length / Foldchange)

summary_data$Section <- as.factor(summary_data$Section)

summary_perprofile <- summary_data


summary_data <- summary_data %>%  filter(Section!="Transition" & Lineprofiletype=="Contact")


Singlecell_values <- summary_data$mean_value[summary_data$Section=="Single-cell segment"]
Contact_values <- summary_data$mean_value[summary_data$Section=="Cell-cell-contact"]
Ftest_variance <- var.test(Singlecell_values, Contact_values)
Ftest_results_df <- data.frame(
  statistic = Ftest_variance$statistic,
  parameter = Ftest_variance$parameter,
  p.value = Ftest_variance$p.value,
  conf.int = unlist(Ftest_variance$conf.int),
  estimate = unlist(Ftest_variance$estimate),
  names = unlist(Ftest_variance$data.name)
)
print(Ftest_results_df)

summary_all <- summary_data %>% 
  group_by(Section) %>% 
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) 


# Create the ggplot object
histo <- ggplot(subset(file_optimallinethickness, Section!="Transition" & Lineprofiletype=="Contact"), aes(x = Value_relative)) +
  geom_histogram(data=subset(file_optimallinethickness, Section =="Cell-cell-contact"& Lineprofiletype=="Contact"), bins=30, 
                 aes(y = ..ncount.., fill=Section), alpha=0.5)+
  geom_histogram(data=subset(file_optimallinethickness, Section =="Single-cell segment"& Lineprofiletype=="Contact"), bins=30, 
                 aes(y = ..ncount.., fill=Section), alpha=0.5)+
  geom_density(aes(Value_relative, after_stat(scaled), color = Section)) +
  facet_wrap(. ~ ID) +
  
  # Add vertical lines for mean +/- SD
  geom_vline(data = summary_data, aes(xintercept = mean_value_relative, linetype=Section)) +
  geom_vline(data = summary_data, aes(xintercept = mean_value_relative - sd_value_relative,linetype=Section), color="gray50") +
  geom_vline(data = summary_data, aes(xintercept = mean_value_relative + sd_value_relative, linetype=Section), color="gray50") 


histo <- histo+custom_theme + theme(legend.position="right")+
  scale_color_manual(values=c("gray","red", "green4"))+
  scale_fill_manual(values=c("red", "green4"))

# Print the plot
print(histo)


summary_all_quantifications <- summary_data %>% 
  group_by(number, ID) %>% 
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) 

ggplot(summary_data, aes(x=Section, y=mean_value_relative))+
  geom_violin(scale="width")+
  geom_boxplot(width=0.2)+
  geom_jitter(width=0.2)

mean_value_relative_all <- mean(summary_data$mean_value_relative[summary_data$Section=="Cell-cell-contact"], na.rm=TRUE)
sd_value_relative_all <- sd(summary_data$mean_value_relative[summary_data$Section=="Cell-cell-contact"], na.rm=TRUE)

transitionzone_mean_all <- mean(data_optimallinethickness$Transitionzone_length[summary_data$Section=="Cell-cell-contact"], na.rm=TRUE)
transitionzone_sd_all <- sd(data_optimallinethickness$Transitionzone_length[summary_data$Section=="Cell-cell-contact"], na.rm=TRUE)

summary_data <- summary_data %>%  group_by(number, ID) %>% 
  mutate(Gcx_density = ifelse(mean_value_relative < 2, "Underrepresented", "Overrepresented"),
         Gcx_density = ifelse(Section=="Single-cell segment" | Section=="Transition", NA, Gcx_density),
         Gcx_contact_enrichment = ifelse(mean_value_relative < mean_value_relative_all - 0.5 * sd_value_relative_all, "Low", "Mid"),
         Gcx_contact_enrichment = ifelse(mean_value_relative > mean_value_relative_all + 0.5 * sd_value_relative_all, "High", Gcx_contact_enrichment),
         Gcx_contact_enrichment = ifelse(Section=="Single-cell segment", NA, Gcx_contact_enrichment),
         Gcx_contact_enrichment = as.factor(Gcx_contact_enrichment),
         Gcx_contact_enrichment = factor(Gcx_contact_enrichment, levels=c("Low", "Mid", "High")))

data_optimallinethickness <- data_optimallinethickness %>% group_by(number, ID) %>% 
  mutate(
         Transitionzone_classification = ifelse(Transitionzone_length < transitionzone_mean_all - 0.5 * transitionzone_sd_all, "Short", "Med"),
         Transitionzone_classification = ifelse(Transitionzone_length > transitionzone_mean_all + 0.5 * transitionzone_sd_all, "Long", Transitionzone_classification),
         Transitionzone_classification = ifelse(Section=="Single-cell segment", NA, Transitionzone_classification),
         Transitionzone_classification = as.factor(Transitionzone_classification),
         Transitionzone_classification = factor(Transitionzone_classification, levels=c("Short", "Med", "Long"))
  )

file_optimallinethickness$Gcx_contact_enrichment <- summary_data$Gcx_contact_enrichment[match(file_optimallinethickness$number, summary_data$number)]
file_optimallinethickness$Gcx_density <- summary_data$Gcx_density[match(file_optimallinethickness$number, summary_data$number)]

file_optimallinethickness$Transitionzone_classification <- data_optimallinethickness$Transitionzone_classification[match(file_optimallinethickness$number, data_optimallinethickness$number)]


summary_perprofile$Gcx_density <- summary_data$Gcx_density[match(summary_perprofile$number, summary_data$number)]
summary_perprofile$Gcx_density <- summary_data$Gcx_density[match(summary_perprofile$number, summary_data$number)]
data_optimallinethickness$Gcx_density <- summary_data$Gcx_density[match(data_optimallinethickness$number, summary_data$number)]
data_optimallinethickness$Gcx_density <- summary_data$Gcx_density[match(data_optimallinethickness$number, summary_data$number)]





## 3c. Plot Example histograms of high/med/low glycocalyx enrichment in cell-cell contact -----
set.seed(123)

selected_numbers <- file_optimallinethickness %>% filter(!is.na(Gcx_contact_enrichment)) %>% 
  group_by(Gcx_contact_enrichment) %>%
  sample_n(1) %>%
  pull(number)

subset_data <- file_optimallinethickness %>%
  filter(number %in% selected_numbers)

set.to.plot <- subset_data

nbins = 30

histo <- ggplot(set.to.plot, 
                aes(x = Value_relative, fill = Section)) +
  geom_histogram(data=subset(set.to.plot, Section =="Cell-cell-contact"& Lineprofiletype=="Contact"), bins=nbins, 
                 aes(y = ..ncount..), alpha=0.5)+
  geom_histogram(data=subset(set.to.plot, Section =="Single-cell segment"& Lineprofiletype=="Contact"), bins=nbins, 
                 aes(y = ..ncount..), alpha=0.5)+
  geom_histogram(data=subset(set.to.plot, Section =="Background"& Lineprofiletype=="Contact"), bins=nbins, 
                 aes(y = ..ncount..), alpha=0.5)+
  geom_density(data=subset(set.to.plot, Section!="Transition"), aes(Value_relative, after_stat(scaled),color=Section), alpha=0)

histo <- histo + facet_nested(.~Gcx_contact_enrichment)

print(histo)  



## 4a. Plot Overview plot with all thresholds per lineplot (X_micron_normalized vs. Gcx enrichment factor)----


#determine  optimal 'span' of loess smoothening method:
# Define a range of span values to test
span_values <- seq(0.1, 1.0, by = 0.1)
data <- data.frame(x=file_optimallinethickness$X_micron_normalized[file_optimallinethickness$Lineprofiletype=="Contact"
                                                   & file_optimallinethickness$Section!="Background"],
                   y=file_optimallinethickness$Value_relative[file_optimallinethickness$Lineprofiletype=="Contact"
                                                         & file_optimallinethickness$Section!="Background"])

data <- na.omit(data)
# Perform cross-validation for each span value
cv_results <- sapply(span_values, function(span) {
  mse <- loess_mse(data,span)
  return(mse)
})

# Find the span value with the lowest cross-validated mean squared error
best_span <- span_values[which.min(cv_results)]
print(paste("Best span for loess smoothening:", best_span))


#file_optimallinethickness_sampled <- sample_n(file_optimallinethickness, size=nrow(file_optimallinethickness)/2)
#file_optimallinethickness_sampled <- subset(file_optimallinethickness_sampled,Lineprofiletype=="Contact" & Section!= "Background" )
  

plotset <- subset(file_optimallinethickness, 
                  Lineprofiletype=="Contact" & 
                    Section!= "Background")
  
k <- ggplot(data=plotset,aes(x=X_micron_normalized))+
    geom_point(data=plotset, aes(x=X_micron_normalized, y=Value_relative, color=Section), size=0.5, stroke=0)+
    geom_line(data=plotset, aes(y=Value_relative, group=number, color=Section), size=1)+
    #geom_smooth(aes(y=Value_relative), color="black", method="gam", span=0.5)+
    #geom_hline(yintercept=2, linetype=2)+
    ylab("MFI of rel. to single-cell- \n segment of all cells")+
    xlab("Plot profile (micron)")
    #ylim(c(0.5,3))+
    #xlim(-5,5)
  
k <- k +
  facet_wrap(.~number, ncol=5)
  

df_forplot <- plotset %>%  group_by(Lineprofiletype, Exposure, Linethickness, number, ImageNR, expnr) %>% 
  distinct(single_cell_mean, cell_contact_mean) %>%  filter(Lineprofiletype=="Contact")

midtransition <- file_optimallinethickness %>% 
  group_by(number, Lineprofiletype,Section) %>% 
  dplyr::summarize(midX_transition = mean(Transition_mid)) %>%  filter(Lineprofiletype=="Contact" & Section=="Transition")


overview <- k + geom_vline(data=midtransition, aes(xintercept=midX_transition))
print(overview)


## 4b. Plot All individual lineplots together with smoothened average -----

stdev_line <- plotset %>% ungroup() %>%  group_by(X_micron_normalized)  %>% 
  dplyr::summarize(sd_perXcoord = sd(Value_relative, na.rm=TRUE))
merged_data <- left_join(plotset, stdev_line, by = "X_micron_normalized")


setforplot <- subset(file_optimallinethickness, Section!="Background")

setforplot <- setforplot %>%  group_by(number) %>%  filter(any(Section=="Transition"))


library(ggplot2)

k <- ggplot(data = subset(setforplot), aes(x = X_micron_normalized, y = Value_relative)) +
  #geom_point(size = 0.5, stroke = 0) +
  geom_line(data=subset(setforplot), aes(y=Value_relative, group=interaction(number), color=Section), size=0.5, alpha=0.5)+
  #geom_ribbon(aes(ymin = predict(loess(Value_relative ~ X_micron_normalized, data = plotset, span = best_span), plotset) - 0.5 *sd(Value_relative),
  #                ymax = predict(loess(Value_relative ~ X_micron_normalized, data = plotset, span = best_span), plotset) + 0.5* sd(Value_relative)),
  #            alpha = 0.5, fill = "gray", stroke=0) +
  ylab("Gcx enrichmentfactor") +
  xlab("Plot profile (micron)") +
  theme(legend.position="none")+
  xlim(-5, 5)+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12))


k <- k + geom_hline(yintercept=2, linetype=2) 

print(k)

lineplot <- k

k <- k + 
  stat_smooth(method = "loess", span = best_span, se = TRUE, color = "black", alpha=0.3, size=0.5, aes(group=Channel))

lineplot.smoothened <- k

facetted <- lineplot.smoothened + facet_wrap(Permeabilization_status~Channel, nrow=1)

print(facetted)

## 4c. Example plot of line trajectory (X_micron_normalized vs. MFI_bgcor)----

example_subset <- subset(file_optimallinethickness, Lineprofiletype=="Contact" & Section!="Background"&
                           CellNR==2 & ImageNR==1 & ROINR==2 & expnr==2)

example_subset_thresholds <- left_join(df_forplot, example_subset)


example <- ggplot(data=example_subset,aes(x=X_micron_normalized))+
  geom_line(aes(y=Value_bgcor, color=Section, group=number), alpha=0.5)+
  geom_point(aes(y=Value_bgcor, color=Section), stroke=0, size=0.5)+
  #geom_smooth(aes(y=Value_relative), color="black", method="gam", span=0.1)+
  ylab("MFI")+
  
  #geom_vline(xintercept=lowerXthr, linetype=2)+
  #geom_vline(xintercept=upperXthr, linetype=2)+
  xlab("Plot profile (micron)")+
  xlim(c(-5,5))
  
  example <- example + 
    
  #  geom_hline(data=example_subset_thresholds, aes(yintercept=single_cell_mfi_mean), linetype="dashed", color="green4")+
  #  geom_hline(data=example_subset_thresholds, aes(yintercept=cell_contact_mfi_mean), linetype="dashed", color="red3")+
  
  
  theme(legend.position="none")
  print(example)



## 4d. Plot Example plot with all thresholds per lineplot (X_micron_normalized vs. Gcx enrichment factor)-----
  
  
  
  
  #determine  optimal 'span' of loess smoothening method:
  # Define a range of span values to test
  span_values <- seq(0.1, 1.0, by = 0.1)
  data <- data.frame(x=file_optimallinethickness$X_micron_normalized[file_optimallinethickness$Lineprofiletype=="Contact"
                                                                     & file_optimallinethickness$Section!="Background"],
                     y=file_optimallinethickness$Value_relative[file_optimallinethickness$Lineprofiletype=="Contact"
                                                                & file_optimallinethickness$Section!="Background"])
  
  data <- na.omit(data)
  # Perform cross-validation for each span value
  cv_results <- sapply(span_values, function(span) {
    mse <- loess_mse(data,span)
    return(mse)
  })
  
  # Find the span value with the lowest cross-validated mean squared error
  best_span <- span_values[which.min(cv_results)]
  print(paste("Best span for loess smoothening:", best_span))
  
  
 # file_optimallinethickness_sampled <- sample_n(file_optimallinethickness, size=nrow(file_optimallinethickness)/2)
 # file_optimallinethickness_sampled <- subset(file_optimallinethickness_sampled,Lineprofiletype=="Contact" & Section!= "Background" )
  
  
  exampleplots <-  file_optimallinethickness %>%  subset(Lineprofiletype=="Contact" & Section!= "Background")
   
  
  identification <- exampleplots %>% group_by(Gcx_density) %>%
    sample_n(1) %>%
    pull(number)
  
  
  exampleplots <- exampleplots %>%
    filter(number %in% identification)
  
  k <- ggplot(data=exampleplots,aes(x=X_micron_normalized))+
    geom_point(data=exampleplots, aes(x=X_micron_normalized, y=Value_relative, shape=Section, color=Section), stroke=0)+
    geom_line(data=exampleplots, aes(y=Value_relative, group=number, color=Section), size=0.5)+
    #geom_smooth(aes(y=Value_relative), color="black", method="gam", span=0.5)+
    #geom_hline(yintercept=2, linetype=2)+
    geom_vline(xintercept=lowerXthr, linetype=2)+
    geom_vline(xintercept=upperXthr, linetype=2)+
    ylab("MFI of rel. to single-cell- \n segment of all cells")+
    xlab("Plot profile (micron)")
  #ylim(c(0.5,3))+
  #xlim(-5,5)
  
  k <- k +
    facet_wrap(.~number, ncol=5)
  
  
  df_forplot <- exampleplots %>%  group_by(Lineprofiletype, Exposure, Linethickness, number, ImageNR, expnr) %>% 
    distinct(single_cell_mean, cell_contact_mean) %>%  filter(Lineprofiletype=="Contact")
  
  k <- k +
    geom_hline(data=df_forplot, aes(yintercept=single_cell_mean), linetype="dashed")+
    geom_hline(data=df_forplot, aes(yintercept=cell_contact_mean), linetype="dashed")
  
  
  overview <- k
  print(overview)
  
  
  
  
## 5a. Violin plot of MFI vs. Section ----
  
  library(rstatix)
  groups <- d
  
  
data_optimallinethickness$Section <- as.factor(data_optimallinethickness$Section)
data_optimallinethickness$Section <- factor(data_optimallinethickness$Section, levels=c("Background", "Single-cell segment", "Transition", "Cell-cell-contact"))

setforplot <- subset(data_optimallinethickness, Section!="Background")  



#setforplot <- setforplot  %>% ungroup() %>% 
#  group_by(Section) %>% 
#  mutate(Value_relative = ifelse(Section=="Single-cell segment", Value_bgcor/ mean(Value_bgcor), Value_relative)
#         ) %>%  group_by(number) %>%  
#  mutate(Value_relative = ifelse(Section!="Single-cell segment", Value_bgcor/ Value_bgcor[Section=="Single-cell segment"], Value_relative)
#  )




setforplot$Section <- as.factor(setforplot$Section)
setforplot$Section <- factor(setforplot$Section, levels=c("Single-cell segment", "Transition", "Cell-cell-contact"))


setforplot <- setforplot %>%  mutate(Permeabilization_status = as.factor(Permeabilization_status),
                                     Permeabilization_status = factor(Permeabilization_status, levels=c("Unpermeabilized", "Saponin"))
                                     )

setforplot <- setforplot %>%  filter(Channel=="ALCAM") 

ggplot(setforplot, aes(x=Permeabilization_status, y=Value_relative))+
  geom_violin(scale="width")+
  geom_boxplot(width=0.5, outlier.shape=NA)+
  geom_point()+
  facet_wrap(.~Section)


ggplot(setforplot, aes(x=Permeabilization_status, y=Value_bgcor))+
  geom_violin(scale="width")+
  geom_boxplot(width=0.5, outlier.shape=NA)+
  geom_point()+
  facet_wrap(.~Section)


shapiro.test(setforplot$Value_bgcor[setforplot$Section=="Single-cell segment"])

violin <- custom_violinplot_with_statistics(subset(setforplot, Section=="Cell-cell-contact"), group_vars=NULL,x_var="Permeabilization_status",
                                            y_var="Value_relative",jitter=F, plot_mean=FALSE)
stats <- violin$statistical_test
print(stats)
violin <- violin$plot + ylab("GCX Enrichment") + theme(legend.position = "none")+
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D"))+geom_point() 

print(violin)


dunn.stats <- setforplot %>% subset(Section=="Cell-cell-contact") %>% ungroup() %>% 
  group_by(Channel) %>%
  dunn_test(Value_relative ~ Permeabilization_status)

print(dunn.stats)



dunn.stats <- setforplot %>% subset(Section=="Cell-cell-contact") %>% ungroup() %>% 
  group_by(Channel) %>%
  wilcox_test(Value_relative ~ Permeabilization_status)

print(dunn.stats)


violin.enrichment <- violin + geom_point()    +  
  scale_y_continuous(breaks=c(0,2,4,6,8,10, 12), limits=(c(0,12))) + ylab("Enrichment factor") + xlab("")

print(violin.enrichment)

median(setforplot$Value_relative[setforplot$Section=="Single-cell segment" & setforplot$Channel=="GCX"])
median(setforplot$Value_relative[setforplot$Section=="Single-cell segment" & setforplot$Channel=="ALCAM"])
median(setforplot$Value_relative[setforplot$Section=="Cell-cell-contact" & setforplot$Channel=="GCX"])
median(setforplot$Value_relative[setforplot$Section=="Cell-cell-contact" & setforplot$Channel=="ALCAM"])


median(setforplot$Value_bgcor[setforplot$Section=="Cell-cell-contact"])


median(setforplot$Value_relative[setforplot$Section=="Single-cell segment"])
median(setforplot$Value_relative[setforplot$Section=="Cell-cell-contact"])


distinct_groups <- setforplot %>% ungroup() %>%  group_by(Channel, Section, CellNR, ImageNR, Subfolder, expnr, Permeabilization_status) %>%  distinct() %>%  
  select(Channel, Section, ImageNR, Subfolder, expnr, Permeabilization_status, CellNR)
print(distinct_groups)

unique_numbers <- setforplot %>%  ungroup() %>%  group_by(expnr,  ImageNR) %>%  summarize(count  = max(row_number()))
print(unique_numbers)



# unpaired statistics
results <- boxplot_stats %>%
  group_by(Azidosugar, Permeabilization, SPAAC_Time) %>%
  summarise(
    wilcox = list(
      wilcox.test(median ~ SPAAC)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(wilcox, function(x) x$p.value),
    statistic = sapply(wilcox, function(x) x$statistic)
  ) %>%
  select(-wilcox)
print(results)

## extract boxplot statistics -----
boxplot_stats <- setforplot %>%
  ungroup() %>%
  group_by(Channel, Section, Permeabilization_status) %>%
  summarise(
    n = n(),
    mean = mean(Value_relative, na.rm = TRUE),
    sd = sd(Value_relative, na.rm = TRUE),
    median = median(Value_relative, na.rm = TRUE),
    Q1 = quantile(Value_relative, 0.25, na.rm = TRUE),
    Q3 = quantile(Value_relative, 0.75, na.rm = TRUE),
    IQR = IQR(Value_relative, na.rm = TRUE),
    lower_whisker = max(min(Value_relative, na.rm = TRUE), Q1 - 1.5 * IQR),
    upper_whisker = min(max(Value_relative, na.rm = TRUE), Q3 + 1.5 * IQR),
    ci_lower = mean - qt(0.975, df = n - 1) * sd / sqrt(n),
    ci_upper = mean + qt(0.975, df = n - 1) * sd / sqrt(n)
  ) %>%
  ungroup()  

print(boxplot_stats)


boxplot_stats <- setforplot %>%
  ungroup() %>%
  group_by(Channel, Section, Permeabilization_status) %>%
  summarise(
    n = n(),
    mean = mean(Value_bgcor, na.rm = TRUE),
    sd = sd(Value_bgcor, na.rm = TRUE),
    median = median(Value_bgcor, na.rm = TRUE),
    Q1 = quantile(Value_bgcor, 0.25, na.rm = TRUE),
    Q3 = quantile(Value_bgcor, 0.75, na.rm = TRUE),
    IQR = IQR(Value_bgcor, na.rm = TRUE),
    lower_whisker = max(min(Value_bgcor, na.rm = TRUE), Q1 - 1.5 * IQR),
    upper_whisker = min(max(Value_bgcor, na.rm = TRUE), Q3 + 1.5 * IQR),
    ci_lower = mean - qt(0.975, df = n - 1) * sd / sqrt(n),
    ci_upper = mean + qt(0.975, df = n - 1) * sd / sqrt(n)
  ) %>%
  ungroup() 

print(boxplot_stats)


## 5b. Plot multipanel image of plot profiles----
dualpanel <- plot_grid(example,lineplot.smoothened, violin, ncol=3, align="vh") + 
  custom_theme
dualpanel <- dualpanel + 
  theme(legend.position="none")
print(dualpanel)

## 6a. Plot Transition zone length vs. Gcx enrichment factor-----

dotplot.subset <- subset(data_optimallinethickness, Section=="Cell-cell-contact") %>%  
  filter(!is.na(Transitionzone_length) & !is.na(Value_relative)) %>% 
  mutate(Transitionzone_length = as.numeric(Transitionzone_length),
         Value_relative = as.numeric(Value_relative))

dotplot.subset$Transitionzone_classification <- as.factor(dotplot.subset$Transitionzone_classification)
dotplot <- custom_lineplot_with_rsquared(data = dotplot.subset, 
                                         regression_type= "linear",
                                         x_var = Transitionzone_length,
                                         y_var = Value_relative,
                                         size=2)

dotplot.stats <- dotplot$rsquared_values
dotplot.graph <- dotplot$plot
dotplot.graph <- dotplot.graph + labs(x="Transition zone length (micron)", y="Gcx enrichment factor")+ ylim(1,2.7) + xlim(0,4.5)
print(dotplot.graph)

avg_transitionlength <- dotplot.subset %>%  ungroup() %>% 
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) 

## 7a. Plot examples of Transition length (long and short) ------

file.to.plot <- subset(file_optimallinethickness, 
                       Lineprofiletype=="Contact" & 
                         expnr==1 & 
                         ImageNR== "4A" &
                         CellNR==1 & 
                         Section!="Background"  | 
                         Lineprofiletype=="Contact" & 
                         expnr==1 & 
                         ImageNR==3 &
                         CellNR==2 & 
                         Section!="Background" )


file.to.plot <- file.to.plot %>% dplyr::rename(maxX = maxX.x)


calcmid <- file.to.plot %>% group_by(number) %>% 
  filter(Section == "Transition") %>%
  dplyr::summarize(minX = min(X_micron_normalized),
            maxX = max(X_micron_normalized))
calcmid <- calcmid %>% 
  mutate(Transition_mid = (maxX + minX)/2) %>% 
  dplyr::select(number, Transition_mid)
file.to.plot <- file.to.plot %>%  select(-Transition_mid, minX, maxX)


file.to.plot <- file.to.plot %>%
  left_join(calcmid, by = "number") 


file.to.plot <- file.to.plot %>% 
  group_by(number) %>% 
  mutate(aligned_X_micron = X_micron_normalized - Transition_mid,
         maxX_aligned = maxX - Transition_mid,
         minX_aligned = minX - Transition_mid)

file.to.plot$consecutive_transitions <- as.factor(file.to.plot$consecutive_transitions)
  

distvssignal <- ggplot(file.to.plot, aes(x=aligned_X_micron, y=Value_relative, color=Section))+
  facet_wrap(.~Transitionzone_classification, ncol=1)

rect_data <- file.to.plot %>%  group_by(number, Section, Transitionzone_classification) %>%  dplyr::summarize(xmin= min(aligned_X_micron),
                                                                        xmax= max(aligned_X_micron),
                                                                        ymin = -Inf,
                                                                        ymax = Inf) %>%  filter(Section=="Transition")

rectangle <- geom_rect( data=rect_data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.15, inherit.aes=FALSE)

distvssignal <- distvssignal+
   rectangle

distvssignal <- distvssignal +
  geom_line(size=0.5, aes(group=number), alpha=1)+
  geom_point(size=1, stroke=0, alpha=1) +
  
  geom_hline(aes(yintercept=single_cell_mean), linetype="dashed", color="green4")+
  geom_hline(aes(yintercept=cell_contact_mean), linetype="dashed", color="red3")+
  xlim(c(-5,5))+
  
  #geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black") +
  
  labs(y="Delta MFI") +
  
  
  labs(x="Distance (micron)")+
  
  
  theme(legend.position="none")

print(distvssignal)


s <- file.to.plot %>%  group_by(number, Transitionzone_classification) %>%
  dplyr::summarize(across(where(is.numeric), mean, na.rm=TRUE)) 
s$mean_value_relative_all <- rep(mean_value_relative_all, nrow(s))

distvssignal <- distvssignal +
  geom_text(data=s,
            aes(x=0, y=mean_value_relative_all, label=Transitionzone_length),
            vjust=0.5,
            hjust=0.5,
            color="black"
  )

print(distvssignal)
## 8a. Plot Intensity of Cell-cell junction vs. Gcx enrichment -----


custom_lineplot_with_rsquared(data = dotplot.subset, 
                              regression_type= "linear",
                              x_var = Value_bgcor,
                              y_var = Value_relative,
                              size=2)

dotplot.enrichment <- ggplot(dotplot.subset, aes(x=Value_bgcor, y=Value_relative))+
  geom_point()+
  geom_smooth(method="lm", formula="y~x", color="black")+
  labs(x="Cell-cell contact MFI", y= "GCX enrichment")
plot(dotplot.enrichment)

## 8. transition plots together ------
transition.plots <- plot_grid(distvssignal, dotplot.graph)
print(transition.plots)

## 9. plot all transition zones -----


allnumbers <- unique(file_optimallinethickness$number[file_optimallinethickness$Section!="Background"])

allnumbers <- na.omit(allnumbers)


for(p in 1:length(allnumbers)){
  print(p)
  subnumber <- allnumbers[p]
  
file.to.plot <- subset(file_optimallinethickness, number==subnumber & Section!="Background")


file.to.plot <- file.to.plot %>% dplyr::rename(maxX = maxX.x)


calcmid <- file.to.plot %>% group_by(number) %>% 
  filter(Section == "Transition") %>%
  dplyr::summarize(minX = min(X_micron_normalized),
                   maxX = max(X_micron_normalized))
calcmid <- calcmid %>% 
  mutate(Transition_mid = (maxX + minX)/2) %>% 
  dplyr::select(number, Transition_mid)
file.to.plot <- file.to.plot %>%  select(-Transition_mid, minX, maxX)


file.to.plot <- file.to.plot %>%
  left_join(calcmid, by = "number") 


file.to.plot <- file.to.plot %>% 
  group_by(number) %>% 
  mutate(aligned_X_micron = X_micron_normalized - Transition_mid,
         maxX_aligned = maxX - Transition_mid,
         minX_aligned = minX - Transition_mid)

file.to.plot$consecutive_transitions <- as.factor(file.to.plot$consecutive_transitions)


distvssignal <- ggplot(file.to.plot, aes(x=X_micron, y=Value_relative, color=Section))+
  facet_wrap(.~Transitionzone_classification, ncol=1)

rect_data <- file.to.plot %>%  group_by(number, Section, Transitionzone_classification) %>%  dplyr::summarize(xmin= min(aligned_X_micron),
                                                                                                              xmax= max(aligned_X_micron),
                                                                                                              ymin = -Inf,
                                                                                                              ymax = Inf) %>%  filter(Section=="Transition")

rectangle <- geom_rect( data=rect_data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.15, inherit.aes=FALSE)

distvssignal <- distvssignal+
  rectangle

distvssignal <- distvssignal +
  geom_line(size=1, aes(group=number), alpha=1)+
  geom_point(size=1.5, stroke=0, alpha=1) +
  
  geom_hline(aes(yintercept=single_cell_mean), linetype="dashed", color="green4", size=1)+
  geom_hline(aes(yintercept=cell_contact_mean), linetype="dashed", color="red3", size=1)+
  #geom_smooth(method="lm", formula = y~log(x), fullrange=TRUE, se=FALSE, color="black") +
  
  labs(y="Delta MFI") +
  
  
  labs(x="Distance (micron)")+
  
  
  theme(legend.position="top")

print(distvssignal)

pdf(paste0("Lineplot_","Expnr_", unique(file.to.plot$expnr), "- Cell_", unique(file.to.plot$CellNR), "- Roi_ ", unique(file.to.plot$ROINR), "- ID_", subnumber, ".pdf"))
print(distvssignal)
dev.off()

}

## save files for superviolin ----

data_optimallinethickness <- data_optimallinethickness %>% filter(Section!="Background")

superviolin <- data_optimallinethickness %>% ungroup() %>%  dplyr::select(Value_bgcor, expnr, Section) %>%  filter(Section!="Background") 
write.csv(x=superviolin, "superviolin.csv")

summary_perexp <- data_optimallinethickness %>% ungroup() %>% group_by(expnr, Section) %>%  summarise(Value_realtive = mean(Value_bgcor)) %>% filter(Section!="Background")

shapiro.test(data_optimallinethickness$Value_bgcor[data_optimallinethickness$Section=="Single-cell segment"])
shapiro.test(file_optimallinethickness$Value_relative[file_optimallinethickness$Section=="Cell-cell-contact"])

wilcox.test(file_optimallinethickness$Value_relative[file_optimallinethickness$Section=="Cell-cell-contact"], mu = 2)

median(data_optimallinethickness$Value_relative[data_optimallinethickness$Section=="Cell-cell-contact"])

wilcox.test(file_optimallinethickness$Value_relative[file_optimallinethickness$Section=="Single-cell segment"], 
            file_optimallinethickness$Value_relative[file_optimallinethickness$Section=="Cell-cell-contact"], 
            paired = F)

kruskal.test(file_optimallinethickness$Value_relative[file_optimallinethickness$Section=="Single-cell segment"],
             file_optimallinethickness$Value_relative[file_optimallinethickness$Section=="Cell-cell-contact"])


# 1. Perform ANOVA
anova_model <- aov(Value_bgcor ~ Section, data = data_optimallinethickness)
summary(anova_model)  # Check if the p-value is significant

# Post hoc pairwise t-tests after significant ANOVA
pairwise_result <- pairwise.t.test(data_optimallinethickness$Value_bgcor, data_optimallinethickness$Section, p.adjust.method = "holm")
print(pairwise_result)


summary_perexp <- data_optimallinethickness %>% ungroup() %>% group_by(expnr, Section) %>%  summarise(Value_bgcor = mean(Value_bgcor)) %>% filter(Section!="Background")



kruskal.test(Value_bgcor ~ Section, data = summary_perexp)
library(FSA)
dunn.results.perexp <- dunnTest(Value_bgcor ~ Section, data = data_optimallinethickness, method = "sidak")

print(dunn.results.perexp$res)


summary_perexp$expnr <- as.factor(summary_perexp$expnr)
summary_perexp$Section <- as.factor(summary_perexp$Section)

# Run Friedman test
friedman_result <- friedman.test(Value_bgcor ~ Section | expnr, data = summary_perexp)
print(friedman_result)

pairwise_results <- pairwise.wilcox.test(summary_perexp$Value_bgcor, 
                                         summary_perexp$Section, 
                                         p.adjust.method = "bonferroni", 
                                         paired = TRUE)
print(pairwise_results)


## dotplot samples intensity ALCAM vs. GCX ------

data_optimallinethickness <- data_optimallinethickness %>%  mutate(Section= as.factor(Section),
                                                                   Section = factor(Section, levels=c("Single-cell segment", "Transition", "Cell-cell-contact")),
                                                                   Permeabilization_status = as.factor(Permeabilization_status),
                                                                   Permeabilization_status = factor(Permeabilization_status, levels=c("Unpermeabilized", "Saponin"))
                                                                   )

ggplot(subset(data_optimallinethickness, Section!="Background" & Channel!="GCX"), aes(x=Section, y=Value_bgcor))+
  geom_point()+
  #geom_violin(scale="width")+
  #geom_boxplot(width=0.5, outlier.shape=NA)+
  geom_line(aes(group=interaction(ImageNR, expnr, CellNR, Permeabilization_status, number), color=Permeabilization_status))+
  facet_wrap(Permeabilization_status~expnr+ImageNR)