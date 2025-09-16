#----packages------
library(ggplot2)
library(gridExtra)
library(tidyr)
library(matrixStats)
library(plyr)
library(dplyr)
library(ggpubr)
library(Rmisc)
library(svDialogs)
library(ggsignif)
library(stats)
library(stringr)
library(remotes)
library(cowplot)
library("report")
library(multcomp)
library(tseries)
library(ggplot2)
library(scales)


##---parameters------


workdir <- choose.dir(getwd(), "Choose folder containing subdirectories with .csv files")

dirlist <- list.dirs(path=workdir, recursive=FALSE)
ndirs <- length(dirlist)

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


highres=TRUE

if(highres==TRUE){
pixelsize=0.042
pxsize=0.042}else{
  pixelsize=0.72
  pxsize=0.72
}

## define theme for plotting------

# Define the custom theme with the "sans" font, black text color, and smaller legend
custom_theme <- theme(
  text = element_text(family = "sans", size = 7, color = "black"),
  legend.key.size = unit(1, "lines"),  # Adjust the size of the legend key
  legend.key.height = unit(0.5, "lines"),  # Adjust the height of the legend key
  # Customize the appearance of geom_text
  axis.text = element_text(size = 6),  # Customize the axis text size
  axis.title = element_text(size = 6),  # Customize the axis title text size
  # Add custom settings for geom_text
  
  
)

theme_set(custom_theme)
##---functions ------


calculate_residuals_and_tests <- function(x, y, filename) {
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
  
  # Save the results to a file
  write.csv(results_df, file = filename, row.names = FALSE)
  
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
  
  # New relationship (y ~ (x/sqrt(x)))
  new_formula <- "y ~ (x/sqrt(x))"
  new_model <- lm(new_formula, data = data.frame(x = x, y = y))
  new_residuals <- residuals(new_model)
  new_jarque <- jarque.bera.test(new_residuals)
  new_sw <- shapiro.test(new_residuals)
  new_r_squared <- summary(new_model)$r.squared
  new_adj_r_squared <- summary(new_model)$adj.r.squared
  new_aic <- AIC(new_model)
  new_bic <- BIC(new_model)
  new_coeffs <- coef(new_model)
  new_intercept <- new_coeffs[1]
  new_slope <- new_coeffs[2]
  new_intercept_p_value <- summary(new_model)$coefficients[1, "Pr(>|t|)"]
  new_slope_p_value <- summary(new_model)$coefficients[2, "Pr(>|t|)"]
  new_intercept_sig <- ifelse(new_intercept_p_value < 0.05, "*", "n.s")
  new_slope_sig <- ifelse(new_slope_p_value < 0.05, "*", "n.s")
  
  # Create a data frame to store the results
  results_df <- data.frame(
    Formula = c(linear_formula, exponential_formula, logarithmic_formula, new_formula),
    JB_p_value = c(linear_jarque$p.value, ifelse(is.null(exponential_jarque), NA, exponential_jarque$p.value), logarithmic_jarque$p.value, new_jarque$p.value),
    SW_p_value = c(linear_sw$p.value, ifelse(is.null(exponential_sw), NA, exponential_sw$p.value), logarithmic_sw$p.value, new_sw$p.value),
    R_squared = c(linear_r_squared, ifelse(is.null(exponential_r_squared), NA, exponential_r_squared), logarithmic_r_squared, new_r_squared),
    Adj_R_squared = c(linear_adj_r_squared, ifelse(is.null(exponential_adj_r_squared), NA, exponential_adj_r_squared), logarithmic_adj_r_squared, new_adj_r_squared),
    AIC = c(linear_aic, ifelse(is.null(exponential_aic), NA, exponential_aic), logarithmic_aic, new_aic),
    BIC = c(linear_bic, ifelse(is.null(exponential_bic), NA, exponential_bic), logarithmic_bic, new_bic),
    Coefficient = c(linear_slope, ifelse(is.null(exponential_slope), NA, exponential_slope), logarithmic_slope, new_slope),
    Intercept = c(linear_intercept, ifelse(is.null(exponential_intercept), NA, exponential_intercept), logarithmic_intercept, new_intercept),
    Coefficient_Significance = c(linear_slope_sig, ifelse(is.null(exponential_slope_sig), NA, exponential_slope_sig), logarithmic_slope_sig, new_slope_sig),
    Intercept_Significance = c(linear_intercept_sig, ifelse(is.null(exponential_intercept_sig), NA, exponential_intercept_sig), logarithmic_intercept_sig, new_intercept_sig),
    Coefficient_p_value = c(linear_slope_p_value, ifelse(is.null(exponential_slope_p_value), NA, exponential_slope_p_value), logarithmic_slope_p_value, new_slope_p_value),
    Intercept_p_value = c(linear_intercept_p_value, ifelse(is.null(exponential_intercept_p_value), NA, exponential_intercept_p_value), logarithmic_intercept_p_value, new_intercept_p_value),
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
  
  if (any(x > log(.Machine$integer.max))) {
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


##---import existing file------
existingfilename <- basename(list.files(workdir, pattern="SNR_summarised.csv"))
if(length(existingfilename==1)){
existingfile <- read.csv2(paste0(workdir,"/",existingfilename), header=T, sep = ",")
analyzedrois <- existingfile$name
}



##---import data------

for(d in 1:ndirs){
  dir <- dirlist[d]
  setwd(dir)
  Filebasename <- basename(list.files(dir, pattern="SNR.csv"))
  for(n in 1:length(Filebasename)){
    filename <- Filebasename[n]
    file <- read.csv2(paste0(dir,"/",filename), header=T, sep = ",", quote="")
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

rm(file)
rm(file2)
file <- separate(file1, (ncol(file1)), into=c("Treatment", "CellNR", "Acquisitiontype", "RoiNR", "Parameter"),sep="_", remove = FALSE)

file$Value <- as.numeric(file$Value)
file$Disctance_micron <- file$X * pixelsize

##---determine thresholds of signal and noise------

uniquenames <- unique(file$name)

if(length(existingfilename==1)){
  newnames <- setdiff(uniquenames,analyzedrois)
}else{
  newnames <- uniquenames
}

if(length(newnames!=0)){
  for(x in 1:length(newnames)){
sub <- subset(file, name==newnames[x])


p <- ggplot(sub, aes(x=X, y=Value))+
        geom_line()+
  ggtitle(newnames[x])

print(p)


  mx <- max(sub$Value)
  Xmax <- sub$X[sub$Value==mx]
  if(highres==TRUE){
    leftthr <- Xmax-40
    rightthr <- Xmax+40
  }else{
    Xmax <- dlgInput("Enter a number:", Xmax)$res
    Xmax <- as.numeric(Xmax)
    leftthr <- Xmax-2
    rightthr <- Xmax+2
  }
  signal <- subset(sub, X>leftthr & X<rightthr)
  noise_extracellular <- subset(sub, X<leftthr)
  noise_intracellular <- subset(sub, X>rightthr)

  snr <- sub %>%  group_by(name) %>%  summarise_at(c("Value"), list(Mean_overall=mean)) %>%  as.data.frame()
  snr$leftthr <- leftthr
  snr$rightthr <- rightthr
  snr$Max_signal <- max(signal$Value)
  snr$Mean_signal <- mean(signal$Value)
  snr$SD_signal <- sd(signal$Value)
  snr$Max_noise_extracellular <- max(noise_extracellular$Value)
  snr$Mean_noise_extracellular <- mean(noise_extracellular$Value)
  snr$SD_noise_extracellular <- sd(noise_extracellular$Value)
  
  snr$Max_noise_intracellular <- max(noise_intracellular$Value)
  snr$Mean_noise_intracellular <- mean(noise_intracellular$Value)
  snr$SD_noise_intracellular <- sd(noise_intracellular$Value)
  
  
  
  plot_thr <- print(p)+
    geom_vline(xintercept=leftthr, linetype=2)+
    geom_vline(xintercept=rightthr, linetype=2)+
    ylab("Pixel Intensity")+
    xlab("Distance (pixels)")
  
  print(plot_thr)
  
  
  ggsave(paste0("Lineplot_threshold_", uniquenames[x],".png"), device="png", path=workdir, width = 5, height = 2,dpi=300)
  
  
  if(x==1){
    snr1 <- snr
  }else{
    snr2 <- snr
    snr1 <- rbind(snr1,snr2)
  }
}


snr1 <- separate(snr1, 1, into=c("Treatment", "CellNR", "Acquisitiontype", "RoiNR", "Parameter"),sep="_", remove = FALSE)

if(highres==TRUE){
snr1$snr_extracellular <- (snr1$Mean_signal-snr1$Mean_noise_extracellular)/(snr1$SD_noise_extracellular)
snr1$snr_intracellular <- (snr1$Mean_signal-snr1$Mean_noise_intracellular)/(snr1$SD_noise_intracellular)
}else{
  
  snr1$snr_extracellular <- (snr1$Max_signal-snr1$Mean_noise_extracellular)/(snr1$SD_noise_extracellular)
  snr1$snr_intracellular <- (snr1$Max_signal-snr1$Mean_noise_intracellular)/(snr1$SD_noise_intracellular)
  
}
}

##---combine previous with new analyses------

if(length(existingfilename==1)){
  if(exists("snr1")){
  existingfile <- existingfile[,2:ncol(existingfile)]
  snr1 <- rbind(existingfile,snr1)
  }else{
    snr1 <- existingfile
  }
}


##---convert all values to numeric of imported file ------
snr1$snr_intracellular <- as.numeric(snr1$snr_intracellular)
snr1$Mean_signal <- as.numeric(snr1$Mean_signal)
snr1$Max_signal <- as.numeric(snr1$Max_signal)

## -- other adjustments ------

file <- file %>% mutate(number= group_indices(.,name, CellNR))

file <- file %>% group_by(number) %>%  mutate(X_micron = X*pxsize)
file <- file %>% group_by(number) %>%  mutate(Xmax = X[Value==max(Value)])

if(highres==TRUE){
file <- file %>%  group_by(number) %>%  mutate(leftthr = Xmax-40)
file <- file %>%  group_by(number) %>%  mutate(rightthr = Xmax+40)

file <- file %>% group_by(number) %>%  mutate(Xmax_micron = X_micron[Value==max(Value)])
file <- file %>%  group_by(number) %>%  mutate(leftthr_micron = Xmax_micron-(40*pxsize))
file <- file %>%  group_by(number) %>%  mutate(rightthr_micron = Xmax_micron+(40*pxsize))
}else{
  file <- file %>%  group_by(number) %>%  mutate(leftthr = Xmax-2)
  file <- file %>%  group_by(number) %>%  mutate(rightthr = Xmax+2)
  
  file <- file %>% group_by(number) %>%  mutate(Xmax_micron = X_micron[Value==max(Value)])
  file <- file %>%  group_by(number) %>%  mutate(leftthr_micron = Xmax_micron-(2*pxsize))
  file <- file %>%  group_by(number) %>%  mutate(rightthr_micron = Xmax_micron+(2*pxsize))
  
}
file <- file %>%
  group_by(number) %>%
  mutate(Section = ifelse(X >= leftthr & X < rightthr, "Signal", "Noise"))
file <- file %>%
  group_by(number) %>%
  mutate(Section = ifelse(X < leftthr, NA, Section))


##---summarize data per cell------


snr1_summarized <- snr1 %>% group_by(Treatment,CellNR) %>%  summarise_if(is.numeric, mean, na.rm = TRUE) %>%  as.data.frame()
##------statistics------


snr.model.statistics <- calculate_residuals_and_tests_limitcorrected(snr1$Mean_signal, snr1$snr_intracellular)
if(highres==TRUE){
  write.csv2(x=snr.model.statistics, "SNR.model.statistics.highres.csv")
}else{
  write.csv2(x=snr.model.statistics, "SNR.model.statistics.lowres.csv")
}



pv <- snr.model.statistics$Coefficient_p_value
pv <- pv[1]

##----plots----

#-----lineplot snr vs fluorescence------



threshSNR <- subset(snr1, CellNR=="554" | CellNR=="115" | CellNR =="424")
threshSNR <- mean(threshSNR$snr_intracellular)

snr1$snr_intracellular <- snr1$Mean_signal/as.numeric(snr1$Mean_noise_intracellular)

snr1 <- snr1 %>%  filter(snr_intracellular < 20)

lineplot <- ggplot(snr1, aes(x=Mean_signal, y=snr_intracellular)) +
  geom_point(aes(color=Treatment), alpha=0.5, size=2, stroke=0) +
  geom_smooth(method="lm", formula= (y ~ (sqrt(x))), se=TRUE, linewidth=0.5, color="black") +
  
  xlab("Mean Signal (A.U.)")+
  ylab("Signal-To-Noise Ratio")+
 # scale_x_continuous(trans = "log10", minor_breaks = NULL,labels = trans_format("log10", function(x) sprintf("%.0f", 10^x)),limits = c(10^1.5, 300)) +
  theme(legend.position = c(0.9, 0.85),
        legend.justification = c(1, 0),
        legend.box.just = "right",
        legend.title=element_blank(),
        #strip.background = element_rect(color = "black", fill = NA, size = 0.5),
        #panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.background = element_rect(fill = "transparent"))+
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)


lineplot <- lineplot+ theme_classic()
print(lineplot)


if(highres==FALSE){
  lineplot <- lineplot + geom_hline(yintercept = threshSNR, linetype=2)
}


print(lineplot)

changePlotSize(lineplot)

saveRDS(lineplot, file="Lineplot_SNR.RData")

## -- save relevant data of snr plot ----

summarised_data <- snr1 %>% dplyr::select(Treatment, Mean_signal, snr_intracellular)


##---- percentage-----

snr1$abovethreshold[snr1$snr_intracellular>threshSNR] <- TRUE
snr1$abovethreshold[snr1$snr_intracellular<threshSNR] <- FALSE


snr1_threshold <- snr1 %>%
  group_by(Treatment, abovethreshold) %>%
  dplyr::summarise(count = n(), .groups = "drop")

snr1_sum <- snr1_threshold %>% group_by(Treatment) %>% summarise(sum=sum(count))
snr1_abovethr <- subset(snr1_threshold, abovethreshold==TRUE)
snr1_abovethr$total_count <- snr1_sum$sum
snr1_abovethr$perc_above_thr <- snr1_abovethr$count/snr1_abovethr$total_count*100

percabovethr_SNR <- ggplot(snr1_abovethr, aes(x=Treatment, y= perc_above_thr))+
  geom_bar(stat="identity", color= "black", fill="white")+
  
  xlab("")+
  ylab("% above threshold")+

  
  theme(axis.text.y = element_text(size=8),
        text = element_text(size = 8),
        legend.text = element_text(size=8),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size=8),
        legend.position = "none") +
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)

print(percabovethr_SNR)


snr1_lowres <- plot_grid(lineplot, percabovethr_SNR, rel_widths=c(0.75, 0.25), ncol=2,align="v")

print(snr1_lowres)

changePlotSize(snr1_lowres)

saveRDS(percabovethr_SNR, file="Perc_above_thr_SNR.RData")
  

## quantile -----

snr1 <- snr1 %>%  ungroup() %>% mutate(conf_interval = quantile(snr_intracellular[Treatment=="DBCO"], .95)) 

snr1$SPAAC <- snr1$Treatment
snr1$SNR <- snr1$snr_intracellular
data_measurement <- snr1


data_measurement <- data_measurement  %>%  mutate(High_intensity = ifelse(SNR>conf_interval, TRUE, FALSE))

smry <- data_measurement %>%  group_by(SPAAC) %>% dplyr::summarise(median_SNR = median(SNR))

ggplot(subset(data_measurement),aes(x=SPAAC, y=SNR))+
  geom_violin()+
  geom_point()

highintensity_counts <- data_measurement %>%  group_by(SPAAC, High_intensity) %>%  dplyr::summarise(Count = n()) 

highintensity_counts <- highintensity_counts %>%  group_by(SPAAC) %>% 
  mutate(Percentage_high = Count[High_intensity==TRUE]/ (Count[High_intensity==TRUE] + Count[High_intensity==FALSE])*100) %>%  filter(High_intensity==TRUE)

highintensity_counts_all <- highintensity_counts %>% group_by(SPAAC) %>% summarize(Percentage_high_sd = sd(Percentage_high),
                                                                                   Percentage_high = mean(Percentage_high),
)

summary_specific_abovetrh <- data_measurement  %>% group_by(SPAAC) %>%  summarise(mean(conf_interval))




summary_all <- data_measurement %>% group_by( Azidosugar, SPAAC) %>% filter(!is.infinite(SNR)) %>% summarize(Median_SNR= median(SNR, na.rm=T),
                                                                                                             Mean_SNR = mean(SNR, na.rm=T),
                                                                                                             Mean_CI = mean(conf_interval)
)




t_test_result <- t.test(highintensity_counts$Percentage_high[highintensity_counts$SPAAC=="THS"], alternative = "greater", mu=5)
print(t_test_result)

ggplot(highintensity_counts, aes(x=SPAAC))+
  geom_col(position="dodge", aes(y=Percentage_high), fill="white", color="black")+
  #geom_errorbar(aes(ymin=Percentage_high, ymax=Percentage_high + Percentage_high_sd))+
  theme_classic()+
  labs(x="", y="% High intensity")



## plot -----
##boxplot/violin plot Treatment vs snr

boxplot <- ggplot(snr1, aes(x=Treatment, y=snr_intracellular))+
  geom_violin()+
  geom_jitter(aes(color=CellNR), size=0.5)+
  
  geom_boxplot(width=0.2, outlier.shape=NA)+
  
  
  theme(axis.text.y = element_text(size=8))+
  theme(text = element_text(size = 8))+
  theme(legend.text = element_text(size=8),legend.key.size = unit(2, "mm"))+
  theme(axis.text=element_text(size=8))+
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)


print(boxplot)

ggsave("Boxplot_SNR.pdf", device="pdf", path=workdir, width = 4, height = 2,dpi=300)


saveRDS(lineplot, file="Boxplot_SNR.RData")


##summarized lineplot

summarized_line <- ggplot(snr1_summarized, aes(x=Mean_signal, y=snr_intracellular))+
  geom_smooth(method="lm",se=FALSE)+
  geom_point(aes(color=Treatment))+
  theme(axis.text.y = element_text(size=8))+
  theme(text = element_text(size = 8))+
  theme(legend.text = element_text(size=8),legend.key.size = unit(2, "mm"))+
  theme(axis.text=element_text(size=8))+
  scale_colour_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)+
  scale_fill_viridis_d(option="turbo", direction=-1, begin=0.2, end=0.8)

print(summarized_line)

ggsave("Lineplot_percell_SNR.pdf", device="pdf", path=workdir, width = 3, height = 2,dpi=300)



## --- plot profile ------

if(highres==TRUE){
examplenr <- unique(file$number[file$CellNR==5 & file$RoiNR==8])
}else{
examplenr <- unique(file$number[file$name=="TMTHSI_5_4_SNR"])
}

profile <- subset(file, number==examplenr)

pp <- ggplot(data=profile, aes(x=X_micron, y=Value))+
  geom_path(aes(color=Section, group=number))+
  labs(y="MFI", x="Profile (micron)")+
  theme(legend.position="none")+
  scale_color_manual(values=c("red", "blue", "black"))

  pp <- pp+
    geom_vline(xintercept = unique(profile$leftthr_micron), linetype=2)+
    geom_vline(xintercept = unique(profile$rightthr_micron), linetype=2)


pp <- pp + custom_theme

print(pp)


pp <- pp + xlim(c(0,40))

print(pp)

changePlotSize(pp)

##------remove unwanted intermediate files-----
rm(signal)
rm(snr)
rm(snr2)
rm(noise)
rm(file1)
rm(dirfiles)
rm(dirfiles2)
rm(sub)

##------save csv files-----------------
write.csv(x=snr1, "SNR_summarised.csv")
write.csv(x=file, "Lineplot_data.csv")


