rm(list = ls()); cat("\014");

manualLinearRegression <- function(formula, data) {
  # Prepare the formula and data
  formula <- as.formula(formula)
  data <- as.data.frame(data)
  
  # Construct the design matrix
  X <- model.matrix(formula, data)
  Y <- model.response(model.frame(formula, data))
  
  # Solve for coefficients (betas)
  betas <- solve(t(X) %*% X) %*% t(X) %*% Y
  
  # Calculate fitted values and residuals
  fitted <- X %*% betas
  residuals <- Y - fitted
  
  # Calculate sum of squares
  SSR <- sum((fitted - mean(Y))^2)
  SSE <- sum(residuals^2)
  SST <- sum((Y - mean(Y))^2)
  
  # Calculate degrees of freedom
  df_total <- length(Y) - 1
  df_model <- length(betas) - 1
  df_residual <- df_total - df_model
  
  # Calculate mean square for regression and residual
  MSR <- SSR / df_model
  MSE <- SSE / df_residual
  
  # Calculate F statistic and p-value
  F_value <- MSR / MSE
  F_p_value <- pf(F_value, df_model, df_residual, lower = FALSE)
  
  # Calculate R-squared and Adjusted R-squared
  R_squared <- SSR / SST
  R_squared_adj <- 1 - (1 - R_squared) * (df_total / df_residual)
  
  # Calculate standard errors, t-values, and p-values for coefficients
  std_errors <- sqrt(diag(MSE * solve(t(X) %*% X)))
  t_values <- betas / std_errors
  p_values <- 2 * pt(-abs(t_values), df_residual)
  
  # Format the results into a list
  results <- list(
    coefficients = data.frame(Estimate = betas, Std.Error = std_errors, t.value = t_values, P.Value = p_values),
    anova = data.frame(SumSq = c(SSR, SSE), Df = c(df_model, df_residual), MeanSq = c(MSR, MSE), F.value = F_value, Pr_F = F_p_value),
    r.squared = R_squared,
    adj.r.squared = R_squared_adj,
    residuals = residuals,
    fitted.values = fitted
  )
  
  # Return the results
  return(results)
}


# 
# # Example 1
# x <- c(151, 174, 138, 186, 128, 136, 179, 163, 152, 131)
# y <- c(63, 81, 56, 91, 47, 57, 76, 72, 62, 48)
# manualLinearRegression(y~x, data=data.frame(x,y))
# 
# # Example 2
# dat <- read.csv("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-September-Bioinformatics-Prerequisites/master/friday/lm_example_data.csv")
# head(dat)
# oneway.model <- manualLinearRegression(expression ~ treatment, data = dat)
# oneway.model
# 
# # Example 3
# # Simulate some data
# set.seed(123) # for reproducibility
# x <- rnorm(100)
# y <- 5 + 3 * x + rnorm(100)
# 
# # Create a data frame
# sim_data <- data.frame(x = x, y = y)
# 
# # Fit the model manually
# fit_sim <- manualLinearRegression(y ~ x, sim_data)
# print(fit_sim)



