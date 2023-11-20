# BIOS625
This R package provides functions for manually fitting linear regression models. It's designed for educational purposes to demonstrate the underlying calculations involved in fitting a linear regression model. The main function in the package, `manualLinearRegression`, takes a formula and data frame as inputs and returns a list of regression outputs including coefficients, residuals, fitted values, an ANOVA table, and R-squared values.


## Overview
The primary function, `manualLinearRegression`, allows users to fit linear models by manually performing calculations typically handled by more abstracted functions like `lm()`. This function is particularly useful for educational purposes and for those who wish to understand the nitty-gritty of linear regression analysis.

## Installation

To install the latest version from GitHub, use:

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install BIOS625 from GitHub
devtools::install_github("Jiayuny54/BIOS625")

# Load BIOS625
library(BIOS625)
```



## Usage
The usage of BIOS625 revolves around the manualLinearRegression function, which performs linear regression analysis manually.

## Syntax:
BIOS625::manualLinearRegression(formula, data)
**formula: A symbolic description of the model to be fitted.**
  **data: The data frame containing the variables specified in formula.**

  ## Example 1: Linear Regression on example dataset
  ```r
library(BIOS625)
dat <- read.csv("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-September-Bioinformatics-Prerequisites/master/friday/lm_example_data.csv")
head(dat)
# sample expression  batch treatment  time temperature
# 1      1     1.2140 Batch1         A time1       11.77
# 2      2     1.4797 Batch1         A time2       12.16
# 3      3     1.0878 Batch1         A time1       10.54
# 4      4     1.4439 Batch1         A time2       10.08
# 5      5     0.6372 Batch1         A time1       12.04
# 6      6     2.1227 Batch1         B time2       13.50
oneway.model <- BIOS625::manualLinearRegression(expression ~ treatment, data = dat)
oneway.model
# $coefficients
# Estimate Std.Error t.value   P.Value
# (Intercept) -38.4551   8.04901  -4.778 1.395e-03
# x             0.6746   0.05191  12.997 1.164e-06
#
# $anova
# SumSq Df  MeanSq F.value      Pr_F
# 1 1787.45  1 1787.45   168.9 1.164e-06
# 2   84.65  8   10.58   168.9 1.164e-06
#
# $r.squared
# [1] 0.9548
#
# $adj.r.squared
# [1] 0.9491
#
# $residuals
# [,1]
# 1  -0.4111
# 2   2.0729
# 3   1.3588
# 4   3.9775
# 5  -0.8951
# 6   3.7081
# 7  -6.3002
# 8   0.4936
# 9  -2.0857
# 10 -1.9189
#
# $fitted.values
# [,1]
# 1  63.41
# 2  78.93
# 3  54.64
# 4  87.02
# 5  47.90
# 6  53.29
# 7  82.30
# 8  71.51
# 9  64.09
# 10 49.92
```



## Example 2: Linear Regression on a mocked Dataset
``` r
# Simulate some data
set.seed(123) # for reproducibility
x <- rnorm(100)
y <- 5 + 3 * x + rnorm(100)

# Create a data frame
sim_data <- data.frame(x = x, y = y)

# Fit the model manually
fit_sim <- BIOS625::manualLinearRegression(y ~ x, sim_data)
print(fit_sim)
# $coefficients
# Estimate Std.Error t.value   P.Value
# (Intercept)    4.897   0.09755   50.20 1.003e-71
# x              2.948   0.10688   27.58 5.575e-48
#
# $anova
# SumSq Df   MeanSq F.value      Pr_F
# 1 716.67  1 716.6673   760.6 5.575e-48
# 2  92.34 98   0.9423   760.6 5.575e-48
#
# $r.squared
# [1] 0.8859
#
# $adj.r.squared
# [1] 0.8847
#
# $residuals
# [,1]
# 1   -0.63701
# 2    0.34761
# 3   -0.06210
# 4   -0.24104
# 5   -0.84203
# 6    0.14777
# 7   -0.65792
# 8   -1.63152
# 9   -0.31346
# 10   0.99842
#
# $fitted.values
# [,1]
# 1    3.24518
# 2    4.21874
# 3    9.49153
# 4    5.10502
# 5    5.27828
# 6    9.95240
# 7    6.25576
# 8    1.16839
# 9    2.87268
# 10   3.58360

```


## Features
+ **Manual Calculation**: Executes all steps in the regression calculation, offering transparency in the model fitting process.
+ **Detailed Output**: Provides a comprehensive output, including coefficients, ANOVA table, and R-squared statistics.
+ **Educational Tool**: Serves as an educational resource for understanding the mathematical foundations of linear regression.
Contributing
+ **Your contributions are welcome! Please feel free to submit pull requests or create issues for bugs and feature requests on the GitHub repository.**

## License
+ BIOS625 is open-source software licensed under the MIT license.
