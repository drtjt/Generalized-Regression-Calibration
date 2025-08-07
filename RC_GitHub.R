# ------------------------------------------------------------------------------
# Estimating Regression Coefficients in a Cox Regression Model, 
#   Where the Effect of Interest has a Non-Linear Relationship
#   and Measured with Error
#
# Date Created: July 30, 2025
#
# Author: Trevor Thomson
#         tthomson@fredhutch.org
# ------------------------------------------------------------------------------


rm(list =ls())
set.seed(5347547)
setwd("C:/Users/tthomson/Dropbox/Fred Hutchinson Cancer Centre_Work/R Code/Women's Health Initiative") # REMOVE

library(MASS)
library(pbapply)
library(nleqslv)
library(Rcpp)
library(RcppArmadillo)
library(stringr)
library(numDeriv) 
library(mvtnorm)
library(cubature) 
library(survival)
library(splines2) 

source("RC_Functions_GitHub.R")

# Specify sample sizes
N1 = 10000
N2 = 500
N3 = 0.2*N2

# Specify coefficients for the measurement error model
a_vec = c(0, 0.8, -0.8) # intercept, effect of Z, effect of V

# Specify the variablity not explained by the (true) biomarker and covariate
sigma2_epsilon = 1

# Specify the measurement error variance for the observed biomarker
sigma2_e = 0.25

# Specify regression parameters in a Cox regression model
# SETTINGS:
#   1: Covariates = (Z, V)
#   2: Covariates = (Z, Z^2, V)
#   3: Covariates = (Z, Z^2, Z^3, V)
#   4: Covariates = (Some Nonlinear function of Z, V)
Covariate_Setting = 1
Beta_Setting = 2 # COMMENT
#Beta_Setting = as.integer(commandArgs(trailingOnly=T)[4]) # UNCOMMENT
if (Covariate_Setting == 1){
  Beta = c(-0.2, 0.1)
  DF_True = 1
}
if (Covariate_Setting == 2){
  Beta = c(-0.2, 0.05, 0.1)
  DF_True = 2
}
if (Covariate_Setting == 3){
  Beta = c(-0.2, 0.05, 0.02, 0.1)
  DF_True = 3
}
if (Covariate_Setting == 4){
  Beta = c(1,0.1) # Note: The "1" is just a placeholder
  DF_True = 1 # Arbitrarily set to 1
}

# Specify the baseline hazard function parameter
lambda0 = 0.02

# Specify the censoring rate
lambda_Cens = 0.2

# Specify the covariance matrix for the covariates
rho = 0.2 # Correlation
Sigma_ZX = matrix(c(1,rho,
                    rho,1), nrow = 2, ncol = 2)

# Obtain "regression parameters" for nuisance parameters
S2 = 1 - rho^2
alpha1 = (a_vec[2] * S2) / (a_vec[2]^2 * S2 + sigma2_epsilon)
alpha0 = - alpha1 * a_vec[1]
alpha2_vec = rho - alpha1 * a_vec[2] * rho - alpha1 * a_vec[3:length(a_vec)]
sigma2 = (S2 * sigma2_epsilon) / (a_vec[2]^2 * S2 + sigma2_epsilon)

# Specify a grid of time points to estimate the survival function
t_grid = seq(0,20, length.out = 101)

# Specify the z-values for estimating the conditional survival function
z_grid = seq(-3, 3, length.out = 11)

seed_initial = 5347 
set.seed(seed_initial)

# -----------------------------------------
#   Analysis Cohort
# -----------------------------------------

# Generate covariates
ZX_1_mat = MASS::mvrnorm(n = N1, 
                         mu = c(0,0),
                         Sigma = Sigma_ZX)

# Generate the rates for the event time
rates = Rates_func(Z_vec = ZX_1_mat[,1],
                   X_mat = as.matrix(ZX_1_mat[,2]),
                   lambda_0 = lambda0,
                   Beta_vec = Beta[-length(Beta)],
                   gamma_vec = Beta[length(Beta)],
                   Covariate_Setting = Covariate_Setting)

# Generate the event time
T_vec = T_HP_vec.m1 = sapply(1:N1,
                             function(i) rexp(1, 
                                              rate = rates[i]))

# Generate the censoring times
CensoringTime_vec = rexp(N1, rate = lambda_Cens)

# Generate right censored event times 
U_vec = pmin(T_vec, CensoringTime_vec) 
delta_vec = I(T_vec <= CensoringTime_vec)*1

# Obtain noisy measurement of biomarker that follows a measurement error model
#   Q_{ij} = a_0 + a_1 Z_i + a_2' X + epsilon_{ij}
#   with j = 1 for all i
epsilon_vec = rnorm(N1, mean = 0, sd = sqrt(sigma2_epsilon))
Q_1_vec = sapply(1:N1,
                 function(i){
                   Q = a_vec[1] + 
                     a_vec[2] * ZX_1_mat[i,1] + 
                     as.numeric(as.matrix(ZX_1_mat[i,2:dim(ZX_1_mat)[2]]) %*% a_vec[3:length(a_vec)]) +
                     epsilon_vec[i]
                   return(Q)
                 })

# Create a dataframe for the analysis cohort
df_1 = data.frame(ID = 1:N1,
                  Time = T_vec,
                  U = U_vec,
                  delta = delta_vec,
                  X = ZX_1_mat[,2:dim(ZX_1_mat)[2]],
                  Z1 = ZX_1_mat[,1],
                  Z2 = ZX_1_mat[,1]^2,
                  Z3 = ZX_1_mat[,1]^3,
                  Z4 = ZX_1_mat[,1]^4,
                  Q1 = Q_1_vec,
                  Q2 = Q_1_vec^2,
                  Q3 = Q_1_vec^3,
                  Q4 = Q_1_vec^4)
df_1_OrderedZ = df_1[order(df_1$Z1),]
df_1_OrderedQ = df_1[order(df_1$Q1),]
rm(df_1)


# -----------------------------------------
#   Biomarker Cohort
# -----------------------------------------

# Generate covariates
ZX_2_mat = MASS::mvrnorm(n = N2, 
                         mu = c(0,0),
                         Sigma = Sigma_ZX)

# Generate measurement errors for each unit
epsilon_mat = MASS::mvrnorm(n = N2,
                            mu = c(0, 0),
                            Sigma = diag(2)*sigma2_epsilon)

e_mat = MASS::mvrnorm(n = N2,
                      mu = c(0, 0),
                      Sigma = diag(2)*sigma2_e)

# Generate Q_{i1} and Q_{i2} for each unit
Q_mat = sapply(1:N2,
               function(i){
                 mu = a_vec[1] + 
                   a_vec[2] * ZX_2_mat[i,1] + 
                   as.numeric(as.matrix(ZX_2_mat[i,2:dim(ZX_2_mat)[2]]) %*% a_vec[3:length(a_vec)]) 
                 Q1 = mu + epsilon_mat[i,1]
                 Q2 = mu + epsilon_mat[i,2]
                 val = c(Q1, Q2)
                 return(val)
               })
Q_mat = t(Q_mat)

# Generate W_{i1} and W_{i2} for each unit
W_mat = sapply(1:N2,
               function(i){
                 W1 = ZX_2_mat[i,1] + e_mat[i,1]
                 W2 = ZX_2_mat[i,1] + e_mat[i,2]
                 val = c(W1, W2)
                 return(val)
               })
W_mat = t(W_mat)

# Generate a dataset for the biomarker cohort
df_2 = data.frame(ID = 1:N2,
                  W1 = W_mat[,1],
                  Q1 = Q_mat[,1],
                  Z = ZX_2_mat[,1],
                  X = ZX_2_mat[,2:dim(ZX_2_mat)[2]])
df_2$ID = df_2$ID + N1

# Generate a dataset for the "replicated biomarker cohort"
df_3 = data.frame(ID = 1:N3,
                  W1 = W_mat[1:N3, 1],
                  W2 = W_mat[1:N3, 2])
df_3$ID = df_3$ID + N1


# ----------------------------------------------------------------------
#         LINEAR REGRESSION CALIBRATION
# ----------------------------------------------------------------------


Results_Polynomial_LRC = EstimatingFunction_Wrapper(df1 = df_1_OrderedQ,
                                                    df2 = df_2,
                                                    df3 = df_3,
                                                    df1_ID_index = 1,
                                                    df1_U_index = 3,
                                                    df1_delta_index = 4,
                                                    df1_Q_index = 10,
                                                    df1_X_index = 5,
                                                    df2_ID_index = 1,
                                                    df2_Response_index = 2,
                                                    df2_Covariates_index = c(3,5),
                                                    df3_ID_index = 1,
                                                    Function_Type = "Polynomial",
                                                    DF = 1,
                                                    Method = "Linear",
                                                    EstimateNuisance_Logi = T,
                                                    mu_vec = NULL,
                                                    Variance_vec = NULL,
                                                    t_grid = t_grid,
                                                    z_grid = z_grid,
                                                    C = 1000,
                                                    seed = 235346,
                                                    x = 0,
                                                    Spline_Z_vec = NULL,
                                                    Boundary.knots = NULL,
                                                    Estimate_CI_Surv = F)

# ----------------------------------------------------------------------
#         QUADRATIC REGRESSION CALIBRATION
# ----------------------------------------------------------------------


Results_Polynomial_QRC = EstimatingFunction_Wrapper(df1 = df_1_OrderedQ,
                                                    df2 = df_2,
                                                    df3 = df_3,
                                                    df1_ID_index = 1,
                                                    df1_U_index = 3,
                                                    df1_delta_index = 4,
                                                    df1_Q_index = 10,
                                                    df1_X_index = 5,
                                                    df2_ID_index = 1,
                                                    df2_Response_index = 2,
                                                    df2_Covariates_index = c(3,5),
                                                    df3_ID_index = 1,
                                                    Function_Type = "Polynomial",
                                                    DF = 2,
                                                    Method = "Quadratic",
                                                    EstimateNuisance_Logi = T,
                                                    mu_vec = NULL,
                                                    Variance_vec = NULL,
                                                    t_grid = t_grid,
                                                    z_grid = z_grid,
                                                    C = 1000,
                                                    seed = 235346,
                                                    x = 0,
                                                    Spline_Z_vec = NULL,
                                                    Boundary.knots = NULL,
                                                    Estimate_CI_Surv = F)


# ----------------------------------------------------------------------
#         CUBIC REGRESSION CALIBRATION
# ----------------------------------------------------------------------


Results_Polynomial_CRC = EstimatingFunction_Wrapper(df1 = df_1_OrderedQ,
                                                    df2 = df_2,
                                                    df3 = df_3,
                                                    df1_ID_index = 1,
                                                    df1_U_index = 3,
                                                    df1_delta_index = 4,
                                                    df1_Q_index = 10,
                                                    df1_X_index = 5,
                                                    df2_ID_index = 1,
                                                    df2_Response_index = 2,
                                                    df2_Covariates_index = c(3,5),
                                                    df3_ID_index = 1,
                                                    Function_Type = "Polynomial",
                                                    DF = 3,
                                                    Method = "General_MC",
                                                    EstimateNuisance_Logi = T,
                                                    mu_vec = NULL,
                                                    Variance_vec = NULL,
                                                    t_grid = t_grid,
                                                    z_grid = z_grid,
                                                    C = 1000,
                                                    seed = 235346,
                                                    x = 0,
                                                    Spline_Z_vec = NULL,
                                                    Boundary.knots = NULL,
                                                    Estimate_CI_Surv = F) # Takes a few minutes
