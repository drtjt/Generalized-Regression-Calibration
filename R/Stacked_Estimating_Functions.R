# Write a function to define the (stacked) estimating function
#   Gamma(Omega) = sum Gamma_i(Omega)


# -- INPUT:
# Beta_vec: Parameters for Effect of Interest (Vector)
# gamma_vec: Parameters for Additional Covariates (Vector)
# Nuisance_Parameters_vec: Parameters for nuisance parameters estimated with additional data (Vector)
# U_vec: Observed Event Times (Vector)
# delta_vec: Observed Event Time Indicators (Vector)
# Q_vec: Covariate for Effect of Interest (Vector)
# X_mat: Additional Covariates (Vector or matrix)
# Function_Type: Functional form for g(Z; beta) (Character)
#                 "Polynomial": Polynomial
#                 "Bernstein": Bernstein Polynomial
#                 "B-Spline": B-Spline
#                 "C-Spline": C-Spline
#                 "I-Spline": I-Spline
#                 "M-Spline": M-Spline
#                 "NC-Spline": Natural Cubic Spline
# DF: Degree for the polynomial or degrees of freedom for the B-Spline (Scalar)
# Method: The method we want to employ (Character)
#           one of {"Cox", "Linear", "Quadratic", "General_1", "General_", "General_MC"}
# C = 1000: Number of Monte Carlo samples
# seed = 235346: seed for "General_MC" (Scalar)
# Nuisance_Regression_EF_func: Nuisance parameter estimating function defined in Function (9) (Function)
# Nuisance_sigma2e_EF_func: Measurement error variance estimating function defined in Function (9) (Function)
# df_1: Dataset to estimate hazard function: (Dataframe)
# df_2: Dataset to estimate nuisance parameters: (Dataframe)
# df_3: Dataset to estimated measurement error: (Dataframe)
#   NOTE: df3 needs to be in "wide format" (all observations for a subject is in a single row)
# Boundary.knots = NULL: Boundary points at which to anchor the splines (Vector)
#
# -- OUTPUT:
# A matrix containing each units contribution to Gamma(Omega) (Matrix)

Stacked_Estimating_Functions = function(Beta_vec,
                                        gamma_vec,
                                        Nuisance_Parameters_vec,
                                        U_vec,
                                        delta_vec,
                                        Q_vec,
                                        X_mat,
                                        Function_Type,
                                        DF = 1,
                                        Method,
                                        C = 1000,
                                        seed = 235346,
                                        Nuisance_Regression_EF_func,
                                        Nuisance_sigma2e_EF_func,
                                        df_1,
                                        df_2,
                                        df_3,
                                        Boundary.knots = NULL){

  n1 = length(U_vec)

  df2.names = names(df_2)[-1]
  df1_match = df_1[,df2.names]
  df1_match_mat = as.matrix(cbind(1, df1_match))

  # Obtain the estimated mean and variance
  mu_Estimated_vec = df1_match_mat %*% Nuisance_Parameters_vec[1:(length(Nuisance_Parameters_vec)-2)]
  mu_Estimated_vec = as.numeric(mu_Estimated_vec)
  Variance_Estimated_vec = rep(Nuisance_Parameters_vec[length(Nuisance_Parameters_vec)-1] - Nuisance_Parameters_vec[length(Nuisance_Parameters_vec)],
                               n1)

  # Compute the estimating function for (beta, gamma)
  vec.1 = EF_GRC_func(Beta_vec = Beta_vec,
                      gamma_vec = gamma_vec,
                      U_vec = U_vec,
                      delta_vec = delta_vec,
                      Q_vec = Q_vec,
                      X_mat = X_mat,
                      mu_vec = mu_Estimated_vec,
                      Variance_vec = Variance_Estimated_vec,
                      Function_Type = Function_Type,
                      DF = DF,
                      Method = Method,
                      C = C,
                      seed = seed,
                      Boundary.knots = Boundary.knots)

  # Compute the estimating function for Nuisance parameters (excluding measurement error variance)
  mat.i2 = Nuisance_Regression_EF_func(par = Nuisance_Parameters_vec[-length(Nuisance_Parameters_vec)],
                                       Sum = F)

  # Compute the estimating function for the measurement error variance
  vec.i3 = Nuisance_sigma2e_EF_func(par = Nuisance_Parameters_vec[length(Nuisance_Parameters_vec)],
                                    Sum = F)

  # Return the stacked estimating function
  val_vec = c(vec.1,
              as.numeric(colSums(mat.i2)),
              sum(vec.i3))

  return(val_vec)
}
