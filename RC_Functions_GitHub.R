# ----- Function file ------
# Source file in "RC_GitHub.R"

# Source C++ functions
sourceCpp('Cpp_Functions.cpp')

# ---------------------------------------------------------------------------
#           TABLE OF CONTENTS
#
# (1) Generate Conditional Event Rate for Cox Regression Model
#
# (2) True Survival Function
#
# (3) Survival Prediction Function at Given Times
#
# (4) H(.) and partial H(.) / partial beta
#
# (5) Generalized Regression Calibration Score Function
#
# (6) Breslow Estimator for Generalized Regression Calibration
#
# (7) Stacked Estimating Functions
#
# (8) Function to compute B(beta, gamma); Lin and Wei (1989)
#
# (9) Estimation Wrapper: "Gen_RC"
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# (1) Generate Conditional Event Rate for Cox Regression Model
# ---------------------------------------------------------------------------

# Write a function that computes
# lambda(t) = lambda_0 * exp(g(Z; Beta) + gamma*X)
#  *** To be used to simulate data under the Cox model

# -- INPUT:
# Z_vec: Covariate that may have a non-linear effect on the hazard (Vector)
# X_mat: Covariates that have a constant effect on the hazard (Matrix)
# lambda_0: (Constant) Baseline hazard function parameter (Scalar)
# Beta_vec: Regression parameters for Z (Vector)
# gamma_vec: Regression parameters for X (Vector)
# Covariate_Setting = 1: Setting for g(.) specification
#     1: Linear Function
#     2: Quadratic Function
#     3: Tertiary Function
#     4: Non-Linear Function ("NonLinear_func")
#
# -- OUTPUT:
# Rates_vec: Rates that satisfy the Cox regression model (Vector)


Rates_func = function(Z_vec,
                      X_mat,
                      lambda_0,
                      Beta_vec,
                      gamma_vec,
                      Covariate_Setting = 1){
  
  # Define a nonlinear function
  NonLinear_func = function(x){
    0.2 * x / sqrt(1+x^2)
  }
  
  N1 =length(Z_vec)
  
  if (Covariate_Setting == 1){
    val = lambda_0 * exp(Z_vec * Beta_vec + 
                           as.numeric(X_mat %*% gamma_vec))
  }
  if (Covariate_Setting == 2){
    val = lambda_0 * exp(Z_vec * Beta_vec[1] + 
                           Z_vec^2 * Beta_vec[2] + 
                           as.numeric(X_mat %*% gamma_vec))
  }
  if (Covariate_Setting == 3){
    val = lambda_0 * exp(Z_vec * Beta_vec[1] + 
                           Z_vec^2 * Beta_vec[2] + 
                           Z_vec^3 * Beta_vec[3] + 
                           as.numeric(X_mat %*% gamma_vec))
  }
  if (Covariate_Setting == 4){
    g = sapply(1:N1,
               function(i){
                 NonLinear_func(Z_vec[i])
               })
    val = lambda_0 * exp(g + as.numeric(X_mat %*% gamma_vec))
  }
  return(val)
}


# ---------------------------------------------------------------------------
# (2) True Survival Function
# ---------------------------------------------------------------------------


# Write a function that will compute 
# exp(-Lambda(t | X,Z)), where
# Lambda(t | X,Z) = Lambda_0(t) * exp(g(Z; beta) + gamma' X)


# -- INPUT:
# Times: Event times (Vector)
# Lambda_0t: Cumulative baseline hazard functions evaluated at time "Times" (Vector)
# Z: Covariate(s) (Scalar)
# X: Additional Covariate(s) (Scalar)
# Beta_vec: Parameters for Z (Vector)
# gamma_Vec: Parameters for X (Vector)
# Spline_Z_vec = NULL: Data for fitting a B-Spline (Vector)
# Function_Type: Functional form for g(Z; beta) (Character)
#                 "Polynomial": Polynomial
#                 "Bernstein": Bernstein Polynomial
#                 "B-Spline": B-Spline
#                 "C-Spline": C-Spline
#                 "I-Spline": I-Spline
#                 "M-Spline": M-Spline
#                 "NC-Spline": Natural Cubic Spline
# DF = 3: Degrees of freedom for the Spline/Polynomial (Scalar) 
# Boundary.knots = NULL: Boundary points at which to anchor the splines (Vector)

# -- OUTPUT:
# A dataframe containing
# (i) Times: Event times 
# (ii) Surv_Probs: Survival probability


True_SF_func = function(Times,
                        Lambda_0t,
                        Z,
                        X,
                        Beta_vec,
                        gamma_vec,
                        Spline_Z_vec = NULL,
                        Function_Type = "Polynomial",
                        DF = 3,
                        Boundary.knots = NULL){
  
  if (length(Beta_vec) != DF) stop("'Beta_vec' needs to be the same as 'DF'!")
  
  if (Function_Type == "Polynomial"){
    pred.Z = outer(Z, 1:DF, "^")
  } else{
    if (Function_Type == "Bernstein"){
      spline_data = splines2::bernsteinPoly(Spline_Z_vec,
                                            degree = DF)
    }
    if (Function_Type == "B-Spline"){
      spline_data = splines2::bSpline(Spline_Z_vec,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "C-Spline"){
      spline_data = splines2::cSpline(Spline_Z_vec,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "I-Spline"){
      spline_data = splines2::iSpline(Spline_Z_vec,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "M-Spline"){
      spline_data = splines2::mSpline(Spline_Z_vec,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "NC-Spline"){
      spline_data = splines2::naturalSpline(Spline_Z_vec,
                                            df = DF,
                                            Boundary.knots = Boundary.knots)
    } 
    pred.Z = suppressWarnings(as.numeric(predict(spline_data,
                                                 newx = Z)))
  }
  Lambda_t = Lambda_0t * exp( as.numeric(pred.Z %*% Beta_vec) + 
                                as.numeric(t(X) %*% gamma_vec) )

  Surv_Probs = exp(-Lambda_t) 
  
  val = data.frame(Times = Times,
                   Surv_Probs = Surv_Probs)
  
  return(val)
}


# ---------------------------------------------------------------------------
# (3) Survival Prediction Function at Given Times
# ---------------------------------------------------------------------------

# Given estimated survival probabilities, write a function that obtains survival
# probabilities at specific times

# -- INPUT:
# Times_Specify: Times that we want to estimate the survival probability (Vector)
# Est_Times: Times that the survival probability is estimated at (Vector)
# Est_Surv: Estimated survival probabilities at time "Est_Times" (Vector)

# -- OUTPUT:
# A dataframe containing
# (i) Times_Specify: Event times 
# (ii) Surv_Probs: Survival probability evaluated at "Times_Specify"


SurvProb_Times_func = function(Times_Specify,
                               Est_Times,
                               Est_Probs){
  
  val = sapply(Times_Specify,
               function(s){
                 
                 val.s = 1 # Initialize

                 ind = which(Est_Times <= s)
                 
                 if (length(ind) > 0){
                   val.s = Est_Probs[max(ind)]
                 }
                 return(val.s)
               })
  return(val)
}




# ---------------------------------------------------------------------------
# (4) H(.) and partial H(.) / partial beta
# ---------------------------------------------------------------------------

# Write a function that computes 
# H(beta, Q, V; theta) = E(exp(g(Z;beta)) | Q, V; theta)
# and its partial derivative wrt beta

# -- INPUT:
# Method: The method we want to employ (Character)
#           one of {"Cox", "Linear", "Quadratic", "General_1", "General_", "General_MC"}
# Function_Type: Functional form for g(Z; beta) (Character)
#                 "Polynomial": Polynomial
#                 "Bernstein": Bernstein Polynomial
#                 "B-Spline": B-Spline
#                 "C-Spline": C-Spline
#                 "I-Spline": I-Spline
#                 "M-Spline": M-Spline
#                 "NC-Spline": Natural Cubic Spline
# Beta_vec: Parameters for Effect of Interest (Vector)
# DF: Degree for the polynomial or degrees of freedom for the Spline (Scalar)
# data: dataframe that contains Q, X, mu, sigma2 (Dataframe)
# seed = 235346: seed for "General_MC" (Scalar)
# C = 1000: Number of Monte Carlo samples
# Boundary.knots = NULL: Boundary points at which to anchor the splines (Vector)
#
# -- OUTPUT:
# List containing
# (i) H_vec: H(.) for each unit (Vector)
# (ii) partialH_mat: Partial derivative for each unit wrt beta (Matrix)


H_partialH_func = function(Method,
                           Function_Type,
                           Beta_vec,
                           DF,
                           data,
                           seed,
                           C,
                           Boundary.knots = NULL){
  
  require(splines2)
  Q_vec = data$Q
  if (any(is.null(Boundary.knots))){ 
    Boundary.knots = range(Q_vec)
  }
  if (Method == "Cox"){
    
    if (length(Beta_vec) != DF) stop("'Beta_vec' needs to be the same as 'DF'!")
    
    # Initialize matrices
    g_mat = matrix(NA, nrow = dim(data)[1], ncol = DF) # Initialize
    colnames(g_mat) = as.character(1:DF)

    if (Function_Type == "Polynomial"){
      spline_data = outer(data$Q, 1:DF, "^")
    }
    if (Function_Type == "Bernstein"){
      spline_data = bernsteinPoly(data$Q,
                                  degree = DF,
                                  Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "B-Spline"){
      spline_data = bSpline(data$Q,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "C-Spline"){
      spline_data = cSpline(data$Q,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "I-Spline"){
      spline_data = iSpline(data$Q,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "M-Spline"){
      spline_data = mSpline(data$Q,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "NC-Spline"){
      spline_data = naturalSpline(data$Q,
                                  df = DF,
                                  Boundary.knots = Boundary.knots)
    }

    # Compute H(Beta) and its partial derivative for each unit
    g_mat = spline_data %*% Beta_vec
    g_vec = apply(g_mat,1,sum)
    H_vec = exp(g_vec)
    partialH_mat = spline_data * H_vec
  }
  if (Method == "Linear"){
    
    if (length(Beta_vec) != 1) stop("'Beta_vec' needs to be a scalar!")
    
    # Compute H(Beta) and its partial derivative for each unit
    g_vec = data$mu + 0.5 * Beta_vec * data$sigma2
    H_vec = exp(Beta_vec * g_vec)
    partialH_mat = as.matrix(data$mu + Beta_vec * data$sigma2) * H_vec
    
  }
  if (Method == "Quadratic"){
    
    if (length(Beta_vec) != 2) stop("'Beta_vec' needs to be a vector of length 2!")

    DF = 2 # Override "DF"
    spline_data = outer(data$mu, 1:DF, "^")
    
    # Initialize matrices
    g_mat = matrix(NA, nrow = dim(data)[1], ncol = DF) # Initialize
    colnames(g_mat) = as.character(1:DF)
    
    # Compute H(Beta) and its partial derivative for each unit
    g_mat = spline_data %*% Beta_vec
    g_vec = apply(g_mat,1,sum)
    H_vec = exp(g_vec)
    partialH_mat = spline_data * H_vec
  }
  if (Method == "General_1"){
    
    if (length(Beta_vec) != DF) stop("length of 'Beta_vec' needs to match 'DF'!")
    
    data = data[order(data$mu),]

    if (Function_Type == "Polynomial"){
      spline_data = outer(data$mu, 1:DF, "^")
    }
    if (Function_Type == "Bernstein"){
      spline_data = bernsteinPoly(data$mu,
                                  degree = DF,
                                  Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "B-Spline"){
      spline_data = bSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "C-Spline"){
      spline_data = cSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "I-Spline"){
      spline_data = iSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "M-Spline"){
      spline_data = mSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "NC-Spline"){
      spline_data = naturalSpline(data$mu,
                                  df = DF,
                                  Boundary.knots = Boundary.knots)
    }
    
    # Initialize matrices
    g_mat = matrix(NA, nrow = dim(data)[1], ncol = DF) # Initialize
    colnames(g_mat) = as.character(1:DF)
    
    # Compute H(Beta) and its partial derivative for each unit
    g_mat = spline_data %*% Beta_vec
    g_vec = apply(g_mat,1,sum)
    H_vec = exp(g_vec)
    partialH_mat = spline_data * H_vec
  }
  if (Method == "General_2"){

    data = data[order(data$mu),]

    if (Function_Type == "Polynomial"){
      spline_data = outer(data$mu, 1:DF, "^")
      spline1_data = outer(data$mu, 1:DF, function(x, j) ifelse(j == 0, 0, j * x^(j - 1))) 
      spline2_data = outer(data$mu, 1:DF, function(x, j) ifelse(j <= 1, 0, j * (j - 1) * x^(j - 2))) 
    }
    if (Function_Type == "Bernstein"){
      spline_data = bernsteinPoly(data$mu,
                                  degree = DF,
                                  Boundary.knots = Boundary.knots)
      spline1_data = deriv(spline_data, derivs = 1)
      spline2_data = deriv(spline_data, derivs = 2)
    }
    if (Function_Type == "B-Spline"){
      spline_data = bSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
      spline1_data = deriv(spline_data, derivs = 1)
      spline2_data = deriv(spline_data, derivs = 2)
    }
    if (Function_Type == "C-Spline"){
      spline_data = cSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
      spline1_data = deriv(spline_data, derivs = 1)
      spline2_data = deriv(spline_data, derivs = 2)
    }
    if (Function_Type == "I-Spline"){
      spline_data = iSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
      spline1_data = deriv(spline_data, derivs = 1)
      spline2_data = deriv(spline_data, derivs = 2)
    }
    if (Function_Type == "M-Spline"){
      spline_data = mSpline(data$mu,
                            df = DF,
                            Boundary.knots = Boundary.knots)
      spline1_data = deriv(spline_data, derivs = 1)
      spline2_data = deriv(spline_data, derivs = 2)
    }
    if (Function_Type == "NC-Spline"){
      spline_data = naturalSpline(data$mu,
                                  df = DF,
                                  Boundary.knots = Boundary.knots)
      spline1_data = deriv(spline_data, derivs = 1)
      spline2_data = deriv(spline_data, derivs = 2)
    }
    
    # Initialize matrices
    g_mat0 = matrix(NA, nrow = dim(data)[1], ncol = DF) # Initialize
    colnames(g_mat0) = as.character(1:DF)
    g_mat1 = g_mat0
    g_mat2 = g_mat0
    
    # Compute H(Beta) and its partial derivative for each unit
    g_mat0 = spline_data %*% Beta_vec
    g_mat1 = spline1_data %*% Beta_vec
    g_mat2 = spline2_data %*% Beta_vec
    g_vec0 = apply(g_mat0,1,sum)
    g_vec1 = apply(g_mat1,1,sum)
    g_vec2 = apply(g_mat2,1,sum)
    H_vec = exp(g_vec0) * (1 + (0.5*data$sigma2 * g_vec1^2) + g_vec2 )
    partialH_mat = (H_vec * spline_data) + (exp(g_vec0) * (data$sigma2 * (g_vec1 + 1) * spline2_data))
  }
  if (Method == "General_MC"){
    
    if (length(Beta_vec) != DF) stop("'Beta_vec' needs to be the same as 'DF'!")

    set.seed(seed = seed)
    
    # Initialize lists to store Monte Carlo versions for H(.) and its partial derivative
    H_MC_list = list()
    partialH_MC_list = list()

    for (mc in 1:C){
      
      # Generate iid N(mu, sigma2) random variables
      Z_mc_vec = rnorm(dim(data)[1], 
                       mean = data$mu, 
                       sd = sqrt(data$sigma2))
      
      # Initialize matrices
      g_mc_mat = matrix(NA, nrow = dim(data)[1], ncol = DF) # Initialize
      colnames(g_mc_mat) = as.character(1:DF)

      if (Function_Type == "Polynomial"){
        spline_mc_data = outer(Z_mc_vec, 1:DF, "^")
      }
      if (Function_Type == "Bernstein"){
        spline_mc_data = bernsteinPoly(Z_mc_vec,
                                       degree = DF,
                                       Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "B-Spline"){
        spline_mc_data = bSpline(Z_mc_vec,
                                 df = DF,
                                 Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "C-Spline"){
        spline_mc_data = cSpline(Z_mc_vec,
                                 df = DF,
                                 Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "I-Spline"){
        spline_mc_data = iSpline(Z_mc_vec,
                                 df = DF,
                                 Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "M-Spline"){
        spline_mc_data = mSpline(Z_mc_vec,
                                 df = DF,
                                 Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "NC-Spline"){
        spline_mc_data = naturalSpline(Z_mc_vec,
                                       df = DF,
                                       Boundary.knots = Boundary.knots)
      }
      
      # Compute H(Beta) and its partial derivative for each unit
      g_mc_mat = spline_mc_data %*% Beta_vec
      g_mc_vec = apply(g_mc_mat,1,sum)
      H_mc_vec = exp(g_mc_vec)
      partialH_mc_mat = spline_mc_data * H_mc_vec

      H_MC_list[[mc]] = H_mc_vec
      partialH_MC_list[[mc]] = partialH_mc_mat
    }
    
    # Take the average of "H_MC_list" and "partialH_MC_list"
    H_vec = Reduce("+", H_MC_list) / length(H_MC_list)
    partialH_mat = Reduce("+", partialH_MC_list) / length(partialH_MC_list)
  }
  val = list(H_vec = H_vec,
             partialH_mat = partialH_mat)
  return(val)
}



# ---------------------------------------------------------------------------
# (5) Generalized Regression Calibration Score Function
# ---------------------------------------------------------------------------


# Write a function that computes the estimating function for the 
#   generalized regression calibration

# -- INPUT:
# Beta_vec: Parameters for Effect of Interest (Vector)
# gamma_vec: Parameters for Additional Covariates (Vector)
# U_vec: Observed Event Times (Vector)
# delta_vec: Observed Event Time Indicators (Vector)
# Q_vec: Covariate for Effect of Interest (Vector)
# X_mat: Additional Covariates (Vector or matrix)
# mu_vec: Mean of [Z | Q, V] for each unit (Vector)
# Variance_vec: Variance of [Z | Q, V] for each unit (Vector)
# Function_Type: Functional form for g(Z; beta) (Character)
#                 "Polynomial": Polynomial
#                 "Bernstein": Bernstein Polynomial
#                 "B-Spline": B-Spline
#                 "C-Spline": C-Spline
#                 "I-Spline": I-Spline
#                 "M-Spline": M-Spline
#                 "NC-Spline": Natural Cubic Spline
# DF: Degree for the polynomial or degrees of freedom for the Spline (Scalar)
# Method: The method we want to employ (Character)
#           one of {"Cox", "Linear", "Quadratic", "General_1", "General_", "General_MC"}
# C = 1000: Number of Monte Carlo samples
# seed = 235346: seed for "General_MC" (Scalar)
# Boundary.knots = NULL: Boundary points at which to anchor the splines (Vector)
#
# -- OUTPUT:
# U(Beta_vec, gamma_vec): Score function evaluated at Beta_vec and gamma_vec (Vector)



EF_GRC_func = function(Beta_vec,
                       gamma_vec,
                       U_vec,
                       delta_vec,
                       Q_vec,
                       X_mat,
                       mu_vec,
                       Variance_vec,
                       Function_Type,
                       DF = 1,
                       Method,
                       C = 1000,
                       seed = 235346,
                       Boundary.knots = NULL){
  
  if (length(Beta_vec) != DF) stop("length of 'Beta_vec' needs to match 'DF'!")

  data = data.frame(U = U_vec,
                    delta = delta_vec,
                    Q = Q_vec,
                    X = X_mat,
                    mu = mu_vec,
                    sigma2 = Variance_vec)
  
  # Get H(.) and its partial derivative
  H_PartialH_list = H_partialH_func(Method = Method,
                                    Function_Type = Function_Type,
                                    Beta_vec = Beta_vec,
                                    DF = DF,
                                    data = data,
                                    seed = seed,
                                    C = C,
                                    Boundary.knots = Boundary.knots)
  H_vec = H_PartialH_list$H_vec
  partialH_mat = H_PartialH_list$partialH_mat
  
  ind1 = which(colnames(data) %in% "Q")
  ind2 = which(colnames(data) %in% "mu")
  X_mat = data[,(ind1+1):(ind2-1)]
  
  # Compute partial_H(Beta) / H(Beta) for each unit
  partialH_Divided_H_mat = partialH_mat / H_vec
  colnames(partialH_Divided_H_mat) = as.character(1:DF)

  data = cbind(data, partialH_Divided_H_mat)
  
  # Create a dataframe for units that experienced the event
  index_delta1 = which(data$delta == 1)
  data_delta1 = data[index_delta1,]
  data_delta1.sorted = data_delta1[order(data_delta1$U),]
  Covariates_delta1.sorted = cbind(data_delta1.sorted[,which(colnames(data) %in% 1:DF)],
                                 data_delta1.sorted[,(ind1+1):(ind2-1)])
  
  # Get indices when each "pseudo"individual enters and leaves the study
  Index_EnterRiskSet = rep(0, dim(data)[1]) # Assuming every unit joins the study at time 0
  Index_LeaveRiskSet = findInterval(data$U,
                                    data_delta1.sorted$U)
  
  # Compute S0 and S1 
  S_01 = S_01_GRC_cpp(H_vec = H_vec,
                      partialH_mat = partialH_mat,
                      gamma_vec = gamma_vec,
                      X_mat = as.matrix(X_mat),
                      m = dim(data_delta1)[1],
                      IndexEnterRiskSet = Index_EnterRiskSet,
                      IndexLeaveRiskSet = Index_LeaveRiskSet)

  S0 = S_01$r0[1:dim(data_delta1)[1]]
  S1 = S_01$r1[1:dim(data_delta1)[1],]
  S1_S0 = S1 / S0 
  W.S1_S0 = Covariates_delta1.sorted - S1_S0
  val = colSums(W.S1_S0)
  
  return(as.numeric(val))
}



# ---------------------------------------------------------------------------
# (6) Breslow Estimator for Generalized Regression Calibration
# ---------------------------------------------------------------------------


# Write a function for the estimator of the baseline hazard function

# -- INPUT:
# Beta_vec: Parameters for Effect of Interest (Vector)
# gamma_vec: Parameters for Additional Covariates (Vector)
# U_vec: Observed Event Times (Vector)
# delta_vec: Observed Event Time Indicators (Vector)
# Q_vec: Covariate for Effect of Interest (Vector)
# X_mat: Additional Covariates (Vector or matrix)
# mu_vec: Mean of [Z | Q, V] for each unit (Vector)
# Variance_vec: Variance of [Z | Q, V] for each unit (Vector)
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
# Boundary.knots = NULL: Boundary points at which to anchor the splines (Vector)
#
# -- OUTPUT:
# A dataframe containing information for the cumulative baseline hazard function


Baseline_Hazard_GRC_func = function(Beta_vec,
                                    gamma_vec,
                                    U_vec,
                                    delta_vec,
                                    Q_vec,
                                    X_mat,
                                    mu_vec,
                                    Variance_vec,
                                    Function_Type,
                                    DF = 1,
                                    Method,
                                    C = 1000,
                                    seed = 235346,
                                    Boundary.knots = NULL){
  
  if (length(Beta_vec) != DF) stop("length of 'Beta_vec' needs to match 'DF'!")
  
  data = data.frame(U = U_vec,
                    delta = delta_vec,
                    Q = Q_vec,
                    X = X_mat,
                    mu = mu_vec,
                    sigma2 = Variance_vec)
  
  # Get H(.) and its partial derivative
  H_PartialH_list = H_partialH_func(Method = Method,
                                    Function_Type = Function_Type,
                                    Beta_vec = Beta_vec,
                                    DF = DF,
                                    data = data,
                                    seed = seed,
                                    C = C,
                                    Boundary.knots = Boundary.knots)
  H_vec = H_PartialH_list$H_vec
  partialH_mat = H_PartialH_list$partialH_mat
  
  ind1 = which(colnames(data) %in% "Q")
  ind2 = which(colnames(data) %in% "mu")
  X_mat = data[,(ind1+1):(ind2-1)]
  
  # Compute partial_H(Beta) / H(Beta) for each unit
  partialH_Divided_H_mat = partialH_mat / H_vec
  colnames(partialH_Divided_H_mat) = as.character(1:DF)
  
  data = cbind(data, partialH_Divided_H_mat)
  
  # Create a dataframe for units that experienced the event
  index_delta1 = which(data$delta == 1)
  data_delta1 = data[index_delta1,]
  data_delta1.sorted = data_delta1[order(data_delta1$U),]
  Covariates_delta1.sorted = cbind(data_delta1.sorted[,which(colnames(data) %in% 1:DF)],
                                   data_delta1.sorted[,(ind1+1):(ind2-1)])
  
  # Get indices when each "pseudo"individual enters and leaves the study
  Index_EnterRiskSet = rep(0, dim(data)[1]) # Assuming every unit joins the study at time 0
  Index_LeaveRiskSet = findInterval(data$U,
                                    data_delta1.sorted$U)
  
  # Compute S0 and S1 
  S_01 = S_01_GRC_cpp(H_vec = H_vec,
                      partialH_mat = partialH_mat,
                      gamma_vec = gamma_vec,
                      X_mat = as.matrix(X_mat),
                      m = dim(data_delta1)[1],
                      IndexEnterRiskSet = Index_EnterRiskSet,
                      IndexLeaveRiskSet = Index_LeaveRiskSet)
  
  S0 = S_01$r0[1:dim(data_delta1)[1]]
  
  # Initialize the numerator of the baseline hazard estimator to be 1
  Numerator = rep(1, dim(data_delta1.sorted)[1])
  
  # Initialize a dataframe that has the information we want
  data_Baseline = data.frame(time = data_delta1.sorted$U,
                             Numerator = Numerator,
                             Denominator = S0)
  
  # Take the cumulative sum of the numerator by "data_Baseline$time"
  data_Baseline$Numerator = ave(data_Baseline$Numerator, data_Baseline$time, FUN=cumsum)
  
  # Retain the indices where the difference in the numerator isn't 1
  data_Baseline = data_Baseline[c(which(diff(data_Baseline$Numerator)!= 1), dim(data_Baseline)[1]),]
  
  # Obtain the baseline hazard estimate
  data_Baseline$lambdahat = data_Baseline$Numerator / data_Baseline$Denominator
  
  # Obtain the cumulative baseline hazard estimate
  data_Baseline$Lambdahat = cumsum(data_Baseline$lambdahat)
  
  # include t=0 - define the estimate to be 0
  data_Baseline = rbind.data.frame(rep(0,dim(data_Baseline)[2]),
                                   data_Baseline)
  
  # include t=max(Data_Cox$time2) - define the estimate to be the last event time
  if ( data_Baseline[dim(data_Baseline)[1],]$time != max(data$U) ){
    last_entry = c(max(data$U), rep(0,dim(data_Baseline)[2]-2), data_Baseline[dim(data_Baseline)[1],]$Lambdahat)
    data_Baseline = rbind.data.frame(data_Baseline,
                                     last_entry)
  }

  return(data_Baseline)
}





# ---------------------------------------------------------------------------
# (7) Stacked Estimating Functions
# ---------------------------------------------------------------------------


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



# ---------------------------------------------------------------------------
# (8) Function to compute B(beta, gamma); Lin and Wei (1989)
# ---------------------------------------------------------------------------

# Write a function that will compute Gamma_i^*(Omega) = Gamma_i(Omega) + o_p(1)
# see Lin and Wei (1989)

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
# DF: Degree for the polynomial or degrees of freedom for the Polynomial or Spline (Scalar)
# Method: The method we want to employ (Character)
#           one of {"Cox", "Linear", "Quadratic", "General_1", "General_", "General_MC"}
# C = 1000: Number of Monte Carlo samples
# seed = 235346: seed for "General_MC" (Scalar)
# Nuisance_Regression_EF_func: Nuisance parameter estimating function defined in Function (9) (Function)
# Nuisance_sigma2e_EF_func: Measurement error variance estimating function defined in Function (9) (Function)
# df_1: Dataset to estimate hazard function: (Dataframe)
# df_2: Dataset to estimate nuisance parameters: (Dataframe)
# df_3: Dataset to estimated measurement error: (Dataframe)
# df1_ID_index = NA: index corresponding to ID variable in "df1" (Scalar)
# df2_ID_index = NA: index corresponding to ID variable in "df2" (Scalar)
# df3_ID_index = NA: index corresponding to ID variable in "df3" (Scalar)
# df2_Response_index: index corresponding to response variablein "df2" (Scalar)
# df2_Covariates_index: index/indices corresponding to additional covariates in "df2" (Vector)
# Boundary.knots = NULL: Boundary points at which to anchor the splines (Vector)
#
# -- OUTPUT:
# A matrix containing each units contribution to Gamma_i^*(Omega) (Matrix)


Stacked_Estimating_Functions_Modified = function(Beta_vec,
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
                                                 df1_ID_index = NA,
                                                 df2_ID_index = NA,
                                                 df3_ID_index = NA,
                                                 df2_Response_index,
                                                 df2_Covariates_index,
                                                 Boundary.knots = NULL){
  
  if (!is.na(df1_ID_index)){
    df1_ID = df_1[,df1_ID_index]
  } else{
    df1_ID = 1:dim(df_1)[1]
  }
  if (!is.na(df2_ID_index)){
    df2_ID = df_2[,df2_ID_index]
  } else{
    df2_ID = (dim(df_1)[1]+1):(dim(df_1)[1]+dim(df_2)[1])
  }
  if (!is.na(df3_ID_index)){
    df3_ID = df_3[,df3_ID_index]
  } else{
    df3_ID = (dim(df_1)[1]+dim(df_2)[1]+1):(dim(df_1)[1]+dim(df_2)[1]+dim(df_3)[1])
  }
  
  df123_ID = unique(c(df1_ID, df2_ID, df3_ID))
  n.row = length(unique(c(df1_ID, df2_ID, df3_ID)))
  
  # Initialize matrices to store contributions for the estimating function
  Gamma_mat = matrix(0, nrow = n.row, ncol = length(Beta_vec) + length(gamma_vec) + length(Nuisance_Parameters_vec))
  
  # Compute the estimating function for Nuisance parameters (excluding measurement error variance)
  mat.i2 = Nuisance_Regression_EF_func(par = Nuisance_Parameters_vec[-length(Nuisance_Parameters_vec)],
                                       Sum = F)
  Gamma_mat[which(df123_ID %in% df2_ID), (length(Beta_vec) + length(gamma_vec)+1):(dim(Gamma_mat)[2]-1) ] = mat.i2
  
  # Compute the estimating function for the measurement error variance
  vec.i3 = Nuisance_sigma2e_EF_func(par = Nuisance_Parameters_vec[length(Nuisance_Parameters_vec)],
                                    Sum = F)
  Gamma_mat[which(df123_ID %in% df3_ID), dim(Gamma_mat)[2]  ] = vec.i3
  
  # ----- The rest of the function is to obtain the "adjusted estimating function"
  #         for (Beta, gamma); see Lin and Wei (1989)

  df2.names = names(df_2)[c(df2_Covariates_index)]
  df1_match = df_1[,df2.names]
  df1_match_mat = as.matrix(cbind(1, df1_match))
  
  # Obtain the estimated mean and variance
  mu_Estimated_vec = df1_match_mat %*% Nuisance_Parameters_vec[1:(length(Nuisance_Parameters_vec)-2)]
  mu_Estimated_vec = as.numeric(mu_Estimated_vec)
  Variance_Estimated_vec = rep(Nuisance_Parameters_vec[length(Nuisance_Parameters_vec)-1] - Nuisance_Parameters_vec[length(Nuisance_Parameters_vec)], 
                               length(mu_Estimated_vec))
  
  if (length(Beta_vec) != DF) stop("length of 'Beta_vec' needs to match 'DF'!")

  data1 = data.frame(U = U_vec,
                     delta = delta_vec,
                     Q = Q_vec,
                     X = X_mat,
                     mu = mu_Estimated_vec,
                     sigma2 = Variance_Estimated_vec)
  
  # Get H(.) and its partial derivative
  H_PartialH_list = H_partialH_func(Method = Method,
                                    Function_Type = Function_Type,
                                    Beta_vec = Beta_vec,
                                    DF = DF,
                                    data = data1,
                                    seed = seed,
                                    C = C,
                                    Boundary.knots = Boundary.knots)
  H_vec = H_PartialH_list$H_vec
  partialH_mat = H_PartialH_list$partialH_mat

  ind1 = which(colnames(data1) %in% "Q")
  ind2 = which(colnames(data1) %in% "mu")
  X_mat = data1[,(ind1+1):(ind2-1)]
  
  # Compute partial_H(Beta) / H(Beta) for each unit
  partialH_Divided_H_mat = partialH_mat / H_vec
  colnames(partialH_Divided_H_mat) = as.character(1:DF)
  data1 = cbind(data1, partialH_Divided_H_mat)

  data1.sorted = data1[order(data1$U),]
  index_original_data1 = order(data1$U)
  Covariates.sorted = cbind(data1.sorted[,which(colnames(data1) %in% 1:DF)],
                            data1.sorted[,(ind1+1):(ind2-1)])
  Covariates.sorted = as.matrix(Covariates.sorted)
  
  # Get indices when each "pseudo"individual enters and leaves the study
  Index_EnterRiskSet_all = rep(0, dim(data1)[1]) # Assuming every unit joins the study at time 0
  Index_LeaveRiskSet_all = findInterval(data1$U,
                                        data1.sorted$U)
  
  # Compute S0 and S1 
  S_01 = S_01_GRC_cpp(H_vec = H_vec,
                      partialH_mat = partialH_mat,
                      gamma_vec = gamma_vec,
                      X_mat = as.matrix(X_mat),
                      m = dim(df_1)[1],
                      IndexEnterRiskSet = Index_EnterRiskSet_all,
                      IndexLeaveRiskSet = Index_LeaveRiskSet_all)
  
  S0 = S_01$r0[1:dim(df_1)[1]]
  S1 = S_01$r1[1:dim(df_1)[1],]
  S1_S0 = S1 / S0 
  
  # Compute "data1.sorted$delta * (Covariates.sorted - S1_S0)"
  part1_mat = data1.sorted$delta * (Covariates.sorted - S1_S0)
  
  # Extract the event time indicators from "data1.sorted" and add in 0 for a time after the final observed time
  delta_vec.sorted = data1.sorted$delta
  
  # Recompute the indices when they leave the risk set
  # (wrt to the ordered data)
  Index_LeaveRiskSet_ordered_all = sapply(1:dim(df_1)[1], 
                                          function(i){
                                            val = max(which(data1.sorted$U %in% data1.sorted$U[i])) # "max" in case of ties
                                            return(val)
                                          })
  
  # Initialize matrices to store results
  part2_mat = matrix(0, nrow = dim(df_1)[1], ncol = dim(part1_mat)[2])
  
  for (i in 1:dim(df_1)[1]){
    
    # Get the index when they enter and leave the risk set
    a1 = Index_EnterRiskSet_all[i] + 1 
    a2 = Index_LeaveRiskSet_ordered_all[i] 
    
    # Compute exp(gamma'X) * H for unit i
    ind.i = index_original_data1[i]
    e_H.i = exp(as.matrix(X_mat)[ind.i,] %*% gamma_vec) * H_vec[ind.i]
    e_H.i = as.numeric(e_H.i)
    
    S1_S0.a1.a2 = S1[a1:a2,] / S0[a1:a2]
    if (a1 == a2) S1_S0.a1.a2 = t(S1_S0.a1.a2)

    d.a1.a2 = -sweep(S1_S0.a1.a2, 
                     2, 
                     Covariates.sorted[i,], 
                     FUN = "-")
    
    val.i = delta_vec.sorted[a1:a2] * e_H.i / S0[a1:a2] * d.a1.a2
    
    # Update "part2_mat[i,]"
    part2_mat[i,] = colSums(val.i) 
  }

  mat.i1 = part1_mat - part2_mat
  Gamma_mat[index_original_data1, 1:(length(Beta_vec) + length(gamma_vec)) ] = mat.i1
  
  return(Gamma_mat)
}


# ---------------------------------------------------------------------------
# (9) Estimation Wrapper: "Gen_RC"
# ---------------------------------------------------------------------------

# Write a function that will estimate the parameters, estimate the variance,
# and estimate survival probabilities

# -- INPUT:
# df1: Dataset to estimate hazard function: (Dataframe)
# df2: Dataset to estimate nuisance parameters: (Dataframe)
# df3: Dataset to estimated measurement error: (Dataframe)
#   NOTE: df3 needs to be in "wide format" (all observations for a subject is in a single row)
# df1_ID_index = NA: index corresponding to ID variable in "df1" (Scalar)
# df1_U_index: Index corresponding to observed event time in "df1" (Scalar)
# df1_delta_index: Index corresponding to observed event time indicator in "df1" (Scalar)
# df1_Q_index: index corresponding to variable of interest in "df1" for the hazard function (Scalar)
# df1_X_index: index/indices corresponding to additional covariates in "df1" for the hazard function (Vector)
# df2_ID_index = NA: index corresponding to ID variable in "df2" (Scalar)
# df2_Response_index: index corresponding to response variablein "df2" (Scalar)
# df2_Covariates_index: index/indices corresponding to additional covariates in "df2" (Vector)
# df3_ID_index = NA: index corresponding to ID variable in "df3" (Scalar)
# Function_Type: Functional form for g(Z; beta) (Character)
#                 "Polynomial": Polynomial
#                 "Bernstein": Bernstein Polynomial
#                 "B-Spline": B-Spline
#                 "C-Spline": C-Spline
#                 "I-Spline": I-Spline
#                 "M-Spline": M-Spline
#                 "NC-Spline": Natural Cubic Spline
# DF: Degree for the polynomial or degrees of freedom for the Polynomial or Spline (Scalar)
# Method: The method we want to employ (Character)
#           one of {"Cox", "Linear", "Quadratic", "General_1", "General_2", "General_MC"}
# EstimateNuisance_Logi = T: Indicator if we want to estimate mean and variance of [Z | Q,V] (Logical)
# mu_vec = NULL: Mean of [Z | Q, V] for each unit (Vector)
# Variance_vec = NULL: Variance of [Z | Q, V] for each unit (Vector)
# t_grid: Times to evaluate survival function (Vector)
# z_grid: Values of Z to evaluate conditional survival function (Vector)
# C = 1000: Number of Monte Carlo samples (Scalar)
# seed = 235346: seed for "General_MC" (Scalar)
# x = 0: Value of X to be evaluated in the survival function (Vector)
# Spline_Z_vec = NULL: Data for fitting a Spline (Vector)
# Boundary.knots = NULL: Boundary points at which to anchor the splines (Vector)
# Estimate_CI_Surv = T: Indicator if we want a Monte Carlo confidence interval for survival functions (Logical)
#
# -- OUTPUT:
# A list containing
#   "Estimates": Estimates of parameters for the hazard function (Vector)
#   "Nuisance_Estimates": Nuisance Parameters Estimates (Vector)
#   "Variance": Variance of estimated parameters in the order of the hazard regression parameters and then nuisance parameters (Matrix)
#   "SurvProb_zgrid_list": Survival probabilities evaluated at each "t_grid" and "z_grid" (List)
#   "fit_EF": Output from "nleqslv::nleqslv" that solves the estimating equation for hazard parameter estimation


Gen_RC = function(df1,
                  df2,
                  df3,
                  df1_ID_index = NA,
                  df1_U_index,
                  df1_delta_index,
                  df1_Q_index,
                  df1_X_index,
                  df2_ID_index = NA,
                  df2_Response_index,
                  df2_Covariates_index,
                  df3_ID_index = NA,
                  Function_Type = "Polynomial",
                  DF = 1,
                  Method = "Cox",
                  EstimateNuisance_Logi,
                  mu_vec = NULL,
                  Variance_vec = NULL,
                  t_grid,
                  z_grid,
                  C = 1000,
                  seed = 235346,
                  x = 0,
                  Spline_Z_vec = NULL,
                  Boundary.knots = NULL,
                  Estimate_CI_Surv = T){
  
  if (length(df1_ID_index) > 1) stop("'df1_ID_index' needs to be a scalar!")
  if (length(df1_U_index) > 1) stop("'df1_U_index' needs to be a scalar!")
  if (length(df1_delta_index) > 1) stop("'df1_delta_index' needs to be a scalar!")
  if (length(df1_Q_index) > 1) stop("'df1_Q_index' needs to be a scalar!")
  if (length(df2_ID_index) > 1) stop("'df2_ID_index' needs to be a scalar!")
  if (length(df2_Response_index) > 1) stop("'df2_Response_index' needs to be a scalar!")
  if (length(df3_ID_index) > 1) stop("'df3_ID_index' needs to be a scalar!")
  
  
  # Store the IDs of "df1" "df2" and "df3"
  df1_ID = df1[,df1_ID_index]
  df2_ID = df2[,df2_ID_index]
  df3_ID = df3[,df3_ID_index]
  
  if (length(df1_ID) != length(unique(df1_ID)) & EstimateNuisance_Logi == T) stop("Function cannot handle common IDs within 'df1' when estimating variance")
  if (length(df2_ID) != length(unique(df2_ID)) & EstimateNuisance_Logi == T) stop("Function cannot handle common IDs within 'df2' when estimating variance")
  if (length(df3_ID) != length(unique(df3_ID)) & EstimateNuisance_Logi == T) stop("Function cannot handle common IDs within 'df3' when estimating variance")
  
  if (class(df1)[1] != "data.frame") df1 = as.data.frame(df1)
  if (length(x) < length(df1_X_index)) x = rep(x[1], length(df1_X_index))
  
  # Initialize "mu_vec" and "Variance_vec", if required
  if (is.null(mu_vec[1])) mu_vec = rep(0, dim(df1)[1])
  if (is.null(Variance_vec[1])) Variance_vec = rep(0, dim(df1)[1])
  
  df2_use = df2 # Initialize
  
  # Initialize estimates for nuisance parameters
  df23_NuisanceEstimates = NA
  
  if (EstimateNuisance_Logi == T){
    
    # ----------------------------------
    #   df2 Nuisance Parameters
    # ----------------------------------
    
    if (class(df2)[1] != "data.frame") df2 = as.data.frame(df2)
    df2_use = df2[,c(df2_Response_index,
                     df2_Covariates_index)]
    colnames(df2_use)[1] = "W1"
    
    #   fit a linear regression model
    fit_lm = lm(W1 ~ .,
                data = df2_use)
    
    # Write estimating functions for the regression parameters
    Nuisance_Regression_Parameter_EF_func = function(alpha_vec,
                                                     df,
                                                     Sum = T){
      X = as.matrix(cbind(1, df[,2:dim(df)[2]]))
      val_mat = as.numeric(df$W1 - X %*% alpha_vec) * X
      if (Sum == T) val_mat = colSums(val_mat)
      return(val_mat)
    }
    Nuisance_Regression_sigma2_EF_func = function(alpha_vec,
                                                  sigma2,
                                                  df,
                                                  Sum = T){
      m2 = dim(df)[1]
      X = as.matrix(cbind(1, df[,2:dim(df)[2]]))
      val = (as.numeric(df$W1 - X %*% alpha_vec))^2 - (m2-length(alpha_vec))/m2 * sigma2
      if (Sum == T) val = sum(val)
      return(val)
    }
    Nuisance_Regression_EF_func = function(par,
                                           Sum = T){
      val1_mat = Nuisance_Regression_Parameter_EF_func(alpha_vec = par[-length(par)],
                                                       df = df2_use,
                                                       Sum = Sum) 
      val2_vec = Nuisance_Regression_sigma2_EF_func(alpha_vec = par[-length(par)],
                                                    sigma2 = par[length(par)],
                                                    df = df2_use,
                                                    Sum = Sum)
      if (Sum == T){
        val = c(val1_mat, val2_vec)
      } else{
        val = cbind(val1_mat, val2_vec)
      }
      
      return(val)
    }
    
    # Solve the estimating equation
    NuisanceParameters_2_EF = nleqslv(x = rep(0, dim(df2_use)[2]+1),
                                      fn = Nuisance_Regression_EF_func,
                                      jacobian = T,
                                      control = list(trace = 0,
                                                     maxit = 500))
    
    # ----------------------------------
    #   df3 Nuisance Parameters
    # ----------------------------------
    
    if (class(df3)[1] != "data.frame") df3 = as.data.frame(df3)
    
    df3_use = df3 # Initialize
    if (!is.na(df3_ID_index)) df3_use = df3_use[,-df3_ID_index]
    
    df3_use$Wbar = rowMeans(df3_use, na.rm = TRUE)  
    
    # The corresponding estimating function is...
    Nuisance_sigma2e_EF_func = function(par,
                                        Sum = T){
      val = rowSums((df3_use[, 1:(dim(df3_use)[2]-1)] - df3_use$Wbar)^2) - par
      if (Sum == T) val = sum(val)
      return(val)
    }
    
    # Solve the estimating equation
    NuisanceParameters_3_EF = nleqslv(x = 0,
                                      fn = Nuisance_sigma2e_EF_func,
                                      jacobian = T,
                                      control = list(trace = 0,
                                                     maxit = 500))
    sigma2_ehat = NuisanceParameters_3_EF$x
    df23_NuisanceEstimates = c(NuisanceParameters_2_EF$x, sigma2_ehat)
    
    # Get the variables in "df1" that were used to fit the linear model
    df2.names = names(df2_use)[-1]
    df1_match = df1[,df2.names]
    if (length(fit_lm$coefficients)-1 != dim(df1_match)[2]) stop("Variables in 'df1' and 'df2' do not match...")
    
    # Define "mu_vec" and "Variance_vec"
    mu_vec = as.numeric(predict(fit_lm,
                                newdata = df1_match))
    Variance_vec = rep(summary(fit_lm)$sigma^2 - sigma2_ehat, dim(df1)[1])
    
  }
  
  # Setup the estimating function
  EF_GRC = function(par){
    EF_GRC_func(Beta_vec = par[1:DF],
                gamma_vec = par[(DF+1):(DF+length(df1_X_index))],
                U_vec = df1[,df1_U_index],
                delta_vec = df1[,df1_delta_index],
                Q_vec = df1[,df1_Q_index],
                X_mat = df1[,df1_X_index],
                mu_vec = mu_vec,
                Variance_vec = Variance_vec,
                Function_Type = Function_Type,
                DF = DF,
                Method = Method,
                Boundary.knots = Boundary.knots)
  }
  
  # Solve the estimating equation
  fit_EF = nleqslv(x = rep(0, DF+length(df1_X_index)),
                   fn = EF_GRC,
                   jacobian = T,
                   control = list(trace = 0,
                                  maxit = 500))
  
  # Make the correction if "Method == 'Quadratic'"
  if (Method == "Quadratic"){
    Betahat = rep(NA, 2)
    sigma2 = Variance_vec[1]
    Betahat[1] = fit_EF$x[1] / (1 + 2*fit_EF$x[2] * sigma2)
    Betahat[2] = fit_EF$x[2] / (1 + 2*fit_EF$x[2] * sigma2)
  } else{
    Betahat = fit_EF$x[1:DF]
  }
  
  # Estimate the baseline hazard function
  p = length(fit_EF$x)
  Baseline_Hazard = Baseline_Hazard_GRC_func(Beta_vec = fit_EF$x[1:DF],
                                             gamma_vec = fit_EF$x[(DF+1):(DF+length(df1_X_index))],
                                             U_vec = df1[,df1_U_index],
                                             delta_vec = df1[,df1_delta_index],
                                             Q_vec = df1[,df1_Q_index],
                                             X_mat = df1[,df1_X_index],
                                             mu_vec = mu_vec,
                                             Variance_vec = Variance_vec,
                                             Function_Type = Function_Type,
                                             DF = DF,
                                             Method = Method,
                                             Boundary.knots = Boundary.knots)
  
  # Make the correction if "Method == 'Quadratic'"
  if (Method == "Quadratic"){
    Correction1 = sqrt(1 - 2*Betahat[2] * sigma2)
    Correction2 = exp(-0.5*sigma2 * Betahat[1]^2 / (1 - 2*Betahat[2]*sigma2))
    Baseline_Hazard$Lambdahat = Baseline_Hazard$Lambdahat * Correction1 * Correction2
  }
  
  # Estimate the variance of parameter estimates
  if (EstimateNuisance_Logi == T){
    
    # Define the stacked estimating function 
    Stacked_EF = function(par){
      BetaGamma_vec = par[1:p]
      Nuisance_par = par[(p+1):length(par)]
      val = Stacked_Estimating_Functions(Beta_vec = BetaGamma_vec[1:DF],
                                         gamma_vec = BetaGamma_vec[(DF+1):(DF+length(df1_X_index))],
                                         Nuisance_Parameters_vec = Nuisance_par,
                                         U_vec = df1[,df1_U_index],
                                         delta_vec = df1[,df1_delta_index],
                                         Q_vec = df1[,df1_Q_index],
                                         X_mat = df1[,df1_X_index],
                                         Function_Type = Function_Type,
                                         DF = DF,
                                         Method = Method,
                                         df_1 = df1,
                                         df_2 = df2_use,
                                         df_3 = df3_use,
                                         Nuisance_Regression_EF_func = Nuisance_Regression_EF_func,
                                         Nuisance_sigma2e_EF_func = Nuisance_sigma2e_EF_func,
                                         Boundary.knots = Boundary.knots)
      return(val)
    }
    
    # Obtain the jacobian
    if (Method %in% "General_MC"){
      A_mat = nleqslv(x = c(fit_EF$x,
                            df23_NuisanceEstimates),
                      fn = Stacked_EF,
                      jacobian = T,
                      control = list(trace = 0,
                                     maxit = 500))$jac
    } else{
      A_mat = numDeriv::jacobian(Stacked_EF,
                                 x = c(fit_EF$x,
                                       df23_NuisanceEstimates))
    }
    
    # Obtain B(.)
    B_mat = crossprod(Stacked_Estimating_Functions_Modified(Beta_vec = fit_EF$x[1:DF],
                                                            gamma_vec = fit_EF$x[(DF+1):(DF+length(df1_X_index))],
                                                            Nuisance_Parameters_vec = df23_NuisanceEstimates,
                                                            U_vec = df1[,df1_U_index],
                                                            delta_vec = df1[,df1_delta_index],
                                                            Q_vec = df1[,df1_Q_index],
                                                            X_mat = df1[,df1_X_index],
                                                            Function_Type = Function_Type,
                                                            DF = DF,
                                                            Method = Method,
                                                            C = C, 
                                                            seed = seed,
                                                            Nuisance_Regression_EF_func = Nuisance_Regression_EF_func,
                                                            Nuisance_sigma2e_EF_func = Nuisance_sigma2e_EF_func,
                                                            df_1 = df1,
                                                            df_2 = df2,
                                                            df_3 = df3,
                                                            df1_ID_index = df1_ID_index,
                                                            df2_ID_index = df2_ID_index,
                                                            df3_ID_index = df3_ID_index,
                                                            df2_Response_index = df2_Response_index,
                                                            df2_Covariates_index = df2_Covariates_index,
                                                            Boundary.knots = Boundary.knots))
    
    # Compute the variance
    Covariance_mat = MASS::ginv(A_mat) %*% B_mat %*% t(MASS::ginv(A_mat))
    Covariance_mat_orig = Covariance_mat
  } else{
    
    # Estimate the variance from the Observed Information matrix
    Covariance_mat = MASS::ginv(-numDeriv::jacobian(EF_GRC,
                                                    x = fit_EF$x))
    Covariance_mat_orig = Covariance_mat
  }
  
  # Obtain Monte Carlo variance estimators if "Method == 'Quadratic'"
  if (Method == "Quadratic"){
    
    # Sample from the estimators asymptotic distribution
    set.seed(seed)
    if (EstimateNuisance_Logi == F){
      mu = fit_EF$x
    } else{
      mu = c(fit_EF$x, df23_NuisanceEstimates)
    }
    Est_C = MASS::mvrnorm(C, 
                          mu = mu, 
                          Sigma = Covariance_mat)
    
    # Initialize a matrix to store results
    Var_Betahat_QuadraticRC = Est_C*0
    
    # Estimates of gamma do not need to be corrected
    Var_Betahat_QuadraticRC[,(DF+1):dim(Var_Betahat_QuadraticRC)[2]] = Est_C[,(DF+1):dim(Var_Betahat_QuadraticRC)[2]]
    
    # Apply the correction to betahat
    Var_Betahat_QuadraticRC[,1] = Est_C[,1] / (1 + 2*Est_C[,2] * sigma2)
    Var_Betahat_QuadraticRC[,2] = Est_C[,2] / (1 + 2*Est_C[,2] * sigma2)
    
    # Update "Covariance_mat"
    Covariance_mat = var(Var_Betahat_QuadraticRC)
  }
  
  # Obtain (baseline) survival predictions
  Baseline_Surv = exp(-Baseline_Hazard$Lambdahat)
  
  # Obtain survival probabilities at times "t_grid"
  Baseline_Surv_tgrid = SurvProb_Times_func(Times_Specify = t_grid,
                                            Est_Times = Baseline_Hazard$time,
                                            Est_Probs = Baseline_Surv)
  
  # Run a for loop to obtain survival probabilities with Z = z_grid[k]
  SurvProb_zgrid_list = list() # Initialize
  SurvProb_zgrid_CI_list = list() # Initialize
  SurvProb_zgrid_Var_list = list() # Initialize
  for (k in 1:length(z_grid)){

    z = z_grid[k]
    
    # Obtain the predicted survival probability
    Est_SurvProb = True_SF_func(Times = t_grid,
                                Lambda_0t = -log(Baseline_Surv_tgrid),
                                Z = z,
                                X = x,
                                Beta_vec = Betahat,
                                gamma_vec = fit_EF$x[(DF+1):(DF+length(df1_X_index))],
                                DF = DF,
                                Function_Type = Function_Type,
                                Spline_Z_vec = Spline_Z_vec,
                                Boundary.knots = Boundary.knots)
    
    SurvProb_zgrid_list[[k]] = Est_SurvProb$Surv_Probs
    
    # Proceed if we want to estimate the survival function confidence interval
    if (Estimate_CI_Surv == T){
      set.seed(seed)
      if (EstimateNuisance_Logi == F){
        mu = fit_EF$x
      } else{
        mu = c(fit_EF$x, df23_NuisanceEstimates)
      }
      Est_C = MASS::mvrnorm(C, 
                            mu = mu, 
                            Sigma = Covariance_mat_orig)
      
      # Initialize a matrix to store results
      SurvProbs.z_mat = matrix(NA, nrow = C, ncol = length(t_grid))
      colnames(SurvProbs.z_mat) = t_grid
      for (cc in 1:C){
        Betahat.c = Est_C[cc, 1:DF]
        gammahat.c = Est_C[cc, (DF+1):(DF+length(df1_X_index))]
        if (EstimateNuisance_Logi == F){
          Nuisance.c = NULL
          sigma2.c = NULL
        } else{
          Nuisance.c = Est_C[cc, (DF+length(df1_X_index)+1):dim(Est_C)[2]]
          sigma2.c = Nuisance.c[length(Nuisance.c)-1]
        }
        
        # Make the correction if "Method == 'Quadratic'"
        if (Method == "Quadratic"){
          Betahat.c_Corrected = rep(NA, 2)
          sigma2.c_Corrected = ifelse(EstimateNuisance_Logi == F,
                                      Variance_vec[1],
                                      sigma2.c)
          Betahat.c_Corrected[1] = Betahat.c[1] / (1 + 2*Betahat.c[2] * sigma2.c_Corrected)
          Betahat.c_Corrected[2] = Betahat.c[2] / (1 + 2*Betahat.c[2] * sigma2.c_Corrected)
        } else{
          Betahat.c_Corrected = Betahat.c
        }
        
        # Estimate the baseline hazard function
        Baseline_Hazard.c = Baseline_Hazard_GRC_func(Beta_vec = Betahat.c_Corrected,
                                                     gamma_vec = gammahat.c,
                                                     U_vec = df1[,df1_U_index],
                                                     delta_vec = df1[,df1_delta_index],
                                                     Q_vec = df1[,df1_Q_index],
                                                     X_mat = df1[,df1_X_index],
                                                     mu_vec = mu_vec,
                                                     Variance_vec = Variance_vec,
                                                     Function_Type = Function_Type,
                                                     DF = DF,
                                                     Method = Method,
                                                     Boundary.knots = Boundary.knots)
        
        # Make the correction if "Method == 'Quadratic'"
        if (Method == "Quadratic"){
          Correction1.c = sqrt(1 - 2*Betahat.c_Corrected[2] * sigma2.c_Corrected)
          Correction2.c = exp(-0.5*sigma2.c_Corrected * Betahat.c_Corrected[1]^2 / (1 - 2*Betahat.c_Corrected[2]*sigma2.c_Corrected))
          Baseline_Hazard.c$Lambdahat = Baseline_Hazard.c$Lambdahat * Correction1.c * Correction2.c
        }
        
        # Obtain (baseline) survival predictions
        Baseline_Surv.c = exp(-Baseline_Hazard.c$Lambdahat)
        
        # Obtain survival probabilities at times "t_grid"
        Baseline_Surv_tgrid.c = SurvProb_Times_func(Times_Specify = t_grid,
                                                    Est_Times = Baseline_Hazard.c$time,
                                                    Est_Probs = Baseline_Surv.c)
        
        # Obtain the predicted survival probability with Z = z, and X = x, with these parameters
        Est_SurvProb.c = True_SF_func(Times = t_grid,
                                      Lambda_0t = -log(Baseline_Surv_tgrid.c),
                                      Z = z,
                                      X = x,
                                      Beta_vec = Betahat.c_Corrected,
                                      gamma_vec = gammahat.c,
                                      DF = DF,
                                      Function_Type = Function_Type,
                                      Spline_Z_vec = Spline_Z_vec,
                                      Boundary.knots = Boundary.knots)
        
        SurvProbs.z_mat[cc,] = Est_SurvProb.c[,2]
      }
      
      
      # Obtain the upper and lower limits
      UL_LL.z = sapply(1:dim(SurvProbs.z_mat)[2],
                       function(tt){
                         data.tt = SurvProbs.z_mat[,tt]
                         val = as.numeric(quantile(data.tt, probs = c(0.025, 0.975)))
                         return(val)
                       })
      rownames(UL_LL.z) = c("2.5%", "97.5%")
      colnames(UL_LL.z) = t_grid
      
      # Obtain pointwise variance of survival probabilities
      Var.z = sapply(1:dim(SurvProbs.z_mat)[2],
                     function(tt){
                       data.tt = SurvProbs.z_mat[,tt]
                       val = as.numeric(var(data.tt))
                       return(val)
                     })
      
      SurvProb_zgrid_Var_list[[k]] = Var.z
      SurvProb_zgrid_CI_list[[k]] = UL_LL.z
    }
  }

  val_list = list(Estimates = c(Betahat, fit_EF$x[(DF+1):(DF+length(df1_X_index))]),
                  Nuisance_Estimates = df23_NuisanceEstimates,
                  Variance = Covariance_mat,
                  Baseline_Hazard = Baseline_Hazard,
                  SurvProb_zgrid_list = SurvProb_zgrid_list,
                  SurvProb_zgrid_CI_list = SurvProb_zgrid_CI_list,
                  SurvProb_zgrid_Var_list = SurvProb_zgrid_Var_list,
                  fit_EF = fit_EF)
  
  return(val_list)
  

}

