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
  data_Baseline$Numerator = stats::ave(data_Baseline$Numerator, data_Baseline$time, FUN=cumsum)

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
