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
