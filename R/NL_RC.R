#' Nonlinear Regression Calibration Under a Cox Regression Model
#'
#' @description
#' Apply regression calibration to estimate parameters in a Cox regression model when the exposure variable is measured with error, and can have a linear or nonlinear effect on the hazard function.
#'
#' @param df1 Data frame corresponding to the analysis cohort that contains right-censored observations of the event time.
#' @param df2 Data frame corresponding to the calibration dataset.
#' @param df3 Data frame corresponding to the replicate dataset. Note that `df3` needs to be in wide format.
#' @param df1_ID_index Column number for the identification variable in `df1`. If not specified (default), it is assumed that no units are shared between `df1` and `df2` / `df3`.
#' @param df1_U_index Column number for the right-censored event time in `df1`.
#' @param df1_delta_index Column number for the event time indicator in `df1`.
#' @param df1_Q_index Column number for the observed exposure variable in `df1`.
#' @param df1_X_index Column number(s) in `df1` for any additional covariates to include in the Cox model. Note that categorical variables need to be transformed to dummy variables.
#' @param df2_ID_index Column number for the identification variable in `df2`. If not specified (default), it is assumed that no units are shared between `df2` and `df1` / `df3`.
#' @param df2_Response_index Column number for the observed biomarker in `df2`.
#' @param df2_Covariates_index Column number(s) in `df2` for any additional covariates in the measurement error model. Note that categorical variables need to be transformed to dummy variables.
#' @param df3_ID_index Column number for the identification variable in `df3`. If not specified (default), it is assumed that no units are shared between `df3` and `df1` / `df2`.
#' @param Function_Type Specifies the functional form of the exposure's relationship with the hazard function. Users can select one of the following:
#' \itemize{
#' \item `"Polynomial"` - Degree `DF` polynomial (default).
#' \item `"Bernstein"` - Linear combination of `DF` Bernstein polynomial basis functions.
#' \item `"B-Spline"` - Linear combination of `DF` B-spline basis functions.
#' \item `"C-Spline"` - Linear combination of `DF` C-spline basis functions.
#' \item `"I-Spline"` - Linear combination of `DF` I-spline basis functions.
#' \item `"M-Spline"` - Linear combination of `DF` M-spline basis functions.
#' \item `"NC-Spline"` - Linear combination of `DF` natural cubic spline basis functions.
#' }
#' @param DF Degree of the polynomial or number of basis functions for splines (default is 1). Note that this parameter is irrelevent whenever `Method = "Linear"` or `Method = "Quadratic"`.
#' @param Method Specifies the estimation procedure. Users can select one of the following:
#' @param EstimateNuisance_Logi
#' \itemize{
#' \item `"Cox"` - Estimate parameters in the Cox model that ignores measurement error (default).
#' \item `"Linear"` - Estimate parameters with regression calibration in a Cox model.
#' \item `"Quadratic"` - Estimate parameters with regression calibration with a quadratic exposure term included in a Cox model.
#' \item `"General_1"` - Estimate parameters with nonlinear regression calibration that applies a first-order Taylor series expansion to simplify the computations.
#' \item `"General_2"` - Estimate parameters with nonlinear regression calibration that applies a second-order Taylor series expansion to simplify the computations.
#' \item `"General_MC"` - Estimate parameters with nonlinear regression calibration that applies a Monte Carlo approximation.
#' }
#' @param mu_vec Vector corresponding to the conditional mean of the exposure variable given available data for each unit in `df1`. If provided, these values are used instead of estimating them (default in NULL).
#' @param Variance_vec Vector corresponding to the conditional variance of the exposure variable given available data for each unit in `df1`. If provided, these values are used instead of estimating them (default is NULL).
#' @param t_grid Vector of time points at which to estimate the conditional survival function.
#' @param z_grid Vector of exposure values at which to evaluate the conditional survival function.
#' @param C Number of Monte Carlo samples for estimating parameters and/or standard errors (default is 1000).
#' @param seed Random number seed (default is 235346).
#' @param x Value of covariates at which the conditional survival function is evaluated (default is 0).
#' @param Spline_Z_vec Vector of exposure values for generating spline basis functions (default is NULL).
#' @param Boundary.knots Boundary points for the spline basis functions (default is NULL).
#' @param Estimate_CI_Surv Logical indicator specifying whether a 95% point-wise Monte Carlo confidence interval is computed (default is TRUE).
#'
#' @return A list containing
#' \itemize{
#' \item `Estimates` - Estimated parameters \eqn{(\beta', \gamma')'} from the specified Cox regression model.
#' \item `Nuisance_Estimates` - Estimated (nuisance) parameters \eqn{\theta} (see Thomson and Huang (2025+)), and measurement error variance.
#' \item `Variance` - Estimated covariance matrix of the estimated parameters, with the first block corresponding to the estimated covariance matrix for the estimated Cox regression parameters.
#' \item `Baseline_Hazard` - Data frame containing the estimated cumulative baseline hazard function at each observed uncensored event time.
#' \item `SurvProb_zgrid_list` - Estimated conditional survival function, with each element of the list corresponding to a fixed value from `z_grid` at times `t_grid`.
#' \item `SurvProb_zgrid_CI_list` - 95% point-wise confidence interval of the conditional survival function for each value in `z_grid` at times `t_grid` (if requested).
#' \item `SurvProb_zgrid_Var_list` - Estimated variance of the conditional survival function for each value in `z_grid` at times `t_grid` (if requested).
#' }
#' @export
#'
#' @details
#' The specified Cox regression model is
#' \deqn{\lambda(t;Z,X) = \lambda_0(t) \exp\{ g(Z; \beta) + \gamma'X \},} where \eqn{Z} denotes an exposure variable (scalar), \eqn{X} are additional covariates, and \eqn{g(Z;\beta)} is a known function of \eqn{Z} with parameter (vector) \eqn{\beta}. Users can specify the functional form of \eqn{g(Z;\beta)} with the `Function_Type` specification.
#'
#' Rather than observing \eqn{Z}, the observed exposure \eqn{Q} follows a measurement error model subject to systematic bias influenced by subject-specific characteristics
#' \deqn{Q = a_0 + a_1 Z + a_2' V + \epsilon,} where \eqn{V = (X',R')'}, and \eqn{R} includes variables that influence the measurement error mechanism.
#'
#' We additionally suppose that we observe a biomarker \eqn{W} that adheres to a classical measurement error model
#' \deqn{W = Z + e.}
#'
#' Under the assumption that the event is "rare", the induced hazard model is
#' \deqn{\lambda(t;Q,V) = \lambda_0(t) \exp\{ \gamma' X \} H(\beta, Q, V; \theta),}
#' where \eqn{H(\beta, Q, V; \theta) = E_Z(\exp\{ g(Z;\beta) \} | Q, V)}, and \eqn{\theta} is a nuisance parameter vector related to \eqn{(a_0, a_1, a_2')'}.
#' \itemize{
#' \item When \eqn{g(Z;\beta) = \beta_1 Z} (`Function_Type = "Polynomial"`, `Method = "Linear"`), \eqn{H(\beta, Q, V; \theta) = \exp\{ \beta_1 E(Z | Q,V) + \beta_1^2 Var(Z | Q,V)/2 \}.}
#' \item When \eqn{g(Z;\beta) = \beta_1 Z + \beta_2 Z^2} (`Function_Type = "Polynomial"`, `Method = "Quadratic"`), \eqn{H(\beta, Q, V; \theta)} has a closed-form solution that leads to a procedure to estimate \eqn{\beta} and \eqn{\gamma} (see Huang and Prentice (2025)).
#' \item When \eqn{g(Z;\beta)} is a linear combination of spline basis functions or a polynomial of degree-3 (or higher), \eqn{H(\beta, Q, V; \theta)} has no analytical solution.
#'  \itemize{
#'  \item We allow users to approximate \eqn{H(\beta, Q, V; \theta)} with a Monte Carlo approximation (recommended) with the option `Method = "General_MC"` specification.
#'  \item We allow users to approximate \eqn{H(\beta, Q, V; \theta)} with a first- or second-order Taylor series expansion around \eqn{E(Z | Q,V)} with the options `Method = "General_1"` and `Method = "General_2"`, respectively.
#'  }
#' }
#'
#' @references
#' Thomson, T.J. and Huang Y. (2025+). Estimating Non-linear Exposure and Event Time Association in the Presence of Exposure Measurement Error. \emph{Submitted}
#'
#' Huang, Y. and Prentice, R.L. (2025). Biomarker-assisted reporting in nutritional epidemiology addressing measurement error in exposure-disease associations. \emph{Biostatistics}. \strong{26}(1), kxaf014.
#'
#' @author
#' Trevor J. Thomson (\email{tthomson@fredhutch.org})
#'
#' @note
#' All of the source code that `NL_RC` uses can be found on GitHub: \url{https://github.com/drtjt/NonLinearRC}
#'
#' @examples
#' data(df_Analysis)
#' data(df_Biomarker)
#' data(df_Replicate)
#'
#' # ----------------------------------------------------------------------
#' #        LINEAR REGRESSION CALIBRATION
#' # ----------------------------------------------------------------------
#'
#' Results_Polynomial_LRC = NL_RC(df1 = df_Analysis,
#'                                df2 = df_Biomarker,
#'                                df3 = df_Replicate,
#'                                df1_ID_index = 1,
#'                                df1_U_index = 3,
#'                                df1_delta_index = 4,
#'                                df1_Q_index = 7,
#'                                df1_X_index = 5,
#'                                df2_ID_index = 1,
#'                                df2_Response_index = 2,
#'                                df2_Covariates_index = c(3,5),
#'                                df3_ID_index = 1,
#'                                Function_Type = "Polynomial",
#'                                DF = 1,
#'                                Method = "Linear",
#'                                EstimateNuisance_Logi = T,
#'                                mu_vec = NULL,
#'                                Variance_vec = NULL,
#'                                t_grid = seq(0,20, length.out = 101),
#'                                z_grid = seq(-3, 3, length.out = 11),
#'                                C = 1000,
#'                                seed = 235346,
#'                                x = 0,
#'                                Spline_Z_vec = NULL,
#'                                Boundary.knots = NULL,
#'                                Estimate_CI_Surv = F) # True beta: -0.2, True gamma: 0.1
#'
#' # ----------------------------------------------------------------------
#' #         QUADRATIC REGRESSION CALIBRATION
#' # ----------------------------------------------------------------------
#'
#' Results_Polynomial_QRC = NL_RC(df1 = df_Analysis,
#'                                df2 = df_Biomarker,
#'                                df3 = df_Replicate,
#'                                df1_ID_index = 1,
#'                                df1_U_index = 3,
#'                                df1_delta_index = 4,
#'                                df1_Q_index = 7,
#'                                df1_X_index = 5,
#'                                df2_ID_index = 1,
#'                                df2_Response_index = 2,
#'                                df2_Covariates_index = c(3,5),
#'                                df3_ID_index = 1,
#'                                Function_Type = "Polynomial",
#'                                DF = 2,
#'                                Method = "Quadratic",
#'                                EstimateNuisance_Logi = T,
#'                                mu_vec = NULL,
#'                                Variance_vec = NULL,
#'                                t_grid = seq(0,20, length.out = 101),
#'                                z_grid = seq(-3, 3, length.out = 11),
#'                                C = 1000,
#'                                seed = 235346,
#'                                x = 0,
#'                                Spline_Z_vec = NULL,
#'                                Boundary.knots = NULL,
#'                                Estimate_CI_Surv = F)
#'
#' # ----------------------------------------------------------------------
#' #         CUBIC REGRESSION CALIBRATION
#' # ----------------------------------------------------------------------
#'
#' Results_Polynomial_CRC = NL_RC(df1 = df_Analysis,
#'                                df2 = df_Biomarker,
#'                                df3 = df_Replicate,
#'                                df1_ID_index = 1,
#'                                df1_U_index = 3,
#'                                df1_delta_index = 4,
#'                                df1_Q_index = 7,
#'                                df1_X_index = 5,
#'                                df2_ID_index = 1,
#'                                df2_Response_index = 2,
#'                                df2_Covariates_index = c(3,5),
#'                                df3_ID_index = 1,
#'                                Function_Type = "Polynomial",
#'                                DF = 3,
#'                                Method = "General_MC",
#'                                EstimateNuisance_Logi = T,
#'                                mu_vec = NULL,
#'                                Variance_vec = NULL,
#'                                t_grid = seq(0,20, length.out = 101),
#'                                z_grid = seq(-3, 3, length.out = 11),
#'                                C = 100,
#'                                seed = 235346,
#'                                x = 0,
#'                                Spline_Z_vec = NULL,
#'                                Boundary.knots = NULL,
#'                                Estimate_CI_Surv = F) # Takes a few minutes

NL_RC = function(df1,
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
    NuisanceParameters_2_EF = nleqslv::nleqslv(x = rep(0, dim(df2_use)[2]+1),
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
    NuisanceParameters_3_EF = nleqslv::nleqslv(x = 0,
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
  fit_EF = nleqslv::nleqslv(x = rep(0, DF+length(df1_X_index)),
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
      A_mat = nleqslv::nleqslv(x = c(fit_EF$x,
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

      # Obtain point-wise variance of survival probabilities
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
                  SurvProb_zgrid_Var_list = SurvProb_zgrid_Var_list)

  return(val_list)
}
