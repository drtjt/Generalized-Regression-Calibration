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

  #require(splines2)
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
    pred.Z = suppressWarnings(as.numeric(stats::predict(spline_data,
                                                        newx = Z)))
  }
  Lambda_t = Lambda_0t * exp( as.numeric(pred.Z %*% Beta_vec) +
                                as.numeric(t(X) %*% gamma_vec) )

  Surv_Probs = exp(-Lambda_t)

  val = data.frame(Times = Times,
                   Surv_Probs = Surv_Probs)

  return(val)
}
