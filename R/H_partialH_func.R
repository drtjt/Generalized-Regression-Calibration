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

  #require(splines2)
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
      spline_data = splines2::bernsteinPoly(data$Q,
                                            degree = DF,
                                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "B-Spline"){
      spline_data = splines2::bSpline(data$Q,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "C-Spline"){
      spline_data = splines2::cSpline(data$Q,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "I-Spline"){
      spline_data = splines2::iSpline(data$Q,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "M-Spline"){
      spline_data = splines2::mSpline(data$Q,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "NC-Spline"){
      spline_data = splines2::naturalSpline(data$Q,
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
      spline_data = splines2::bernsteinPoly(data$mu,
                                            degree = DF,
                                            Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "B-Spline"){
      spline_data = splines2::bSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "C-Spline"){
      spline_data = splines2::cSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "I-Spline"){
      spline_data = splines2::iSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "M-Spline"){
      spline_data = splines2::mSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "NC-Spline"){
      spline_data = splines2::naturalSpline(data$mu,
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
      spline_data = splines2::bernsteinPoly(data$mu,
                                            degree = DF,
                                            Boundary.knots = Boundary.knots)
      # spline1_data = splines2::deriv(spline_data, derivs = 1)
      # spline2_data = splines2::deriv(spline_data, derivs = 2)
      spline1_data = splines2::bernsteinPoly(data$mu,
                                             degree = DF,
                                             Boundary.knots = Boundary.knots,
                                             derivs = 1)
      spline2_data = splines2::bernsteinPoly(data$mu,
                                             degree = DF,
                                             Boundary.knots = Boundary.knots,
                                             derivs = 2)
    }
    if (Function_Type == "B-Spline"){
      spline_data = splines2::bSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
      # spline1_data = splines2::deriv(spline_data, derivs = 1)
      # spline2_data = splines2::deriv(spline_data, derivs = 2)
      spline1_data = splines2::bSpline(data$mu,
                                       df = DF,
                                       Boundary.knots = Boundary.knots,
                                       derivs = 1)
      spline2_data = splines2::bSpline(data$mu,
                                       df = DF,
                                       Boundary.knots = Boundary.knots,
                                       derivs = 2)
    }
    if (Function_Type == "C-Spline"){
      spline_data = splines2::cSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
      # spline1_data = splines2::deriv(spline_data, derivs = 1)
      # spline2_data = splines2::deriv(spline_data, derivs = 2)
      spline_data = splines2::cSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
    }
    if (Function_Type == "I-Spline"){
      spline_data = splines2::iSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
      # spline1_data = splines2::deriv(spline_data, derivs = 1)
      # spline2_data = splines2::deriv(spline_data, derivs = 2)
      spline1_data = splines2::iSpline(data$mu,
                                       df = DF,
                                       Boundary.knots = Boundary.knots,
                                       derivs = 1)
      spline2_data = splines2::iSpline(data$mu,
                                       df = DF,
                                       Boundary.knots = Boundary.knots,
                                       derivs = 2)
    }
    if (Function_Type == "M-Spline"){
      spline_data = splines2::mSpline(data$mu,
                                      df = DF,
                                      Boundary.knots = Boundary.knots)
      # spline1_data = splines2::deriv(spline_data, derivs = 1)
      # spline2_data = splines2::deriv(spline_data, derivs = 2)
      spline1_data = splines2::mSpline(data$mu,
                                       df = DF,
                                       Boundary.knots = Boundary.knots,
                                       derivs = 1)
      spline2_data = splines2::mSpline(data$mu,
                                       df = DF,
                                       Boundary.knots = Boundary.knots,
                                       derivs = 2)
    }
    if (Function_Type == "NC-Spline"){
      spline_data = splines2::naturalSpline(data$mu,
                                            df = DF,
                                            Boundary.knots = Boundary.knots)
      # spline1_data = splines2::deriv(spline_data, derivs = 1)
      # spline2_data = splines2::deriv(spline_data, derivs = 2)
      spline1_data = splines2::naturalSpline(data$mu,
                                             df = DF,
                                             Boundary.knots = Boundary.knots,
                                             derivs = 1)
      spline2_data = splines2::naturalSpline(data$mu,
                                             df = DF,
                                             Boundary.knots = Boundary.knots,
                                             derivs = 2)
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
      Z_mc_vec = stats::rnorm(dim(data)[1],
                              mean = data$mu,
                              sd = sqrt(data$sigma2))

      # Initialize matrices
      g_mc_mat = matrix(NA, nrow = dim(data)[1], ncol = DF) # Initialize
      colnames(g_mc_mat) = as.character(1:DF)

      if (Function_Type == "Polynomial"){
        spline_mc_data = outer(Z_mc_vec, 1:DF, "^")
      }
      if (Function_Type == "Bernstein"){
        spline_mc_data = splines2::bernsteinPoly(Z_mc_vec,
                                                 degree = DF,
                                                 Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "B-Spline"){
        spline_mc_data = splines2::bSpline(Z_mc_vec,
                                           df = DF,
                                           Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "C-Spline"){
        spline_mc_data = splines2::cSpline(Z_mc_vec,
                                           df = DF,
                                           Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "I-Spline"){
        spline_mc_data = splines2::iSpline(Z_mc_vec,
                                           df = DF,
                                           Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "M-Spline"){
        spline_mc_data = splines2::mSpline(Z_mc_vec,
                                           df = DF,
                                           Boundary.knots = Boundary.knots)
      }
      if (Function_Type == "NC-Spline"){
        spline_mc_data = splines2::naturalSpline(Z_mc_vec,
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
