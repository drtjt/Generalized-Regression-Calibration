#' Simulated Analysis Cohort Dataset
#'
#' @description
#' A simulated \emph{analysis cohort} dataset to illustrate the use of `NL_RC`.
#'
#' @format A data frame with 10,000 rows and 7 variables, ordered by variable "Q1":
#'  \describe{
#'  \item{ID}{A variable that uniquely identifies each unit and (potentially) used to link units in `df_Biomarker` and/or `df_Replicate`.}
#'  \item{Time}{Simulated event time for each unit.}
#'  \item{U}{Observed event time that is subject to random right-censoring.}
#'  \item{delta}{Observed event time indicator: \eqn{I(T \leq C)}.}
#'  \item{X}{Additional covariate used to generate the event time.}
#'  \item{Z1}{True exposure variable for each unit. Note that this variable would not be observed in practice.}
#'  \item{Q1}{The observed (noisy) exposure variable for each unit.}
#'  }
#'
#' @source {R script to regenerate simulated data can be found at \url{https://github.com/drtjt/NonLinearRC}}
#'
#' @examples
#' data(df_Analysis)
#'
"df_Analysis"



#' Simulated Biomarker Dataset
#'
#' @description
#' A simulated \emph{biomarker} dataset to illustrate the use of `NL_RC`.
#'
#' @format A data frame with 500 rows and 5 variables:
#'  \describe{
#'  \item{ID}{A variable that uniquely identifies each unit and (potentially) used to link units in `df_Analysis` and/or `df_Replicate`.}
#'  \item{W1}{The observed biomarker for each unit.}
#'  \item{Q1}{The observed (noisy) exposure variable for each unit.}
#'  \item{Z}{True exposure variable for each unit. Note that this variable would not be observed in practice.}
#'  \item{X}{Additional covariate used in the measurement error model.}
#'  }
#'
#' @source {R script to regenerate simulated data can be found at \url{https://github.com/drtjt/NonLinearRC}}
#'
#' @examples
#' data(df_Biomarker)
#'
"df_Biomarker"


#' Simulated Replicate Dataset
#'
#' @description
#' A simulated \emph{replicate} dataset (in wide format) to illustrate the use of `NL_RC`.
#'
#' @format A data frame with 100 rows and 3 variables (two replicated observations of the biomarker for each unit):
#'  \describe{
#'  \item{ID}{A variable that uniquely identifies each unit and (potentially) used to link units in `df_Analysis` and/or `df_Biomarker`.}
#'  \item{W1}{The first measurement of the biomarker for each unit.}
#'  \item{W2}{The second measurement of the biomarker for each unit.}
#'  }
#'
#'
#' @source {R script to regenerate simulated data can be found at \url{https://github.com/drtjt/NonLinearRC}}
#'
#' @examples
#' data(df_Replicate)
#'
"df_Replicate"
