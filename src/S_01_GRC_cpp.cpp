#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
List S_01_GRC_cpp(NumericVector H_vec,
                  NumericMatrix partialH_mat,
                  NumericVector gamma_vec,
                  NumericMatrix X_mat,
                  int m,
                  IntegerVector IndexEnterRiskSet,
                  IntegerVector IndexLeaveRiskSet) {

  int n = H_vec.size();
  int p1 = partialH_mat.ncol();
  int p2 = X_mat.ncol();
  int M1 = m + 1;

  // Initialize vectors and matrices
  NumericVector R0(M1); // initialize
  NumericVector R0A(M1); // initialize
  NumericVector R0B(M1); // initialize
  NumericMatrix R1(M1, p1 + p2); // initialize
  NumericMatrix R1A(M1, p1 + p2); // initialize
  NumericMatrix R1B(M1, p1 + p2); // initialize

  // Run a for loop over n to update R0A, R0B, R1A, R1B
  for (int i = 0; i < n; ++i) {

    int a = IndexEnterRiskSet[i];
    int b = IndexLeaveRiskSet[i];

    // Compute gamma'X
    double gammaX = 0.0;
    for (int j = 0; j < p2; ++j) {
      gammaX += gamma_vec[j] * X_mat(i, j);
    }

    double e_gammaX_i = std::exp(gammaX);
    double R0_i = e_gammaX_i * H_vec[i];

    // Compute R1beta_i
    // NumericVector R1beta_i(p1);
    // for (int j = 0; j < p1; ++j) {
    //   R1beta_i[j] = e_gammaX_i * partialH_mat(i, j);
    // }
    NumericVector R1beta_i = e_gammaX_i * partialH_mat(i,_); // can this replace the above?

    // Compute R1gamma_i
    // NumericVector R1gamma_i(p2);
    // for (int j = 0; j < p2; ++j) {
    //   R1gamma_i[j] = e_gammaX_i * H_vec[i] * X_mat(i, j);
    // }
    NumericVector R1gamma_i = e_gammaX_i * H_vec[i] * X_mat(i,_);

    // Combine R1beta_i and R1gamma_i into R1_i
    NumericVector R1_i(p1 + p2);
    for (int j = 0; j < p1; ++j) {
      R1_i[j] = R1beta_i[j];
    }
    for (int j = 0; j < p2; ++j) {
      R1_i[p1 + j] = R1gamma_i[j];
    }

    // Update R0A and R0B
    R0A[a] += R0_i;
    R0B[b] += R0_i;

    // Update R1A[a,] and R1B[b,]
    R1A(a,_) =  R1A(a,_) + R1_i;
    R1B(b,_) =  R1B(b,_) + R1_i;
  }

  // Update "R0[1]" and "R1[1,]"
  R0(0) = R0A(0) - R0B(0);
  R1(0,_) = R1A(0,_) - R1B(0,_);

  // Run a for loop over 2:M to update "R0[2:M1]" and "R1[2:M1,]"
  for (int j = 1 ; j < M1 ; ++j) {

    // Update "R0[j+1]" and "R1[j+1,]"
    R0(j) = R0(j-1) + R0A(j) - R0B(j);
    R1(j,_) = R1(j-1,_) + R1A(j,_) - R1B(j,_);
  }

  return List::create(Named("r0") = R0,
                      Named("r1") = R1);
}
