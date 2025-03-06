// mcem_single_dye.cpp
// Compile with: sourceCpp("mcem_single_dye.cpp")
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double estimate_single_dye_q_cpp(NumericVector N_obs, int Tsteps) {
  double num = 0.0, den = 0.0;
  // Loop from t = 0 to Tsteps-2 because we compare N_obs[t] with N_obs[t+1]
  for (int t = 0; t < Tsteps - 1; t++) {
    num += N_obs[t] - N_obs[t+1];
    den += N_obs[t];
  }
  // Avoid division by zero
  if (den < 1e-12)
    return 0.0;
  return num / den;
}
