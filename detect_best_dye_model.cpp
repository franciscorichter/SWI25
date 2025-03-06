// detect_best_dye_model.cpp
// Compile with: sourceCpp("detect_best_dye_model.cpp")

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

//------------------------------------------------------------
// Helper function: allZeros
// Checks if all elements of a NumericVector are zero.
//------------------------------------------------------------
bool allZeros(const NumericVector& vec) {
  for (int i = 0; i < vec.size(); i++) {
    if (vec[i] != 0) return false;
  }
  return true;
}

//------------------------------------------------------------
// Function: sample_latent_once_cpp
// Purpose: For a given observed total vector N_obs, number of dyes K, 
//          number of time steps Tsteps, and current guess of q (q_guess),
//          simulate one latent trajectory.
//------------------------------------------------------------
// [[Rcpp::export]]
List sample_latent_once_cpp(const NumericVector& N_obs, int K, int Tsteps, const NumericVector& q_guess) {
  // subpop: (Tsteps+1) x K matrix storing latent subpopulation counts.
  // D: Tsteps x K matrix storing decays at each time step.
  NumericMatrix subpop(Tsteps + 1, K);
  NumericMatrix D(Tsteps, K);
  
  // Initialize: at t = 0, assume equal partition of N_obs[0]
  double init_val = N_obs[0] / K;
  for (int k = 0; k < K; k++) {
    subpop(0, k) = init_val;
  }
  
  // Loop over time steps
  for (int t = 0; t < Tsteps; t++) {
    // Determine total decays at time t: Dt = N_obs[t] - N_obs[t+1]
    int Dt = N_obs[t] - N_obs[t + 1];
    
    // Current subpopulation counts (row t of subpop)
    NumericVector Nk_minus(K);
    for (int k = 0; k < K; k++) {
      Nk_minus[k] = subpop(t, k);
    }
    
    // Compute "pressure" for each dye: press[k] = q_guess[k] * Nk_minus[k]
    NumericVector press(K);
    double s = 0.0;
    for (int k = 0; k < K; k++) {
      press[k] = q_guess[k] * Nk_minus[k];
      s += press[k];
    }
    if (s < 1e-12) {
      for (int k = 0; k < K; k++) {
        press[k] = 1.0 / K;
      }
    } else {
      for (int k = 0; k < K; k++) {
        press[k] /= s;
      }
    }
    
    // Sample decays across dyes from a multinomial distribution.
    std::vector<int> dtemp(K);
    std::vector<double> press_arr(K);
    for (int k = 0; k < K; k++) {
      press_arr[k] = press[k];
    }
    rmultinom(Dt, &press_arr[0], K, &dtemp[0]);
    
    // Clamp decays if necessary and update subpop for next time step.
    for (int k = 0; k < K; k++) {
      double dec = dtemp[k];
      if (dec > Nk_minus[k])
        dec = Nk_minus[k];
      D(t, k) = dec;
      subpop(t + 1, k) = Nk_minus[k] - dec;
    }
  }
  
  return List::create(Named("D") = D,
                      Named("subpop") = subpop);
}

//------------------------------------------------------------
// Function: log_f_function_cpp
// Purpose: Compute the log-density of the latent sample under the 
//          "target" model (log f).
//------------------------------------------------------------
double log_f_function_cpp(const NumericMatrix& D, const NumericMatrix& subpop, 
                          const NumericVector& q_guess, int Tsteps, int K) {
  double total = 0.0;
  for (int t = 0; t < Tsteps; t++) {
    for (int k = 0; k < K; k++) {
      double n_ktm1 = subpop(t, k);
      double d_kt = D(t, k);
      if (d_kt > n_ktm1)
        return -INFINITY;
      total += std::lgamma(n_ktm1 + 1) - std::lgamma(d_kt + 1) - std::lgamma(n_ktm1 - d_kt + 1);
      total += d_kt * log(std::max((double)q_guess[k], 1e-15));
      total += (n_ktm1 - d_kt) * log(std::max(1.0 - q_guess[k], 1e-15));
    }
  }
  return total;
}

//------------------------------------------------------------
// Function: log_g_function_cpp
// Purpose: Compute the log-density of the latent sample under the 
//          proposal distribution (log g).
//------------------------------------------------------------
double log_g_function_cpp(const NumericMatrix& D, const NumericMatrix& subpop, 
                          const NumericVector& q_guess, const NumericVector& N_obs, 
                          int Tsteps, int K) {
  double total = 0.0;
  for (int t = 0; t < Tsteps; t++) {
    int Dt = N_obs[t] - N_obs[t + 1];
    NumericVector Nk_minus(K);
    for (int k = 0; k < K; k++) {
      Nk_minus[k] = subpop(t, k);
    }
    NumericVector press(K);
    double s = 0.0;
    for (int k = 0; k < K; k++) {
      press[k] = q_guess[k] * Nk_minus[k];
      s += press[k];
    }
    if (s < 1e-12) {
      for (int k = 0; k < K; k++) {
        press[k] = 1.0 / K;
      }
    } else {
      for (int k = 0; k < K; k++) {
        press[k] /= s;
      }
    }
    total += std::lgamma(Dt + 1);
    for (int k = 0; k < K; k++) {
      total -= std::lgamma(D(t, k) + 1);
      total += D(t, k) * log(std::max((double)press[k], 1e-15));
    }
  }
  return total;
}

//------------------------------------------------------------
// Function: mcem_multi_dye_cpp
// Purpose: The main MCEM algorithm that repeatedly generates M latent samples,
//          computes weights, updates the q parameters, and records iteration history.
//------------------------------------------------------------
// [[Rcpp::export]]
List mcem_multi_dye_cpp(NumericVector N_obs, int K, int Tsteps, NumericVector q_init, 
                        int M = 200, int max_iter = 20, double tol = 1e-4) {
  NumericVector q_curr = clone(q_init);
  
  // History: store iteration number and corresponding q values.
  std::vector<NumericVector> history;
  NumericVector initHist(K + 1);
  initHist[0] = 0;
  for (int k = 0; k < K; k++) {
    initHist[k + 1] = q_curr[k];
  }
  history.push_back(initHist);
  
  for (int iter = 1; iter <= max_iter; iter++) {
    std::vector<List> sample_list(M);
    NumericVector logf_vec(M), logg_vec(M);
    
    // Generate M latent samples and compute log f and log g for each.
    for (int m = 0; m < M; m++) {
      List latent = sample_latent_once_cpp(N_obs, K, Tsteps, q_curr);
      NumericMatrix D = latent["D"];
      NumericMatrix subpop = latent["subpop"];
      double lf = log_f_function_cpp(D, subpop, q_curr, Tsteps, K);
      double lg = log_g_function_cpp(D, subpop, q_curr, N_obs, Tsteps, K);
      logf_vec[m] = lf;
      logg_vec[m] = lg;
      sample_list[m] = latent;
    }
    
    // Compute importance weights: w = exp(logf - logg)
    NumericVector w_raw(M);
    double w_sum = 0.0;
    for (int m = 0; m < M; m++) {
      w_raw[m] = exp(logf_vec[m] - logg_vec[m]);
      w_sum += w_raw[m];
    }
    if (w_sum < 1e-12) {
      for (int m = 0; m < M; m++) {
        w_raw[m] = 1.0 / M;
      }
      w_sum = 1.0;
    }
    
    // Update q: accumulate numerator and denominator over samples and time steps t = 1,...,Tsteps-1.
    NumericVector num(K, 0.0), den(K, 0.0);
    for (int m = 0; m < M; m++) {
      List latent = sample_list[m];
      NumericMatrix D = latent["D"];
      NumericMatrix subpop = latent["subpop"];
      for (int t = 1; t < Tsteps; t++) {
        for (int k = 0; k < K; k++) {
          num[k] += w_raw[m] * D(t, k);
          den[k] += w_raw[m] * subpop(t, k);
        }
      }
    }
    NumericVector q_new(K);
    for (int k = 0; k < K; k++) {
      q_new[k] = (den[k] < 1e-12) ? 0.0 : num[k] / den[k];
    }
    
    // Check convergence: maximum absolute difference between q_new and q_curr.
    double diff_q = 0.0;
    for (int k = 0; k < K; k++) {
      diff_q = std::max(diff_q, fabs(q_new[k] - q_curr[k]));
    }
    q_curr = q_new;
    
    // Record iteration history.
    NumericVector histRow(K + 1);
    histRow[0] = iter;
    for (int k = 0; k < K; k++) {
      histRow[k + 1] = q_curr[k];
    }
    history.push_back(histRow);
    
    Rcpp::Rcout << "Iter " << iter << ": q = ";
    for (int k = 0; k < K; k++) {
      Rcpp::Rcout << q_curr[k] << " ";
    }
    Rcpp::Rcout << "(diff = " << diff_q << ")\n";
    
    if (diff_q < tol) {
      Rcpp::Rcout << "MCEM converged.\n";
      break;
    }
  }
  
  // Convert history vector to an R list.
  List historyList(history.size());
  for (size_t i = 0; i < history.size(); i++) {
    historyList[i] = history[i];
  }
  
  return List::create(Named("q_hat") = q_curr,
                      Named("history") = historyList);
}

//------------------------------------------------------------
// Function: detect_best_dye_model_cpp
// Purpose: Run the MCEM algorithm for K = 1 and K = 2 models,
//          and choose the best model based on a simple heuristic.
//          If the estimated decay rate for the second dye is below a threshold,
//          then the one-dye model is preferred.
//------------------------------------------------------------
// [[Rcpp::export]]
List detect_best_dye_model_cpp(NumericVector N_obs, int Tsteps, 
                               double q_init1, NumericVector q_init2,
                               int M = 200, int max_iter = 20, double tol = 1e-4,
                               double threshold = 1e-3) {
  // If the observed data is all zeros, return a one-dye model with zero rate.
  if (allZeros(N_obs)) {
    NumericVector q_zero(1, 0.0);
    return List::create(Named("best_K") = 1,
                        Named("res1") = List::create(Named("q_hat") = q_zero),
                        Named("res2") = R_NilValue);
  }
  
  // Run MCEM for one-dye model.
  NumericVector q_init_1(1);
  q_init_1[0] = q_init1;
  List res1 = mcem_multi_dye_cpp(N_obs, 1, Tsteps, q_init_1, M, max_iter, tol);
  
  // Run MCEM for two-dye model.
  List res2 = mcem_multi_dye_cpp(N_obs, 2, Tsteps, q_init2, M, max_iter, tol);
  NumericVector q_hat_2 = res2["q_hat"]; // q_hat for K = 2
  
  // Simple heuristic: if the second dye's estimated rate is below threshold, choose K = 1.
  int best_K = (q_hat_2[1] < threshold) ? 1 : 2;
  
  return List::create(Named("best_K") = best_K,
                      Named("res1") = res1,
                      Named("res2") = res2);
}

