/*==================================================================
THIS IS THE GLOBAL FUNCTION FILE FOR BOTH FE AND VF SCHEME 
===================================================================*/
#ifndef GLOBAL_FUNCTION_H
#define GLOBAL_FUNCTION_H

// --- Implementation of Exact Solution--
double exact_solution(double x, double t, double eta) {
    double lambda = M_PI * M_PI + (eta * eta) / 4.0;
    double steady = std::exp(eta * x);
    double transient = std::exp(0.5 * eta * x) * std::exp(-lambda * t) * std::sin(M_PI * x);
    return steady + transient;
}

// ------Implementation stationnaire Exact Solution  for (Q13)----
double exact_solution_stationary(double x, double eta) {
    return std::exp(eta * x);
}
#endif