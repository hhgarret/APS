#ifndef ACOUSTICPOSITIONING_MLPART_V2_2D_HPP
#define ACOUSTICPOSITIONING_MLPART_V2_2D_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <span>
#include <armadillo>
#include <nlopt.hpp>

// Define Speed of sound in m/s
//const double c = 343;

// Structure to track historical values of x and gradients
struct OptimizationData {
    std::vector<std::vector<double>> x_history;  // Stores x at each iteration
    std::vector<std::vector<double>> grad_history; // Stores gradients at each iteration
};

// Function declarations
double objective_function(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
void run_LMOptimization(OptimizationData& opt_data);
//arma::mat calculate_absolute_positions(const arma::mat& d_calculated);

#endif // ACOUSTICPOSITIONING_LMPART_HPP