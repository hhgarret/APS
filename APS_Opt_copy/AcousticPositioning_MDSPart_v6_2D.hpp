#ifndef ACOUSTICPOSITIONING_MDSPART_V6_2D_HPP
#define ACOUSTICPOSITIONING_MDSPART_V6_2D_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <span>
#include <armadillo>
#include <nlopt.hpp>

// Structure to track historical values of x and gradients
struct OptimizationData_MDS
{
    std::vector<std::vector<double>> x_history;    // Stores x at each iteration
    std::vector<std::vector<double>> grad_history; // Stores gradients at each iteration
};

// Function declarations
arma::mat calculate_absolute_positions(const arma::mat &d_calculated, const int &rotation_type);
double objective_function_MDS(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
arma::mat run_MDSOptimization(OptimizationData_MDS &opt_data_MDS);

#endif