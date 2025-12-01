o#include "AcousticPositioning_MDSPart_v6_2D.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <span>
#include <armadillo>
#include <nlopt.hpp>

// 2D version, no z at all.

// v6 update: p0,p1 and p2 are not fixed anymore. In previous version that they were fixed, we had
// big errors on the z-axis because we fixed p0,p1 and p2 to be on z=0 plane.

// v3 update: This version uses nlopt optimization to find the closest coordinates to distances
// used MDS geometric method with reflected points on z=0 plane as initial points

// MDS using geometric method
arma::mat calculate_absolute_positions(const arma::mat &d_calculated, const int &rotation_type)
{ // load calculated distance

  // Number of sensors
  int sensor_num = d_calculated.n_rows;

  // initialize calculated coordination matrix, two possiilities where the p3
  // has either positive or negative value on the z-axis.
  arma::mat coordinates_calculated_positive_2D = arma::zeros(sensor_num, 2);
  arma::mat coordinates_calculated_negative_2D = arma::zeros(sensor_num, 2);

  // Add the coordinates for our fixed 3 points
  // P0=[0,0], Pa=[A,0]; Pb=[B,C]
  // Pa and Pb are the nodes closest to P0.
  arma::uvec closest_neighbors_idx = sort_index(d_calculated.row(0));
  // The first point is P0, so Pa is the second node with two possibilities
  arma::uword Pa_idx = closest_neighbors_idx(1);
  arma::uword Pb_idx = closest_neighbors_idx(2);
  coordinates_calculated_positive_2D.row(0) = {0, 0};
  coordinates_calculated_negative_2D.row(0) = {0, 0};
  coordinates_calculated_positive_2D.row(Pa_idx) = {d_calculated(0, Pa_idx), 0};
  coordinates_calculated_negative_2D.row(Pa_idx) = {d_calculated(0, Pa_idx), 0};
  // Calculate B and C in Pb=[Pb_x,Pb_y]
  double Pb_x = (pow(d_calculated(0, Pb_idx), 2) - pow(d_calculated(Pa_idx, Pb_idx), 2) + pow(d_calculated(Pa_idx, 0), 2)) / (2 * d_calculated(Pa_idx, 0));
  double Pb_y_squared = pow(d_calculated(0, Pb_idx), 2) - pow(Pb_x, 2);
  double Pb_y_positive = sqrt(Pb_y_squared);
  double Pb_y_negative = -Pb_y_positive;
  coordinates_calculated_positive_2D.row(Pb_idx) = {Pb_x, Pb_y_positive};
  coordinates_calculated_negative_2D.row(Pb_idx) = {Pb_x, Pb_y_negative};

  // Next, we start form the next node and consider two possibilities: first
  // for Pb_y_positive and then for Pb_y_negative.At the
  // end we will have two possibilities for final coordinates.
  // Positive
  arma::mat A_positive(3, 2);
  arma::colvec b_positive(3);
  // arma::mat A_positive(2, 2);
  // arma::colvec b_positive(2);
  arma::rowvec r1_positive = coordinates_calculated_positive_2D.row(0);
  arma::rowvec r2_positive = coordinates_calculated_positive_2D.row(Pa_idx);
  arma::rowvec r3_positive = coordinates_calculated_positive_2D.row(Pb_idx);
  arma::colvec r4_positive(2);
  // Negative
  arma::mat A_negative(3, 2);
  arma::colvec b_negative(3);
  arma::rowvec r1_negative = coordinates_calculated_negative_2D.row(0);
  arma::rowvec r2_negative = coordinates_calculated_negative_2D.row(Pa_idx);
  arma::rowvec r3_negative = coordinates_calculated_negative_2D.row(Pb_idx);
  arma::colvec r4_negative(2);

  for (int i = 1; i < sensor_num; i++)
  {
    // Making sure the other two nodes, Pa and Pb, are not considered again.
    if (i != Pa_idx && i != Pb_idx)
    {
      //  Positive
      A_positive.row(0) = 2 * (r1_positive - r2_positive);
      A_positive.row(1) = 2 * (r1_positive - r3_positive);
      A_positive.row(2) = 2 * (r2_positive - r3_positive);
      b_positive(0) = dot(r1_positive, r1_positive) - dot(r2_positive, r2_positive) + pow(d_calculated(i, Pa_idx), 2) - pow(d_calculated(i, 0), 2);
      b_positive(1) = dot(r1_positive, r1_positive) - dot(r3_positive, r3_positive) + pow(d_calculated(i, Pb_idx), 2) - pow(d_calculated(i, 0), 2);
      b_positive(2) = dot(r2_positive, r2_positive) - dot(r3_positive, r3_positive) + pow(d_calculated(i, Pb_idx), 2) - pow(d_calculated(i, Pa_idx), 2);
      // r4_positive = solve(A_positive,b_positive);
      r4_positive = pinv(A_positive) * b_positive;
      // Negative
      A_negative.row(0) = 2 * (r1_negative - r2_negative);
      A_negative.row(1) = 2 * (r1_negative - r3_negative);
      A_negative.row(2) = 2 * (r2_negative - r3_negative);
      b_negative(0) = dot(r1_negative, r1_negative) - dot(r2_negative, r2_negative) + pow(d_calculated(i, Pa_idx), 2) - pow(d_calculated(i, 0), 2);
      b_negative(1) = dot(r1_negative, r1_negative) - dot(r3_negative, r3_negative) + pow(d_calculated(i, Pb_idx), 2) - pow(d_calculated(i, 0), 2);
      b_negative(2) = dot(r2_negative, r2_negative) - dot(r3_negative, r3_negative) + pow(d_calculated(i, Pb_idx), 2) - pow(d_calculated(i, Pa_idx), 2);
      // r4_negative = solve(A_negative,b_negative);
      r4_negative = pinv(A_negative) * b_negative;
      // Update the coordination matrix
      coordinates_calculated_positive_2D.row(i) = r4_positive.t();
      coordinates_calculated_negative_2D.row(i) = r4_negative.t();
    }
  }

  arma::mat calculated_absolute_positions;

  switch (rotation_type)
  {
  case 1:
    calculated_absolute_positions = coordinates_calculated_positive_2D;
    break;
  case 2:
    calculated_absolute_positions = coordinates_calculated_negative_2D;
    break;
  default:
    // Handle unexpected rotation_type values, if necessary
    calculated_absolute_positions = coordinates_calculated_positive_2D;
    break;
  }




  // Code made by Noah to perform MDS. Seems to provide a more accurate first guess??
  int numdimensions = 2;
  arma::mat distanceMatrix = d_calculated;
  distanceMatrix = pow(distanceMatrix, numdimensions);
  int N_channels = sensor_num;
  arma::mat H = arma::eye(N_channels, N_channels) - (1.0 / N_channels) * arma::ones(N_channels) * arma::trans(arma::ones(N_channels));
  arma::mat gramianMatrix = -0.5 * H * distanceMatrix * H;
  arma::vec eigenValues;
  arma::mat eigenVectors;
  arma::eig_sym(eigenValues, eigenVectors, gramianMatrix);
  eigenValues = eigenValues.tail(numdimensions);
  eigenVectors = eigenVectors.tail_cols(numdimensions);
  eigenVectors = arma::normalise(eigenVectors, 2, 0);

  arma::mat coordinateMatrix = arma::trans(eigenVectors * sqrt(arma::diagmat(eigenValues)));
  //arma::mat coordinateMatrix = arma::trans(eigenVectors * (arma::diagmat(eigenValues)));
  coordinateMatrix = arma::trans(coordinateMatrix);


  return coordinateMatrix;
  return calculated_absolute_positions;
}

// Objective function
// Note: If any changes are made to x and grad inside the objective_function, those changes will be reflected outside
// the function. This is because they are passed by reference.
// having const for x variable ensures that the function does not unintentionally alter x.
double objective_function_MDS(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{

  auto *data = static_cast<OptimizationData_MDS *>(my_func_data);

  // Track the current x values
  data->x_history.push_back(x);

  // Convert from std::vector to Armadillo vector
  arma::rowvec x_arma = arma::conv_to<arma::rowvec>::from(x);

  // Total number of optimization variables
  int opt_var_num = x_arma.n_elem;
  int sensor_num = (opt_var_num / 2);
  int n_rows = sqrt(sensor_num);
  arma::uword n_rows_u = static_cast<arma::uword>(n_rows);

  // Convert vector of variables to matrix
  arma::mat x_mat = x_arma.as_row();
  // Extract sensor positions
  arma::mat coordination_matrix = trans(reshape(x_mat, 2, sensor_num));

  // Load calculated relative distance from file
  arma::mat d_calculated;
  d_calculated.load("APS_relative_distance_calculated.txt");

  // Compute the objective function value
  // Define t_hat elements
  arma::mat d_hat = arma::zeros(sensor_num, sensor_num);
  for (int i = 0; i < sensor_num; ++i)
  {
    for (int j = 0; j < sensor_num; ++j)
    {
      if (i != j)
      {
        arma::rowvec Pij = coordination_matrix.row(i) - coordination_matrix.row(j);
        d_hat(i, j) = arma::norm(Pij);
      }
    }
  }
  //arma::mat residual_mat = pow(d_hat, 2) - pow(d_calculated, 1);
  //CHANGE BY HARLEY: WHY POW 2 INSTEAD OF POW 1?
  arma::mat residual_mat = pow(d_hat, 2) - pow(d_calculated, 2);
  double objective = accu(square(residual_mat));

  // If gradients are requested, compute them
  if (!grad.empty())
  {
    // Calculate gradient matrix
    // Define gradient vector, they are g[0] to g[3*(N-3)-1] each related
    // to the gradient with respect to each optimization variable.
    // N-3 sensor node positions (3 sensors are fixed)
    // First, we define a matrix form for each gradient and call it element
    // by element. Then we create the gradient vector from the matrix.
    // Note that to simplify the calculations, we calculate gradient for all
    // N sensor nodes. We then remove the values related to the 3 fixed nodes.
    // ---------- Gradient Matrix Calculations
    // ---- Pi
    arma::mat gradient_mat_Pi_x = arma::zeros(sensor_num, sensor_num);
    arma::mat gradient_mat_Pi_y = arma::zeros(sensor_num, sensor_num);
    for (int i = 0; i < sensor_num; ++i)
    {
      for (int j = 0; j < sensor_num; ++j)
      {
        if (i != j)
        {
          arma::rowvec Pij = coordination_matrix.row(i) - coordination_matrix.row(j);
          double dij = arma::norm(Pij);
          // ---- Pi
          gradient_mat_Pi_x(i, j) = Pij(0) * (pow(dij, 2) - pow(d_calculated(i, j), 2));
          gradient_mat_Pi_y(i, j) = Pij(1) * (pow(dij, 2) - pow(d_calculated(i, j), 2));
        }
      }
    }

    // ---------- Gradient Vector
    arma::rowvec grad_arma = arma::zeros(1, 2 * (sensor_num));
    for (int k = 1; k <= sensor_num; ++k)
    {
      // Gradient with respect to x-axis of sensor k position
      grad_arma(2 * k - 2) = 4 * sum(gradient_mat_Pi_x.row(k - 1));
      // Gradient with respect to y-axis of sensor k position
      grad_arma(2 * k - 1) = 4 * sum(gradient_mat_Pi_y.row(k - 1));
    }

    // convert from armadillo vector to std vector
    grad = arma::conv_to<std::vector<double>>::from(grad_arma);
    // Save the gradient for this iteration
    data->grad_history.push_back(grad);
  }
  else
  {
    // If gradient is not requested, push a zero-filled vector
    data->grad_history.push_back(std::vector<double>(x.size(), 0.0));
  }
  return objective;
}

arma::mat run_MDSOptimization(OptimizationData_MDS &opt_data_MDS)
{
  double minf;

  // Load relative distance matrix
  arma::mat d_calculated;
  d_calculated.load("APS_relative_distance_calculated.txt");
  int sensor_num = d_calculated.n_rows;
  int opt_var_num = 2 * (sensor_num);
  int n_rows = sqrt(sensor_num);
  arma::uword n_rows_u = static_cast<arma::uword>(n_rows);

  // Use MDS as initial guess
  arma::mat coormat = calculate_absolute_positions(d_calculated, 2);

  arma::rowvec opt_var_vec = coormat.as_row();

  // Convert from Armadillo vector to std vector
  std::vector<double> x = arma::conv_to<std::vector<double>>::from(opt_var_vec);

  // Create the optimizer
  // LN_xx: gradient free
  // LD_xx: with gradient  //LD_LBFGS   LD_SLSQP
  //nlopt::opt optimizer(nlopt::LD_LBFGS, opt_var_num);
  //HARLEY CHANGE: trying different optimizers
  nlopt::opt optimizer(nlopt::LD_MMA, opt_var_num);
  // Set the objective function
  optimizer.set_min_objective(objective_function_MDS, &opt_data_MDS);
  // Set optimization parameters
  /*
  optimizer.set_stopval(1e-30);
  optimizer.set_xtol_rel(1e-30);
  optimizer.set_ftol_rel(1e-30);
  optimizer.set_maxeval(10000);
  optimizer.set_initial_step(100);
  */

  //optimizer.set_stopval(-std::numeric_limits<double>::infinity());
  optimizer.set_stopval(.00001);
  optimizer.set_xtol_rel(1e-20);
  optimizer.set_ftol_rel(1e-20);
  optimizer.set_maxeval(10000);
  optimizer.set_initial_step(10);

  arma::mat calculated_coordination_matrix;

  try
  {
    nlopt::result result = optimizer.optimize(x, minf);
    // Extract optimal solution
    // Convert the solution from std::vector to Armadillo vector
    arma::rowvec opt_var_vec_result = arma::conv_to<arma::rowvec>::from(x);
    // Convert vector of variables to matrix
    // Total number of optimization variables
    int opt_var_num = opt_var_vec_result.n_elem;
    int sensor_num = (opt_var_num / 2);
    // Convert vector of variables to matrix
    arma::mat x_mat = opt_var_vec_result.as_row();
    // Extract sensor positions
    arma::mat optimized_coordination_matrix = trans(reshape(x_mat, 2, sensor_num));

    calculated_coordination_matrix = optimized_coordination_matrix;
    //calculated_coordination_matrix = coormat; //test by HARLEY to experiment w/o MDS optimization

    // Print results to terminal
    std::cout << "MDS Min objective function: " << minf << std::endl;
    std::cout << "MDS Total number of iterations: " << optimizer.get_numevals() << std::endl;
    std::cout << "MDS Optimization result: " << result << std::endl;
  }
  catch (const nlopt::roundoff_limited &e)
  {
    std::cerr << "MDS Optimization stopped due to roundoff errors: " << e.what() << std::endl;
  }
  catch (const nlopt::forced_stop &e)
  {
    std::cerr << "MDS Optimization was forcibly stopped: " << e.what() << std::endl;
  }
  catch (const std::invalid_argument &e)
  {
    std::cerr << "MDS Invalid argument: " << e.what() << std::endl;
  }
  catch (const std::bad_alloc &e)
  {
    std::cerr << "MDS Memory allocation error: " << e.what() << std::endl;
  }
  catch (const std::runtime_error &e)
  {
    std::cerr << "MDS Runtime error: " << e.what() << std::endl;

    // Diagnostics: Get and print last known optimum
    double last_opt_f = optimizer.last_optimum_value();
    nlopt::result last_result = optimizer.last_optimize_result();
    // Get the last x vector
    std::vector<double> last_x = x; // x contains the last values used by optimizer

    // Convert the last solution to get coordinates and wind vector
    arma::rowvec opt_var_vec_result = arma::conv_to<arma::rowvec>::from(last_x);
    arma::mat x_mat = opt_var_vec_result.as_row();
    arma::mat last_coordination_matrix = trans(reshape(x_mat, 2, sensor_num));

    // Print diagnostics
    std::cerr << "MDS Last known objective value: " << last_opt_f << std::endl;
    std::cerr << "MDS Last optimization result code: " << last_result << std::endl;
    std::cout << "MDS Total number of iterations performed: " << optimizer.get_numevals() << std::endl;
    

    calculated_coordination_matrix = last_coordination_matrix;
  }
  catch (std::exception &e)
  {
    std::cerr << "MDS An unknown error occurred: " << e.what() << std::endl;
    // Diagnostics: Get and print last known optimum
    double last_opt_f = optimizer.last_optimum_value();
    nlopt::result last_result = optimizer.last_optimize_result();
    std::cerr << "MDS Last known objective value: " << last_opt_f << std::endl;
    std::cerr << "MDS Last optimization result code: " << last_result << std::endl;
    std::cout << "MDS Total number of iterations performed: " << optimizer.get_numevals() << std::endl;
  }

  return calculated_coordination_matrix;
}