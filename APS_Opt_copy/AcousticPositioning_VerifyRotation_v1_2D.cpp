#include <iostream>
#include <vector>
#include <cmath>
#include <span>
#include <armadillo>
#include <nlopt.hpp> 

// Structure to track historical values of x and gradients
struct OptimizationData {
    std::vector<std::vector<double>> x_history;  // Stores x at each iteration
    std::vector<std::vector<double>> grad_history; // Stores gradients at each iteration
};

// Calculate rotated coordinates from rotation angle and the translation point
arma::mat calculate_rotated_coordinates(const arma::mat& input_coordination_matrix, double rotation_angle, const arma::rowvec& translation_point) {

  // Rotation angle in radians  
  double theta = rotation_angle;

  
  // Define 2D rotation matrix R
  arma::mat R = arma::zeros(2,2);
  R(0,0) = cos(theta);
  R(0,1) = -sin(theta);
  R(1,0) = sin(theta);
  R(1,1) = cos(theta);

  // Calculate the rotated coordinates
  int sensor_num = input_coordination_matrix.n_rows;
  arma::mat rotated_coordination_matrix = arma::zeros(sensor_num,2);
  for (int i = 0; i < sensor_num; ++i) {
    arma::colvec temp_ij = R*(input_coordination_matrix.row(i).t()-translation_point.t());
    rotated_coordination_matrix.row(i) = temp_ij.t();
  }
  
  return rotated_coordination_matrix;
}

// Objective function 
// Note: If any changes are made to x and grad inside the objective_function, those changes will be reflected outside 
// the function. This is because they are passed by reference. 
// having const for x variable ensures that the function does not unintentionally alter x.

double objective_function(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
  auto* data = static_cast<OptimizationData*>(my_func_data);

  // Track the current x values
  data->x_history.push_back(x);

  // Convert from std::vector to Armadillo vector
  arma::rowvec x_arma = arma::conv_to<arma::rowvec>::from(x);
  
  // Rename the variables
  double theta = std::fmod(x_arma(0), 360.0);  // Rotation angle
  double x_a = x_arma(1);    // Translation in x
  double y_a = x_arma(2);    // Translation in y
  arma::colvec a = {{x_a}, {y_a}};

  // Load calculated (optimized) coordination matrix
  arma::mat calculated_coordination_matrix;
  calculated_coordination_matrix.load("APS_coordinates_optimized.txt");

  // Load actual (original) coordination matrix
  arma::mat actual_coordination_matrix;
  actual_coordination_matrix.load("APS_coordinates_original.txt");

  // Total number of sensors
  int sensor_num = actual_coordination_matrix.n_rows;

  
  // Define the 2D rotation matrix R
  arma::mat R = arma::zeros(2,2);
  R(0,0) = cos(theta);
  R(0,1) = -sin(theta);
  R(1,0) = sin(theta);
  R(1,1) = cos(theta);

  // Compute the objective function value
  arma::rowvec residual_vect = arma::zeros(1, sensor_num);
  for (int i = 0; i < sensor_num; ++i) {
    arma::colvec eij = actual_coordination_matrix.row(i).t() - R * (calculated_coordination_matrix.row(i).t() - a);
    residual_vect(i) = arma::dot(eij, eij);
  }
  double objective = 0.5 * arma::accu(residual_vect);
  // std::cout << "Current objective function value = " << objective << std::endl;

  // If gradients are requested, compute them
  if (!grad.empty()) {
    // Calculate gradient matrix
    // Define gradient vector, they are g[0] to g[5] each related
    // to the gradient with respect to each optimization variable.
    // First, we define a vector form for each gradient and call it element
    // by element for each sensor. Then we accumulate all of them to get one
    // single value for each gradient. We wil then create the main gradient vector.
  
    // ---------- Individual Gradients
    // ---- Rotation Matrix Gradients
      arma::mat R_gradient_theta = arma::zeros(2,2);
      R_gradient_theta(0,0) = -sin(theta);
      R_gradient_theta(0,1) = -cos(theta);
      R_gradient_theta(1,0) = cos(theta);
      R_gradient_theta(1,1) = -sin(theta);

      arma::rowvec gradient_vect_theta = arma::zeros(1, sensor_num);
      arma::rowvec gradient_vect_xa = arma::zeros(1, sensor_num);
      arma::rowvec gradient_vect_ya = arma::zeros(1, sensor_num);

      for (int i = 0; i < sensor_num; ++i) {
          arma::colvec eij = actual_coordination_matrix.row(i).t() - R * (calculated_coordination_matrix.row(i).t() - a);
          arma::colvec grad_eij_theta = -R_gradient_theta * (calculated_coordination_matrix.row(i).t() - a);
          arma::colvec grad_eij_xa = R.col(0);
          arma::colvec grad_eij_ya = R.col(1);
          gradient_vect_theta(i) = arma::dot(grad_eij_theta, eij) + arma::dot(eij, grad_eij_theta);
          gradient_vect_xa(i) = arma::dot(grad_eij_xa, eij) + arma::dot(eij, grad_eij_xa);
          gradient_vect_ya(i) = arma::dot(grad_eij_ya, eij) + arma::dot(eij, grad_eij_ya);
      }
      double gradient_theta = 0.5 * arma::accu(gradient_vect_theta);
      double gradient_xa = 0.5 * arma::accu(gradient_vect_xa);
      double gradient_ya = 0.5 * arma::accu(gradient_vect_ya);

      // Main gradient vector
      arma::rowvec grad_arma = {gradient_theta, gradient_xa, gradient_ya};

      // Convert from Armadillo vector to std vector
      grad = arma::conv_to<std::vector<double>>::from(grad_arma);

      // Save the gradient for this iteration
      data->grad_history.push_back(grad);
  } else {
      // If gradient is not requested, push a zero-filled vector
      data->grad_history.push_back(std::vector<double>(x.size(), 0.0));
  }
  return objective;
}

int main() {
    OptimizationData opt_data;  // Initialize tracking structure
    double minf;

    // Initialize optimizer by considering all zero values
    // Load calculated (optimized) coordination matrix
    arma::mat calculated_coordination_matrix;
    calculated_coordination_matrix.load("APS_coordinates_optimized.txt");
    // Set initial translation point to be P0
    arma::rowvec opt_var_vec = {0,{calculated_coordination_matrix(0,0)},{calculated_coordination_matrix(0,1)}};
    // arma::ones(1,3);
    opt_var_vec.print("Initial starting points: ");
 
    // Convert from Armadillo vector to std vector
    std::vector<double> x = arma::conv_to< std::vector<double> >::from(opt_var_vec);

    // Create the optimizer
    // LN_xx: gradient free
    // LD_xx: with gradient  //LD_LBFGS   LD_SLSQP
    //nlopt::opt optimizer(nlopt::LD_LBFGS, 3); 
    nlopt::opt optimizer(nlopt::LD_MMA, 3);
    // Set the objective function
    optimizer.set_min_objective(objective_function, &opt_data);
    // Set optimization parameters
    optimizer.set_stopval(1e-30);
    optimizer.set_xtol_rel(1e-30);
    optimizer.set_ftol_rel(1e-30);
    optimizer.set_maxeval(10000);
    optimizer.set_initial_step(1000);
  

    try {
        nlopt::result result = optimizer.optimize(x, minf);
        // Extract optimal solution
        // Convert the solution from std::vector to Armadillo vector
        arma::rowvec opt_var_vec_result = arma::conv_to< arma::rowvec >::from(x);
        // Extract rotation angles
        double optimized_rotation_angles = std::fmod(opt_var_vec_result(0), 360.0);
        // Extract translation point
        arma::rowvec optimized_translation_point = opt_var_vec_result.subvec(1,2);
        
        // Print results to terminal
        std::cout << "Min objective function: " << minf << std::endl;
        std::cout << "Total number of iterations: " << optimizer.get_numevals() << std::endl;
        std::cout << "Optimization result: " << result << std::endl;
        std::cout << "Optimized Rotation Angle (deg): " << optimized_rotation_angles* (180.0/3.141592653589793238463) << std::endl;
        optimized_translation_point.print("Optimized Translation Point:");
        opt_var_vec_result.save("APS_RotationVerification_optimized_angle_translation_point.txt", arma::raw_ascii);

        // Calculate the rotated coordination
        arma::mat calculated_coordination_matrix;
        calculated_coordination_matrix.load("APS_coordinates_optimized.txt");
        //calculated_coordination_matrix.load("APS_coordinates_last.txt");
        arma::mat  rotated_coordination_matrix = calculate_rotated_coordinates(calculated_coordination_matrix,optimized_rotation_angles,optimized_translation_point);
        rotated_coordination_matrix.save("APS_RotationVerification_optimized_rotated_coordinates.txt", arma::raw_ascii);
        //arma::mat  test_matrix = calculate_rotated_coordinates(calculated_coordination_matrix,{0.7,0.7,0.7},{1,1,1});
        //test_matrix.save("APS_RotationVerification_optimized_rotated_coordinates_test.txt", arma::raw_ascii);
        
    } 
    catch (const nlopt::roundoff_limited& e) {
    std::cerr << "Optimization stopped due to roundoff errors: " << e.what() << std::endl;
  } 
  catch (const nlopt::forced_stop& e) {
    std::cerr << "Optimization was forcibly stopped: " << e.what() << std::endl;
  } 
  catch (const std::invalid_argument& e) {
    std::cerr << "Invalid argument: " << e.what() << std::endl;
  } 
  catch (const std::bad_alloc& e) {
    std::cerr << "Memory allocation error: " << e.what() << std::endl;
  } 
  catch (const std::runtime_error& e) {
    std::cerr << "Runtime error: " << e.what() << std::endl;
    // Diagnostics: Get and print last known optimum
    double last_opt_f = optimizer.last_optimum_value();
    nlopt::result last_result = optimizer.last_optimize_result();
    // Get the last x vector
    std::vector<double> last_x = x;  // x contains the last values used by optimizer
    
    // Convert the last solution to get coordinates and wind vector
    arma::rowvec opt_var_vec_result = arma::conv_to< arma::rowvec >::from(last_x);
    // Extract rotation angles
    double last_rotation_angles = std::fmod(opt_var_vec_result(0), 360);
    // Extract translation point
    arma::mat last_translation_point = opt_var_vec_result.subvec(1,2);

    // Print diagnostics
    std::cerr << "Last known objective value: " << last_opt_f << std::endl;
    std::cerr << "Last optimization result code: " << last_result << std::endl;
    std::cout << "Total number of iterations performed: " << optimizer.get_numevals() << std::endl;
    
    // Print last values
    std::cout << "Optimized Rotation Angles: " << last_rotation_angles << std::endl;
    last_translation_point.print("Last Translation Point:");

    // Save the last result to text file as the output
    // Calculate the rotated coordination
    arma::mat calculated_coordination_matrix;
    calculated_coordination_matrix.load("APS_coordinates_optimized.txt");
    arma::mat  rotated_coordination_matrix = calculate_rotated_coordinates(calculated_coordination_matrix,last_rotation_angles,last_translation_point);
    rotated_coordination_matrix.save("APS_RotationVerification_optimized_rotated_coordinates.txt", arma::raw_ascii);
    
  }
    catch (std::exception& e) {
      std::cerr << "An unknown error occurred: " << e.what() << std::endl;
    // Diagnostics: Get and print last known optimum
    double last_opt_f = optimizer.last_optimum_value();
    nlopt::result last_result = optimizer.last_optimize_result();
    std::cerr << "Last known objective value: " << last_opt_f << std::endl;
    std::cerr << "Last optimization result code: " << last_result << std::endl;
    std::cout << "Total number of iterations performed: " << optimizer.get_numevals() << std::endl;
    }
}