// LM optimization part
#include "AcousticPositioning_LMPart_v2_2D.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <span>
#include <armadillo>
#include <nlopt.hpp> 

// Define Speed of sound in m/s
const double c = 343;

// Objective function 
// Note: If any changes are made to x and grad inside the objective_function, those changes will be reflected outside 
// the function. This is because they are passed by reference. 
// having const for x variable ensures that the function does not unintentionally alter x.
double objective_function(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
    
    auto* data = static_cast<OptimizationData*>(my_func_data);

    // Track the current x values
    data->x_history.push_back(x);

    // Convert from std::vector to Armadillo vector
    arma::rowvec x_arma = arma::conv_to< arma::rowvec >::from(x);

    // Total number of optimization variables
    int opt_var_num = x_arma.n_elem;
    int sensor_num = (opt_var_num/2)-1;

    // Convert vector of variables to matrix
    arma::mat x_mat = x_arma.as_row();
    arma::mat opt_var_mat = trans(reshape(x_mat,2,sensor_num+1));
    // Extract wind vector
    arma::rowvec wind_vector = opt_var_mat.row(sensor_num);
    // Extract sensor positions
    arma::mat coordination_matrix = opt_var_mat.rows(0,sensor_num-1);
    // Load calculated time of flight from file
    arma::mat t_calculated;
    t_calculated.load("t_measured.txt");
   

    // Compute the objective function value
    // Define t_hat elements
    arma::mat t_hat = arma::zeros(sensor_num,sensor_num);
    for (int i = 0; i < sensor_num; ++i) {
      for (int j = 0; j < sensor_num; ++j) {
        if (i != j) {
          arma::rowvec Pij = coordination_matrix.row(i)-coordination_matrix.row(j); 
          double t_hat_denum = arma::norm(Pij)*c + dot(wind_vector,Pij)+ 1e-12;
          t_hat(i,j) = pow(arma::norm(Pij),2) / t_hat_denum;
        }
      }
    }
    arma::mat residual_mat = t_calculated-t_hat; 
    double objective = accu(square(residual_mat));

    // If gradients are requested, compute them
    if (!grad.empty()) {
      // Calculate gradient matrix
      // Define gradient vector, they are g[0] to g[3*(N+1)] each related
      // to the gradient with respect to each optimization variable.
      // N sensor node positions and 1 wind vector
      // First, we define a matrix form for each gradient and call it element
      // by element. Then we create the gradient vector from the matrix.
      // ---------- Gradient Matrix Calculations
      // ---- Wind
      arma::mat gradient_mat_wx = arma::zeros(sensor_num,sensor_num);
      arma::mat gradient_mat_wy = arma::zeros(sensor_num,sensor_num);
      //arma::mat gradient_mat_wz = arma::zeros(sensor_num,sensor_num);
      // ---- Pi
      arma::mat gradient_mat_Pi_x = arma::zeros(sensor_num,sensor_num);
      arma::mat gradient_mat_Pi_y = arma::zeros(sensor_num,sensor_num);
      //arma::mat gradient_mat_Pi_z = arma::zeros(sensor_num,sensor_num);
      // ---- Pj
      arma::mat gradient_mat_Pj_x = arma::zeros(sensor_num,sensor_num);
      arma::mat gradient_mat_Pj_y = arma::zeros(sensor_num,sensor_num);
      //arma::mat gradient_mat_Pj_z = arma::zeros(sensor_num,sensor_num);
      for (int i = 0; i < sensor_num; ++i) {
        for (int j = 0; j < sensor_num; ++j) {
          if (i != j) {
            arma::rowvec Pij = coordination_matrix.row(i)-coordination_matrix.row(j); 
            double grad_denum = pow(arma::norm(Pij)*c + dot(wind_vector,Pij),2)+ 1e-12;
            // ---- Wind
            gradient_mat_wx(i,j) = pow(arma::norm(Pij),2)*(coordination_matrix(i,0)-coordination_matrix(j,0))/grad_denum;
            gradient_mat_wy(i,j) = pow(arma::norm(Pij),2)*(coordination_matrix(i,1)-coordination_matrix(j,1))/grad_denum;
            //gradient_mat_wz(i,j) = pow(arma::norm(Pij),2)*(coordination_matrix(i,2)-coordination_matrix(j,2))/grad_denum;
            // ---- Pi
            gradient_mat_Pi_x(i,j) = -(c*arma::norm(Pij)*Pij(0)+pow(arma::norm(Pij),2)*wind_vector(0))/grad_denum;
            gradient_mat_Pi_y(i,j) = -(c*arma::norm(Pij)*Pij(1)+pow(arma::norm(Pij),2)*wind_vector(1))/grad_denum;        
            //gradient_mat_Pi_z(i,j) = -(c*arma::norm(Pij)*Pij(2)+pow(arma::norm(Pij),2)*wind_vector(2))/grad_denum;
            // ---- Pj
            gradient_mat_Pj_x(i,j) = -gradient_mat_Pi_x(i,j);
            gradient_mat_Pj_y(i,j) = -gradient_mat_Pi_y(i,j);
            //gradient_mat_Pj_z(i,j) = -gradient_mat_Pi_z(i,j);
          }
        }
      }
      // ---------- Gradient Vector Calculations
      // Calculate the dot product of the residual and gradient matrices
      arma::mat gradient_dot_prod_wx = gradient_mat_wx % residual_mat;
      arma::mat gradient_dot_prod_wy = gradient_mat_wy % residual_mat;
      //arma::mat gradient_dot_prod_wz = gradient_mat_wz % residual_mat;
      arma::mat gradient_dot_prod_Pi_x = gradient_mat_Pi_x % residual_mat;
      arma::mat gradient_dot_prod_Pi_y = gradient_mat_Pi_y % residual_mat;
      //arma::mat gradient_dot_prod_Pi_z = gradient_mat_Pi_z % residual_mat;
      arma::mat gradient_dot_prod_Pj_x = gradient_mat_Pj_x % residual_mat;
      arma::mat gradient_dot_prod_Pj_y = gradient_mat_Pj_y % residual_mat;
      //arma::mat gradient_dot_prod_Pj_z = gradient_mat_Pj_z % residual_mat;

      // ---------- Gradient Vector
      arma::rowvec grad_arma = arma::zeros(1,2*(sensor_num+1));
      for (int k = 1; k <= sensor_num; ++k) {
        // Gradient with respect to x-axis of sensor k position
        grad_arma(2*k-2) = 2*(sum(gradient_dot_prod_Pi_x.row(k-1))+sum(gradient_dot_prod_Pj_x.col(k-1)));
        // Gradient with respect to y-axis of sensor k position
        grad_arma(2*k-1) = 2*(sum(gradient_dot_prod_Pi_y.row(k-1))+sum(gradient_dot_prod_Pj_y.col(k-1)));
        // Gradient with respect to z-axis of sensor k position
        //grad_arma(3*k-1) = 2*(sum(gradient_dot_prod_Pi_z.row(k-1))+sum(gradient_dot_prod_Pj_z.col(k-1)));
      }
      // Gradient with respect to x-axis of wind vector
      grad_arma(2*sensor_num) = 2*accu(gradient_dot_prod_wx);
      // Gradient with respect to y-axis of wind vector
      grad_arma(2*sensor_num+1) = 2*accu(gradient_dot_prod_wy);
      // Gradient with respect to z-axis of wind vector
     //grad_arma(3*sensor_num+2) = 2*accu(gradient_dot_prod_wz);
      
      // convert from armadillo vector to std vector
      grad = arma::conv_to<std::vector<double>>::from(grad_arma);
      // Save the gradient for this iteration
      data->grad_history.push_back(grad);  
    } 
    else {
        // If gradient is not requested, push a zero-filled vector
        data->grad_history.push_back(std::vector<double>(x.size(), 0.0));
    }
    return objective;
}

void run_LMOptimization(OptimizationData& opt_data) {
    double minf;

    // Initialize optimizaer by values calculated from previous problems
    // load calculated time of flight from file
    arma::mat t_calculated;
    t_calculated.load("t_measured.txt");
    // load calculated coordination from file
    arma::mat p_calculated;
    p_calculated.load("APS_MDS_coordinates_calculated.txt");
    // load calculated wind velocity vector from file
    arma::rowvec w_calculated;
    w_calculated.load("APS_wind_calculated.txt");  
    // Set w_z=0, consider wind to be 2 dimensional
    //w_calculated(2) = 0;
    // Save number of sensor nodes and optimization variables
    int sensor_num = p_calculated.n_rows;
    int opt_var_num = 2*sensor_num+2;
    // Create optimization variables in matrix and row format
    arma::mat opt_var_mat = join_vert(p_calculated,w_calculated);
    arma::rowvec opt_var_vec = opt_var_mat.as_row();
 
    // Convert from Armadillo vector to std vector
    std::vector<double> x = arma::conv_to< std::vector<double> >::from(opt_var_vec);

    // Create the optimizer
    // LN_xx: gradient free
    // LD_xx: with gradient  //LD_LBFGS   LD_SLSQP
    //nlopt::opt optimizer(nlopt::LD_LBFGS, opt_var_num); 
    nlopt::opt optimizer(nlopt::LD_MMA, opt_var_num);
    // Set the objective function
    optimizer.set_min_objective(objective_function, &opt_data);
    // Set optimization parameters
    optimizer.set_stopval(1e-30);
    optimizer.set_xtol_rel(1e-30);
    optimizer.set_ftol_rel(1e-30);
    optimizer.set_maxeval(10000);
    //optimizer.set_initial_step(100);
  

    try {
        nlopt::result result = optimizer.optimize(x, minf);
        // Extract optimal solution
        // Convert the solution from std::vector to Armadillo vector
        arma::rowvec opt_var_vec_result = arma::conv_to< arma::rowvec >::from(x);
        // Convert vector of variables to matrix
        arma::mat x_mat = opt_var_vec_result.as_row();
        arma::mat opt_var_mat = trans(reshape(x_mat,2,sensor_num+1));
        // Extract wind vector
        arma::rowvec optimized_wind_vector = opt_var_mat.row(sensor_num);
        // Extract sensor positions
        arma::mat optimized_coordination_matrix = opt_var_mat.rows(0,sensor_num-1);
        // Save optimal solution in a text file
        optimized_wind_vector.save("APS_wind_optimized.txt", arma::raw_ascii);
        optimized_coordination_matrix.save("APS_coordinates_optimized.txt", arma::raw_ascii);
        
        // Print results to terminal
        std::cout << "Min objective function: " << minf << std::endl;
        std::cout << "Total number of iterations: " << optimizer.get_numevals() << std::endl;
        std::cout << "Optimization result: " << result << std::endl;
        //optimized_wind_vector.print("Optimized Wind Vector:");
        //optimized_coordination_matrix.print("Optimized coordination matrix:");
        
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
    arma::rowvec opt_var_vec_result = arma::conv_to<arma::rowvec>::from(last_x);
    arma::mat x_mat = opt_var_vec_result.as_row();
    arma::mat opt_var_mat = trans(reshape(x_mat, 2, sensor_num + 1));
    
    // Extract last wind vector and coordination matrix
    arma::rowvec last_wind_vector = opt_var_mat.row(sensor_num);
    arma::mat last_coordination_matrix = opt_var_mat.rows(0, sensor_num-1);
    
    // Save the last values to files with different names, and same name
    last_wind_vector.save("APS_wind_last.txt", arma::raw_ascii);
    last_coordination_matrix.save("APS_coordinates_last.txt", arma::raw_ascii);
    last_wind_vector.save("APS_wind_optimized.txt", arma::raw_ascii);
    last_coordination_matrix.save("APS_coordinates_optimized.txt", arma::raw_ascii);
    
    // Print diagnostics
    std::cerr << "Last known objective value: " << last_opt_f << std::endl;
    std::cerr << "Last optimization result code: " << last_result << std::endl;
    std::cout << "Total number of iterations performed: " << optimizer.get_numevals() << std::endl;
    
    // Print last values
    last_wind_vector.print("Last Wind Vector:");
    
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