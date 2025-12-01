#include "AcousticPositioning_LMPart_v2_2D.hpp"
#include "AcousticPositioning_MDSPart_v6_2D.hpp"

// Define Speed of sound in m/s
const double c = 343;

int main()
{
  // 1 - Load measured time of flights from text file
  std::cout << "Step 1 - Load measured time of flights from text file" << std::endl;
  arma::mat t_measured;
  t_measured.load("t_measured.txt");
  int sensor_num = t_measured.n_rows;

  // 2 - Calculate relative distances from Equation (10)
  std::cout << "Step 2 - Calculate relative distances from Equation (10)" << std::endl;
  arma::mat d_calculated = arma::zeros(sensor_num, sensor_num);
  for (int i = 0; i < sensor_num; ++i)
  {
    for (int j = 0; j < sensor_num; ++j)
    {
      if (i != j)
      {
        //d_calculated(i, j) = 2 * c / (1 / (t_measured(i, j) + 1e-12) + 1 / (t_measured(j, i) + 1e-12));
        d_calculated(i, j) = 0.5 * c * (t_measured(i, j) + t_measured(j, i));
      }
    }
  }
  d_calculated.save("APS_relative_distance_calculated.txt", arma::raw_ascii);

  // 3 - Run Multi Dimensional Scaling (MDS) algorithm to calculate sensor coordinates
  std::cout << "Step 3 - Run Multi Dimensional Scaling (MDS) algorithm to calculate sensor coordinates" << std::endl;
  OptimizationData_MDS opt_data_MDS;
  arma::mat coordination_matrix = run_MDSOptimization(opt_data_MDS);
  coordination_matrix.save("APS_MDS_coordinates_calculated.txt", arma::raw_ascii);

  // 4 - Calculate wind vector using Linear Least Squares (LLS) method
  std::cout << "Step 4 - Calculate wind vector using Linear Least Squares (LLS) method" << std::endl;
  // std::cout << "Note: z element of the wind vector is considered 0." << std::endl;
  arma::mat A = arma::zeros(sensor_num * (sensor_num - 1) / 2, 2);
  arma::mat b = arma::zeros(sensor_num * (sensor_num - 1) / 2, 1);
  double ij = 0;
  for (int i = 0; i < sensor_num - 1; ++i)
  {
    for (int j = i + 1; j < sensor_num; ++j)
    {
      arma::mat temp_A = (coordination_matrix.row(i) - coordination_matrix.row(j)) / d_calculated(i, j);
      // Remove third column, related to z axis values
      // temp_A.shed_col(2);
      A.row(ij) = temp_A;
      b.row(ij) = 0.5 * (1 / t_measured(i, j) - 1 / t_measured(j, i)) * d_calculated(i, j);
      ij++;
    }
  }
  arma::vec wind_calculated_col = solve(A, b);
  arma::mat wind_mat = trans(arma::mat(wind_calculated_col));
  arma::rowvec wind_calculated = wind_mat.as_row();
  wind_calculated.print("Calculated Wind vector 2D: ");
  wind_calculated.save("APS_wind_calculated.txt", arma::raw_ascii);

  /*
  // Not doing LM optimization
  std::cout << "The results without final optimization is saved." << std::endl;
  coordination_matrix.save("APS_coordinates_optimized.txt", arma::raw_ascii);
  */

  // 5 - Calculated optimized wind and coordination matrices using Levenberg-Marquardt (LM) algorithm

  // Run the optimization
  std::cout << "Step 5 - Calculated optimized wind and coordination matrices using Levenberg-Marquardt (LM) algorithm" << std::endl;

  // Using nlopt library
  OptimizationData opt_data;
  run_LMOptimization(opt_data);

  return 0;
}