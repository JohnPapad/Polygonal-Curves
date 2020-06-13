#ifndef DISTANCES_H
#define DISTANCES_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

using Eigen::MatrixXd;
using Eigen::JacobiSVD;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;
using Eigen::Map;


vector<vector<double>> convert_matrix_to_vector(MatrixXd &mtrx);





#endif // DISTANCES_H
