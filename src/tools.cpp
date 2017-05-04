#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

Eigen::IOFormat Tools::NumpyArrayFormat(5, 0, ",", ",", "[", "]", "[", "]");

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse = VectorXd(4);
    rmse.setZero();

    // compute sum of squared residuals
    for (int i = 0; i < estimations.size(); i++) {
        rmse = rmse.array() + (estimations[i] - ground_truth[i]).array().pow(2);
    }

    // calculate mean of squared residuals
    rmse = rmse / estimations.size();

    // calculate the square root
    rmse = rmse.array().sqrt();

    return rmse;
}
