//
// Created by charlie on 4/19/17.
//

#include "measurement_models.h"
#include "Eigen/Dense"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

VectorXd MeasurementModels::LidarCRTV(VectorXd x) {
    VectorXd z_pred = VectorXd(2);
    z_pred(0) = x(0);
    z_pred(1) = x(1);
    return z_pred;
}

MatrixXd MeasurementModels::LidarCRTV_vec(MatrixXd X) {
    MatrixXd Z_pred = MatrixXd(2, X.row(0).size());
    for (int i; i < X.row(0).size(); i++) {
        Z_pred.col(i) = LidarCRTV(X.col(i));
    }
    return Z_pred;
}

VectorXd MeasurementModels::RadarCRTV(VectorXd x) {
    VectorXd z_pred = VectorXd(3);
    z_pred(0) = sqrt(pow(x(0), 2) + pow(x(1), 2));
    z_pred(1) = x(4);
    z_pred(2) = x(3);
}

MatrixXd MeasurementModels::RadarCRTV_vec(MatrixXd X) {
    MatrixXd Z_pred = MatrixXd(2, X.row(0).size());
    for (int i; i < X.row(0).size(); i++) {
        Z_pred.col(i) = RadarCRTV(X.col(i));
    }
    return Z_pred;
}
