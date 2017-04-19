//
// Created by Charles Schultz on 4/19/17.
//

#include "process_models.h"
#include "Eigen/Dense"
#include <cmath>

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

VectorXd ProcessModels::CRTV(VectorXd x, double dt) {

    double dt2 = dt*dt;

    // create predicted state vector
    VectorXd x_pred = x;

    // Deterministic state
    double v = x(2);
    double y = x(3);
    double yd = x(4);

    // Stochastic state
    double a = x(5);
    double ydd = x(6);

    // zero yaw rate equations
    if (x(4) < 0.00001) {
        x_pred(0) += v*cos(y)*dt;
        x_pred(1) += v*sin(y)*dt;
    }
    // default equations
    else {
        x_pred(0) += (v/yd)*(sin(y+yd*dt)-sin(y));
        x_pred(1) += (v/yd)*(-cos(y+yd*dt)+cos(y));
        x_pred(3) += yd*dt;
    }

    // add stochastic components
    x_pred(0) += 0.5*dt2*cos(y)*a;
    x_pred(1) += 0.5*dt2*sin(y)*a;
    x_pred(2) += a*dt;
    x_pred(3) += 0.5*dt2*ydd;
    x_pred(4) += ydd*dt;

    return x_pred;
}

MatrixXd ProcessModels::CRTV_vec(MatrixXd X, double dt) {

    MatrixXd X_pred = X;

    for (int i=0; i < X_pred.row(0).size(); i++) {
        X_pred.col(i) = CRTV(X.col(i), dt);
    }

    return X_pred;
}



