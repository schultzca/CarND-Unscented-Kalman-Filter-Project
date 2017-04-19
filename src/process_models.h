//
// Created by Charles Schultz on 4/19/17.
//

#ifndef UNSCENTEDKF_STATE_TRANSITION_MODELS_H
#define UNSCENTEDKF_STATE_TRANSITION_MODELS_H

#endif //UNSCENTEDKF_STATE_TRANSITION_MODELS_H

#include "Eigen/Dense"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class ProcessModels {
public:

    /**
     * Constant turn rate and velocity model. This model expects a state
     * vector composed of 7 states. The states are x and y position, velocity
     * magnitude, yaw angle, yaw rate, longitudinal acceleration, and yaw
     * acceleration. This function returns the predicted state given time
     * delta as input.
     * @param (VectorXd) x State vector (dim 7)
     * @param (double) dt Time step in seconds
     * @return (VectorXd) Predicted state vector
     */
    static VectorXd CRTV(VectorXd x, double dt);

    /**
     * Vectorized version of function. Takes matrix composed of state vectors
     * and returns matrix of predicted state vectors.
     * @param (MatrixXd) x State vectors (dim 7, n)
     * @param (double) dt Time step in seconds
     * @return (MatrixXd) Predicted state vectors
     */
    static MatrixXd CRTV_vec(MatrixXd X, double dt);

};