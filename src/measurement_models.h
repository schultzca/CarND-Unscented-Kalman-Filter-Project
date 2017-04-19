//
// Created by Charles Schultz on 4/19/17.
//

#ifndef UNSCENTEDKF_MEASUREMENT_MODELS_H
#define UNSCENTEDKF_MEASUREMENT_MODELS_H

#endif //UNSCENTEDKF_MEASUREMENT_MODELS_H

#include "Eigen/Dense"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class MeasurementModels {
    /**
     * CRTV radar state vector to measurement space model. This function
     * converts the CRTV state vector to the corresponding measurement
     * vector. More specifically it computes longitudinal distance, heading
     * angle, and longitudinal velocity from the state vector.
     * @param (VectorXd) x State vector (dim 5)
     * @return (VectorXd) radar measurement vector (dim 3)
     */
    static VectorXd RadarCRTV(VectorXd x);

    /**
     * Vectorized version of RadarCRTV.
     * @param (MatrixXd) X State vector matrix (dim 5, n)
     * @return (MatrixXd) radar measurement matrix (dim 3, n)
     */
    static MatrixXd RadarCRTV_vec(MatrixXd X);

    /**
     * CRTV lidar state vector to measurement space model.
     * @param (VectorXd) x State vector
     * @return (VectorXd) lidar measurement vector (dim 2)
     */
    static VectorXd LidarCRTV(VectorXd x);

    /**
     * Vectorized version of LidarCRTV.
     * @param (MatrixXd) X State vector matrix (dim 5, n)
     * @return (MatrixXd) lidar measurement matrix (dim 2, n)
     */
    static MatrixXd LidarCRTV_vec(MatrixXd X);
};