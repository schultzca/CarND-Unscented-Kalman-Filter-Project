#include "ukf.h"
#include "tools.h"
#include "process_models.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    // Dimension of state vector
    n_x_ = x_.size();

    // Dimension of augmented state vector
    n_aug_ = n_x_ + 2;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /***************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {

        // first measurement
        cout << "UKF: " << endl;

        // initialize state vector
        x_ = VectorXd(5);
        x_.fill(0.0);

        // initialize covariance matrix
        P_ = MatrixXd(5,5);
        P_ <<   1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1000, 0, 0,
                0, 0, 0, 1000, 0,
                0, 0, 0, 0, 1000;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            double rho = meas_package.raw_measurements_(0);   // radial distance
            double phi = meas_package.raw_measurements_(1);   // heading
            double rhod = meas_package.raw_measurements_(2);  // radial velocity

            x_(0) = rho * cos(phi);    // px
            x_(1) = rho * sin(phi);    // py
            x_(2) = rhod;
            x_(3) = 0;
            x_(4) = 0;

        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_(0) = meas_package.raw_measurements_[0];   // px
            x_(1) = meas_package.raw_measurements_[1];   // py
            x_(2) = 0;
            x_(3) = 0;
            x_(4) = 0;
        }

        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;

    UKF::Prediction(dt);

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    // create augmented state vector
    VectorXd xa = VectorXd(7);
    xa.head(5) = x_;
    xa(5) = 0;
    xa(6) = 0;

    // create augmented state covariance matrix
    MatrixXd Pa = MatrixXd(xa.size(), xa.size());
    Pa.fill(0);
    Pa.topLeftCorner(x_.size(), x_.size()) = P_;
    Pa(5,5) = std_a_ * std_a_;
    Pa(6,6) = std_yawdd_ * std_yawdd_;

    // sigma point spread parameter
    double lambda = 3 - xa.size();

    // generate sigma points
    MatrixXd Xsig = GenerateSigmaPoints(xa, Pa, lambda);

    // compute predicted sigma points
    MatrixXd Xsig_pred = ProcessModels::CRTV_vec(Xsig, delta_t);

    // create mean predicted state
    VectorXd xp = VectorXd(xa.size());

    // create predicted covariance matrix
    MatrixXd Pp = MatrixXd(xa.size(), xa.size());

    // compute predicted mean and covariance matrix
    PredictMeanAndCovariance(Xsig_pred, lambda, xp, Pp);

    // Set state to predicted state
    x_ = xp.head(5);

    // Set state covariance matrix
    P_ = Pp.topLeftCorner(5, 5);

    // Set predicted sigma points for measurement update
    Xsig_pred_ = Xsig_pred.topRows(5);

}

/**
 * Generate sigma points given mean vector and corresponding
 * covariance matrix. This does not have any implicit assumptions
 * regarding the model. Augment x and P prior to calling this
 * function.
 * @param (VectorXd) x State vector
 * @param (MatrixXd) P State covariance matrix
 * @return (MatrixXd) Augmented sigma point matrix
 */
MatrixXd UKF::GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda) {

    int n_x = x.size();

    // create augmented sigma point matrix
    MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

    // compute square root of state covariance matrix
    MatrixXd Psq = P.llt().matrixL();

    // set first sigma point to mean value
    Xsig.col(0) = x;

    // set remaining sigma points
    for (int i=0; i < n_x; i++){
        Xsig.col(i+1) = x + sqrt(lambda + n_x) * Psq.col(i);
        Xsig.col(i + n_x + 1) = x - sqrt(lambda + n_x) * Psq.col(i);
    }

    return Xsig;
}

/**
 * Compute predicted mean and covariance matrix.
 * @param {MatrixXd} Xsig Sigma points
 * @param {double} lambda spreading parameter
 * @param {VectorXd} xp Predicted state vector (output)
 * @param {MatrixXd} Pp Predicted covariance matrix (output)
 */
void UKF::PredictMeanAndCovariance(MatrixXd Xsig_pred, double lambda, VectorXd &xp, MatrixXd &Pp) {

    int n_x = Xsig_pred.col(0).size();

    // compute weights
    VectorXd weights = VectorXd(2 * n_x + 1);
    double weight0 = lambda/(lambda + n_x);
    weights(0) = weight0;
    for (int i = 1; i < 2 * n_x + 1; i++) {
        double weight = 0.5/(lambda + n_x);
        weights(i) = weight;
    }

    // predicted state mean
    xp.fill(0);
    for (int i = 0; i < 2 * n_x + 1; i++) {
        xp += weights(i) * Xsig_pred.col(i);
    }

    //predicted state covariance matrix
    Pp.fill(0.0);
    for (int i = 0; i < 2 * n_x + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - xp;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Pp = Pp + weights(i) * x_diff * x_diff.transpose() ;
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
}
