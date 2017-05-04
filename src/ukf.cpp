#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * State prediction function.
 * @param x augmented state vector (dim 7)
 * @param dt time step in seconds
 * @return predicted state vector (dim 5)
 */
VectorXd StatePrediction(VectorXd x, double dt) {

    double px = x(0);
    double py = x(1);
    double vm = x(2);
    double yaw = x(3);
    double yawd = x(4);
    double va = x(5);
    double vy = x(6);

    VectorXd x_pred = VectorXd(5);
    x_pred.setZero();

    if (fabs(yawd) < 0.001) {
        x_pred << vm * cos(yaw) * dt + 0.5 * pow(dt, 2) * cos(yaw) * va,
                  vm * sin(yaw) * dt * 0.5 * pow(dt, 2) * sin(yaw) * va,
                  va * dt,
                  0.5 * pow(dt, 2) * vy,
                  vy * dt;
    } else {
        x_pred << (vm / yawd) * (sin(yaw + yawd * dt) - sin(yaw)) + 0.5 * pow(dt, 2) * cos(yaw) * va,
                  (vm / yawd) * (cos(yaw) - cos(yaw + yawd * dt)) * 0.5 * pow(dt, 2) * sin(yaw) * va,
                  va * dt,
                  yawd * dt + 0.5 * pow(dt, 2) * vy,
                  vy * dt;
    }

    x_pred = x_pred + x.head(5);

    return x_pred;
}

/**
 * Radar measurement function.
 * @param x state vector (dim 5)
 * @return measurement vector
 */
VectorXd RadarMeasurement(VectorXd x) {
    double px = x(0);
    double py = x(1);
    double vm = x(2);
    double yaw = x(3);
    double yawd = x(4);

    VectorXd z = VectorXd(3);

    const double th = 0.001;

    // avoid division by zero
    if (sqrt(pow(px, 2) + pow(py, 2)) > th) {
        z(0) = sqrt(pow(px, 2) + pow(py, 2));
        z(1) = atan2(py, px);
        z(2) = (px * vm * cos(yaw) + py * vm * sin(yaw)) / z(0);
    } else {
        z(0) = th;
        z(1) = atan2(th, th);
        z(2) = (px * vm * cos(yaw) + py * vm * sin(yaw)) / z(0);
    }

    return z;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = false;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = false;

    // initial state vector
    x_.setZero(5);

    // initial covariance matrix
    P_.setIdentity(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2.5;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.5;

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
    n_x_ = 5;

    // Dimension of augmented state vector
    n_aug_ = n_x_ + 2;

    // Sigma point spread parameter
    lambda_ = 3 - n_aug_;
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

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            double rho = meas_package.raw_measurements_(0);   // radial distance
            double psi = meas_package.raw_measurements_(1);   // heading
            double rhod = meas_package.raw_measurements_(2);  // radial velocity

            x_(0) = rho * cos(psi);    // px
            x_(1) = rho * sin(psi);    // py
            x_(2) = rhod;
            x_(3) = psi;
            x_(4) = 0;

            if (fabs(rho) > 0.001) {
                is_initialized_ = true;
            }
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

            if (sqrt(pow(x_(0), 2) + pow(x_(1), 2)) > 0.001) {
                is_initialized_ = true;
            }
        }

        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;

    const double dt = 0.05;

    while (delta_t > 0.1) {
        Prediction(dt);
        delta_t -= dt;
    }
    UKF::Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        double px = meas_package.raw_measurements_(0);
        double py = meas_package.raw_measurements_(1);
        if (sqrt(pow(px, 2) + pow(py, 2)) > 0.001) {
            UpdateLidar(meas_package);
        }
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    int n_sig = 2 * n_aug_ + 1;

    MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig);
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    MatrixXd A_aug = MatrixXd(n_aug_, n_aug_);
    MatrixXd ones_aug = MatrixXd(1, n_aug_); // dynamic ones vector
    ones_aug = ones_aug.setOnes();
    MatrixXd ones_sig = MatrixXd(1, n_sig);
    ones_sig.setOnes();

    x_aug.setZero();
    x_aug.head(5) << x_;

    P_aug.setZero();
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = pow(std_a_, 2);
    P_aug(6, 6) = pow(std_yawdd_, 2);

    A_aug = P_aug.llt().matrixL();

    // row vector to broad cast mean vector as matrix


    Xsig_aug << x_aug,
                x_aug * ones_aug + sqrt(lambda_ + n_aug_) * A_aug,
                x_aug * ones_aug - sqrt(lambda_ + n_aug_) * A_aug;

    MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig);
    Xsig_pred.col(0) = StatePrediction(Xsig_aug.col(0), delta_t);
    for (int i = 1; i < n_sig; i++) {
        Xsig_pred.col(i) = StatePrediction(Xsig_aug.col(i), delta_t);
    }

    VectorXd weights = VectorXd(n_sig);
    MatrixXd W = MatrixXd(n_sig, n_sig);

    weights.fill(0.5 / (lambda_ + n_aug_));
    weights(0) = lambda_ / (lambda_ + n_aug_);

    W = MatrixXd(weights.asDiagonal());

    x_ = Xsig_pred * weights;
    P_ = (Xsig_pred - x_ * ones_sig) * W * (Xsig_pred - x_ * ones_sig).transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    // measurement vector
    VectorXd z = meas_package.raw_measurements_;

    // measurement function
    MatrixXd H = MatrixXd(2, 5);
    H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;

    // measurement noise covariance
    MatrixXd R = MatrixXd(2, 2);
    R(0, 0) = pow(std_laspx_, 2);
    R(1, 1) = pow(std_laspy_, 2);

    VectorXd z_pred = H * x_;
    VectorXd y = z - z_pred;
    MatrixXd S = H * P_ * H.transpose() + R;

    // calculate kalman gain
    MatrixXd K = P_ * H.transpose() * S.inverse();

    // update state estimate and estimate covariance
    MatrixXd I = MatrixXd(n_x_, n_x_);
    x_ = x_ + K * y;
    P_ = (I - K * H) * P_;

    NIS_laser_ = (z - z_pred).transpose() * S * (z - z_pred);
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
