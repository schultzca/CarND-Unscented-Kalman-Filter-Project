#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

VectorXd ProcessModel(VectorXd x, double dt) {

    double   px = x(0);
    double   py = x(1);
    double    v = x(2);
    double  psi = x(3);
    double dpsi = x(4);
    double   va = x(5);
    double   vp = x(6);

    VectorXd xp(5);
    xp.setZero();

    VectorXd dx(5);
    dx.setZero();

    double dt2 = dt*dt;

    // handle division by zero
    if (fabs(dpsi) < 0.001) {
        dx << v*cos(psi)*dt + 0.5*dt2*cos(psi)*va,
              v*sin(psi)*dt + 0.5*dt2*sin(psi)*va,
              dt*va,
              0.5*dt2*vp,
              dt*vp;
    }
    else {
        dx << (v/dpsi)*(sin(psi + dpsi*dt) - sin(psi)) + 0.5*dt2*cos(psi)*va,
              (v/dpsi)*(cos(psi) - cos(psi + dpsi*dt)) + 0.5*dt2*sin(psi)*va,
              dt*va,
              dt*dpsi + 0.5*dt2*vp,
              dt*vp;
    }

    xp = x.head(5) + dx;

    return xp;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // initial state vector
    x_.setZero(5);

    // initial covariance matrix
    P_.setIdentity(5, 5);

    // Dimension of state vector
    n_x_ = 5;

    // Dimension of augmented state vector
    n_aug_ = n_x_ + 2;

    // Sigma point spread parameter
    lambda_ = 3 - n_aug_;

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
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if (!is_initialized_) {
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            double px = meas_package.raw_measurements_(0);
            double py = meas_package.raw_measurements_(1);
            x_(0) = px;
            x_(1) = py;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double   r = meas_package.raw_measurements_(0);
            double psi = meas_package.raw_measurements_(1);
            double  dr = meas_package.raw_measurements_(2);
            x_(0) = r*cos(psi);
            x_(1) = r*sin(psi);
            x_(2) = dr;
            x_(3) = psi;
        }

        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;

        return;
    }

    // compute time step in seconds since last measurement
    double delta_t = (meas_package.timestamp_ - time_us_)/1000000.;
    time_us_ = meas_package.timestamp_;

    double const dt = 0.05;
    while (delta_t > dt) {
        Prediction(dt);
        delta_t -= dt;
    }
    Prediction(delta_t);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    int n_sig = 2 * n_aug_ + 1;

    // create augmented state vector
    VectorXd xa(n_aug_);
    xa.setZero();
    xa.head(n_x_) = x_;

    // create process noise covariance matrix
    VectorXd Q(2);
    Q << std_a_ * std_a_, std_yawdd_ * std_yawdd_;

    // create augmented state covariance matrix
    MatrixXd Pa(n_aug_, n_aug_);
    Pa.fill(0.0);
    Pa.topLeftCorner(n_x_, n_x_) = P_;
    Pa.bottomRightCorner(2, 2) = MatrixXd(Q.asDiagonal());

    // compute cholesky decomposition of Pa
    MatrixXd Aa(n_aug_, n_aug_);
    Aa = Pa.llt().matrixL();

    // compute augmented sigma point matrix
    MatrixXd Xa(n_aug_, n_sig);
    MatrixXd ora = MatrixXd::Ones(1, n_aug_);   // ones row vector (1, n_aug_)
    Xa << xa,
          xa * ora + sqrt(lambda_ + n_aug_) * Aa,
          xa * ora - sqrt(lambda_ + n_aug_) * Aa;

    // compute predicted sigma points
    MatrixXd Xp(n_x_, n_sig);
    for (int i = 0; i < n_sig; i++) {
        Xp.col(i) = ProcessModel(Xa.col(i), delta_t);
    }

    // calculate sigma point weights
    VectorXd weights(n_sig);
    weights(0) = lambda_ / (lambda_ + n_aug_);
    weights.tail(n_sig - 1).fill(0.5 / (lambda_ + n_aug_));

    // calculate predicted mean
    x_ = Xp * weights;

    MatrixXd ors = MatrixXd::Ones(1, n_sig);
    P_ = (Xp - x_*ors) * MatrixXd(weights.asDiagonal()) * (Xp - x_*ors).transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

}
