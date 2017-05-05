#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

double ConstrainAngle(double x){
    x = fmod(x + M_PI, M_2_PI);
    if (x < 0)
        x += M_2_PI;
    return x - M_PI;
}

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
VectorXd RadarMeasurement(VectorXd x){
    double  px = x(0);
    double  py = x(1);
    double  vm = x(2);
    double  yaw = x(3);
    double  yawd = x(4);

    VectorXd Z = VectorXd(3);

    if (sqrt(pow(px*px + py*py,2))>0.001){
        Z(0) = sqrt(px*px+py*py);
        Z(1) = atan2(py,px);
        Z(2) = (px*vm*cos(yaw) + py*vm*sin(yaw))/Z(0);
        return Z;
    } else {
        Z(0) = sqrt(0.001);
        Z(1) = atan2(0.001,0.001);
        Z(2) = (px*vm*cos(yaw) + py*vm*sin(yaw))/Z(0);
        return Z;
    }
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

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
        if (fabs(meas_package.raw_measurements_(0)) > 0.001) {
            UpdateRadar(meas_package);
        }
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

    Xsig_pred_ = Xsig_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    VectorXd z = meas_package.raw_measurements_;

    MatrixXd H = MatrixXd(2, n_x_);
    H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;

    MatrixXd R = MatrixXd(2, 2);
    R << pow(std_laspx_, 2), 0,
         0, pow(std_laspy_, 2);

    VectorXd z_pred = H * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P_ * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H) * P_;

    NIS_laser_ = (z - z_pred).transpose() * S * (z - z_pred);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    int n_z = 3;
    int n_sig = 2 * n_aug_ + 1;

    VectorXd weights = VectorXd(n_sig);
    weights(0) = lambda_ / (lambda_ + n_aug_);
    weights.tail(n_sig - 1).fill(0.5 / (lambda_ + n_aug_));

    // convert sigma points to measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sig);
    for (int i = 0; i < n_sig; i++) {
        Zsig.col(i) = RadarMeasurement(Xsig_pred_.col(i));
    }

    // predict state mean
    VectorXd z_pred = VectorXd(3);
    z_pred = Zsig * weights;

//    cout << z_pred << endl;

    MatrixXd ones = MatrixXd(1, n_sig);
    ones.setOnes();

    MatrixXd Z_diff = MatrixXd(n_z, n_sig);
    Z_diff = Zsig - z_pred * ones;

//    cout << Z_diff << endl;

    // calculate measurement covariance
    MatrixXd S = MatrixXd(n_z, n_z);
    S = Z_diff * MatrixXd(weights.asDiagonal()) * Z_diff.transpose();

    VectorXd r_radar = VectorXd(n_z);
    r_radar << std_radr_, std_radphi_, std_radrd_;

    S = S + MatrixXd(r_radar.asDiagonal());

    cout << S << endl;

    MatrixXd X_diff = MatrixXd(n_x_, n_sig);
    X_diff = Xsig_pred_ - x_ * ones;

    cout << X_diff << endl;

    // calculate cross correlation matrix
    MatrixXd T = MatrixXd(n_x_, n_z);
    T = X_diff * MatrixXd(weights.asDiagonal()) * Z_diff.transpose();

    cout << T << endl;

    // calculate kalman gain
    MatrixXd K = T * S.inverse();

    cout << K << endl;

    VectorXd z = VectorXd(3);
    z = meas_package.raw_measurements_;

    x_ = x_ + K * (z - z_pred);
    P_ = P_ - K * S * K.transpose();

    NIS_radar_ = (z - z_pred).transpose() * S * (z - z_pred);
}
