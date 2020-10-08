#include "tracking.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

Tracking::Tracking() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // create a 4D state vector, we don't know yet the value of the x state
  kf_.x_ = VectorXd(4);

  // state covariance matrix
  kf_.P_ = MatrixXd(4, 4);
  kf_.P_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;

  // measurement covariance
  kf_.R_ = MatrixXd(2, 2);
  kf_.R_ << 0.0225, 0, 0, 0.0225;

  // measurement matrix
  kf_.H_ = MatrixXd(2, 4);
  kf_.H_ << 1, 0, 0, 0, 0, 1, 0, 0;

  // state transition matrix
  kf_.F_ = MatrixXd(4, 4);
  kf_.F_ << 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  // set the acceleration noise component
  noise_ax = 5;
  noise_ay = 5;
}

Tracking::~Tracking() {}

void Tracking::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  if (!is_initialized_) {
    cout << "Kalman filter init" << endl;

    // set the state with the initial location and zero velocity
    kf_.x_ << measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1], 0, 0;
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  // compute the time elapsed between the current and previous measurements
  // dt: expressed in second
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // TODO: YOUR CODE HERE
  // 1. Modify the F matrix so that the time is integrated
  kf_.F_(0, 2) = dt;
  kf_.F_(0, 3) = dt;
  // 2. Set the process covariance matrix Q
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  kf_.Q_ = MatrixXd(4, 4);
  kf_.Q_ <<
    dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
    0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
    dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
    0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  // 3. Call the Kalman Filter predict() function
  kf_.Predict();
  // 4. Call the Kalman Filter update() function
  kf_.Update(measurement_pack.raw_measurements_);
  //      with the most recent raw measurements_

  cout << "x= " << kf_.x_ << endl;
  cout << "P_= " << kf_.P_ << endl;
}
