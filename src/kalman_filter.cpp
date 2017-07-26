#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // Predict, same for laser and radar
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  // update function for laser

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  //new estimate
  int x_size = x_.size();
  MatrixXd I  = MatrixXd::Identity(x_size, x_size);
  x_ = x_ + (K*y);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  //update function for radar
  // rho_pred
  float rho_pred = pow(pow(x_[0], 2) + pow(x_[1], 2), .5);
  if(rho_pred<0.00001)
  {
      rho_pred=0.00001;
  }

  // theta_pred
  float theta_pred = atan2(x_[1], x_[0]);  
  if(x_[0] == 0) theta_pred = 0; // C++ returns undefined for atan(0,0). This is to return 0 instead of undefined

  // rho_dot_pred
  float rho_dot_pred = (x_[0] * x_[2] + x_[1] * x_[3]) / rho_pred;

  VectorXd z_pred(3);
  z_pred  << rho_pred, theta_pred, rho_dot_pred;

  VectorXd y = z - z_pred;
  y(1) = atan2(sin(y(1)), cos(y(1)));

  MatrixXd PHt =  P_ * H_.transpose();
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();

  // new estimates
  x_ = x_ + K*y;
  int x_size = x_.size();
  MatrixXd I  = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
