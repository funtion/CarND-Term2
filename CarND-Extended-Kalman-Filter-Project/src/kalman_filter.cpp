#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  this->x_ = x_in;
  this->P_ = P_in;
  this->F_ = F_in;
  this->H_ = H_in;
  this->R_ = R_in;
  this->Q_ = Q_in;
}

void KalmanFilter::Predict() {
  this->x_ = this->F_ * this->x_;
  MatrixXd Ft = this->F_.transpose();
  this->P_ = this->F_ * this->P_ * Ft + this->Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  const int n = this->x_.size();
  MatrixXd Ht = H_.transpose();
  VectorXd y = z - this->H_ * this->x_;
  MatrixXd S = this->H_ * this->P_ * Ht + this->R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = this->P_ * Ht * Si;

  this->x_ = this->x_ + K * y;
  MatrixXd I = MatrixXd::Identity(n, n);
  this->P_ = (I - K * this->H_) * this->P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  const double PI = acos(-1);
  if (z(1) > PI || z(1) < -PI) {
    // some time we will have invalid input
    return;
  }

  const int n = this->x_.size();
  double p = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  double phi = atan2(x_(1), x_(0));

  double v = (x_(0)*x_(2) + x_(1)*x_(3)) / p;

  VectorXd h = VectorXd(3);
  h << p, phi, v;
  
  VectorXd y = z - h;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = this->H_ * this->P_ * Ht + this->R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = this->P_ * Ht * Si;

  this->x_ = this->x_ + K * y;
  MatrixXd I = MatrixXd::Identity(n, n);
  this->P_ = (I - K * this->H_) * this->P_;
}
