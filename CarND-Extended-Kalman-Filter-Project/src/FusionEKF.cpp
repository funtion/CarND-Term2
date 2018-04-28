#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  this->R_laser_ = MatrixXd(2, 2);
  this->R_radar_ = MatrixXd(3, 3);
  this->H_laser_ = MatrixXd(2, 4);
  this->Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  this->R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  this->R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  this->H_laser_ << 1, 0, 0, 0,
                    0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  cout<<measurement_pack.sensor_type_<<endl;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;

    MatrixXd P = MatrixXd(4, 4);
    P << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 10, 0,
         0, 0, 0, 10;
    MatrixXd F = MatrixXd(4, 4);
    F << 1, 0, 1, 0,
         0, 1, 0, 1,
         0, 0, 1, 0,
         0, 0, 0, 1;
    MatrixXd Q = MatrixXd::Zero(2, 2);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        VectorXd x(4);
        double p = measurement_pack.raw_measurements_(0);
        double phi = measurement_pack.raw_measurements_(1);
        double v = measurement_pack.raw_measurements_(2);

        x << p * cos(phi), p * sin(p), v * cos(phi), v * sin(p);

        this->ekf_.Init(x, P, F, this->H_laser_, this->R_radar_, Q);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      VectorXd x(4);
      x << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
      this->ekf_.Init(x, P, F, this->H_laser_, this->R_laser_, Q);
    }

    this->previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  cout<<measurement_pack.timestamp_<<endl;
  double dt = (measurement_pack.timestamp_ - this->previous_timestamp_)/1000000.0;
  MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, dt, 0,
        0, 1, 0, dt,
        0, 0, 1, 0,
        0, 0, 0, 1;
  this->previous_timestamp_ = measurement_pack.timestamp_;
  MatrixXd Qv = MatrixXd(2, 2);
  Qv << 9, 0,
        0, 9;
  MatrixXd G = MatrixXd(4, 2);
  double dtt2 = dt * dt / 2;
  G << dtt2, 0,
       0, dtt2,
       dt, 0,
       0, dt;
  MatrixXd Gt = G.transpose();
  MatrixXd Q = G * Qv * Gt;

  this->ekf_.F_ = F;
  this->ekf_.Q_ = Q;

  this->ekf_.Predict();
  cout <<"Predict\n";
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  cout<<"raw_measurements_\n";
  cout<<measurement_pack.raw_measurements_<<endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    this->Hj_ = this->tools.CalculateJacobian(this->ekf_.x_);
    this->ekf_.H_ = this->Hj_;
    this->ekf_.R_ = this->R_radar_;
    this->ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    this->ekf_.H_ = this->H_laser_;
    this->ekf_.R_ = this->R_laser_;
    this->ekf_.Update(measurement_pack.raw_measurements_);
  }
  cout<<"Update\n";
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << "------------------------\n";
}
