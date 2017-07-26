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

  // dit is waar jy al die initial values insit, André

  // begin state, dit gaan verander na stap 1
  is_initialized_ = false;

  // timestamp vir begin
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);  //noise vir laser
  R_radar_ = MatrixXd(3, 3);  // noise vir radar
  H_laser_ = MatrixXd(2, 4);  // H vir Laser


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) 
  {
    // initial all variables 
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ <<  1, 0, 0,    0,
                0, 1, 0,    0,
                0, 0, 1000,  0,
                0, 0, 0,  1000;

    // these are sensor specific
    ekf_.H_ = MatrixXd::Zero(2,4);
    ekf_.F_ = MatrixXd::Zero(4,4);
    ekf_.Q_ = MatrixXd::Zero(4,4);


    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    {
      // for when the first raw is radar
      float rho = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float rhodot = measurement_pack.raw_measurements_[2]; 

      float px = rho * cos(theta);
      float py = rho * sin(theta);
      float vx = 0.0;
      float vy = 0.0;

      ekf_.x_ << px,py,vx,vy;


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
    {
      // when it's laser (lasar with an A sounds cooler)
      ekf_.x_ <<  measurement_pack.raw_measurements_[0], 
                  measurement_pack.raw_measurements_[1],      
                  0.0, 
                  0.0;
    }

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/ 1000000.0;
  float dt2 = pow(dt,2);
  float dt3 = pow(dt,3);
  float dt4 = pow(dt,4);

  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ <<    1, 0, dt, 0,
                0, 1, 0, dt,
                0, 0, 1, 0,
                0, 0, 0, 1;

  float noise_ax = 9.0;
  float noise_ay = 9.0;

  ekf_.Q_ <<  dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
        0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
        dt3/2*noise_ax, 0, dt2*noise_ax, 0,
        0, dt3/2*noise_ay, 0, dt2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    if(true) // to check performance of laser/ radar alone. better accuracy with both
    {


    Tools tools;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
  } else {

    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);


  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
