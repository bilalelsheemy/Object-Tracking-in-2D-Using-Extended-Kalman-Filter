#include "FusionEKF.h"
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
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  //Hj_ = MatrixXd(3, 4);

  // measurement space matrix - LASER
  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

 ekf_.P_ = MatrixXd(4, 4);
 ekf_.P_ << 1, 0, 0, 0,
     0, 1, 0, 0,
     0, 0, 1000, 0,
     0, 0, 0, 1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

FusionEKF fusionekf_;
Tools tools_;

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
   float dt = (measurement_pack.timestamp_ - fusionekf_.previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
   float dt_2 = dt * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;
   fusionekf_.previous_timestamp_ = measurement_pack.timestamp_ ;

   cout<<"dt = "<<dt<<endl;
   cout<<"previous_timestamp_ = "<<fusionekf_.previous_timestamp_<<endl;

   int noise_ax = 9;
   int noise_ay = 9;

  fusionekf_.Q = MatrixXd(4, 4);
  fusionekf_.Q <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  fusionekf_.F = MatrixXd(4, 4);
  fusionekf_.F << 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;




  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // first measurement
    cout << "EKF: " << endl;

    ekf_.x_ = VectorXd(4);
    ekf_.Hj_ = MatrixXd(3, 4);
    ekf_.H_ = MatrixXd(2, 4);

    ekf_.Q_ = MatrixXd(4, 4);

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ = fusionekf_.F;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float x_pos, y_pos, x_vel, y_vel;
      float rho, phi, rho_vel;


      rho = measurement_pack.raw_measurements_[0];
      phi = measurement_pack.raw_measurements_[1];
      rho_vel = measurement_pack.raw_measurements_[2];

      float radar_px = rho * cos(phi);
      float radar_py = rho * sin(phi);


      Eigen::MatrixXd Jacob = tools_.CalculateJacobian(ekf_.x_);
      if (Jacob != Eigen::MatrixXd::Identity(3, 4))
      {
        ekf_.Hj_ = tools_.CalculateJacobian(ekf_.x_);
      }
      cout<<"One.1"<<endl;
      ekf_.x_ << radar_px, radar_py, 0, 0;
      ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_,
                              ekf_.Hj_, fusionekf_.R_radar_, ekf_.Q_);

    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      ekf_.H_ = fusionekf_.H_laser_;
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0; // 3.122427e-01, 5.803398e-01, 0, 0;

      ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_,
                               ekf_.H_, fusionekf_.R_laser_, ekf_.Q_);

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    //fusionekf_.previous_timestamp_ = measurement_pack.timestamp_ ;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds. >> Done
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix. >> Done.
   */
   ekf_.Q_ = fusionekf_.Q;
   ekf_.F_ = fusionekf_.F;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices. >> Done.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    //cout<< measurement_pack.raw_measurements_ <<endl;
    Eigen::MatrixXd Jacob = tools_.CalculateJacobian(ekf_.x_);
    if (Jacob != Eigen::MatrixXd::Identity(3, 4))
    {
      ekf_.Hj_ = tools_.CalculateJacobian(ekf_.x_);
    }
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, fusionekf_.R_radar_);
  }

  else {
    // Laser updates
    ekf_.H_ = fusionekf_.H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_, fusionekf_.R_laser_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
