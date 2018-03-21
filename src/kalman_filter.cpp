#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
#define PI 3.14159265
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
  /**
  TODO:
    * predict the state
  */
  this->x_ = this->F_*this->x_;
  this->P_ = this->F_*this->P_*this->F_.transpose() + this->Q_;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &R_laser) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  Eigen::VectorXd err = z - this->H_*this->x_;
  Eigen::MatrixXd S = this->H_*this->P_*this->H_.transpose() + R_laser;
  Eigen::MatrixXd Ka_gain = this->P_*this->H_.transpose()*S.inverse();

  // Update the state as per the new measurements
  this->x_ = this->x_+ (Ka_gain * err);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(4, 4);
  this->P_ = (I - (Ka_gain * this->H_))*this->P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &R_radar) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  Eigen::VectorXd H_xprime = VectorXd(3);
  H_xprime << sqrt(this->x_(0)*this->x_(0)+ this->x_(1)*this->x_(1)),
            atan2(this->x_(1),this->x_(0)),
            ((this->x_(0)*this->x_(2)) + (this->x_(1)*this->x_(3)))/(sqrt(this->x_(0)*this->x_(0)+ this->x_(1)*this->x_(1)));



  Eigen::VectorXd err = z - H_xprime;

  Eigen::MatrixXd S = this->Hj_*this->P_*this->Hj_.transpose() + R_radar;

  Eigen::MatrixXd Ka_gain = this->P_*this->Hj_.transpose()*S.inverse();

  // Normalizing Phi angle mesaurement
  if (err(1) > (PI))
  {
    err(1) = err(1) - (2*PI);
  }
  else if(err(1) < (-PI))
  {
    err(1) = err(1) + (2*PI);
  }

  // Update the state as per the new measurements

  this->x_ = this->x_+ (Ka_gain * err);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(4, 4);
  this->P_ = (I - (Ka_gain * this->Hj_))*this->P_;

}
