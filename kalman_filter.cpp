#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
#include <math.h>
/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  UpdateState(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  if (x_[0] == 0 && x_[1] == 0)
    return;
  double rho = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  double phi = atan2(x_[1], x_[0]);
  double rhodot;
  if (fabs(rho) < 0.001){
      rhodot = 0;
  }
  else{
      rhodot = (x_[0] * x_[2] + x_[1] * x_[3])/rho;
  }
  VectorXd h = VectorXd(3);
  h << rho, phi, rhodot;
  VectorXd y = z - h;
  while ( y(1) > M_PI || y(1) < -M_PI ) {
    if ( y(1) > M_PI ) {
      y(1) -= M_PI;
    }
    else {
      y(1) += M_PI;
    }
  }
  UpdateState(y);
  
}

void KalmanFilter::UpdateState(const VectorXd &y){
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
