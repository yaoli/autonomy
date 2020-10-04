#include "Eigen/Dense"
#include "ukf.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
int main() {

  // Create a UKF instance
  UKF ukf;

  /**
   * Programming assignment calls
   */
  // The state vector is 5D of [P_x, P_y, v, phi, phi_dot].
  // Therefore, the number of sigma points is 2*5+1=11
  //MatrixXd Xsig = MatrixXd(5, 11);
  //ukf.GenerateSigmaPoints(&Xsig);
  
  //MatrixXd Xsig_aug = MatrixXd(7, 15);
  //ukf.AugmentedSigmaPoints(&Xsig_aug);

  //MatrixXd Xsig_pred = MatrixXd(15, 5);
  //ukf.SigmaPointPrediction(&Xsig_pred);

  //VectorXd x_pred = VectorXd(5);
  //MatrixXd P_pred = MatrixXd(5, 5);
  //ukf.PredictMeanAndCovariance(&x_pred, &P_pred);

  VectorXd z_out = VectorXd(3);
  MatrixXd S_out = MatrixXd(3, 3);
  ukf.PredictRadarMeasurement(&z_out, &S_out);
  // print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  return 0;
}
