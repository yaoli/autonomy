#include "Eigen/Dense"
#include "ukf.h"
#include <iostream>

using Eigen::MatrixXd;

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

  MatrixXd Xsig_pred = MatrixXd(15, 5);
  ukf.SigmaPointPrediction(&Xsig_pred);
       
  // print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  return 0;
}
