#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  const int n = estimations.size();
  VectorXd result = VectorXd::Zero(estimations[0].size());

  for(int i = 0; i < n; i++) {
    const VectorXd& est = estimations.at(i);
    const VectorXd& truth = ground_truth.at(i);
    VectorXd diff = est - truth;
    diff = diff.array() * diff.array();
    result += diff;
  }
  
  result = result / n;
  return result.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj = MatrixXd::Zero(3, 4);

  //recover state parameters
  const double px = x_state(0);
  const double py = x_state(1);
  const double vx = x_state(2);
  const double vy = x_state(3);

  const double constA = sqrt(px * px + py * py);
  if( fabs(constA) > 1e-5) {
    const double constB = constA * constA;
    const double constC = constA * constB;
    Hj << px/constA, py/constA, 0, 0,
          -py/constB, px/constB, 0, 0,
          py * (vx * py - vy * px) / constC, px * (vy * px - vx * py) / constC, px/constA, py/constA;
  }

  return Hj;
}
