#include <iostream>
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