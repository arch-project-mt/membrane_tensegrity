#ifndef rmsd
#define rmsd
#include "rmsd_io.h"
#include <cmath>

RMSDResult CalcRMSD(Eigen::Matrix3Xd P, Eigen::Matrix3Xd Q,
                    std::vector<double> W, bool output_each_distance = false) {
  if (P.cols() != Q.cols())
    throw "CalcRMSD(): input data mis-match";
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 3);
  int p_cols_num = P.cols();
  for (int i = 0; i < p_cols_num; i++) {
    H += W[i] * Q.col(i) * P.col(i).transpose();
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU |
                                               Eigen::ComputeThinV);
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  (d > 0.0) ? d = 1.0 : d = -1.0;
  I(2, 2) = d;
  Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();
  if (output_each_distance) {
    Eigen::Matrix3Xd RQ = R * Q;
    std::ofstream each_distance;
    std::string main_path = "./";
    each_distance.open(main_path + "each_distance_result.csv");
    for (int i = 0; i < p_cols_num; i++) {
      auto dist = (P.col(i) - RQ.col(i));
      each_distance << dist.norm() << ",";
    }
    each_distance.close();
  }
  double sd = 0;
  for (int i = 0; i < p_cols_num; i++) {
    sd += W[i] * P.col(i).transpose() * P.col(i);
    sd += W[i] * Q.col(i).transpose() * Q.col(i);
  }
  Eigen::MatrixXd RH = R * H;
  sd -= 2 * RH.trace();
  if (sd < 0)
    sd = 0;
  RMSDResult rmsd_result;
  rmsd_result.rmsd_result = std::sqrt(sd / p_cols_num);
  rmsd_result.optR = R;
  return rmsd_result;
}

ConformationPair MoveToOrigin(Eigen::Matrix3Xd P, Eigen::Matrix3Xd Q,
                              std::vector<double> W) {
  Eigen::Vector3d p{0, 0, 0}, q{0, 0, 0};
  ConformationPair PQ_pair;
  double w_sum = 0.0;
  int p_cols_num = P.cols();
  for (int i = 0; i < p_cols_num; i++) {
    w_sum += W[i];
  }
  for (int i = 0; i < p_cols_num; i++) {
    p += (W[i] * P.col(i)) / w_sum;
  }
  for (int i = 0; i < p_cols_num; i++) {
    q += (W[i] * Q.col(i)) / w_sum;
  }
  Eigen::Matrix3Xd X = P.colwise() - p;
  Eigen::Matrix3Xd Y = Q.colwise() - q;
  PQ_pair.P = X;
  PQ_pair.Q = Y;
  return PQ_pair;
}
#endif
