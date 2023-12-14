#ifndef rmsd_struct
#define rmsd_struct
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include <cmath>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct ConformationPair {
  Eigen::Matrix3Xd P, Q;
};

struct AcNumberPQpdb {
  std::string ac_number;
  std::string pdb_id;
};
#endif
