#pragma once

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/SparseLU"

#include <iostream>
#include <fstream>
#include <string>


class DataIO {
public:
  DataIO() {}
  ~DataIO() {}

  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXi F;

  void readObj(std::string filename);
  void writeObj(std::string filename);

  // VEF to cube edges
  void fitInputDataStructure(
    std::vector<Eigen::Vector3d>& cube,
    std::vector<std::pair<int, int>>& edges);

  // cube edges to VEF
  void fitOutputDataStructure(
    std::vector<Eigen::Vector3d>& cube,
    std::vector<std::pair<int, int>>& edges);

private:
  void initVEF();

  Eigen::MatrixXi face2edge(Eigen::MatrixXi& F);
};
