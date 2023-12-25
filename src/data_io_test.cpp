#include <iostream>
#include "data_io.hpp"

int main() {
  std::cout << "Hello, World!" << std::endl;

  Eigen::Matrix V = Eigen::MatrixXd();
  Eigen::Matrix E = Eigen::MatrixXd();
  Eigen::Matrix F = Eigen::MatrixXi();
  dataIO::readObj("../bunny.obj", V, E, F);

  auto adj = dataIO::makeAdjacencyMatrix(V, E, F);

  dataIO::writeObj("bunny2.obj", V, E, F);

  return 0;
}
