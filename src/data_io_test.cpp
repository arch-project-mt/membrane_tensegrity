#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "geometry.h"
#include "rmsd.h"
#include "data_io.hpp"


int main() {
  std::cout << "Hello, World!" << std::endl;

  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXi F;

  dataIO::readObj("../iofiles/bunny.obj", V, E, F);
  F = dataIO::edge2face(E);
  E = Eigen::MatrixXi();
  dataIO::writeObj("../iofiles/bunny_e2f.obj", V, E, F);

  dataIO::readObj("../iofiles/bunny.obj", V, E, F);
  E = dataIO::face2edge(F);
  F = Eigen::MatrixXi();
  dataIO::writeObj("../iofiles/bunny_f2e.obj", V, E, F);

  // auto adj = dataIO::makeAdjacencyMatrix(V, E, F);

  // dataIO::writeObj("../iofiles/bunny2.obj", V, E, F);

  return 0;
}
