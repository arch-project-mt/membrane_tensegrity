#include <iostream>
#include <fstream>
#include <string>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>

#include "rmsd.h"
#include "geometry.h"


namespace dataIO{
    void readObj(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &E, Eigen::MatrixXi &F);

    void writeObj(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &E, Eigen::MatrixXi &F);

    Eigen::SparseMatrix<int> makeAdjacencyMatrix(Eigen::MatrixXd &V, Eigen::MatrixXi &E, Eigen::MatrixXi &F);
}
