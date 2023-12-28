#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/SparseLU"

// std
#include <iostream>
#include <fstream>
#include <string>


namespace dataIO{

    void initVEF(Eigen::MatrixXd &V, Eigen::MatrixXd &E, Eigen::MatrixXi &F);

    void readObj(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXd &E, Eigen::MatrixXi &F);
    
    void writeObj(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXd &E, Eigen::MatrixXi &F);

    Eigen::SparseMatrix<int> makeAdjacencyMatrix(Eigen::MatrixXd &V, Eigen::MatrixXd &E, Eigen::MatrixXi &F);

    void fitInputDataStructure(Eigen::MatrixXd &V,
                      Eigen::MatrixXd &E,
                      Eigen::MatrixXi &F,
                      std::vector<Eigen::Vector3d> &cube,
                      std::vector<std::pair<int, int>> &edges);
    
    void fitOutputDataStructure(Eigen::MatrixXd &V,
                      Eigen::MatrixXd &E,
                      Eigen::MatrixXi &F,
                      std::vector<Eigen::Vector3d> &cube,
                      std::vector<std::pair<int, int>> &edges);
}
