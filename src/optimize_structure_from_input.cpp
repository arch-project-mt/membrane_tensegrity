#include "geometry.h"
#include "rmsd.h"
#include "data_io.hpp"
const int NUM_VERTICES = 8;

int main() {
//  std::vector <Eigen::Vector3d> cube;
//  cube.push_back(Eigen::Vector3d(0, 0, 0));
//  cube.push_back(Eigen::Vector3d(0, 1, 0));
//  cube.push_back(Eigen::Vector3d(0, 0, 1));
//  cube.push_back(Eigen::Vector3d(1, 0, 0));
//  cube.push_back(Eigen::Vector3d(1, 1, 0));
//  cube.push_back(Eigen::Vector3d(1, 0, 1));
//  cube.push_back(Eigen::Vector3d(0, 1, 1));
//  cube.push_back(Eigen::Vector3d(1, 1, 1));
//  WriteToCSV(cube, "original_coordinates.csv");
//  std::vector <std::pair<int, int>> edges = {
//          {0, 3},
//          {0, 2},
//          {0, 1},
//          {1, 0},
//          {1, 4},
//          {1, 6},
//          {2, 0},
//          {2, 5},
//          {2, 6},
//          {3, 0},
//          {3, 4},
//          {3, 5},
//          {4, 1},
//          {4, 3},
//          {4, 7},
//          {5, 2},
//          {5, 3},
//          {5, 7},
//          {6, 1},
//          {6, 2},
//          {6, 7},
//          {7, 4},
//          {7, 5},
//          {7, 6}};

  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXi F;
  dataIO::readObj("cube_test.obj", V, E, F);
  std::vector <Eigen::Vector3d> cube;

  for (int i = 0; i < V.rows(); i++) {
    cube.push_back(Eigen::Vector3d(V.row(i)));
  }
  std::vector <std::pair<int, int>> edges;
  for (int i = 0; i < E.rows(); i++) {
    edges.push_back(std::make_pair(E(i, 0), E(i, 1)));
  }

  AdjacencyType adjacencyList;
  MakeAdjacencyList(cube, edges, adjacencyList);
  writeIndicesToCSV(adjacencyList, "original_edges.csv");

  for (auto &vertex: adjacencyList) {
    const VertexInfo &p = vertex.first;
    SortVerticesCounterClockwise(p.position, vertex.second);
  }
  Eigen::MatrixXd A_x(NUM_VERTICES, NUM_VERTICES);
  Eigen::MatrixXd A_y(NUM_VERTICES, NUM_VERTICES);
  Eigen::MatrixXd A_z(NUM_VERTICES, NUM_VERTICES);
  Eigen::VectorXd b(NUM_VERTICES);
  ComputeMatrixAndVector(A_x, A_y, A_z, b, adjacencyList);
  Eigen::SparseMatrix<double> AAT_x = calculateAAT(A_x);
  Eigen::SparseMatrix<double> AAT_y = calculateAAT(A_y);
  Eigen::SparseMatrix<double> AAT_z = calculateAAT(A_z);
  Eigen::MatrixXd d = computeD(AAT_x, AAT_y, AAT_z, A_x, A_y, A_z, b);
  std::vector <Eigen::Vector3d> updatedCube;
  for (int i = 0; i < NUM_VERTICES; i++) {
    Eigen::Vector3d di = d.row(i);
    updatedCube.push_back(cube[i] + di);
  }

  WriteToCSV(updatedCube, "updated_coordinates1.csv");
  bool is_improved = true;
  int max_iteration = 10;
  int iter_num = 1;
  while (iter_num < max_iteration && is_improved) {
    Eigen::MatrixXd target_structure =
            openMatrixData2("original_coordinates.csv");
    Eigen::MatrixXd updated_structure = openMatrixData2(
            "updated_coordinates" + std::to_string(iter_num) + ".csv");
    std::vector<double> default_weights;
    for (int i = 0; i < target_structure.cols(); i++) {
      default_weights.push_back(1.0);
    }
    MoveToOrigin(target_structure, updated_structure, default_weights);
    RMSDResult res =
            CalcRMSD(target_structure, updated_structure, default_weights);
    std::cout << res.rmsd_result << std::endl;
    for (int i = 0; i < updated_structure.cols(); i++) {
      updatedCube[i] = res.optR * updated_structure.col(i);
    }
    AdjacencyType updatedAdjacencyList;
    MakeAdjacencyList(updatedCube, edges, updatedAdjacencyList);
    ComputeMatrixAndVector(A_x, A_y, A_z, b, adjacencyList);
    AAT_x = calculateAAT(A_x);
    AAT_y = calculateAAT(A_y);
    AAT_z = calculateAAT(A_z);
    Eigen::VectorXd c(NUM_VERTICES);
    Eigen::VectorXd cube_x(NUM_VERTICES), cube_y(NUM_VERTICES),
            cube_z(NUM_VERTICES), updatedCube_x(NUM_VERTICES),
            updatedCube_y(NUM_VERTICES), updatedCube_z(NUM_VERTICES);
    for (int i = 0; i < NUM_VERTICES; i++) {
      cube_x(i) = cube[i][0];
      cube_y(i) = cube[i][1];
      cube_z(i) = cube[i][2];
      updatedCube_x(i) = updatedCube[i][0];
      updatedCube_y(i) = updatedCube[i][1];
      updatedCube_z(i) = updatedCube[i][2];
    }
    for (int i = 0; i < NUM_VERTICES; i++) {
      c(i) = b(i) + A_x.row(i) * (cube_x - updatedCube_x);
      c(i) = b(i) + A_y.row(i) * (cube_y - updatedCube_y);
      c(i) = b(i) + A_z.row(i) * (cube_z - updatedCube_z);
    }
    d = computeD(AAT_x, AAT_y, AAT_z, A_x, A_y, A_z, c);
    updated_structure.transpose() = target_structure.transpose() + d;
    MoveToOrigin(target_structure, updated_structure, default_weights);
    RMSDResult new_res =
            CalcRMSD(target_structure, updated_structure, default_weights);
    std::cout << new_res.rmsd_result << std::endl;
    for (int i = 0; i < NUM_VERTICES; i++) {
      updatedCube[i] = updated_structure.col(i);
    }
    is_improved = res.rmsd_result - new_res.rmsd_result > 0;
    iter_num++;
    WriteToCSV(updatedCube,
               "updated_coordinates" + std::to_string(iter_num) + ".csv");
  }
  std::cout << iter_num << std::endl;
  return 0;
}
