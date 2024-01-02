#include "geometry.h"
#include "rmsd.h"
#include "data_io.hpp"


int main()
{
  Eigen::MatrixXd V = Eigen::MatrixXd();
  Eigen::MatrixXi E = Eigen::MatrixXi();
  Eigen::MatrixXi F = Eigen::MatrixXi();
  std::vector<Eigen::Vector3d> cube;
  std::vector<std::pair<int, int>> edges;
  std::string in_filename = "../iofiles/target_cube.obj";
  std::string out_filename = "../iofiles/updated_cube";

  dataIO::readObj(in_filename, V, E, F);
  dataIO::fitInputDataStructure(V, E, F, cube, edges);
  const int NUM_VERTICES = cube.size();

  AdjacencyType adjacencyList;
  MakeAdjacencyList(cube, edges, adjacencyList);
  writeIndicesToCSV(adjacencyList, "../iofiles/original_edges.csv");
  for (auto &vertex : adjacencyList)
  {
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
  std::vector<Eigen::Vector3d> updatedCube;
  for (int i = 0; i < NUM_VERTICES; i++)
  {
    Eigen::Vector3d di = d.row(i);
    updatedCube.push_back(cube[i] + di);
  }
  
  dataIO::fitOutputDataStructure(V, E, F, updatedCube, edges);
  dataIO::writeObj(out_filename + "1.obj", V, E, F);

  dataIO::readObj(in_filename, V, E, F);
  Eigen::MatrixXd target_structure = V.transpose();

  bool is_improved = true;
  int max_iteration = 10;
  int iter_num = 1;
  while (iter_num < max_iteration && is_improved)
  {
    dataIO::readObj(
            out_filename + std::to_string(iter_num) + ".obj",
            V, E, F);
    Eigen::MatrixXd updated_structure = V.transpose();
    
  dataIO::readObj(in_filename, V, E, F);
  Eigen::MatrixXd target_structure = V.transpose();

    std::vector<double> default_weights;
    for (int i = 0; i < target_structure.cols(); i++)
    {
      default_weights.push_back(1.0);
    }
    MoveToOrigin(target_structure, updated_structure, default_weights);
    RMSDResult res =
        CalcRMSD(target_structure, updated_structure, default_weights);
    for (int i = 0; i < updated_structure.cols(); i++)
    {
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
    for (int i = 0; i < NUM_VERTICES; i++)
    {
      cube_x(i) = cube[i][0];
      cube_y(i) = cube[i][1];
      cube_z(i) = cube[i][2];
      updatedCube_x(i) = updatedCube[i][0];
      updatedCube_y(i) = updatedCube[i][1];
      updatedCube_z(i) = updatedCube[i][2];
    }
    for (int i = 0; i < NUM_VERTICES; i++)
    {
      c(i) = b(i) + A_x.row(i) * (cube_x - updatedCube_x);
      c(i) = b(i) + A_y.row(i) * (cube_y - updatedCube_y);
      c(i) = b(i) + A_z.row(i) * (cube_z - updatedCube_z);
    }
    d = computeD(AAT_x, AAT_y, AAT_z, A_x, A_y, A_z, c);
    updated_structure.transpose() = target_structure.transpose() + d;
    MoveToOrigin(target_structure, updated_structure, default_weights);
    RMSDResult new_res =
        CalcRMSD(target_structure, updated_structure, default_weights);

    for (int i = 0; i < NUM_VERTICES; i++)
    {
      updatedCube[i] = updated_structure.col(i);
    }
    float diff = res.rmsd_result - new_res.rmsd_result;
    is_improved = diff > 0;
    iter_num++;

    // output
    std::cout << "iter" << iter_num << " : rmsd = ";
    std::cout << new_res.rmsd_result << " | diff = " << diff << std::endl;

    dataIO::fitOutputDataStructure(V, E, F, updatedCube, edges);
    dataIO::writeObj(
      out_filename + std::to_string(iter_num) + ".obj",
       V, E, F);
  }
  std::cout << "Optimization finished in " << iter_num << " iterations" << std::endl;
  return 0;
}
