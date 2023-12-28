#include "geometry.h"
#include "rmsd.h"
#include "data_io.hpp"

const int NUM_VERTICES = 8;


int main()
{
  // initialize data structures
  Eigen::MatrixXd V = Eigen::MatrixXd();
  Eigen::MatrixXd E = Eigen::MatrixXd();
  Eigen::MatrixXi F = Eigen::MatrixXi();  
  std::vector<Eigen::Vector3d> cube;
  std::vector<std::pair<int, int>> edges;

  dataIO::readObj("../iofiles/target_cube.obj", V, E, F);
  dataIO::fitInputDataStructure(V, E, F, cube, edges);

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
  
  // WriteToCSV(updatedCube, "../iofiles/updated_coordinates1.csv");
  dataIO::fitOutputDataStructure(V, E, F, updatedCube, edges);
  dataIO::writeObj("../iofiles/test_cube_updated1.obj", V, E, F);

  bool is_improved = true;
  int max_iteration = 10;
  int iter_num = 1;
  while (iter_num < max_iteration && is_improved)
  {
    // Eigen::MatrixXd target_structure =
    //     openMatrixData2("../iofiles/original_coordinates.csv");
    // Eigen::MatrixXd updated_structure = openMatrixData2(
    //     "../iofiles/updated_coordinates" + 
    //     std::to_string(iter_num) + ".csv");

    // Eigen::MatrixXd target_structure;
    // Eigen::MatrixXd updated_structure;
    Eigen::MatrixXd V;
    Eigen::MatrixXd E;
    Eigen::MatrixXi F;
    dataIO::readObj("../iofiles/target_cube.obj", V, E, F);  
    dataIO::fitInputDataStructure(V, E, F, cube, edges);
    Eigen::MatrixXd target_structure = V.transpose();

    std::cout << target_structure << std::endl;
    dataIO::readObj(
      "../iofiles/test_cube_updated" + std::to_string(iter_num) + ".obj", 
      V, E, F);
    dataIO::fitInputDataStructure(V, E, F, updatedCube, edges);
    Eigen::MatrixXd updatedCube = V.transpose();
    
    

    std::vector<double> default_weights;
    for (int i = 0; i < target_structure.cols(); i++)
    {
      default_weights.push_back(1.0);
    }
    MoveToOrigin(target_structure, updated_structure, default_weights);
    RMSDResult res =
        CalcRMSD(target_structure, updated_structure, default_weights);
    std::cout << res.rmsd_result << std::endl;
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
    std::cout << new_res.rmsd_result << std::endl;
    for (int i = 0; i < NUM_VERTICES; i++)
    {
      updatedCube[i] = updated_structure.col(i);
    }
    is_improved = res.rmsd_result - new_res.rmsd_result > 0;
    iter_num++;

    // WriteToCSV(updatedCube,
    //            "../iofiles/updated_coordinates" + std::to_string(iter_num) + ".csv");
    dataIO::fitOutputDataStructure(V, E, F, updatedCube, edges);
    dataIO::writeObj(
      "../iofiles/test_cube_updated" + std::to_string(iter_num) + ".obj",
       V, E, F);
  }
  std::cout << iter_num << std::endl;
  return 0;
}
