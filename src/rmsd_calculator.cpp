#include "rmsd.h"
#include <random>

int main(int argc, char **argv) {
  std::string main_path = "./";
  std::ofstream myfile;

  myfile.open(main_path + "rmsd_res.csv");
  myfile << "RMSD_1,"
         << "RMSD_2,"
         << "RMSD_3" << std::endl;
  Eigen::MatrixXd p, q, p_2, q_2, p_3, q_3;
  // openMatrixData assumes the
  // csv is 3 rows n columns, where n is the number
  // of vertices
  p = openMatrixData(main_path +
                     "original_points.csv"); // Specify the path of coordinates
                                             // of one structure.
  q = openMatrixData(main_path +
                     "moved_points.csv"); // Specify the path of coordinates of
                                          // the other structure
  std::vector<double> default_weights;
  for (int i = 0; i < p.cols(); i++) {
    default_weights.push_back(1.0);
  }
  ConformationPair wPQ_pair = MoveToOrigin(p, q, default_weights);
  RMSDResult rmsd_result =
      CalcRMSD(wPQ_pair.P, wPQ_pair.Q, default_weights, false);
  myfile << rmsd_result.rmsd_result << ",";
  // openMatrixData2 assumes the
  // csv is n rows 3 columns, where n is the number
  // of vertices. You can choose either way!
  p_2 = openMatrixData2(
      main_path +
      "coord_ex_20231216-155855.csv"); // Specify the path of coordinates
                                       // of one structure.
  q_2 = openMatrixData2(
      main_path +
      "coord_ex_20231216-155618.csv"); // Specify the path of coordinates
                                       // of the other structure.
  std::vector<double> default_weights2;
  for (int i = 0; i < p_2.cols(); i++) {
    default_weights2.push_back(1.0);
  }
  ConformationPair wPQ_pair_2 = MoveToOrigin(p_2, q_2, default_weights2);
  RMSDResult rmsd_result_2 =
      CalcRMSD(wPQ_pair_2.P, wPQ_pair_2.Q, default_weights2,
               true); // Set true when you want to output the
                      // distance between each vertex
  myfile << rmsd_result_2.rmsd_result << ",";
  p_3 = openMatrixData2(
      main_path + "original_coordinates.csv"); // Specify the path of
                                               // coordinates of one structure.
  q_3 = openMatrixData2(
      main_path + "updated_coordinates.csv"); // Specify the path of coordinates
                                              // of the other structure.
  std::vector<double> default_weights3;
  for (int i = 0; i < p_3.cols(); i++) {
    default_weights3.push_back(1.0);
  }
  ConformationPair wPQ_pair_3 = MoveToOrigin(p_3, q_3, default_weights3);
  RMSDResult rmsd_result_3 =
      CalcRMSD(wPQ_pair_3.P, wPQ_pair_3.Q, default_weights3,
               true); // Set true when you want to output the
                      // distance between each vertex
  myfile << rmsd_result_3.rmsd_result << std::endl;
  myfile.close();
}
