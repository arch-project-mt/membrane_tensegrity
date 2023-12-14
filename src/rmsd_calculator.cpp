#include "rmsd.h"
#include "rmsd_io.h"
#include "rmsd_struct.h"
#include <random>

int main(int argc, char **argv) {
  int total_residue_length = 10;
  std::string main_path =
      "/Users/koyanobunsho/Desktop/architecture/membrane_tensegrity/src/";
  std::ofstream myfile;

  myfile.open(main_path + "rmsd_res.csv");
  myfile << "RMSD" << std::endl;
  Eigen::MatrixXd p, q;
  p = openMatrixData(main_path + "original_points.csv");
  q = openMatrixData(main_path + "moved_points.csv");
  std::vector<double> default_weights;
  for (int i = 0; i < total_residue_length; i++) {
    default_weights.push_back(1.0);
  }
  ConformationPair wPQ_pair = MoveToOrigin(p, q, default_weights);
  double rmsd_result = CalcRMSD(wPQ_pair.P, wPQ_pair.Q, default_weights);
  myfile << rmsd_result << std::endl;
  myfile.close();
}
