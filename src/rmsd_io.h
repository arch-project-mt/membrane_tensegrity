#ifndef rmsd_io
#define rmsd_io
#include "rmsd_struct.h"

template <typename T>
bool getFileContent(std::string fileName, std::vector<T> &Flexibility) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fileName << std::endl;
    return false;
  }
  T val;
  while (in >> val) {
    Flexibility.push_back(val);
  }
  in.close();
  return true;
}

Eigen::MatrixXd openMatrixData(std::string fileToOpen) {
  std::vector<double> matrixEntries;
  std::ifstream matrixDataFile(fileToOpen);
  std::string matrixRowString;
  std::string matrixEntry;
  int matrixRowNumber = 0;
  while (getline(matrixDataFile, matrixRowString)) {
    std::stringstream matrixRowStringStream(matrixRowString);
    while (getline(matrixRowStringStream, matrixEntry, ',')) {
      matrixEntries.push_back(stod(matrixEntry));
    }
    matrixRowNumber++;
  }
  return Eigen::Map<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
      matrixEntries.data(), matrixRowNumber,
      matrixEntries.size() / matrixRowNumber);
}
bool getFileStrContent(
    std::string fileName,
    std::vector<std::pair<std::string, std::string>> &pq_pair) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fileName << std::endl;
    return false;
  }
  std::string val_p, val_q;
  while (in >> val_p >> val_q) {
    pq_pair.push_back(make_pair(val_p, val_q));
  }
  in.close();
  return true;
}
#endif
