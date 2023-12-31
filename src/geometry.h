#ifndef geometry
#define geometry
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/SparseLU"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

struct VertexInfo {
  Eigen::Vector3d position;
  int index;
};

struct VertexInfoHash {
  size_t operator()(const VertexInfo &v) const {
    size_t h1 = std::hash<double>()(v.position.x());
    size_t h2 = std::hash<double>()(v.position.y());
    size_t h3 = std::hash<double>()(v.position.z());
    size_t h4 = std::hash<int>()(v.index);
    return h1 ^ h2 ^ h3 ^ h4; // or use another combining method
  }
};

struct VertexInfoEqual {
  bool operator()(const VertexInfo &a, const VertexInfo &b) const {
    return a.index == b.index && a.position.isApprox(b.position);
  }
};
using AdjacencyType = std::unordered_map<VertexInfo, std::vector<VertexInfo>,
                                         VertexInfoHash, VertexInfoEqual>;
double CalculateAngle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
  /*
  Since cos(θ) = (v1 \cdot v2) / (|v1| |v2|),
    θ = cos^{-1}((v1 \cdot v2) / (|v1| |v2|))
  */
  double v1_dot_v2 = v1.dot(v2);
  double magnitudeProduct = v1.norm() * v2.norm();
  if (magnitudeProduct == 0) {
    return 0;
  }
  return std::acos(v1_dot_v2 / magnitudeProduct);
}

struct VectorHash {
  size_t operator()(const Eigen::Vector3d &v) const {
    return std::hash<double>()(v.x()) ^ std::hash<double>()(v.y()) ^
           std::hash<double>()(v.z());
  }
};

struct VectorEqual {
  bool operator()(const Eigen::Vector3d &a, const Eigen::Vector3d &b) const {
    return a.isApprox(b);
  }
};

double cotangent(double angle) { return 1.0 / tan(angle); }

Eigen::Vector3d computeThetaPvP(const Eigen::Vector3d &v_p,
                                const std::vector<VertexInfo> &neighbors,
                                const std::vector<double> &gamma_qp,
                                const std::vector<double> &xi_qp) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (size_t i = 0; i < neighbors.size(); ++i) {
    Eigen::Vector3d v_q = neighbors[i].position;
    double distance =
        (v_p - v_q).norm(); // distance系統の実装が正しいのかを後で確認する
    res +=
        (cotangent(gamma_qp[i]) + cotangent(xi_qp[i])) / distance * (v_q - v_p);
  }
  return res;
}

// Function to compute the second expression θ_p(v_q)
Eigen::Vector3d computeThetaPvQ(const Eigen::Vector3d &v_p,
                                const Eigen::Vector3d &v_q,
                                const Eigen::Vector3d &v_q_plus,
                                const Eigen::Vector3d &v_q_minus,
                                double gamma_qp, double xi_qp) {
  double distance_pq = (v_p - v_q).squaredNorm();
  double distance_q_plus_q = (v_q_plus - v_q).norm();
  double distance_q_minus_q = (v_q_minus - v_q).norm();
  Eigen::Vector3d term1 =
      (cotangent(xi_qp) + cotangent(gamma_qp)) / distance_pq * (v_p - v_q);
  Eigen::Vector3d term2 =
      (v_q_plus - v_q) / (distance_pq * distance_q_plus_q * sin(xi_qp));
  Eigen::Vector3d term3 =
      (v_q_minus - v_q) / (distance_pq * distance_q_minus_q * sin(gamma_qp));
  return term1 - term2 - term3;
}

Eigen::SparseMatrix<double> calculateAAT(Eigen::MatrixXd &A) {
  Eigen::SparseMatrix<double> sparseA = A.sparseView();
  Eigen::SparseMatrix<double> AAT = sparseA * sparseA.transpose();
  return AAT;
}

void replaceNaNsWithZeros(Eigen::SparseMatrix<double> &matrix) {
  for (int k = 0; k < matrix.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
      if (std::isnan(it.valueRef())) {
        it.valueRef() = 0.0;
      }
    }
  }
}

Eigen::MatrixXd computeD(Eigen::SparseMatrix<double> &AAT_x,
                         Eigen::SparseMatrix<double> &AAT_y,
                         Eigen::SparseMatrix<double> &AAT_z,
                         Eigen::MatrixXd &A_x, Eigen::MatrixXd &A_y,
                         Eigen::MatrixXd &A_z, const Eigen::MatrixXd &b) {
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  Eigen::SparseMatrix<double> AAT = AAT_x + AAT_y + AAT_z;
  replaceNaNsWithZeros(AAT);
  solver.compute(AAT);

  if (solver.info() != Eigen::Success) {
    std::cout << "Factorization failed" << std::endl;
  }
  // Solve (AA^T)x = b
  Eigen::VectorXd x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    std::cout << "Solving failed" << std::endl;
  }
  // Compute d = A^T x
  Eigen::MatrixXd d(A_x.rows(), 3);
  d.col(0) = A_x.transpose() * x;
  d.col(1) = A_y.transpose() * x;
  d.col(2) = A_z.transpose() * x;
  return d;
}
void SortVerticesCounterClockwise(const Eigen::Vector3d &center,
                                  std::vector<VertexInfo> &neighbors) {
  std::sort(neighbors.begin(), neighbors.end(),
            [&center](const VertexInfo &a, const VertexInfo &b) -> bool {
              Eigen::Vector3d dirA = a.position - center;
              Eigen::Vector3d dirB = b.position - center;
              double angleA = atan2(dirA.y(), dirA.x());
              double angleB = atan2(dirB.y(), dirB.x());
              return angleA < angleB;
            });
}

double SumOfAnglesAtVertex(const Eigen::Vector3d &vertex,
                           const std::vector<VertexInfo> &neighbors) {
  double sumOfAngles = 0.0;
  int numNeighbors = neighbors.size();
  for (int i = 0; i < numNeighbors; ++i) {
    Eigen::Vector3d edge1 = neighbors[i].position - vertex;
    Eigen::Vector3d edge2 = neighbors[(i + 1) % numNeighbors].position - vertex;
    sumOfAngles += CalculateAngle(edge1, edge2);
  }
  return sumOfAngles;
}
void WriteToCSV(const std::vector<Eigen::Vector3d> &data,
                const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }
  for (const auto &vec : data) {
    file << vec.x() << "," << vec.y() << "," << vec.z() << std::endl;
  }

  file.close();
}

void writeIndicesToCSV(const AdjacencyType &map, const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  for (const auto &pair : map) {
    file << pair.first.index << ",";
    for (const VertexInfo &v : pair.second) {
      file << v.index << ",";
    }
    file << "\n";
  }

  file.close();
}
void ComputeMatrixAndVector(Eigen::MatrixXd &A_x, Eigen::MatrixXd &A_y,
                            Eigen::MatrixXd &A_z, Eigen::VectorXd &b,
                            AdjacencyType &adjacencyList) {
  A_x.setZero();
  A_y.setZero();
  A_z.setZero();
  b.setZero();
  for (const auto &vertex : adjacencyList) {
    const VertexInfo &p = vertex.first;
    std::vector<double> gamma_qp_vec, xi_qp_vec;
    for (const VertexInfo &q : vertex.second) {
      if (p.position.isApprox(q.position)) {
        continue;
      }
      auto &connectedVerticesq = adjacencyList[q];
      if (connectedVerticesq.size() < 2) {
        continue;
      }
      const Eigen::Vector3d &q_minus = connectedVerticesq[0].position;
      const Eigen::Vector3d &q_plus = connectedVerticesq[1].position;
      Eigen::Vector3d v1 = q_minus - q.position;
      Eigen::Vector3d v2 = q_plus - q.position;
      double gamma_qp = CalculateAngle(v2, p.position - q.position);
      double xi_qp = CalculateAngle(v1, p.position - q.position);
      gamma_qp_vec.push_back(gamma_qp);
      xi_qp_vec.push_back(xi_qp);
      auto &connectedVerticesp = adjacencyList[p];
      if (connectedVerticesp.size() < 2) {
        continue;
      }
      const Eigen::Vector3d &p_minus = connectedVerticesp[0].position;
      const Eigen::Vector3d &p_plus = connectedVerticesp[1].position;
      Eigen::Vector3d v1p = p_minus - p.position;
      Eigen::Vector3d v2p = p_plus - p.position;
      double gamma_pq = CalculateAngle(v2p, q.position - p.position);
      double xi_pq = CalculateAngle(v1p, q.position - p.position);
      if (gamma_pq != 0 && xi_pq != 0) {
        Eigen::Vector3d thetaQvP = computeThetaPvQ(
            q.position, p.position, p_plus, p_minus, gamma_pq, xi_pq);
        A_x(p.index, q.index) = thetaQvP[0];
        A_y(p.index, q.index) = thetaQvP[1];
        A_z(p.index, q.index) = thetaQvP[2];
      }
    }
    if (gamma_qp_vec.size() > 0) {
      Eigen::Vector3d thetaPvP = computeThetaPvP(
          vertex.first.position, vertex.second, gamma_qp_vec, xi_qp_vec);
      A_x(p.index, p.index) = thetaPvP[0];
      A_y(p.index, p.index) = thetaPvP[1];
      A_z(p.index, p.index) = thetaPvP[2];
    }
    b(p.index) = 2 * M_PI - SumOfAnglesAtVertex(p.position, vertex.second);
  }
}

void MakeAdjacencyList(const std::vector<Eigen::Vector3d> &cube,
                       const std::vector<std::pair<int, int>> &edges,
                       AdjacencyType &adjacencyList) {
  for (const auto &edge : edges) {
    VertexInfo connectedVertexInfoFrom, connectedVertexInfoTo;
    connectedVertexInfoFrom.position = cube[edge.first];
    connectedVertexInfoFrom.index = edge.first;
    connectedVertexInfoTo.position = cube[edge.second];
    connectedVertexInfoTo.index = edge.second;
    adjacencyList[connectedVertexInfoFrom].push_back(connectedVertexInfoTo);
  }
}

#endif
