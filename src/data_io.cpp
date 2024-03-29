#include "data_io.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <utility>

void DataIO::initVEF()
{
  V = Eigen::MatrixXd();
  E = Eigen::MatrixXi();
  F = Eigen::MatrixXi();
}

void DataIO::readObj(std::string filename)
{
  initVEF();

  std::ifstream infile(filename + ".obj");
  std::string line, v, f, l;

  while (std::getline(infile, line))
  {
    if (line[0] == 'v' && line[1] == ' ')
    {
      std::istringstream iss(line);
      double x, y, z;
      iss >> v >> x >> y >> z;
      V.conservativeResize(V.rows() + 1, 3);
      V.row(V.rows() - 1) << x, y, z;
    }
    if (line[0] == 'l' && line[1] == ' ')
    {
      std::istringstream iss(line);
      int x, y;
      iss >> l >> x >> y;
      E.conservativeResize(E.rows() + 1, 2);
      E.row(E.rows() - 1) << x - 1, y - 1;
    }
    else if (line[0] == 'f' && line[1] == ' ')
    {
      std::istringstream iss(line);
      int x, y, z;
      iss >> f >> x >> y >> z;
      F.conservativeResize(F.rows() + 1, 3);
      F.row(F.rows() - 1) << x - 1, y - 1, z - 1;
    }
  }
  infile.close();
}

void DataIO::writeObj(std::string filename)
{
  std::ofstream outfile(filename + ".obj");

  for (int i = 0; i < V.rows(); i++)
  {
    outfile << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
  }
  std::vector<std::pair<int, int>> output_edges;
  for (int i = 0; i < E.rows(); i++)
  {
    int edge_start = E(i, 0);
    int edge_end = E(i, 1);
    if (std::find(output_edges.begin(), output_edges.end(),
      std::make_pair(edge_start, edge_end)) != output_edges.end()) {
      continue;
    }
    outfile << "l " << E(i, 0) + 1 << " " << E(i, 1) + 1 << std::endl;
    output_edges.push_back(std::make_pair(edge_start, edge_end));
    output_edges.push_back(std::make_pair(edge_end, edge_start));
  }

  for (int i = 0; i < F.rows(); i++)
  {
    outfile << "f " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << std::endl;
  }

  outfile.close();
}

Eigen::MatrixXi DataIO::face2edge(Eigen::MatrixXi& F)
{
  Eigen::MatrixXi E;
  for (int i = 0; i < F.rows(); i++)
  {
    E.conservativeResize(E.rows() + 3, 2);
    E.row(E.rows() - 3) << F(i, 0), F(i, 1);
    E.row(E.rows() - 2) << F(i, 1), F(i, 2);
    E.row(E.rows() - 1) << F(i, 2), F(i, 0);
  }

  std::set<std::pair<int, int>> E_set;
  for (int i = 0; i < E.rows(); i++)
  {
    E_set.insert({ E(i, 0), E(i, 1) });
    E_set.insert({ E(i, 1), E(i, 0) });
  }

  E = Eigen::MatrixXi::Zero(E_set.size(), 2);
  int i = 0;
  for (auto e = E_set.begin(); e != E_set.end(); e++)
  {
    E.row(i) << e->first, e->second;
    i++;
  }

  return E;
}

void DataIO::fitInputDataStructure(
  std::vector<Eigen::Vector3d>& cube,
  std::vector<std::pair<int, int>>& edges)
{
  for (int i = 0; i < V.rows(); i++)
  {
    cube.push_back(Eigen::Vector3d(V.row(i)));
  }
  if (F.rows() > 0 && E.rows() == 0)
  {
    E = face2edge(F);
  }
  for (int i = 0; i < E.rows(); i++)
  {
    edges.push_back({ E(i, 0), E(i, 1) });
  }
}

void DataIO::fitOutputDataStructure(
  std::vector<Eigen::Vector3d>& cube,
  std::vector<std::pair<int, int>>& edges)
{
  initVEF();
  for (int i = 0; i < cube.size(); i++)
  {
    V.conservativeResize(V.rows() + 1, 3);
    V.row(V.rows() - 1) << cube[i].x(), cube[i].y(), cube[i].z();
  }

  for (int i = 0; i < edges.size(); i++)
  {
    E.conservativeResize(E.rows() + 1, 2);
    E.row(E.rows() - 1) << edges[i].first, edges[i].second;
  }
}
