#include "miniproject_helpers.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace miniproject {

namespace fs = std::filesystem;

Grid loadDetectionFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("loadDetectionFile: cannot open " + path);
  }

  Grid grid;
  std::string line;
  size_t expectedCols = 0;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::vector<int> row;
    int v;
    while (ls >> v) row.push_back(v);
    if (row.empty()) continue;
    if (grid.empty()) {
      expectedCols = row.size();
    } else if (row.size() != expectedCols) {
      throw std::runtime_error("loadDetectionFile: ragged row in " + path);
    }
    grid.push_back(std::move(row));
  }
  if (grid.empty()) {
    throw std::runtime_error("loadDetectionFile: empty file " + path);
  }
  return grid;
}

std::vector<Grid> loadDataset(const std::string& datasetDir) {
  std::vector<fs::path> files;
  for (const auto& entry : fs::directory_iterator(datasetDir)) {
    if (!entry.is_regular_file()) continue;
    const auto& p = entry.path();
    const std::string name = p.filename().string();
    // Match data_file*.txt; skip GIFs and anything else.
    if (name.rfind("data_file", 0) == 0 && p.extension() == ".txt") {
      files.push_back(p);
    }
  }
  if (files.empty()) {
    throw std::runtime_error("loadDataset: no data_file*.txt in " + datasetDir);
  }
  std::sort(files.begin(), files.end());

  std::vector<Grid> frames;
  frames.reserve(files.size());
  std::pair<int, int> shape{-1, -1};
  for (const auto& f : files) {
    Grid g = loadDetectionFile(f.string());
    auto s = gridShape(g);
    if (shape.first < 0) {
      shape = s;
    } else if (s != shape) {
      throw std::runtime_error("loadDataset: inconsistent grid shape at " + f.string());
    }
    frames.push_back(std::move(g));
  }
  return frames;
}

Trajectory loadCoords(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("loadCoords: cannot open " + path);
  }
  Trajectory coords;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    int x, y;
    if (ls >> x >> y) coords.emplace_back(x, y);
  }
  if (coords.empty()) {
    throw std::runtime_error("loadCoords: no coordinates parsed from " + path);
  }
  return coords;
}

std::pair<int, int> gridShape(const Grid& g) {
  if (g.empty()) return {0, 0};
  return {static_cast<int>(g.size()), static_cast<int>(g.front().size())};
}

void printGrid(const Grid& g) {
  for (const auto& row : g) {
    for (size_t x = 0; x < row.size(); ++x) {
      std::cout << row[x] << (x + 1 < row.size() ? ' ' : '\n');
    }
  }
}

void writeMarginal(const std::string& path, const std::vector<std::vector<double>>& posterior) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("writeMarginal: cannot open " + path);
  }
  out << std::fixed << std::setprecision(6);
  for (const auto& row : posterior) {
    for (size_t i = 0; i < row.size(); ++i) {
      out << row[i] << (i + 1 < row.size() ? ' ' : '\n');
    }
  }
}

void writeTrajectory(const std::string& path, const Trajectory& trajectory) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("writeTrajectory: cannot open " + path);
  }
  for (const auto& [x, y] : trajectory) {
    out << x << ' ' << y << '\n';
  }
}

} // namespace miniproject
