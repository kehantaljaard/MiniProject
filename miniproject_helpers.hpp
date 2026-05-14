#ifndef MINIPROJECT_HELPERS_HPP
#define MINIPROJECT_HELPERS_HPP

#include <string>
#include <vector>
#include <utility>

namespace miniproject {

// Observation grid for a single time step.
// grid[y][x] == 1 means cell (x, y) generated a detection.
using Grid = std::vector<std::vector<int>>;

// (x, y) coordinate pair. Used for wumpus positions and occupied cells.
using Coord = std::pair<int, int>;
using Trajectory = std::vector<Coord>;

// Load one detection file (e.g. data_file007.txt).
// File is a whitespace-separated matrix of 0/1 ints; row index is y, col index is x.
Grid loadDetectionFile(const std::string& path);

// Load every data_fileNNN.txt in datasetDir in sorted (lexicographic) order.
// Returns one Grid per time step. All grids are validated to share the same shape.
std::vector<Grid> loadDataset(const std::string& datasetDir);

// Load a coordinate file (wumpus_trajectory.txt or occupied-cells file).
// Each line is "x y".
Trajectory loadCoords(const std::string& path);

// Grid dimensions: {rows (height, y-extent), cols (width, x-extent)}.
std::pair<int, int> gridShape(const Grid& g);

// Print a grid to stdout, one row per line, space-separated.
void printGrid(const Grid& g);

void writeMarginal(const std::string& path, const std::vector<std::vector<double>>& posterior);

void writeTrajectory(const std::string& path, const Trajectory& trajectory);


} // namespace miniproject

#endif // MINIPROJECT_HELPERS_HPP
