/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_logging.hpp"  // initLogging
#include "prlite_testing.hpp"

// emdw headers
#include "emdw.hpp"
#include "discretetable.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"

// project helpers
#include "miniproject_helpers.hpp"

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <cctype>  // toupper
#include <string>  // string
#include <memory>
#include <set>
#include <map>
#include <algorithm>
#include <limits>
#include <random>
#include <iomanip>
#include <filesystem>
#include <sstream>

using namespace std;
using namespace emdw;
using namespace miniproject;
namespace fs = std::filesystem;

struct RunConfig {
  string datasetDir;
  string outputDir;
  double pw;
  double pc;
  string cgType;  // LTRIP, JTREE, or BETHE
};

int main(int argc, char *argv[]) {

  initLogging(argv[0]);
  prlite::TestCase::runAllTests();
  
  if (argc < 5) {
    cerr << "Usage: " << argv[0]
         << " <datasetDir> <outputDir> <pw> <pc> [LTRIP|JTREE|BETHE]" << endl;
    return 1;
  }

  RunConfig cfg;
  cfg.datasetDir = argv[1];
  cfg.outputDir  = argv[2];
  cfg.pw         = stod(argv[3]);
  cfg.pc         = stod(argv[4]);
  cfg.cgType     = (argc > 5) ? argv[5] : "LTRIP";

  unsigned seedVal = emdw::randomEngine.getSeedVal();
  // cout << seedVal << endl;
  emdw::randomEngine.setSeedVal(seedVal);

  ClusterGraph* cgPtr = nullptr;
  try {

    //*********************************************************
    // Predefine types and constants
    //*********************************************************

    typedef int T;
    typedef DiscreteTable<T> DT;
    double defProb = 0.0;


    fs::create_directories(cfg.outputDir); 
    //*********************************************************
    // Load detection data
    //*********************************************************
    cout << "Loading dataset from " << cfg.datasetDir << " ... " << flush;
    vector<Grid> detections = loadDataset(cfg.datasetDir);
    auto [H, W] = gridShape(detections.front());
    const int lastT = static_cast<int>(detections.size()) - 1;
    const int G     = W * H;
    cout << "loaded " << (lastT + 1) << " time steps, grid "
         << W << " x " << H << " (G = " << G << "), parameters pw = "
         << cfg.pw << ", pc = " << cfg.pc << endl;


    //*********************************************************
    // Define the RV ids for X_0, ..., X_T
    //*********************************************************
    vector<RVIdType> X(lastT + 1);
    for (int t = 0; t <= lastT; ++t) {
      X[t] = static_cast<RVIdType>(t);
    }
    //*********************************************************
    // Build the shared cell domain {0, 1, ..., G-1}
    //*********************************************************
    rcptr<vector<T>> cellDom(new vector<T>(G));
    for (int i = 0; i < G; ++i) (*cellDom)[i] = i;
 
    //*********************************************************
    // Prior P(X_0)
    //*********************************************************
    map<vector<T>, FProb> emptyProbs;
    rcptr<Factor> prior = uniqptr<DT> (new DT({X[0]}, {cellDom}, 1.0 / static_cast<double>(G), emptyProbs));    

    //*********************************************************
    // Transition factor P(X_{t+1} | X_t) with clipping
    //*********************************************************
    const int directions[4][2] = {
      {-1, 0}, //up
      {1, 0},  //down
      {0, -1}, //left
      {0, 1},  //right
    };

    auto rowOf  = [W](int idx) { return idx / W; };
    auto colOf  = [W](int idx) { return idx % W; };
    auto cellOf = [W](int row, int col) { return row * W + col; };
    auto clipMove =
        [W, H, &cellOf](int row, int col, int drow, int dcol) {
          int nr = row + drow;
          int nc = col + dcol;
          if (nr < 0 || nr >= H) nr = row;
          if (nc < 0 || nc >= W) nc = col;
          return cellOf(nr, nc);
        };
      
    map<vector<T>, FProb> transProbs;
    for (int i = 0; i < G; ++i) {
      const int r = rowOf(i);
      const int c = colOf(i);
      for (int d = 0; d < 4; ++d) {
        const int j = clipMove(r, c, directions[d][0], directions[d][1]);
        transProbs[{T(i), T(j)}].prob += 0.25;
      }
    }

    //*********************************************************
    // Observation factors: phi_t(X_t = j)
    //*********************************************************
    const double phiOn  = cfg.pw / cfg.pc;
    const double phiOff = (1.0 - cfg.pw) / (1.0 - cfg.pc);

    //*********************************************************
    // Assemble factors
    //*********************************************************

    vector<rcptr<Factor>> factors;
    factors.reserve(1 + lastT + (lastT + 1));
    factors.push_back(prior);

    // Transition factors; one per step. Sparse map is shared.
    for (int t = 0; t < lastT; ++t) {
      rcptr<Factor> trans = uniqptr<DT>(new DT(
          {X[t], X[t + 1]}, {cellDom, cellDom}, defProb, transProbs));
      factors.push_back(trans);
    }

    // Observation factors; one per time step, populated from detections.
    for (int t = 0; t <= lastT; ++t) {
      const Grid& det = detections[t];
      map<vector<T>, FProb> phiProbs;
      for (int row = 0; row < H; ++row) {
        for (int col = 0; col < W; ++col) {
          const int j = cellOf(row, col);
          const int d = det[row][col];
          phiProbs[{T(j)}] = (d == 1) ? phiOn : phiOff;
        }
      }
      rcptr<Factor> phi = uniqptr<DT>(new DT(
          {X[t]}, {cellDom}, defProb, phiProbs));
      factors.push_back(phi);
    }
     
    cout << "Built " << factors.size() << " factors: 1 prior + "
         << lastT << " transitions + " << (lastT + 1) << " observations."
         << endl;


    //*********************************************************
    // Build cluster graph
    //*********************************************************
    cout << "Constructing cluster graph (" << cfg.cgType << ") ... " << flush;
    map<RVIdType, AnyType> obsv;  // no RV-level observations; evidence baked into phi_t
 
    if (cfg.cgType == "JTREE") {
      cgPtr = new ClusterGraph(ClusterGraph::JTREE, factors, obsv);
    } else if (cfg.cgType == "BETHE") {
      cgPtr = new ClusterGraph(ClusterGraph::BETHE, factors, obsv);
    } else {
      cgPtr = new ClusterGraph(ClusterGraph::LTRIP, factors, obsv);
    }
    ClusterGraph& cg = *cgPtr;
    cout << "done." << endl;

    //*********************************************************
    // Run belief propagation
    //*********************************************************
    cout << "Running loopyBP_CG ... " << flush;
    map<Idx2, rcptr<Factor>> msgs;
    MessageQueue msgQ;
    unsigned nMsgs = loopyBP_CG(cg, msgs, msgQ);
    cout << "done. " << nMsgs << " messages sent." << endl;    


    //*********************************************************
    // Extract time step marginals and argmax trajectory
    //*********************************************************
    cout << "Extracting marginals and writing output ... " << flush;
    Trajectory mapTrajectory;
    mapTrajectory.reserve(lastT + 1);
 
    for (int t = 0; t <= lastT; ++t) {
      rcptr<Factor> q = queryLBP_CG(cg, msgs, {X[t]})->normalize();
 
      vector<vector<double>> posterior(H, vector<double>(W, 0.0));
      int bestRow = 0, bestCol = 0;
      double bestProb = -1.0;
 
      for (int row = 0; row < H; ++row) {
        for (int col = 0; col < W; ++col) {
          const int j = cellOf(row, col);
          rcptr<Factor> reduced = q->observeAndReduce({X[t]}, {T(j)});
          const double p = reduced->potentialAt(
              emdw::RVIds{}, emdw::RVVals{}, true);
          posterior[row][col] = p;
          if (p > bestProb) {
            bestProb = p;
            bestRow  = row;
            bestCol  = col;
          }
        }
      }
 
      ostringstream nameStream;
      nameStream << cfg.outputDir << "/marginal_"
                 << setw(3) << setfill('0') << t << ".txt";
      writeMarginal(nameStream.str(), posterior);
 
      // Trajectory uses (x, y) = (col, row), matching the ground-truth format.
      mapTrajectory.emplace_back(bestCol, bestRow);
    }
 
    writeTrajectory(cfg.outputDir + "/trajectory.txt", mapTrajectory);
    cout << "done." << endl;
    cout << "Output written to " << cfg.outputDir << endl;

    delete cgPtr;
    return 0;

  } // try

  catch (char msg[]) {
    cerr << msg << endl;
    delete cgPtr;
    return 1;
  } // catch

  catch (const string& msg) {
    cerr << msg << endl;
    delete cgPtr;
    throw;
  } // catch

  catch (const exception& e) {
    cerr << "Unhandled exception: " << e.what() << endl;
    delete cgPtr;
    throw;
  } // catch

  catch(...) {
    cerr << "An unknown exception / error occurred\n";
    delete cgPtr;
    throw;
  } // catch

} // main
