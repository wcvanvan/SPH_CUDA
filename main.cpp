#include <iostream>
#include <string>
#include "compute.h"
#include "visualizer.h"

int main(int argc, char *argv[]) {
  bool runCompute = true;
  bool runVisualize = true;

  if (argc == 2) {
    std::string arg = argv[1];

    if (arg == "-c") {
      runVisualize = false;
    } else if (arg == "-v") {
      runCompute = false;
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      std::cerr << "Usage: " << argv[0] << " [-c | -v]\n";
      return 1;
    }
  } else if (argc > 2) {
    std::cerr << "Too many arguments.\n";
    std::cerr << "Usage: " << argv[0] << " [-c | -v]\n";
    return 1;
  }

  if (runCompute) compute();
  if (runVisualize) visualize();

  return 0;
}
