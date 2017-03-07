#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include "headers/Coords.h"
#include "headers/Star.h"
namespace chr = std::chrono;
using std::cout;

int main(const int argc, const char* argv[]) {
  std::string filePath;
  unsigned plateID;
  Star::takeInput(argc, argv, filePath, plateID);

  chr::time_point<chr::system_clock> start = chr::system_clock::now();

  std::vector<Star> stars;
  Star::readStars(stars, filePath);
  Star::findProjectionCoordinates(stars);
  double stddev;
  Star::solveLinearEquation(stars, stddev);
    
  //cout << "stddev = " << stddev << " arcseconds\n";
  unsigned counts[] = {0,0,0,0,0,0,0,0};
  printf("%5s %8s %8s %6s %9s %9s %9s %9s  %6s\n",
         "catid","x","y","scid","xi","xifit","eta","etafit","distance");
  for (auto& s : stars) {
    double distance;
    s.printStar(distance);//, stddev);
    double delta = distance / stddev;
    counts[int(delta)]++;
  }
  //cout << "1σ = " << 100*double(counts[0])/stars.size() << "\%\n";
  //cout << "2σ = " << 100*double(counts[0]+counts[1])/stars.size() << "\%\n";


  chr::duration<double> elapsed_seconds = chr::system_clock::now() - start;
  printf("Elapsed time: %.4fs\n", elapsed_seconds.count());
}

// TO DO
  // redo stddev to take sum of each individual xi, eta coord

  // rerun program, remeasure the ones that are furthest off
    // ask Nigel if stddev calculation is ok, I've no idea what I'm doing

  // go and measure refstars for each plate, using catalog2 for easier reading (kind of)