#include <iostream>
#include <vector>
#include <string>
#include <array>
#include "headers/Coords.h"
#include "headers/Star.h"
#include "headers/Definitions.h"
using std::cout;

int main(const int argc, const char* argv[]) {
  std::string filePath;
  unsigned plateID;
  Star::takeInput(argc, argv, filePath, plateID);
  std::vector<Star> stars;
  Cartesian mid;
  Star::readStars(stars, mid, filePath);
  Coords tangentPoint;
  Star::findProjectionCoordinates(stars, tangentPoint);
  std::array<double,6> coeff;
  Star::solveLinearEquation(stars, coeff);

  std::array<double,6> asteroidXi  = {mid.x, mid.y, 1, 0, 0, 0};
  std::array<double,6> asteroidEta = {0, 0, 0, mid.x, mid.y, 1};
  double xi = 0.0, eta = 0.0;
  for (int i = 0; i < 6; i++) {
    xi  += asteroidXi[i]  * coeff[i];
    eta += asteroidEta[i] * coeff[i];
  }
  double xiRMS, etaRMS;
  double rms = Star::rms(stars, xiRMS, etaRMS);
  printf("%5s %8s %8s %6s %9s %9s %9s %9s %6s %7s\n", "catid","x","y","scid","xi","xifit","eta","etafit","arcsec","σ");
  for (auto& s : stars) {
    s.printStar(rms);
  }


  Coords asteroidCoords;
  Coords::inverseGnomonic(xi, eta, tangentPoint, asteroidCoords);
  cout << asteroidCoords.toString() << '\n';
}

// TO DO
  // regularise coords before passing to matrix?