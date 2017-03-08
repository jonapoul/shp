#include <iostream>
#include <vector>
#include <string>
#include <array>
#include "headers/Coords.h"
#include "headers/Star.h"
#include "headers/Definitions.h"
using std::cout;

int main(const int argc, 
         const char* argv[]) {
  std::string filePath, asteroidName;
  unsigned plateID;
  Star::takeInput(argc, argv, filePath, plateID, asteroidName);
  std::vector<Star> stars;
  double x, y;    // pixel coords of the asteroid trail's midpoint
  Star::readStars(stars, x, y, filePath);
  Coords tangentPoint;
  Star::findProjectionCoordinates(stars, tangentPoint);
  Star::matchStars(stars, filePath, plateID, asteroidName, tangentPoint);
  std::array<double,6> coeff;
  Star::solveLinearEquation(stars, coeff);

  std::array<double,6> asteroidXi  = {x, y, 1, 0, 0, 0};
  std::array<double,6> asteroidEta = {0, 0, 0, x, y, 1};
  double xi = 0.0, eta = 0.0;
  for (int i = 0; i < 6; i++) {
    xi  += asteroidXi[i]  * coeff[i];
    eta += asteroidEta[i] * coeff[i];
  }
  double xiRMS, etaRMS;
  double rms = Star::rms(stars, xiRMS, etaRMS);
  printf("%5s %8s %8s %5s %9s %9s %9s %9s %6s %7s\n", "catid","x","y","scid","xi","xifit","eta","etafit","arcsec","σ");
  for (auto& s : stars) s.printStar(rms);

  Coords asteroidCoords;
  Coords::inverseGnomonic(xi, eta, tangentPoint, asteroidCoords);
  cout << asteroidCoords.toString() << '\n';
}

// TO DO
  // Read refstars file with 3 or 4 stars on it
    // calculate coefficients with those 3 or 4
    // open catalog file and grab the brightest 1000 or so as Catalog structs
    // open supercosmos file and grab 1000 brightest as Supercosmos structs
    // calculate expected xi/eta then B1950 RA/DEC of Catalog objects
    // find distance to all the Supercosmos objects
    // if less than whatever, match up and add to Stars vector
    // carry on as normal
    // temporarily convert to 

  // regularise coords before passing to matrix?