#include <iostream>
#include <vector>
#include <string>
#include <array>
#include "headers/Coords.h"
#include "headers/Star.h"
#include "headers/Definitions.h"
#include "headers/Plate.h"
using std::cout;

int main(const int argc, 
         const char* argv[]) {
  std::string filePath, asteroidName;
  unsigned plateID;
  Star::takeInput(argc, argv, filePath, plateID, asteroidName);
  std::vector<Star> stars;
  double x, y;    // pixel coords of the asteroid trail's midpoint
  Star::readStars(stars, x, y, filePath);
  Coords tangentPoint = Star::findProjectionCoordinates(stars);
  Star::matchStars(stars, filePath, plateID, asteroidName, tangentPoint);
  std::array<double,6> coeff = Star::solveLinearEquation(stars);
  double rms = Star::rms(stars, false);
  Coords asteroidCoords = Star::xyToCoords(x, y, coeff, tangentPoint);

  //printf("%5s %8s %8s %5s %9s %9s %9s %9s %7s %7s\n", "catid","x","y","scid","xi","xifit","eta","etafit","arcsec","σ");
  //for (auto& s : stars) s.printStar(rms);
  rms = Star::rms(stars, true);
  for (auto& c : asteroidName) c = toupper(c);
  cout << stars.size() << " stars matched for asteroid " << asteroidName << " and plate " << plateID << '\n';
  cout <<  "Fitted RA/DEC coordinates:\n" << asteroidCoords << '\n'; 
  double xMM, yMM;
  Coords::coordsToPlatePosition(asteroidCoords, plateID, xMM, yMM);
  cout << "Fitted x/y plate position:\n" << xMM << ' ' << yMM << '\n';
  
              // USE THIS BIT FOR FINDING RA/DEC OF MEASURED COORDINATES
  //cout << "Measured RA/DEC:\n";
  //cout << Coords::mmToCoords(302.535, 194.441, Coords(63.25000, -30.00000)) << '\n';
}

// TO DO
  // Take a picture of 1996FG3, starting with 9163 snce the battery died
  // Retake picture of 2000ev70? sc_reduce doesnt like it

  // intensity centroid, automatic asteroid coords? Sextractor?
  // put precovery coords/time into neo checker to see if its actually something else