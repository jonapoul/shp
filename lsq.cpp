#include <iostream>
#include <vector>
#include <string>
#include "headers/Coords.h"
#include "headers/Star.h"
using std::cout;

int main(const int argc, const char* argv[]) {

  std::string filePath;
  unsigned plateID;
  Star::takeInput(argc, argv, filePath, plateID);
  std::vector<Star> stars;
  Star::readStars(stars, filePath);
  for (auto& s : stars) s.printStar();
}

// TO DO
  // rob that least squares algorithm from http://www.vilipetek.com/2013/10/07/polynomial-fitting-in-c-using-boost/

  // make a file with a line for each star
    // put together one with 15+ stars for each plate GOD DAMNIT