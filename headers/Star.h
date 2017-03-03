#ifndef STAR_H
#define STAR_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "Coords.h"
using std::cout;

class Star {
private:
  unsigned catid_;  // ref number in catalog file
  double x_;        // x pixel
  double y_;        // y pixel
  unsigned scid_;   // ref number in supercosmos txt file
  Coords j2000_;     
  Coords b1950_;
  double xi_;       // xi coordinate in radians
  double eta_;      // eta coordinate in radians

public:
  Star(const unsigned catid = 0,
       const double x = 0.0,
       const double y = 0.0,
       const unsigned scid = 0,
       const Coords& j2000 = {},
       const Coords& b1950 = {},
       const double xi = 0.0,
       const double eta = 0.0) 
      : catid_(catid), x_(x), y_(y), scid_(scid), j2000_(j2000), b1950_(b1950), 
        xi_(xi), eta_(eta) { }
  Star(const Star& s)
      : catid_(s.catid_), x_(s.x_), y_(s.y_), scid_(s.scid_), j2000_(s.j2000_),
        b1950_(s.b1950_), xi_(s.xi_), eta_(s.eta_) { }

  inline unsigned catID() const { return catid_; }
  inline double x() const { return x_; }
  inline double y() const { return y_; }
  inline unsigned scID() const { return scid_; }
  inline Coords j2000() const { return j2000_; }
  inline Coords b1950() const { return b1950_; }
  inline double xi() const { return xi_; }
  inline double eta() const { return eta_; }
  inline void setJ2000(const Coords& j2000) { j2000_ = j2000; }
  inline void setB1950(const Coords& b1950) { b1950_ = b1950; }
  inline void setXi(const double xi) { xi_ = xi; }
  inline void setEta(const double eta) { eta_ = eta; }

  void printStar() const {
      printf("%5d %8.2f %8.2f %5d %10.7f %10.7f\n", catid_, x_, y_, scid_, j2000_.ra(DEG), j2000_.dec(DEG));
  }

  static void takeInput(const int argc,
                        const char* argv[],
                        std::string& path,
                        unsigned& plateNumber) {
    std::string filename;
    if (argc < 3) {
      cout << "Asteroid name: ";
      getline(std::cin, filename);
      cout << "Plate number: ";
      std::string numStr;
      getline(std::cin, numStr);
      /* ADD A FILESYSTEM CHECKER HERE SOMETIME */
      plateNumber = stoi(numStr);
    } else {
      filename = std::string(argv[1]);
      plateNumber = stoi(std::string(argv[2]));
    }
    path = "./images/" + filename + "/" + std::to_string(plateNumber) + "/refstars.txt";
  }

  static void readStars(std::vector<Star>& stars, const std::string& path) {
    std::ifstream file(path);
    if (file.is_open()) {
      while (!file.eof()) {
        std::string buffer;
        getline(file, buffer);
        if (buffer[0] == '#' || buffer.length() == 0)
          continue;
        std::stringstream ss(buffer);
        unsigned catalogID, supercosmosID;
        double xPix, yPix, ra2k, dec2k;
        ss >> catalogID >> xPix >> yPix >> supercosmosID >> ra2k >> dec2k;
        stars.push_back( Star(catalogID, xPix, yPix, supercosmosID, Coords(ra2k,dec2k)) );
      }
      file.close();
    } else {
      cout << "Couldn't open " << path << "\n";
      exit(1);
    }
  }

  /*
    Converts between two RA/DEC systems based on time epochs
      t0 = year basis of coordinates c0
      t  = year of output coordinates 

      From Duffett-Smith
  */
  static Coords convertEpoch(const double t0,
                             const double t,
                             const Coords& c0) {
    double ra0 = c0.ra(RAD);
    double dec0 = c0.dec(RAD);

    double N = t - t0;
    double S1 = (3.07327 + 1.33617*sin(ra0)*tan(dec0)) * N;
    S1 /= 3600.0;

    double ra = S1 + (c0.ra(DEG) / 15.0);
    ra *= 15.0;
    while (ra < 0)   ra += 360.0;
    while (ra > 360) ra -= 360.0;

    double S2 = 20.0426 * cos(ra0) * N;
    S2 /= 3600.0;
    double dec = c0.dec(DEG) + S2;

    if (dec < -90) {
      dec = -90 + abs(dec+90);
      ra += 180.0;
      if (ra > 360) ra -= 360;
    } else if (dec > 90) {
      dec = 90 - abs(dec-90);
      ra += 180.0;
      if (ra > 360) ra -= 360;
    }
    return Coords(ra, dec);
  }

  static Coords convertEpoch2(const double t0,
                             const double t,
                             const Coords& c0) {
    return Coords();
  }
};

#endif