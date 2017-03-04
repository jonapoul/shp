#ifndef STAR_H
#define STAR_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <experimental/filesystem>
#include "Coords.h"
using std::cout;
namespace fs = std::experimental::filesystem;

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
      printf("%5d %8.2f %8.2f  %06d  %s  %s\n", catid_, x_, y_, scid_, j2000_.toString(true).c_str(), b1950_.toString(true).c_str());
  }

  static void takeInput(const int argc,
                        const char* argv[],
                        std::string& path,
                        unsigned& plateNumber) {
    fs::path filePath = "./images";
    if (argc == 3) { // if all options where given as arguments
      filePath /= fs::path(argv[1]) / fs::path(argv[2]);
      if (!fs::exists(filePath)) {
        cout << filePath.string() << " doesn't exist\n";
        exit(1);
      }
      for (auto& itr : fs::directory_iterator(filePath)) {
        if (itr.path().filename().string() == "refstars.txt") {
          path = itr.path().string();
          return;
        }
      }
    } else if (argc == 1) { // if no arguments were given ask for asteroid name
      cout << "  Asteroids:\n";
      for (auto& itr : fs::directory_iterator(filePath)) {
        if (is_directory(itr.path())) {
          cout << '\t' << itr.path().stem().string() << '\n';
        }
      }
      cout << "Option: ";
      while (true) {
        std::string buf;
        getline(std::cin, buf);
        std::transform(buf.begin(), buf.end(), buf.begin(), ::tolower);
        buf.erase( remove_if(buf.begin(), buf.end(), isspace), buf.end() );
        if (fs::exists(filePath / buf)) {
          filePath /= buf; 
          path = filePath.string();
          break;
        } else {
          cout << "Try again: ";
        }
      }
    } else if (argc == 2) {
      filePath /= fs::path(argv[1]);
    }
    if (argc <= 2) {  // if the plate number wasnt given
      cout << "  Plates:\n";
      for (auto& itr : fs::directory_iterator(filePath)) {
        if (is_directory(itr.path())) {
          cout << '\t' << itr.path().stem().string() << '\n';
        }
      }
      cout << "Option: ";
      while (true) {
        std::string buf;
        getline(std::cin, buf);
        std::transform(buf.begin(), buf.end(), buf.begin(), ::tolower);
        buf.erase( remove_if(buf.begin(), buf.end(), isspace), buf.end() );
        if (fs::exists(filePath / buf / "refstars.txt")) {
          filePath /= buf / fs::path("refstars.txt"); 
          path = filePath.string() ;
          plateNumber = stoi(buf);
          break;
        } else {
          cout << "Try again: ";
        }
      }
    }
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
        Coords j2000 = {ra2k, dec2k};
        Coords b1950 = Coords::convertEpoch(2000, 1950, j2000);
        stars.push_back( Star(catalogID, xPix, yPix, supercosmosID, j2000, b1950) );
      }
      file.close();
    } else {
      cout << "Couldn't open " << path << "\n";
      exit(1);
    }
  }
};

#endif