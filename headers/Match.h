#ifndef MATCH_H
#define MATCH_H

#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "Plate.h"
#include "Coords.h"
using std::cout;
class Plate;
class Coords;

class Match {
private:
  unsigned count_;
  Plate p_;
  Coords c_;
  double mag_;
  std::pair<double,double> start_;
  std::pair<double,double> mid_;
  std::pair<double,double> end_;
  std::string uncer_;

public:
  Match(const unsigned count = 0,
        const Plate& p = {},
        const Coords& c = {},
        const double mag = 0.0,
        const std::pair<double,double> start = {0.0, 0.0},
        const std::pair<double,double> mid = {0.0, 0.0},
        const std::pair<double,double> end = {0.0, 0.0},
        const std::string& uncer = "")
      : count_(count), p_(p), c_(c), mag_(mag), start_(start), 
        mid_(mid), end_(end), uncer_(uncer) { }
  Match(const Match& m)
      : count_(m.count_), p_(m.p_), c_(m.c_), mag_(m.mag_), start_(m.start_), 
        mid_(m.mid_), end_(m.end_), uncer_(m.uncer_) { }

  inline unsigned count() const { return count_; }
  inline Plate p() const { return p_; }
  inline Coords c() const { return c_; }
  inline double mag() const { return mag_; }
  inline std::pair<double,double> start() const { return start_; }
  inline std::pair<double,double> mid() const { return mid_; }
  inline std::pair<double,double> end() const { return end_; }
  inline std::string uncer() const { return uncer_; }


  /*
    Goes through all matched plates and prints the relevant info about them all
    This is a LITTLE BIT OF A MESS but it works
  */
  static void printMatches(const std::vector<Match>& m) {
    const int SIZE = 48;  // defines width of each printed data box
    size_t length;

    for (int i = 0; i < int(m.size()); i += 2) {
      bool canPrint = !(i == m.size() - 1 && m.size() % 2 == 1);
      std::stringstream ss;
      char buffer[50];

      // Match number and Plate ID
      sprintf(buffer, "%03d\tplate ID       = %d", m[i].count(), m[i].p().id());
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length+3, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "%03d\tplate ID       = %d", m[i+1].count(), m[i+1].p().id());
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // UT Date
      ss << "\tUT date/time   = " << Plate::gregorianToString(m[i].p().gregorian());
      ss << ", " << Plate::julianToUT(m[i].p().julian());
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        ss << "\tUT date        = " << Plate::gregorianToString(m[i+1].p().gregorian());
        ss << ", " << Plate::julianToUT(m[i+1].p().julian());
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Julian Date
      sprintf(buffer, "\tJulian date    = %.5f", m[i].p().julian());
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "\tJulian date    = %.5f", m[i+1].p().julian());
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Plate centre's RA/DEC coordinates in degrees
      sprintf(buffer, "\tPlate Centre   = (%.5f, %.5f)°", m[i].p().coords().ra(DEG), m[i].p().coords().dec(DEG));
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length+1, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "\tPlate Centre   = (%.5f, %.5f)°", m[i+1].p().coords().ra(DEG), m[i+1].p().coords().dec(DEG));
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Object's RA/DEC coordinates in degrees
      sprintf(buffer, "\tObject Coords  = (%.5f, %.5f)°", m[i].c().ra(DEG), m[i].c().dec(DEG));
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length+1, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "\tObject Coords  = (%.5f, %.5f)°", m[i+1].c().ra(DEG), m[i+1].c().dec(DEG));
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Object's positional coordinates on the plate, in millimetres from the bottom left corner
      sprintf(buffer, "\tPlate Position = (%.3f, %.3f) mm", m[i].mid().first, m[i].mid().second);
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "\tPlate Position = (%.3f, %.3f) mm", m[i+1].mid().first, m[i+1].mid().second);
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // x/y positional uncertainties in mm
      sprintf(buffer, "\tUncertainty    = %s mm", m[i].uncer().c_str());
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length+2, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "\tUncertainty    = %s mm", m[i+1].uncer().c_str());
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Length of the object's exposure trail on the plate in millimetres, plus the directional vector
      char buffer2[50];
      double dx = m[i].end().first  - m[i].start().first;
      double dy = m[i].end().second - m[i].start().second;
      double drift = sqrt(dx*dx + dy*dy);
      sprintf(buffer, "\tTrail length   = %.2f mm", drift);
      ss << buffer;
      sprintf(buffer, "\tTrail vector   = (%.3f, %.3f)", dx, dy);
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        dx = m[i+1].end().first  - m[i+1].start().first;
        dy = m[i+1].end().second - m[i+1].start().second;
        drift = sqrt(dx*dx + dy*dy);
        sprintf(buffer2, "\tTrail length   = %.2f mm", drift);
        ss << buffer2;
        sprintf(buffer2, "\tTrail vector   = (%.3f, %.3f)", dx, dy);
      }
      cout << ss.str() << '\n';
      ss.str("");
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        ss << buffer2;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Object's apparent magnitude on the plate
      if ( abs(m[i].mag() - UNKNOWN_MAGNITUDE) < 1e-6 ) {
        sprintf(buffer, "\tMagnitude      = UNKNOWN");
      } else {
        sprintf(buffer, "\tMagnitude      = %.2f", m[i].mag());
      }
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        if ( abs(m[i+1].mag() - UNKNOWN_MAGNITUDE) < 1e-6 ) {
          sprintf(buffer, "\tMagnitude      = UNKNOWN");
        } else {
          sprintf(buffer, "\tMagnitude      = %.2f", m[i+1].mag());
        }
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Plate quality grade
      ss << "\tPlate Grade    = " << m[i].p().grade();
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        ss << "\tPlate Grade    = " << m[i+1].p().grade();
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Plate exposure time in minutes
      sprintf(buffer, "\tExposure       = %.1f mins", m[i].p().exposure() / 60.0);
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "\tExposure       = %.1f mins", m[i+1].p().exposure() / 60.0);
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Sort-of signal to noise ratio
      double snr1 = Ephemeris::counts(m[i].p().exposure(), m[i].mag()) / m[i].p().countLimit();
      double snr2 = Ephemeris::counts(m[i+1].p().exposure(), m[i+1].mag()) / m[i+1].p().countLimit();
      if ( abs(m[i].mag() - UNKNOWN_MAGNITUDE) < 1e-6 ) {
        sprintf(buffer, "\tSNR            = UNKNOWN");
      } else if (snr1 > 1e5){
        sprintf(buffer, "\tSNR            = %.2e", snr1);
      } else {
        sprintf(buffer, "\tSNR            = %.4f", snr1);
      }
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length, ' ') << "| ";
      if (canPrint) {
        if ( abs(m[i+1].mag() - UNKNOWN_MAGNITUDE) < 1e-6 ) {
          sprintf(buffer, "\tSNR            = UNKNOWN");
        } else if (snr2 > 1e5) {
          sprintf(buffer, "\tSNR            = %.2e", snr2);
        } else {
          sprintf(buffer, "\tSNR            = %.4f", snr2);
        }
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      // Plate filter peak wavelength
      sprintf(buffer, "\tFilter λ peak  = %dÅ", m[i].p().wavelength());
      ss << buffer;
      length = SIZE-ss.str().length();
      ss << std::string(length+2, ' ') << "| ";
      if (canPrint) {
        sprintf(buffer, "\tFilter λ peak  = %dÅ", m[i+1].p().wavelength());
        ss << buffer;
      }
      cout << ss.str() << '\n';
      ss.str("");

      if (canPrint)
        cout << std::string(SIZE*2.4, '-') << '\n';
      else
        cout << std::string(SIZE*1.2, '-') << '\n';
    }
    for (int i = 0; PRINT_FOR_SPREADSHEET && i < m.size(); i++) {
      double snr = Ephemeris::counts(m[i].p().exposure(), m[i].mag()) / m[i].p().countLimit();
      printf("%d %.3f %.3f %.3f %.3f\n", m[i].p().id(), m[i].p().julian(), snr, m[i].mid().first, m[i].mid().second);
    }
  }
};

#endif