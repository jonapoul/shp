#ifndef STAR_H
#define STAR_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <array>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <experimental/filesystem>
#include "Coords.h"
#include "Definitions.h"
using std::cout;
namespace fs = std::experimental::filesystem;
namespace ublas = boost::numeric::ublas;

class Star {
private:
  unsigned catid_;  // ref number in catalog file
  double x_;        // x pixel
  double y_;        // y pixel
  unsigned scid_;   // ref number in supercosmos txt file
  Coords j2000_;    // J2000 RA/DEC
  Coords b1950_;    // B1950 RA/DEC
  double xi_;       // xi coordinate in radians
  double eta_;      // eta coordinate in radians
  double xifit_;    // least squares fitted xi value
  double etafit_;   // least squares fitted eta value

public:
  Star(const unsigned catid = 0,
       const double x = 0.0,
       const double y = 0.0,
       const unsigned scid = 0,
       const Coords& j2000 = {},
       const Coords& b1950 = {},
       const double xi = 0.0,
       const double eta = 0.0,
       const double xifit = 0.0,
       const double etafit = 0.0) 
      : catid_(catid), x_(x), y_(y), scid_(scid), j2000_(j2000), b1950_(b1950), 
        xi_(xi), eta_(eta), xifit_(xifit), etafit_(etafit) { }
  Star(const Star& s)
      : catid_(s.catid_), x_(s.x_), y_(s.y_), scid_(s.scid_), j2000_(s.j2000_),
        b1950_(s.b1950_), xi_(s.xi_), eta_(s.eta_), xifit_(s.xifit_), 
        etafit_(s.etafit_) { }

  inline unsigned catID() const { return catid_; }
  inline double x() const { return x_; }
  inline double y() const { return y_; }
  inline unsigned scID() const { return scid_; }
  inline Coords j2000() const { return j2000_; }
  inline Coords b1950() const { return b1950_; }
  inline double xi() const { return xi_; }
  inline double eta() const { return eta_; }
  inline double xiFit() const { return xifit_; }
  inline double etaFit() const { return etafit_; }
  inline void setJ2000(const Coords& j2000) { j2000_ = j2000; }
  inline void setB1950(const Coords& b1950) { b1950_ = b1950; }
  inline void setXi(const double xi) { xi_ = xi; }
  inline void setEta(const double eta) { eta_ = eta; }
  inline void setXiFit(const double xifit) { xifit_ = xifit; }
  inline void setEtaFit(const double etafit) { etafit_ = etafit; }

  void printStar(const double rms) const {
    printf("%5d %8.2f %8.2f ", catid_, x_, y_);
    //printf("%s  %s  ", j2000_.toString(true).c_str(), b1950_.toString(true).c_str());
    printf("%06d %9.6f %9.6f ", scid_, xi_, xifit_);
    printf("%9.6f %9.6f ", eta_, etafit_);
    //printf("%9.6f  %9.6f  ", 100*(xi_-xifit_)/xi_, 100*(eta_-etafit_)/eta_);
    double dx = xi_-xifit_;
    double dy = eta_-etafit_;
    /* angular distance between actual and fitted positions, in arcsecs */
    double distance = sqrt(dx*dx + dy*dy) * RAD_TO_DEG * 3600;
    printf("%6.4f ", distance);
    printf("%6.3f ", distance/rms);

    printf("\n");
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
      std::vector<std::string> asteroids;
      for (auto& itr : fs::directory_iterator(filePath)) {
        if (is_directory(itr.path()))
          asteroids.push_back(itr.path().stem().string());
      }
      for (auto& a : asteroids) cout << '\t' << a << '\n';
      std::sort(asteroids.begin(), asteroids.end(), std::less<std::string>());
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
      std::vector<std::string> plates;
      for (auto& itr : fs::directory_iterator(filePath)) {
        if (is_directory(itr.path()))
          plates.push_back(itr.path().stem().string());
      }
      auto plateSorting = [](const std::string& a, const std::string& b) {
        int ai = stoi(a), bi = stoi(b);
        return ai < bi;
      };
      std::sort(plates.begin(), plates.end(), plateSorting);
      for (auto& p : plates) cout << '\t' << p << '\n';
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

  static void readStars(std::vector<Star>& stars, 
                        Cartesian& midpoint, 
                        const std::string& path) {
    std::ifstream file(path);
    if (file.is_open()) {
      while (!file.eof()) {
        std::string buffer;
        getline(file, buffer);
        std::stringstream ss(buffer);
        if (buffer[0] == '#' || buffer.length() == 0)
          continue;
        if (buffer[0] == '@') {
          char temp;
          Cartesian c1, c2;
          ss >> temp >> c1.x >> c1.y >> c2.x >> c2.y;
          midpoint.x = c1.x + abs(c1.x - c2.x)/2.0;
          midpoint.y = c1.y + abs(c1.y - c2.y)/2.0;
          continue;
        }
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

  static void findProjectionCoordinates(std::vector<Star>& stars, 
                                        Coords& tangentPoint) {
    // first finding the star closest to the image centre
    double xMax = 0, xMin = 10000, yMax = 0, yMin = 10000;
    for (auto& s : stars) {
      if      (s.x() > xMax) xMax = s.x();
      else if (s.x() < xMin) xMin = s.x();
      if      (s.y() > yMax) yMax = s.y();
      else if (s.y() < yMin) yMin = s.y();
    }
    double xMid = xMin + (xMax - xMin) / 2.0;
    double yMid = yMin + (yMax - yMin) / 2.0;
    double minDistance = 10000;
    // finding the centremost reference star
    for (auto& s : stars) {
      double dx = s.x() - xMid;
      double dy = s.y() - yMid;
      double distance = sqrt(dx*dx + dy*dy);
      if (distance < minDistance) {
        minDistance = distance;
        // setting that centre star as the tangent point
        tangentPoint = s.b1950();
      }
    }
    for (auto& s : stars) {
      double xi, eta;
      int status;
      Coords::gnomonic(s.b1950(), tangentPoint, xi, eta, status);
      s.setXi(xi);
      s.setEta(eta);
    }
  }

  /*
    Solves the linear equation AX = Y
      A = matrix of x, y coordinates
      X = vector of 6 coefficients
      Y = vector of xi, eta coordinates
    Then adds the fitted xi,eta values to each Star array object
  */
  static void solveLinearEquation(std::vector<Star>& stars,
                                  std::array<double,6>& coefficients) {
    ublas::matrix<double> A(stars.size()*2, 6);
    ublas::vector<double> Y(stars.size()*2);
    for (size_t i = 0; i < 2*stars.size(); i += 2) {
      // matrix A = (x1 y1 1  0  0  0)
      //            (0  0  0  x1 y1 1)
      //            (x2 y2 1  0  0  0)
      //            (0  0  0  x2 y2 1) 
      //            (................) etc, for each set of x, y coordinates 
      A(i,0)   = A(i+1,3) = stars[i/2].x();
      A(i,1)   = A(i+1,4) = stars[i/2].y();
      A(i,2)   = A(i+1,5) = 1;
      A(i,3)   = A(i,4)   = A(i,5)   = 0;
      A(i+1,0) = A(i+1,1) = A(i+1,2) = 0;

      // vector Y = (xi1, eta1, xi2, eta2, ...);
      Y(i)   = stars[i/2].xi();
      Y(i+1) = stars[i/2].eta();
    }
    ublas::matrix<double> At(ublas::trans(A));
    ublas::matrix<double> At_A = prod(At, A);

    // inverting (A^T * A)
    ublas::permutation_matrix<std::size_t> pm(At_A.size1());
    int status = lu_factorize(At_A, pm);
    if (status != 0) {
      cout << "Problem with LU factorisation\n";
      exit(1);
    }
    ublas::matrix<double> inverse = ublas::identity_matrix<double>(At_A.size1());
    try { 
      lu_substitute(At_A, pm, inverse); 
    } catch (...) { 
      cout << "Problem with LU substitution\n"; 
      exit(1);
    }
    ublas::vector<double> At_Y = prod(At, Y);

    // final vector of 6 least squares coefficients
    // X = (A^T * A)^-1 * A^T * Y
    ublas::vector<double> X = prod(inverse, At_Y);
    std::copy(X.begin(), X.end(), coefficients.begin());

    double sumSq = 0.0;
    for (auto& s : stars) {
      double xifit = X(0)*s.x() + X(1)*s.y() + X(2);
      s.setXiFit(xifit);
      double etafit = X(3)*s.x() + X(4)*s.y() + X(5);
      s.setEtaFit(etafit);

      double dxi  = (s.xi()  - xifit)  * RAD_TO_DEG * 3600;  // units of arcsec
      double deta = (s.eta() - etafit) * RAD_TO_DEG * 3600;
      double distance = dxi*dxi + deta*deta;
      sumSq += distance;
    }
  }

  /*
    Outputs rms of xi and eta differences, both in units of arcsecs
  */
  static double rms(const std::vector<Star>& stars, double& xiRMS, double& etaRMS) {
    xiRMS = etaRMS = 0.0;
    for (const auto& s : stars) {
      xiRMS  += pow(s.xi() - s.xiFit(),  2);
      etaRMS += pow(s.eta()- s.etaFit(), 2);
    }
    xiRMS  = sqrt(xiRMS  / stars.size()) * RAD_TO_DEG * 3600;
    etaRMS = sqrt(etaRMS / stars.size()) * RAD_TO_DEG * 3600;
    cout << "dξ = " << xiRMS << ' ' << "dη = " << etaRMS << '\n';
    double rms = sqrt(xiRMS*xiRMS + etaRMS*etaRMS);
    cout << "rms = " << rms << '\n';
    return rms;
  }
};

#endif