#ifndef STAR_H
#define STAR_H

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
public:
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

  struct SuperCosmosStar {
    unsigned id;
    Coords b1950;
    Coords j2000;
    double magnitude;
    double xi;
    double eta;
  };
  struct CatalogStar {
    unsigned id;
    double x;
    double y;
    double magnitude;
    double xi;
    double eta;
  };

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

  inline double x() const { return x_; }
  inline double y() const { return y_; }
  inline Coords j2000() const { return j2000_; }
  inline Coords b1950() const { return b1950_; }
  inline double xi() const { return xi_; }
  inline double eta() const { return eta_; }
  inline double xifit() const { return xifit_; }
  inline double etafit() const { return etafit_; }

  /*
    Prints all info about the star
  */
  void printStar(const double rms, const bool oneSpaceOnly = false) const {
    double dx = xi_-xifit_;
    double dy = eta_-etafit_;
    /* angular distance between actual and fitted positions, in arcsecs */
    double distance = sqrt(dx*dx + dy*dy) * RAD_TO_DEG * 3600;

    if (oneSpaceOnly) {
      cout << catid_ << ' ' << x_ << ' ' << y_ << ' ' << scid_ << ' ' << xi_;
      cout << ' ' << xifit_ << ' ' << eta_ << ' ' << etafit_ << ' ' << distance;
      cout << ' ' << distance/rms << '\n';
    } else {
      printf("%5d %8.2f %8.2f ", catid_, x_, y_);
      //printf("%s  %s  ", j2000_.toString(true).c_str(), b1950_.toString(true).c_str());
      //printf("%s ", b1950_.toString(true).c_str());
      printf("%5d %9.6f %9.6f ", scid_, xi_, xifit_);
      printf("%9.6f %9.6f ", eta_, etafit_);
      //printf("%9.6f  %9.6f  ", 100*(xi_-xifit_)/xi_, 100*(eta_-etafit_)/eta_);
      printf("%7.4f ", distance);
      printf("%6.3f ", distance/rms);

      printf("\n");
    }
  }

  /*
    Uses the arguments passed intoto the program to find the catalog/supercosmos
    files. If not enough arguments were passed, it'll ask the user for input
  */
  static void takeInput(const int argc,
                        const char* argv[],
                        std::string& path,
                        unsigned& plateNumber,
                        std::string& asteroidName) {
    fs::path filePath = "./images";
    if (argc == 3) { // if all options where given as arguments
      filePath /= fs::path(argv[1]) / fs::path(argv[2]);
      asteroidName = std::string(argv[1]);
      plateNumber = stoi(std::string(argv[2]));
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
      cout << "Couldn't find refstars.txt\n";
      exit(1);
    } else if (argc == 1) { // if no arguments were given, ask for asteroid name
      cout << "  Asteroids:\n\t";
      std::vector<std::string> asteroids;
      for (auto& itr : fs::directory_iterator(filePath)) {
        if (is_directory(itr.path()))
          asteroids.push_back(itr.path().stem().string());
      }
      std::sort(asteroids.begin(), asteroids.end(), std::less<std::string>());
      unsigned longest = 0;
      for (auto a : asteroids) if (a.length() > longest) longest = a.length();
      for (int i = 1; i <= asteroids.size(); i++) {
        unsigned spaces = longest-asteroids[i-1].length();
        cout << asteroids[i-1] << std::string(spaces+4, ' ');
        if (i % 5 == 0 && i != asteroids.size()) cout << "\n\t";
      } 
      cout << "\n  Option: "; 
      while (true) {
        std::string buf;
        getline(std::cin, buf);
        std::transform(buf.begin(), buf.end(), buf.begin(), ::tolower);
        buf.erase( remove_if(buf.begin(), buf.end(), isspace), buf.end() );
        if (fs::exists(filePath / buf)) {
          filePath /= buf; 
          asteroidName = buf;
          path = filePath.string();
          break;
        } else {
          cout << "Try again: ";
        }
      }
    } else if (argc == 2) {
      filePath /= fs::path(argv[1]);
      asteroidName = std::string(argv[1]);
    }
    if (argc <= 2) {  // if the plate number wasnt given
      std::vector<std::string> plates;
      for (auto& itr : fs::directory_iterator(filePath)) {
        if (is_directory(itr.path()))
          plates.push_back(itr.path().stem().string());
      }
      if (plates.size() == 1) {
        plateNumber = stoi(plates[0]);
        path = filePath / plates[0] / fs::path("refstars.txt");
        return;
      }
      auto plateSorting = [](const std::string& a, const std::string& b) {
        int ai = stoi(a), bi = stoi(b);
        return ai < bi;
      };
      std::sort(plates.begin(), plates.end(), plateSorting);
      cout << "  Plates:\n\t";
      unsigned longest = 0;
      for (auto p : plates) if (p.length() > longest) longest = p.length();
      for (int i = 1; i <= plates.size(); i++) {
        unsigned spaces = longest-plates[i-1].length();
        cout << plates[i-1] << std::string(spaces+4, ' ');
        if (i % 5 == 0 && i != plates.size()) cout << "\n\t";
      } 
      cout << "\n  Option: ";
      while (true) {
        std::string buf;
        getline(std::cin, buf);
        std::transform(buf.begin(), buf.end(), buf.begin(), ::tolower);
        buf.erase( remove_if(buf.begin(), buf.end(), isspace), buf.end() );
        if (fs::exists(filePath / buf / "refstars.txt")) {
          filePath /= buf / fs::path("refstars.txt"); 
          path = filePath.string() ;
          plateNumber = stoi(buf);
          return;
        } else {
          cout << "Try again: ";
        }
      }
    }
  }

  /*
    Reads the initial matched stars in as an array of Star objects
  */
  static void readStars(std::vector<Star>& stars, 
                        double& x, 
                        double& y, 
                        const std::string& path) {
    std::ifstream file(path);
    if (file.is_open()) {
      for (std::string buffer; getline(file, buffer); ) {
        std::stringstream ss(buffer);
        if (buffer[0] == '#' || buffer.length() == 0) {
          continue;
        } else if (buffer[0] == '@') {
          char temp;
          double x1, x2, y1, y2;
          ss >> temp >> x1 >> y1 >> x2 >> y2;
          x = (x1 + x2) / 2.0;
          y = (y1 + y2) / 2.0;
        } else {
          unsigned catalogID, supercosmosID;
          double xPix, yPix, ra2k, dec2k;
          ss >> catalogID >> xPix >> yPix >> supercosmosID >> ra2k >> dec2k;
          Coords j2000 = {ra2k, dec2k};
          Coords b1950 = Coords::convertEpoch(2000, 1950, j2000);
          stars.push_back( {catalogID, xPix, yPix, supercosmosID, j2000, b1950} );
        }
      }
      file.close();
    } else {
      cout << "Couldn't open " << path << "\n";
      exit(1);
    }
  }

  /*
    Finds the centremost star in the array of Star objects, and calculates the
    xi/eta coords of each other star with respect to that centre star
  */
  static Coords findProjectionCoordinates(std::vector<Star>& stars) {
    // first finding the star closest to the image centre
    double xMax = 0, xMin = 10000, yMax = 0, yMin = 10000;
    for (auto& s : stars) {
      if      (s.x_ > xMax) xMax = s.x_;
      else if (s.x_ < xMin) xMin = s.x_;
      if      (s.y_ > yMax) yMax = s.y_;
      else if (s.y_ < yMin) yMin = s.y_;
    }
    double xMid = xMin + (xMax - xMin) / 2.0;
    double yMid = yMin + (yMax - yMin) / 2.0;
    double minDistance = 10000;
    // finding the centremost reference star
    Coords tangentPoint;
    for (auto& s : stars) {
      double dx = s.x_ - xMid;
      double dy = s.y_ - yMid;
      double distance = sqrt(dx*dx + dy*dy);
      if (distance < minDistance) {
        minDistance = distance;
        // setting that centre star as the tangent point
        tangentPoint = s.b1950_;
      }
    }
    for (auto& s : stars) {
      double xi, eta;
      Coords::gnomonic(s.b1950_, tangentPoint, xi, eta);
      s.xi_ = xi;
      s.eta_ = eta;
    }
    return tangentPoint;
  }

  /*
    Solves the linear equation AX = Y
      A = matrix of x, y coordinates
      X = vector of 6 coefficients
      Y = vector of xi, eta coordinates
    Then adds the fitted xi,eta values to each Star array object
  */
  static std::array<double,6> solveLinearEquation(std::vector<Star>& stars) {
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
      Y(i)   = stars[i/2].xi_;
      Y(i+1) = stars[i/2].eta_;
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
    } catch (std::exception& e) { 
      cout << "Problem with LU substitution: " << e.what() << "\n"; 
      exit(1);
    }
    ublas::vector<double> At_Y = prod(At, Y);

    // final vector of 6 least squares coefficients
    // X = (A^T * A)^-1 * A^T * Y
    ublas::vector<double> X = prod(inverse, At_Y);
    std::array<double,6> coeff;
    std::copy(X.begin(), X.end(), coeff.begin());
    for (auto& s : stars) {
      double xi  = coeff[0]*s.x_ + coeff[1]*s.y_ + coeff[2];
      double eta = coeff[3]*s.x_ + coeff[4]*s.y_ + coeff[5];
      s.xifit_ = xi;
      s.etafit_ = eta;
    }
    return coeff;
  }

  /*
    Outputs rms of xi and eta differences, both in units of arcsecs
  */
  static double rms(const std::vector<Star>& stars, 
                    const bool verbose = false) {
    double xiRMS = 0.0;
    double etaRMS = 0.0;
    for (const auto& s : stars) {
      xiRMS  += pow(s.xi_ - s.xifit_,  2);
      etaRMS += pow(s.eta_- s.etafit_, 2);
    }
    xiRMS  = sqrt(xiRMS  / stars.size()) * RAD_TO_DEG * 3600;
    etaRMS = sqrt(etaRMS / stars.size()) * RAD_TO_DEG * 3600;
    double rms = sqrt(xiRMS*xiRMS + etaRMS*etaRMS);
    if (verbose)
      cout << "σξ = " << xiRMS << " ση = " << etaRMS << " σ = " << rms << " arcsecs\n";
    return rms;
  }

  /*
    Takes an initial set of coefficients calculated from the refstars.txt file,
    then uses them to match up stars from catalog.txt and supercosmos.txt,
    adding them to the array if they aren't too far from the rms
  */
  static void matchStars(std::vector<Star>& stars,
                         const std::string& filePath,
                         const unsigned& plateNumber,
                         const std::string& asteroidName, 
                         Coords& tangentPoint) {
    if (stars.size() == 0) {
      cout << "Not enough stars in refstars.txt\n";
      exit(1);
    }
    // getting initial coefficients from the file stars
    std::array<double,6> coeff = solveLinearEquation(stars);
    
    // reading the catalog file
    fs::path parent = fs::path(filePath).parent_path();
    std::string catalogfilename = asteroidName + "_" + std::to_string(plateNumber) + "_10_catalog.dat";
    fs::path catalogpath = parent / catalogfilename;
    std::ifstream catalogfile(catalogpath.string());
    std::vector<CatalogStar> catalog;
    if (catalogfile.is_open()) {
      for (std::string buffer; getline(catalogfile, buffer); ) {
        std::stringstream ss(buffer);
        if (buffer.empty() || buffer[0] == '#') 
          continue;
        unsigned id;
        double x, y, magnitude;
        ss >> id >> x >> y >> magnitude;
        double xi  = coeff[0]*x + coeff[1]*y + coeff[2];
        double eta = coeff[3]*x + coeff[4]*y + coeff[5];
        CatalogStar cs = {id, x, y, magnitude, xi, eta};
        catalog.push_back(cs);
      }
      catalogfile.close();
      auto catalogsort = [](CatalogStar& a,CatalogStar& b) { return a.magnitude < b.magnitude; };
      std::sort(catalog.begin(), catalog.end(), catalogsort);
      // remove all but the brightest objects
      for (int i = catalog.size()-1; i >= 2000; i--)
        catalog.erase(catalog.begin()+i);
    } else {
      cout << "Couldn't open " << catalogfilename << "\n";
      exit(1);
    }

    // reading the supercosmos file
    std::string scfilename = std::to_string(plateNumber) + "supercosmos.txt";
    fs::path scpath = parent / scfilename;
    std::ifstream supercosmosfile(scpath.string());
    std::vector<SuperCosmosStar> supercosmos;
    if (supercosmosfile.is_open()) {
      for (std::string buffer; getline(supercosmosfile, buffer); ) {
        std::stringstream ss(buffer);
        if (buffer.length() == 0 || buffer[0] == '#') 
          continue;
        unsigned id;
        double ra2k, dec2k, magnitude, temp;
        double mags[4];
        ss >> id >> ra2k >> dec2k;
        ss >> temp >> temp >> temp >> temp;
        ss >> mags[0] >> mags[1] >> mags[2] >> mags[3];
        magnitude = (mags[0] + mags[1] + mags[2] + mags[3]) / 4.0;
        Coords j2000 = {ra2k, dec2k};
        Coords b1950 = Coords::convertEpoch(2000, 1950, j2000);
        double xi, eta;
        int status = Coords::gnomonic(b1950, tangentPoint, xi, eta);
        if (status == 0) {
          SuperCosmosStar scs = {id, b1950, j2000, magnitude, xi, eta};
          supercosmos.push_back(scs);
        }
      }
      supercosmosfile.close();
      auto supercosmossort = [](SuperCosmosStar& a,SuperCosmosStar& b) { return a.magnitude < b.magnitude; };
      std::sort(supercosmos.begin(), supercosmos.end(), supercosmossort);
      // remove all but the brightest objects
      size_t numStars = supercosmos.size();
      for (int i = numStars; i >= numStars/2; i--)
        supercosmos.erase(supercosmos.begin()+i);
    } else {
      cout << "Couldn't open " << scfilename << "\n";
      exit(1);
    }


    // match them up if the two objects are within 1 arcsec of each other
    double threshold = 1.0;  // in arcseconds
    for (const auto& c : catalog) {
      Star s;
      for (const auto& sc : supercosmos) {
        double dxi = (sc.xi - c.xi) * RAD_TO_DEG * 3600;
        double deta = (sc.eta - c.eta) * RAD_TO_DEG * 3600;
        double distance = sqrt(dxi*dxi + deta*deta);
        if (distance < threshold)
          s = Star(c.id, c.x, c.y, sc.id, sc.j2000, sc.b1950, (sc.xi+c.xi)/2.0, (sc.eta+c.eta)/2.0);
      }
      bool hasAlreadybeenAdded = false;
      for (const auto& star : stars) {
        if (star.catid_ == s.catid_) {
          hasAlreadybeenAdded = true;
          break;
        }
      }
      if (!hasAlreadybeenAdded && s.catid_ != 0) {
        stars.push_back(s);
      }
    }
    // sort in order of catalog id number for ease of reading
    auto starsort = [](Star& a,Star& b) { return a.catid_ < b.catid_; };
    std::sort(stars.begin(), stars.end(), starsort);

    // refitting and removing outliers 5 times, to improve rms
    for (int i = 0; i < 5; i++) {
      tangentPoint = findProjectionCoordinates(stars);
      coeff = solveLinearEquation(stars);
      double rms = Star::rms(stars, false);
      double furthestDistance = 0;
      for (int i = stars.size()-1; i >= 0; i--) {
        double xi, eta;
        Coords::gnomonic(stars[i].b1950(), tangentPoint, xi, eta);
        stars[i].xi_ = xi;
        stars[i].eta_ = eta;
        double dx = stars[i].xi_ - stars[i].xifit_;
        double dy = stars[i].eta_ - stars[i].etafit_;
        double distance = sqrt(dx*dx + dy*dy) * RAD_TO_DEG * 3600;
        if (distance/rms > 1.7) stars.erase(stars.begin()+i);
      }
    }
  }

  /*
    Converts x/y pixel position to xi/eta image positions wrt the tangent point
  */
  static Coords xyToCoords(const double x, 
                           const double y, 
                           const std::array<double,6>& coeff, 
                           const Coords& tangentPoint) {
    double asteroidXi[]  = {x, y, 1, 0, 0, 0};
    double asteroidEta[] = {0, 0, 0, x, y, 1};
    double xi = 0.0, eta = 0.0;
    for (int i = 0; i < 6; i++) {
      xi  += asteroidXi[i]  * coeff[i];
      eta += asteroidEta[i] * coeff[i];
    }
    Coords asteroidCoords;
    Coords::inverseGnomonic(xi, eta, tangentPoint, asteroidCoords);
    return asteroidCoords;
  }

  /*
    Prints all info about the Star array members, followed by the RMS of the fit    
  */
  static void printAllStars(const std::vector<Star>& stars) {
    double rms = Star::rms(stars, false);
    printf("%5s %8s %8s %5s %9s %9s %9s %9s %7s %7s\n", 
           "catid","x","y","scid","xi","xifit","eta","etafit","arcsec","σ");
    for (auto& s : stars) 
      s.printStar(rms);
    rms = Star::rms(stars, true);
  }

  /*
  **  Decompose an [x,y] linear fit into its constituent parameters:
  **  zero points, scales, nonperpendicularity and orientation.
  **  Given:
  **     coeffs    double[6]     transformation coefficients (see note)
  **  Returned:
  **     *xz       double        x zero point
  **     *yz       double        y zero point
  **     *xs       double        x scale
  **     *ys       double        y scale
  **     *perp     double        nonperpendicularity (radians)
  **     *orient   double        orientation (radians)
  **  The model relates two sets of [x,y] coordinates as follows.
  **  Naming the elements of coeffs:
  **     coeffs[0] = a
  **     coeffs[1] = b
  **     coeffs[2] = c
  **     coeffs[3] = d
  **     coeffs[4] = e
  **     coeffs[5] = f
  **  The model transforms coordinates [x1,y1] into coordinates
  **  [x2,y2] as follows:
  **     x2 = a + b*x1 + c*y1
  **     y2 = d + e*x1 + f*y1
  **  The transformation can be decomposed into four steps:
  **     1)  Zero points:
  **             x' = xz + x1
  **             y' = yz + y1
  **     2)  Scales:
  **             x'' = xs*x'
  **             y'' = ys*y'
  **     3)  Nonperpendicularity:
  **             x''' = cos(perp/2)*x'' + sin(perp/2)*y''
  **             y''' = sin(perp/2)*x'' + cos(perp/2)*y''
  **     4)  Orientation:
  **             x2 = cos(orient)*x''' + sin(orient)*y'''
  **             y2 =-sin(orient)*y''' + cos(orient)*y'''
  **
  **  Copyright P.T.Wallace.  All rights reserved.
  **  Last revision:   22 September 1995
  */
  static void coeffsDecomposition(const std::array<double,6>& coeffs, 
                                  double& x0, 
                                  double& y0, 
                                  double& xscale, 
                                  double& yscale, 
                                  double& nonperpendicularity, 
                                  double& orientation,
                                  const bool verbose = false) {
  /* Copy the six coefficients */
    double a = coeffs[0];
    double b = coeffs[1];
    double c = coeffs[2];
    double d = coeffs[3];
    double e = coeffs[4];
    double f = coeffs[5];

  /* Scales */
    double rb2e2 = sqrt(b*b + e*e);
    double rc2f2 = sqrt(c*c + f*f);
    if ((b*f - c*e ) >= 0.0)
      xscale = rb2e2;
    else {
      b = -b;
      c = -c;
      xscale = -rb2e2;
    }
    yscale = rc2f2;

  /* Non-perpendicularity */
    nonperpendicularity = ((c != 0.0 || f != 0.0) ? atan2(c,f) : 0.0) +
                          ((e != 0.0 || b != 0.0) ? atan2(e,b) : 0.0);

  /* Orientation */
    double ws = (c*rb2e2) - (e*rc2f2);
    double wc = (b*rc2f2) + (f*rb2e2);
    orientation = (ws != 0.0 || wc != 0.0) ? atan2 (ws,wc) : 0.0;

  /* Zero corrections */
    double hp = nonperpendicularity / 2.0;
    double shp = sin(hp);
    double chp = cos(hp);
    double sor = sin(orientation);
    double cor = cos(orientation);
    double det = xscale*yscale*(chp + shp)*(chp - shp);
    if (fabs(det) > 0.0) {
      x0 = yscale * (a*((chp*cor) - (shp*sor)) - d*((chp*sor) + (shp*cor)))/det;
      y0 = xscale * (a*((chp*sor) - (shp*cor)) + d*((chp*cor) + (shp*sor)))/det;
    }
    else {
      x0 = 0.0;
      y0 = 0.0;
    }
    if (verbose) {
      printf("Coefficients:\n");
      for (auto c : coeffs) printf("%.4e ", c);
      printf("\n   Coefficient Decomposition:\n");
      printf("\t                 x0 = %.3f\n", x0);
      printf("\t                 y0 = %.3f\n", y0);
      printf("\t             xScale = %.3f arcsecs/pixel\n", xscale*648000/M_PI);
      printf("\t             yScale = %.3f arcsecs/pixel\n", yscale*648000/M_PI);
      printf("\tnonperpendicularity = %.3f\n", nonperpendicularity);
      printf("\t        orientation = %.3f\n", orientation);
    }
  }
};

#endif