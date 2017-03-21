#ifndef COORDS_H
#define COORDS_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Definitions.h"
using std::cout;
namespace ublas = boost::numeric::ublas;

class Coords {
private:
  double radRA_;
  double radDEC_;
  double degRA_;
  double degDEC_; 

public:
  Coords(const double degRA = 0.0, 
         const double degDEC = 0.0)
      : degRA_(degRA), degDEC_(degDEC) { 
    radRA_  = degRA  * DEG_TO_RAD; 
    radDEC_ = degDEC * DEG_TO_RAD; 
  }
  Coords(const Coords& c)
      : radRA_(c.radRA_), radDEC_(c.radDEC_), 
        degRA_(c.degRA_), degDEC_(c.degDEC_) { 
  }
  /*
    Assumes RA is in format XXhXXmXX.XXXs and DEC in format +XX:XX:XX.XX
  */
  Coords(const std::string& ra,
         const std::string& dec) {
    int hrs = stoi(ra.substr(0,2));
    int min = stoi(ra.substr(3,2));
    double sec = stod(ra.substr(6,6));
    degRA_ = hrs*15.0 + min/4.0 + sec/240.0;
    radRA_ = degRA_ * DEG_TO_RAD;

    bool isPositive = (dec[0] == '+');
    int degrees = stoi(dec.substr(1,2));
    int arcmins = stoi(dec.substr(4,2));
    double arcsecs = stod(dec.substr(7,5));
    degDEC_ = degrees + arcmins/60.0 + arcsecs/3600.0;
    degDEC_ *= (isPositive ? 1.0 : -1.0);
    radDEC_ = degDEC_ * DEG_TO_RAD;
  }

  inline double ra(const AngleUnit& specifier) const { 
    return (specifier == DEG) ? degRA_ : radRA_;
  }

  inline double dec(const AngleUnit& specifier) const {
    return (specifier == DEG) ? degDEC_ : radDEC_;
  }

  inline void setRA(const double ra, 
                    const AngleUnit& specifier) {
    if (specifier == DEG) {
      degRA_ = ra;
      radRA_ = ra * DEG_TO_RAD;
    } else {
      radRA_ = ra;
      degRA_ = ra * RAD_TO_DEG;
    }
  }

  inline void setDEC(const double dec,
                     const AngleUnit& specifier) {
    if (specifier == DEG) {
      degDEC_ = dec;
      radDEC_ = dec * DEG_TO_RAD;
    } else {
      radDEC_ = dec;
      degDEC_ = dec * RAD_TO_DEG;
    }
  }

  friend std::ostream& operator<<(std::ostream& stream, 
                                  const Coords& c) {
    stream << c.toString();
    return stream;
  }

  /*
    Takes a plate catalog record string and pulls the RA/DEC coordinates from it
  */
  void parseCoordsFromPlate(const std::string& record) {
    double hour = stod(record.substr(20,2));
    double mins = stod(record.substr(22,2));
    double secs = stod(record.substr(24,1)) * 6.0;
    degRA_ = hour*15.0 + mins/4.0 + secs/240.0;
    radRA_ = degRA_ * DEG_TO_RAD;

    double degrees = stod(record.substr(26,2));
    double arcmins = stod(record.substr(28,2));
    double arcsecs = 0.0;
    bool isPositive = (record[25] != '-');
    degDEC_ = degrees + arcmins/60.0 + arcsecs/3600.0;
    degDEC_ *= (isPositive ? 1.0 : -1.0);
    radDEC_ = degDEC_ * DEG_TO_RAD;
  }

  /*
    Converts a RA coordinate to sexagesimal format:
      "XXhXXmXXs"
  */
  std::string RAtoString() const {
    std::stringstream output;
    double decimal = degRA_ / 15.0;
    int h = int(decimal);
    if (h < 10) output << '0';
    output << h << 'h';
    decimal = (decimal - h) * 60.0;
    int m = int(decimal);
    if (m < 10) output << '0';
    output << m << 'm';
    decimal -= m;
    double s = decimal * 60.0;
    if (s < 10) output << '0';
    output << std::fixed << std::setprecision(3) << s << 's';
    return output.str();
  }

  /*
    Converts a DEC coordinate to sexagesimal format:
      " XX^XX'XX" "
  */
  std::string DECtoString() const {
    std::stringstream output;
    double decimal = degDEC_;
    if (decimal < 0) {
      output << '-';
      decimal *= -1;
    } else output << '+';
    int d = int(decimal);
    if (d < 10) output << '0';
    output << d << 'd';
    decimal = (decimal - d) * 60.0;
    int m = int(decimal);
    if (m < 10) output << '0';
    output << m << '\'';
    decimal -= m;
    double s = decimal * 60.0;
    if (s < 10) output << '0';
    output << std::fixed << std::setprecision(3) << s << '\"';
    return output.str();
  }

  /*
    Returns the ra/dec of the coordinates as a string (in decimal degrees)
  */
  std::string toString(const bool isSexagesimal = false) const {
    if (isSexagesimal) 
      return RAtoString() + ", " + DECtoString();
    else
      return std::to_string(degRA_) + "\t" + std::to_string(degDEC_);
  }

  /*
    Calculates the cosine of the angular distance between two spherical coordinates
  */
  static double cosAngularDistance(const Coords& a, 
                                   const Coords& b) {
    return sin(a.radDEC_)*sin(b.radDEC_) + cos(a.radDEC_)*cos(b.radDEC_)*cos(a.radRA_ - b.radRA_);
  }

  /*
    Calculating the angular distance between two spherical points.
    Defaults to output in degrees
  */
  static double angularDistance(const Coords& a, 
                                const Coords& b, 
                                const bool inDegrees = true) {
    return acos(cosAngularDistance(a, b)) * (inDegrees ? RAD_TO_DEG : 1.0);
  }

  /*
    Gnomonic projection of spherical coordinates onto tangent plane
    Given:
      c       Coords    RA/DEC of point to be projected
      c0      Coords    RA/DEC of tangent point
    Returned:
      xi,eta    double    rectangular coordinates on tangent plane
      status    int     status: 0 = OK, star on tangent plane
                      1 = error, star too far from axis
                      2 = error, antistar on tangent plane
                      3 = error, antistar too far from axis
  */
  static int gnomonic(const Coords& c, 
                      const Coords& c0, 
                      double& xi, 
                      double& eta) {
    double ra   = c.radRA_;
    double dec  = c.radDEC_;
    double ra0  = c0.radRA_;
    double dec0 = c0.radDEC_;
    
    double sin_dec0 = sin(dec0);
    double cos_dec0 = cos(dec0);
    double sin_dec  = sin(dec);
    double cos_dec  = cos(dec);
    double dRA      = ra - ra0;
    double sin_dRA  = sin(dRA);
    double cos_dRA  = cos(dRA);

    // Cosine of the angular distance between the object and the tangent point
    double denom = cosAngularDistance(c, c0);

    // Handle vectors too far from axis
    double TINY = 1e-6;
    int status;
    if ( denom > TINY ) { 
      status = 0;
    } else if ( denom >= 0.0 ) {
        status = 1;
        denom = TINY;
    } else if ( denom > -TINY ) {
        status = 2;
        denom = -TINY;
    } else {
        status = 3;
    }
    // Compute tangent plane coordinates (even in dubious cases)
    xi  = ( cos_dec*sin_dRA ) / denom;
    eta = ( sin_dec*cos_dec0 - cos_dec*sin_dec0*cos_dRA ) / denom;
    return status;
  }

  /*
    Transform tangent plane coordinates into spherical.
    Given:
      xi,eta  double   tangent plane rectangular coordinates
      c0      Coords   spherical coordinates of tangent point
    Returned:
      c       Coords   spherical coordinates of the point at (xi,eta)
  */
  static void inverseGnomonic(const double xi, 
                              const double eta, 
                              const Coords& c0, 
                              Coords& c) {
    double dec0 = c0.radDEC_;
    double ra0  = c0.radRA_;
    double sin_dec0 = sin(dec0);
    double cos_dec0 = cos(dec0);
    double denom = cos_dec0 - eta*sin_dec0;
    double ra  = atan2(xi, denom) + ra0;
    while (ra > 2.0*M_PI) ra -= 2.0*M_PI;
    while (ra < 0.0)      ra += 2.0*M_PI;
    c.setRA(ra, RAD);
    double dec = atan2(sin_dec0 + eta*cos_dec0, sqrt(xi*xi + denom*denom));
    c.setDEC(dec, RAD);
  }

  /*
    Converts spherical RA/DEC coordinates into 3D cartesian coordinates, assuming unit sphere
  */
  static void toCartesian(const Coords &c, 
                          double &x, 
                          double &y, 
                          double &z) {
    double ra  = c.radRA_;
    double dec = c.radDEC_;
    x = cos(dec) * cos(ra);
    y = cos(dec) * sin(ra);
    z = sin(dec);
  }

  /*
    Converts 3D cartesian coordinates to spherical RA/DEC, assuming unit sphere
  */
  static void toSpherical(const double x, 
                          const double y, 
                          const double z, 
                          Coords& c) {
    double ra = atan2(y, x);
    while (ra < 0.0)    ra += 2.0*M_PI;
    while (ra > 2.0*M_PI)   ra -= 2.0*M_PI;
    c.setRA(ra, RAD);
    c.setDEC(M_PI/2.0 - atan2(sqrt(x*x + y*y), z), RAD);
  }

  /*
    Linearly interpolates between two spherical polar coordinates.
    1. Transform both RA/DEC coordinates into 3D cartesian coordinates
    2. Linearly interpolate between each x, y z value using the time values t0, t and t1
    3. Normalise each value
    4. Transform back to spherical polars
  */
  static Coords linInterp(const Coords& c0, 
                          const double t0, 
                          const Coords& c1, 
                          const double t1, 
                          const double t) {
    double x0, y0, z0, x1, y1, z1;
    toCartesian(c0, x0, y0, z0);
    toCartesian(c1, x1, y1, z1);
    double x = x0 + (t-t0)*(x1-x0)/(t1-t0);
    double y = y0 + (t-t0)*(y1-y0)/(t1-t0);
    double z = z0 + (t-t0)*(z1-z0)/(t1-t0);
    double mag = sqrt(x*x + y*y + z*z);
    x /= mag;
    y /= mag;
    z /= mag;
    Coords output;
    toSpherical(x, y, z, output);
    return output;
  }

  /*
    Converts the output of a gnomonic transformation to mm
      1) convert radians to arcseconds
      2) convert arcsecs to millimetres, using the given rate of 67.12"/mm
      3) add a shift of half the plate size
    This gives the output as relative to the bottom left corner of the plate
  */
  static double radsToMM(const double rads) {
    return (rads * 3600.0 * 180.0) / (ARCSECS_PER_MM * M_PI) + (PLATE_SIZE / 2.0);
  }

  /*
    Inverse of the above
  */
  static double mmToRads(const double mm) {
    return (mm - (PLATE_SIZE / 2.0)) * (ARCSECS_PER_MM * M_PI) / (3600.0 * 180.0);
  }

  static Coords mmToCoords(const double x,
                           const double y,
                           const Coords& tangentPoint) {
    double xi  = mmToRads(x);
    double eta = mmToRads(y);
    Coords output;
    inverseGnomonic(xi, eta, tangentPoint, output);
    return output;
  }

  /*
    Converts between two RA/DEC systems based on time epochs
      year1 = basis of coordinates c1
      year2 = year of output coordinates 

    HIGH PRECISION VERSION
    From Duffett-Smith p73-75
  */
  static Coords convertEpoch(const unsigned year1,
                             const unsigned year2,
                             const Coords& c1) {
    // i would use this from Plate.h, but circular includes are a pain in the arse
    auto gregorianToJulian = [](double day, int month, int year) -> double {
      if (month < 3) { 
        year--; 
        month += 12; 
      }
      int B;
      if (year > 1582) {
        int A = int(year / 100.0);
        B = 2 - A + int(A/4);
      } else B = 0;
      int C = (year < 0) ? int(365.25 * year - 0.75) : int(365.25 * year);
      int D = int(30.6001 * (month+1));
      return B + C + D + day + 1720994.5;
    };
    auto createMatrix = [](double T){
      double xi    = (0.640616*T + 0.0000839*T*T + 0.0000050*T*T*T) * DEG_TO_RAD;
      double z     = (0.640616*T + 0.0003041*T*T + 0.0000051*T*T*T) * DEG_TO_RAD;
      double theta = (0.556753*T - 0.0001185*T*T - 0.0000116*T*T*T) * DEG_TO_RAD;
      double CX = cos(xi),    SX = sin(xi);
      double CZ = cos(z),     SZ = sin(z);
      double CT = cos(theta), ST = sin(theta);
      ublas::matrix<double> m(3,3);
      m(0,0) =  CX*CT*CZ-SX*SZ; m(0,1) =  CX*CT*SZ+SX*CZ; m(0,2) = CX*ST;
      m(1,0) = -SX*CT*CZ-CX*SZ; m(1,1) = -SX*CT*SZ+CX*CZ; m(1,2) = -SX*ST;
      m(2,0) = -ST*CZ;          m(2,1) = -ST*SZ;          m(2,2) = CT;
      return m;
    };

    double JD1 = gregorianToJulian(1, 1, year1);
    double JD2 = gregorianToJulian(1, 1, year2);
    double T1 = (JD1 - 2451545.0) / 36525.0;
    double T2 = (JD2 - 2451545.0) / 36525.0;

    ublas::matrix<double> P_prime = createMatrix(T1);
    ublas::matrix<double> P = ublas::trans(createMatrix(T2));
    double ra1 = c1.ra(RAD), dec1 = c1.dec(RAD);
    ublas::vector<double> v(3);
    v(0) = cos(ra1)*cos(dec1);
    v(1) = sin(ra1)*cos(dec1);
    v(2) = sin(dec1);
    ublas::vector<double> s = prod(P_prime, v);
    ublas::vector<double> w = prod(P, s);
    double m = w(0), n = w(1), p = w(2);

    double ra2  = atan2(n, m) * RAD_TO_DEG;
    while (ra2 > 360) ra2 -= 360;
    while (ra2 < 0)   ra2 += 360;
    double dec2 = asin(p) * RAD_TO_DEG;
    return Coords(ra2, dec2);
  }

  static void coordsToPlatePosition(const Coords& asteroidCoords,
                                    const unsigned plateID,
                                    double& x,
                                    double& y) {
    std::string plateNumStr = std::to_string(plateID);
    Coords plateCentre;
    std::ifstream catalogfile("catalog.txt");
    if (catalogfile.is_open()) {
      while (plateNumStr.length() < 5)
        plateNumStr = ' ' + plateNumStr;
      for (std::string buffer; getline(catalogfile, buffer); ) {
        if (buffer.empty()) 
          continue;
        if (buffer.substr(2, 5) == plateNumStr) {
          plateCentre.parseCoordsFromPlate(buffer);
          break;
        }
      }
      catalogfile.close();
      double xi, eta;
      Coords::gnomonic(asteroidCoords, plateCentre, xi, eta);
      x = radsToMM(xi);
      y = radsToMM(eta);
    } else {
      cout << "catalog.txt couldn't be opened\n";
      exit(1);
    }
  }
};

#endif