#ifndef COORDS_H
#define COORDS_H

#include <iostream>
#include <string>
#include <math.h>
using namespace std;

class Coords {
private:
	double m_radRA;
	double m_radDEC;
	double m_degRA;
	double m_degDEC;

	static constexpr double degToRad = M_PI/180.0;
	static constexpr double radToDeg = 180.0/M_PI;

public:
	Coords(const double degRA = 0.0, const double degDEC = 0.0)
		: m_degRA(degRA), m_degDEC(degDEC) { m_radRA = degRA*degToRad; m_radDEC = degDEC*degToRad; }
	Coords(const Coords& c)
		: m_radRA(c.m_radRA), m_radDEC(c.m_radDEC), m_degRA(c.m_degRA), m_degDEC(c.m_degDEC) { }

	double getDegRA()  const { return m_degRA; }
	double getDegDEC() const { return m_degDEC; }
	double getRadRA()  const { return m_radRA; }
	double getRadDEC() const { return m_radDEC; }

	void setDegRA (const double r) { m_degRA = r;  m_radRA = r*degToRad; }
	void setDegDEC(const double d) { m_degDEC = d; m_radDEC = d*degToRad; }
	void setRadRA (const double r) { m_radRA = r;  m_degRA = r*radToDeg; }
	void setRadDEC(const double d) { m_radDEC = d; m_degDEC = d*radToDeg; }

	void parseCoordsFromPlate(const string& record);
	string RAtoString() const;
	string DECtoString() const;
	static double cosAngularDistance(const Coords& a, const Coords& b);
	static double angularDistance(const Coords& a, const Coords& b, const bool inDegrees = true);
	static void gnomonic(const Coords& c, const Coords& c0, double& x, double& y, int& status);
	static void inverseGnomonic(const double x, const double y, const Coords& c0, Coords& c);
	static void toCartesian(const Coords& c, double& x, double& y, double& z);
	static void toSpherical(const double x, const double y, const double z, Coords& c);
	static Coords interpolate(const Coords& c0, const double x0, const Coords& c1, const double x1, const double x);

	friend Coords operator+(const Coords& c1, const Coords& c2);
};


void Coords::parseCoordsFromPlate(const string& record) {
	int hour =	stoi(record.substr(20,2));
	int mins =	stoi(record.substr(22,2));
	int secs =	stoi(record.substr(24,1)) * 6.0;
	m_degRA = hour*15.0 + mins/4.0 + secs/240.0;
	m_radRA = m_degRA * degToRad;

	int degrees = stoi(record.substr(26,2));
	int arcmins = stoi(record.substr(28,2));
	int arcsecs = 0;
	bool isPositive = (record[25] != '-');
	m_degDEC = degrees + arcmins/600.0;
	m_degDEC *= (isPositive ? 1.0 : -1.0);
	m_radDEC = m_degDEC * degToRad;
}

string Coords::RAtoString() const {
	string output = "";
	double decimal = m_degRA / 15.0;
	int h = int(decimal);
	if (h < 10) output += '0';
	output += to_string(h) + 'h';
	decimal = (decimal - h) * 60.0;
	int m = int(decimal);
	if (m < 10) output += '0';
	output += to_string(m) + 'm';
	decimal -= m;
	int s = round(decimal * 60.0);
	if (s < 10) output += '0';
	output += to_string(s) + 's';
	return output;
}

string Coords::DECtoString() const {
	string output = "";
	double decimal = m_degDEC;
	if (decimal < 0) {
		output += '-';
		decimal *= -1;
	} else output += '+';
	int d = int(decimal);
	if (d < 10) output += '0';
	output += to_string(d) + '\370';
	decimal = (decimal - d) * 60.0;
	int m = int(decimal);
	if (m < 10) output += '0';
	output += to_string(m) + '\'';
	decimal -= m;
	int s = int(decimal * 60.0);
	if (s < 10) output += '0';
	output += to_string(s) + '\"';
	return output;
}

double Coords::cosAngularDistance(const Coords& a, const Coords& b) {
	double a_ra  = a.m_radRA;
	double a_dec = a.m_radDEC;
	double b_ra  = b.m_radRA;
	double b_dec = b.m_radDEC;
	return sin(a_dec)*sin(b_dec) + cos(a_dec)*cos(b_dec)*cos(a_ra-b_ra);
}

double Coords::angularDistance(const Coords& a, const Coords& b, const bool inDegrees) {
	return acos(cosAngularDistance(a, b)) * (inDegrees ? radToDeg : 1.0);
}

/*
**  Gnomonic projection of spherical coordinates onto tangent plane
**  Given:
**     c      	   Coords   RA/DEC of point to be projected in radians
**     c0          Coords   RA/DEC of tangent point in radians
**  Returned:
**     x,y    	   double   rectangular coordinates on tangent plane
**     *j          int      status:   0 = OK, star on tangent plane
**                                    1 = error, star too far from axis
**                                    2 = error, antistar on tangent plane
**                                    3 = error, antistar too far from axis
*/
void Coords::gnomonic(const Coords& c, const Coords& c0, double& x, double& y, int& status) {
	double ra   = c.m_radRA;
	double dec  = c.m_radDEC;
	double ra0  = c0.m_radRA;
	double dec0 = c0.m_radDEC;
	
	double sin_dec0 = sin(dec0);
	double cos_dec0 = cos(dec0);
	double sin_dec  = sin(dec);
    double cos_dec  = cos(dec);
    double dRA      = ra - ra0;
    double sin_dRA  = sin(dRA);
    double cos_dRA  = cos(dRA);

	/* Reciprocal of star vector length to tangent plane */
    double denom = cosAngularDistance(c, c0);

    double TINY = 1e-6;
	/* Handle vectors too far from axis */
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

	/* Compute tangent plane coordinates (even in dubious cases) */
   	x = cos_dec*sin_dRA / denom;
   	y = ( sin_dec*cos_dec0 - cos_dec*sin_dec0*cos_dRA ) / denom;
}

/*
**  Transform tangent plane coordinates into spherical.
**  Given:
**     xi,eta      double   tangent plane rectangular coordinates
**     raz,decz    double   spherical coordinates of tangent point
**  Returned:
**     *ra,*dec    double   spherical coordinates (0-2pi,+/-pi/2)
*/
void Coords::inverseGnomonic(const double x, const double y, const Coords& c0, Coords& c) {
	double dec0 = c0.m_radDEC;
	double ra0  = c0.m_radRA;
	double sin_dec0 = sin(dec0);
	double cos_dec0 = cos(dec0);
	double denom = cos_dec0 - y*sin_dec0;
	double atan = atan2(x, denom) + ra0;
	while (atan > 2.0*M_PI) atan -= 2.0*M_PI;
	while (atan < 0.0)	  	atan += 2.0*M_PI;
	double dec = atan2(sin_dec0 + y*cos_dec0, sqrt(x*x + denom*denom)) * radToDeg;
	c.setRadRA(atan);
	c.setRadDEC( atan2(sin_dec0 + y*cos_dec0, sqrt(x*x + denom*denom)) );
}

void Coords::toCartesian(const Coords &c, double &x, double &y, double &z) {
	double ra = c.m_radRA;
	double dec = c.m_radDEC;
	x = cos(dec) * cos(ra);
	y = cos(dec) * sin(ra);
	z = sin(dec);
}

void Coords::toSpherical(const double x, const double y, const double z, Coords& c) {
	double ra = atan2(y, x);
	while (ra < 0.0) 		ra += 2.0*M_PI;
	while (ra > 2.0*M_PI) 	ra -= 2.0*M_PI;
	c.setRadRA(ra);
	c.setRadDEC( M_PI/2.0 - atan2(sqrt(x*x + y*y), z) );
}

/*
	Linearly interpolates between two spherical polar coordinates.
	1. Transform both RA/DEC coordinates into 3D cartesian coordinates
	2. Linearly interpolate between each x, y z value using the time values t0, t and t1
	3. Normalise each value
	4. Transform back to spherical polars
*/
Coords Coords::interpolate(const Coords& c0, const double t0, const Coords& c1, const double t1, const double t) {
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

Coords operator+(const Coords& c1, const Coords& c2) {
	double ra  = c1.m_degRA  + c2.m_degRA;
	double dec = c1.m_degDEC + c2.m_degDEC;
	
	if (ra > 360.0)    ra -= 360.0;
	else if (ra < 0.0) ra += 360.0;

	if (dec > 90.0) {
		dec = 180.0 - dec;
		ra  = 360.0 - ra;
	} else if (dec < -90.0) {
		dec = -180.0 - dec;
		ra  = 360.0 - ra;
	}
	return Coords(ra, dec);
}

#endif