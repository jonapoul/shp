#ifndef COORDS_H
#define COORDS_H

#include <iostream>
#include <string>
#include <math.h>
#include "RA.h"
#include "DEC.h"
using namespace std;

class Coords {
private:
	RA  m_ra;
	DEC m_dec;

public:
	Coords(const RA& r = {}, const DEC& d = {})
		: m_ra(r), m_dec(d) { }
	Coords(const Coords& c)
		: m_ra(c.m_ra), m_dec(c.m_dec) { }
	Coords(const double ra, const double dec)   // both of these in DEGREES, not radians
		: m_ra(ra), m_dec(dec) { m_ra.setRadians(ra*M_PI/180); m_dec.setRadians(dec*M_PI/180); }

	RA  getRA()  const { return m_ra; };
	DEC getDEC() const { return m_dec; };

	void setRA (const RA& r)  { m_ra = r; };
	void setDEC(const DEC& d) { m_dec = d; };

	void parseCoordsFromPlate(const string& record);
	static double cosAngularDistance(const Coords& a, const Coords& b);
	static double angularDistance(const Coords& a, const Coords& b);
	static void gnomonic(const Coords& c, const Coords& c0, double& x, double& y, int& status);
	static void inverseGnomonic(const double x, const double y, const Coords& c0, Coords& c);
	static void toCartesian(const Coords& c, double& x, double& y, double& z);
	static void toSpherical(const double x, const double y, const double z, Coords& c);
	static Coords interpolate(const Coords& c0, const double x0, const Coords& c1, const double x1, const double x);

	friend Coords operator+(const Coords& c1, const Coords& c2);
};


void Coords::parseCoordsFromPlate(const string& record) {
	if (record.length() < 31) {
		cerr << "string passed to Coords::parseFromPlateRecord() is too short\n";
		return;
	}
	m_ra = {
		stoi(record.substr(20,2)),
		stoi(record.substr(22,2)),
		stoi(record.substr(24,1)) * 6.f
	};
	m_dec = {
		!(record[25] == '-'),		// if the sign column is '-', set isPositive to true
		stoi(record.substr(26,2)),
		stoi(record.substr(28,2)),
		0.f
	};
}

double Coords::cosAngularDistance(const Coords& a, const Coords& b) {
	double a_ra  = a.m_ra.getRadians();
	double a_dec = a.m_dec.getRadians() * (a.m_dec.isPositive() ? 1 : -1);
	double b_ra  = b.m_ra.getRadians();
	double b_dec = b.m_dec.getRadians() * (b.m_dec.isPositive() ? 1 : -1);

	double cosx   = sin(a_dec)*sin(b_dec) + cos(a_dec)*cos(b_dec)*cos(a_ra-b_ra);
	return cosx;
}

double Coords::angularDistance(const Coords& a, const Coords& b) {
	return acos(cosAngularDistance(a, b)) * 180 / M_PI;			// in degrees
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
	double ra   = c.getRA().getRadians();
	double dec  = c.getDEC().getRadians();
	double ra0  = c0.getRA().getRadians();
	double dec0 = c0.getDEC().getRadians();
	
	/* Trig functions */
	double sin_dec0 = sin(dec0);
	double cos_dec0 = cos(dec0);
	double sin_dec  = sin(dec);
    double cos_dec  = cos(dec);
    double dRA      = ra - ra0;
    double sin_dRA  = sin(dRA);
    double cos_dRA  = cos(dRA);

	/* Reciprocal of star vector length to tangent plane */
    double denom = cosAngularDistance(c, c0);

    const double TINY = 1e-6;
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
	double dec0 = c0.getDEC().getRadians();
	double ra0  = c0.getRA().getRadians();
	double sin_dec0 = sin(dec0);
	double cos_dec0 = cos(dec0);
	double denom = cos_dec0 - y*sin_dec0;
	double atan = atan2(x, denom) + ra0;
	while (atan > 2*M_PI) atan -= 2*M_PI;
	while (atan < 0)	  atan += 2*M_PI;
	double ra = atan * 180/M_PI;
	double dec = atan2(sin_dec0 + y*cos_dec0, sqrt(x*x + denom*denom)) * 180/M_PI;
	c.setRA(ra);
	c.setDEC(dec);
}

void Coords::toCartesian(const Coords &c, double &x, double &y, double &z) {
	double ra = c.getRA().getRadians();
	double dec = c.getDEC().getRadians();
	x = cos(dec) * cos(ra);
	y = cos(dec) * sin(ra);
	z = sin(dec);
}

void Coords::toSpherical(const double x, const double y, const double z, Coords& c) {
	double ra = atan2(y, x) * 180/M_PI;
	double dec = (M_PI/2 - atan2(sqrt(x*x + y*y), z)) * 180/M_PI;
	c.setRA(ra);
	c.setDEC(dec);
}

/*
	Linearly interpolates between two spherical polar coordinates.
	1. Transform both coordinates into 3D cartesian coordinates
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
	double ra  = c1.getRA().getDegrees()  + c2.getRA().getDegrees();
	double dec = c1.getDEC().getDegrees() + c2.getDEC().getDegrees();
	
	if (ra > 360) 	 ra -= 360;
	else if (ra < 0) ra += 360;

	if (dec > 90) {
		dec = 180 - dec;
		ra  = 360 - ra;
	} else if (dec < -90) {
		dec = -180 - dec;
		ra  = 360 - ra;
	}
	RA  r(ra);
	DEC d(dec);
	r.fixRA();
	d.fixDEC();
	return Coords(r, d);
}

#endif