#ifndef COORDS_H
#define COORDS_H
#include <iostream>
#include <vector>
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
	Coords(const double degRA = 0.0, 
	       const double degDEC = 0.0)
			: m_degRA(degRA), m_degDEC(degDEC) { 
		m_radRA = degRA*degToRad; 
		m_radDEC = degDEC*degToRad; 
	}
	Coords(const Coords& c)
			: m_radRA(c.m_radRA), m_radDEC(c.m_radDEC), m_degRA(c.m_degRA), m_degDEC(c.m_degDEC) { 
	}

	inline double getDegRA() const { return m_degRA; }
	inline double getDegDEC() const { return m_degDEC; }
	inline double getRadRA() const { return m_radRA; }
	inline double getRadDEC() const { return m_radDEC; }
	inline void setDegRA(const double r) { m_degRA = r; m_radRA = r*degToRad; }
	inline void setDegDEC(const double d) { m_degDEC = d; m_radDEC = d*degToRad; }
	inline void setRadRA(const double r) { m_radRA = r;  m_degRA = r*radToDeg; }
	inline void setRadDEC(const double d) { m_radDEC = d; m_degDEC = d*radToDeg; }

	/*
		Takes a plate catalog record string and pulls the RA/DEC coordinates from it
	*/
	void parseCoordsFromPlate(const string& record) {
		double hour = stod(record.substr(20,2));
		double mins = stod(record.substr(22,2));
		double secs = stod(record.substr(24,1)) * 6.0;
		m_degRA = hour*15.0 + mins/4.0 + secs/240.0;
		m_radRA = m_degRA * degToRad;

		double degrees = stod(record.substr(26,2));
		double arcmins = stod(record.substr(28,2));
		double arcsecs = 0.0;
		bool isPositive = (record[25] != '-');
		m_degDEC = degrees + arcmins/60.0 + arcsecs/3600.0;
		m_degDEC *= (isPositive ? 1.0 : -1.0);
		m_radDEC = m_degDEC * degToRad;
	}

	/*
		Converts a RA coordinate to sexagesimal format:
			" XXhXXmXXs "
	*/
	string RAtoString() const {
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

	/*
		Converts a DEC coordinate to sexagesimal format:
			" XX^XX'XX" "
	*/
	string DECtoString() const {
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

	/*
		Returns the ra/dec of the coordinates as a string (in decimal degrees)
	*/
	string toString() const {
		string output = to_string(m_degRA);
		output += ", " + to_string(m_degDEC);
		return output;
	}

	/*
		Calculates the cosine of the angular distance between two spherical coordinates
	*/
	static double cosAngularDistance(const Coords& a, 
	                                 const Coords& b) {
		return sin(a.m_radDEC)*sin(b.m_radDEC) + cos(a.m_radDEC)*cos(b.m_radDEC)*cos(a.m_radRA - b.m_radRA);
	}

	/*
		Calculating the angular distance between two spherical points.
		Defaults to output in degrees
	*/
	static double angularDistance(const Coords& a, 
	                              const Coords& b, 
	                              const bool inDegrees = true) {
		return acos(cosAngularDistance(a, b)) * (inDegrees ? radToDeg : 1.0);
	}

	/*
		Gnomonic projection of spherical coordinates onto tangent plane
		Given:
			c 			Coords 		RA/DEC of point to be projected
			c0 			Coords 		RA/DEC of tangent point
		Returned:
			xi,eta		double 		rectangular coordinates on tangent plane
			status		int 		status: 0 = OK, star on tangent plane
											1 = error, star too far from axis
											2 = error, antistar on tangent plane
											3 = error, antistar too far from axis
	*/
	static void gnomonic(const Coords& c, 
	                     const Coords& c0, 
	                     double& xi, 
	                     double& eta, 
	                     int& status) {
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

		// Cosine of the angular distance between the object and the tangent point
	    double denom = cosAngularDistance(c, c0);

	    // Handle vectors too far from axis
	    double TINY = 1e-6;
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
	}

	/*
		Transform tangent plane coordinates into spherical.
		Given:
			xi,eta 		double   tangent plane rectangular coordinates
			c0 			Coords   spherical coordinates of tangent point
		Returned:
			c 			Coords   spherical coordinates of the point at (xi,eta)
	*/
	static void inverseGnomonic(const double xi, 
	                            const double eta, 
	                            const Coords& c0, 
	                            Coords& c) {
		double dec0 	= c0.m_radDEC;
		double ra0  	= c0.m_radRA;
		double sin_dec0 = sin(dec0);
		double cos_dec0 = cos(dec0);
		double denom 	= cos_dec0 - eta*sin_dec0;
		double atan 	= atan2(xi, denom) + ra0;
		while (atan > 2.0*M_PI) atan -= 2.0*M_PI;
		while (atan < 0.0)	  	atan += 2.0*M_PI;
		c.setRadRA(atan);
		c.setRadDEC( atan2(sin_dec0 + eta*cos_dec0, sqrt(xi*xi + denom*denom)) );
	}

	/*
		Converts spherical RA/DEC coordinates into 3D cartesian coordinates, assuming unit sphere
	*/
	static void toCartesian(const Coords &c, 
	                        double &x, 
	                        double &y, 
	                        double &z) {
		double ra  = c.m_radRA;
		double dec = c.m_radDEC;
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
};

#endif