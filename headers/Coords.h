#ifndef COORDS_H
#define COORDS_H

/*
	This ifdef bit is to silence some Boost errors
	Robbed from http://www.vilipetek.com/2013/10/07/polynomial-fitting-in-c-using-boost/
*/
#ifdef BOOST_UBLAS_TYPE_CHECK
#undef BOOST_UBLAS_TYPE_CHECK
#endif
#define BOOST_UBLAS_TYPE_CHECK 0
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
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
	static double cosAngularDistance(const Coords& a, const Coords& b) {
		double a_ra  = a.m_radRA;
		double a_dec = a.m_radDEC;
		double b_ra  = b.m_radRA;
		double b_dec = b.m_radDEC;
		return sin(a_dec)*sin(b_dec) + cos(a_dec)*cos(b_dec)*cos(a_ra-b_ra);
	}

	/*
		Calculating the angular distance between two spherical points.
		Defaults to output in degrees
	*/
	static double angularDistance(const Coords& a, const Coords& b, const bool inDegrees = true) {
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
	static void gnomonic(const Coords& c, const Coords& c0, double& xi, double& eta, int& status) {
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
	static void inverseGnomonic(const double xi, const double eta, const Coords& c0, Coords& c) {
		double dec0 	= c0.m_radDEC;
		double ra0  	= c0.m_radRA;
		double sin_dec0 = sin(dec0);
		double cos_dec0 = cos(dec0);
		double denom 	= cos_dec0 - eta*sin_dec0;
		double atan 	= atan2(xi, denom) + ra0;
		while (atan > 2.0*M_PI) atan -= 2.0*M_PI;
		while (atan < 0.0)	  	atan += 2.0*M_PI;
		double dec = atan2(sin_dec0 + eta*cos_dec0, sqrt(xi*xi + denom*denom)) * radToDeg;
		c.setRadRA(atan);
		c.setRadDEC( atan2(sin_dec0 + eta*cos_dec0, sqrt(xi*xi + denom*denom)) );
	}

	/*
		Converts spherical RA/DEC coordinates into 3D cartesian coordinates, assuming unit sphere
	*/
	static void toCartesian(const Coords &c, double &x, double &y, double &z) {
		double ra  = c.m_radRA;
		double dec = c.m_radDEC;
		x = cos(dec) * cos(ra);
		y = cos(dec) * sin(ra);
		z = sin(dec);
	}

	/*
		Converts 3D cartesian coordinates to spherical RA/DEC, assuming unit sphere
	*/
	static void toSpherical(const double x, const double y, const double z, Coords& c) {
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
	static Coords linInterp(const Coords& c0, const double t0, const Coords& c1, 
							const double t1, const double t) {
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
		Interpolates between an arbitrary number of coordinate/time points using a
		least-squares fitting process.
			1) Convert all spherical coordinates (RA, DEC) to 3D cartesian (x, y, z)
			2) Perform a polynomial fit to a given degree
			3) Calculate the value of each (x, y, z) value at a time t
			4) Normalise the cartesian vector so the magnitude is exactly 1
			5) Convert each set back to RA/DEC
	*/
	static Coords polyInterp(const vector<Coords>& coords, const vector<double>& time, 
							 const double t, const int degree) {
		size_t size = coords.size();
		if (size != time.size()) {
			cout << "The two vectors passed to Coords::polyInterp() are different sizes: coords is ";
			cout << coords.size() << " and time is " << time.size() << ".\nExiting...\n\n";
			exit(1);
		}
		vector<double> x(size), y(size), z(size);
		for (size_t i = 0; i < size; i++)
			toCartesian(coords[i], x[i], y[i], z[i]);

		// quadratic fit, I would maybe do more but I'd lose precision in the julian dates
		vector<double> xCoeff = polyfit(time, x, degree);
		vector<double> yCoeff = polyfit(time, y, degree);
		vector<double> zCoeff = polyfit(time, z, degree);

		// calculating the interpolated x,y,z values based on the 3 fitted polynomials
		double xt = polyvalue(xCoeff, t);
		double yt = polyvalue(yCoeff, t);
		double zt = polyvalue(zCoeff, t);

		// normalising
		double mag = sqrt(xt*xt + yt*yt + zt*zt);
		xt /= mag;
		yt /= mag;
		zt /= mag;
		// sending it back to spherical
		Coords output;
		toSpherical(xt, yt, zt, output);
		return output;
	}

	/*
		Fits a set of X values to their respective Y values using a polynomial equation
		to a given degree of accuracy, using a least-squares fitting process via the 
		boost::numeric::ublas functions

		Shamelessly robbed from http://www.vilipetek.com/2013/10/07/polynomial-fitting-in-c-using-boost/
	*/
	template<typename T>
	static vector<T> polyfit(const vector<T>& oX, const vector<T>& oY, int nDegree) {
		using namespace boost::numeric::ublas;
		if (oX.size() != oY.size()) 
			throw invalid_argument( "X and Y vector sizes do not match" );

		// more intuative to increment nDegree
		nDegree++;
		size_t nCount = oX.size();
		matrix<T> oXMatrix(nCount, nDegree);
		matrix<T> oYMatrix(nCount, 1);
		
		// copy y matrix
		for (size_t i = 0; i < nCount; i++)
			oYMatrix(i, 0) = oY[i];

		// create the X matrix
		for (size_t nRow = 0; nRow < nCount; nRow++) {
			T nVal = 1.0;
			for (int nCol = 0; nCol < nDegree; nCol++) {
				oXMatrix(nRow, nCol) = nVal;
				nVal *= oX[nRow];
			}
		}
		// transpose X matrix
		matrix<T> oXtMatrix( trans(oXMatrix) );
		// multiply transposed X matrix with X matrix
		matrix<T> oXtXMatrix( prec_prod(oXtMatrix, oXMatrix) );
		// multiply transposed X matrix with Y matrix
		matrix<T> oXtYMatrix( prec_prod(oXtMatrix, oYMatrix) );
		// lu decomposition
		permutation_matrix<int> pert(oXtXMatrix.size1());
		const size_t singular = lu_factorize(oXtXMatrix, pert);
		// must be singular
		BOOST_ASSERT( singular == 0 );
		// backsubstitution
		lu_substitute(oXtXMatrix, pert, oXtYMatrix);
		// copy the result to coeff
		return std::vector<T>( oXtYMatrix.data().begin(), oXtYMatrix.data().end() );
	}

	/*
		Calculates the y-value of a point y = a0 + a1*x + a2*x^2 ... ai*x^i
		Given a vector of coefficients 'a' and an x coordinate to calculate for
	*/
	template<typename T>
	static T polyvalue(const vector<T>& coeff, const T x) {
		T y = 0.0;
		for (int i = 0; i < (int)coeff.size(); i++) {
			y += coeff[i] * pow(x, i);
		}
		return y;
	}
};

#endif