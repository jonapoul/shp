#ifndef PLATE_H
#define PLATE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <utility>
#include <boost/algorithm/string.hpp>
#include "Coords.h"
using namespace std;

class Plate {
private:
	int m_id;				// plate identification number
	Coords m_coords;		// combines RA and DEC to one object
	string m_gregorian;		// UT gregorian date string, formatted as yymmdd
	double m_julian;		// julian date in the middle of exposure
	double m_exp;			// exposure time in seconds
	char m_grade;			// graded plate quality, A being best and C being worst
	double m_magLimit;		// faintest object that can be seen on the plate, decided by the emulsion/filter/etc

public:
	Plate(const int num=0, const Coords& c={}, const string& g="", const double j=0.0, const double exp=0.0, const char grade=' ', const double lim=0.0)
		: m_id(num), m_coords(c), m_gregorian(g), m_julian(j), m_exp(exp), m_grade(grade), m_magLimit(lim) { }
	Plate(const Plate& p)
		: m_id(p.m_id), m_coords(p.m_coords), m_gregorian(p.m_gregorian), m_julian(p.m_julian), m_exp(p.m_exp), m_grade(p.m_grade), m_magLimit(p.m_magLimit) { }

	inline int id() const { return m_id; };
	inline Coords coords() const { return m_coords; };
	inline string gregorian() const { return m_gregorian; }
	inline double julian() const { return m_julian; };
	inline double exposure() const { return m_exp; }
	inline char grade() const { return m_grade; };
	inline double magLimit() const { return m_magLimit; }
	inline void setID(const int id) { m_id = id; };
	inline void setCoords(const Coords& c) { m_coords = c; };
	inline void setGregorian(const string& g) { m_gregorian = g; }
	inline void setJulian(const double j) { m_julian = j; };
	inline void setExposure(const double e) { m_exp = e; }
	inline void setGrade(const char grade) { m_grade = grade; };
	inline void setMagLimit(const double lim) { m_magLimit = lim; }

	/*
		Reads a line from the plate catalog file and fills in the relevant fields of the Plate object
		buffer = full string representing all of one plate record
		Outputs a boolean flag to indicate whether the plate is valid for our purposes
	*/
	bool parsePlateString(const string& buffer) {
		// plate suffix column
		// if this isn't blank it indicates shenanigans while recording the image
		// T = tracked shot, M = multiple shots, P = full-aperture prism (distorted)
		if (buffer[7] == 'T' || buffer[7] == 'M' || buffer[7] == 'P') 
			return false;
		// checking for the word "TEST", so we can throw it out
		if (buffer.substr(16,4) == "TEST" || buffer.substr(15,4) == "TEST")
			return false;
		// plate number, for later reference
		m_id = stoi(buffer.substr(2, 5));
		// reading the RA/DEC coordinates in sexagesimal format
		m_coords.parseCoordsFromPlate(buffer);
		// julian date calculated from gregorian
		m_gregorian = buffer.substr(30, 6);
		string lst = buffer.substr(36, 4);
		// if any of the digits of the lst are blank spaces/letters, reject it
		if (!isdigit(lst[0]) || 
		    !isdigit(lst[1]) || 
		    !isdigit(lst[2]) || 
		    !isdigit(lst[3]))
			return false;
		m_julian = convertDate(m_gregorian, lst);

		// adding half of the exposure time to the julian date
		// this means that the outputted position is the position of the object halfway through the exposure
		// the substring is in the format "mmmt" where t = tenth of a minute
		m_exp = stod(buffer.substr(52,4))/14400.0;
		m_julian += m_exp/2.0;
		// plate quality grade, from A to C
		m_grade = buffer[56];
		// calculating the faintest apparent magnitude that would be visible on the plate
		m_magLimit = limitingMagnitude(buffer);

		return true;
	}

	/*
		Prints out a plate object's info. Used for testing
	*/
	void printPlate() const {
		printf("ID = %6d ", m_id);
		printf("RA = %s ", m_coords.RAtoString().c_str());
		printf("DEC = %s ", m_coords.DECtoString().c_str());
		printf("GD = %s ", m_gregorian.c_str());
		printf("JD = %.2f ", m_julian);
		printf("Grade = \'%c\'\n", m_grade);
	}	

	/* 
		Takes in a 6-char date string, in the format "yymmdd", outputs a Julian date as a double.
		This assumes that all dates are between 1st Jan 1917 and 31st Dec 2016, which is fine for
		these plate catalogs since they don't go beyond ~2003

		Shamelessly pilfered from "Practical Astronomy with your Calculator or Spreadsheet"
	*/
	static double convertDate(const string& gregorianDate, const string& lst) {
		int year = stoi(gregorianDate.substr(0, 2));
		if (year < 100)
			year += (year < 17) ? 2000 : 1900;
		int month = stoi(gregorianDate.substr(2, 2));
		int day   = stoi(gregorianDate.substr(4, 2));
		int hour  = stoi(lst.substr(0, 2));
		int mins  = stoi(lst.substr(2, 2));

		double julian = gregorianToJulian(double(day), month, year);
		double gst = LSTtoGST(hour, mins, 0.0, 149.07);
		double ut = GSTtoUT(gst, julian);
		double frac = (ut - 12.0) / 24.0;
		julian += frac;

		return julian;
	}

	/*
		Converts a UT Gregorian date to a floating point Julian date
	*/
	static double gregorianToJulian(double day, int month, int year) {
		if (month < 3) { 
			year--; 
			month += 12; 
		}
		int B;
		if (year > 1582) {
			int A = int(year / 100.0);
			B = 2 - A + int(A/4);
		} else B = 0;
		int C;
		if (year < 0) {
			C = int(365.25 * year - 0.75) ;
		} else C = int(365.25 * year);
		int D = int(30.6001 * (month+1));
		return B + C + D + day + 1720994.5;
	}

	/*
		Converting a floating poit Julian Date to UT Gregorian date string as "dd/mm/yyyy"
	*/
	static string julianToGregorian(double jd) {
		jd += 0.5;
		int I = int(jd);
		double F = jd - I;
		int B;
		if (I > 2299160) {
			int A = int( (I-1867216.25)/36524.25 );
			B = I + A - int(A/4.0) + 1;
		}
		else B = I;
		int C = B + 1524;
		int D = int( (C-122.1)/365.25 );
		int E = int( 365.25*D );
		int G = int( (C-E)/30.6001 );
		int d = int(C - E + F - int(30.6001*G));
		int m = (G < 13.5) ? G-1 : G-13;
		int y = (m > 2.5) ? D-4716 : D-4715;
		
		string output;
		if (d < 10) output += '0';
		output += to_string(d) + '/';
		if (m < 10) output += '0';
		output += to_string(m) + '/';
		output += to_string(y);
		return output;
	}

	/*
		Takes the LST at Siding Springs observatory and returns the Greenwich Sidereal Time
	*/
	static double LSTtoGST(const int hour, const int min, const double sec, const double longitude) {
		double lstF = hour + (min/60.0) + (sec/3600.0);
		double longHour = longitude/15.0; 	// in hours
		double gst = lstF - longHour;
		while (gst > 24) gst -= 24;
		while (gst < 0) gst += 24;
		return gst;
	}

	/*
		Converts decimal Greenwich Sidereal Time to decimal Universal Time
		Also uses the precalculated julian date
	*/
	static double GSTtoUT(const double gst, const double JD) {
		double S = JD - 2451545;
		double T = S / 36525.0;
		double T0 = 6.697374558 + (2400.051336 * T) + (0.000025862 * T * T);
		while (T0 > 24) T0 -= 24;
		while (T0 < 0)  T0 += 24;
		double A = gst - T0;
		while (A > 24) A -= 24;
		while (A < 0)  A += 24;
		double UT = A * 0.9972695663;
		return UT;
	}

	/*
		Takes a Plate array reference and a filename
		Opens the filename and reads all valid plate records into the array
	*/
	static void readPlateCatalog(vector<Plate>& plates, const string& filename) {
		ifstream platesFile(filename);
		if (platesFile.is_open()) {
			while (!platesFile.eof()) {
				string buffer;
				getline(platesFile, buffer);
				if (buffer.length() > 0) {
					Plate p;
					// only adds it to the array if the plate is valid
					if (p.parsePlateString(buffer))
						plates.push_back(p);
				}
			}
			platesFile.close();
		} else {
			cout << "Plate file \"" << filename << "\" is not valid\n";
		}
	}

	/*
		Prints all relevant info about a plate/ephemeris match
	*/
	static void printMatch(const Plate& p, const Coords& interp, const int count, const double mag, const pair<double,double>& start, const pair<double,double>& mid, const pair<double,double>& end) {

		// coordinate conversion from xi/eta to x/y coordinates from bottom left of plate (in mm)
		double degreesToMillimetres = 3600.0 / 67.12;
		double radsToDegs = 180.0 / M_PI;
		double x = (mid.first  * degreesToMillimetres * radsToDegs) + (354.5/2.0);
		double y = (mid.second * degreesToMillimetres * radsToDegs) + (354.5/2.0);

		// calculating the drift distance between the start and end of exposure (in mm)
		double dx = end.first  - start.first;
		double dy = end.second - start.second;
		double drift = sqrt(dx*dx + dy*dy);

		printf("%03d", count);
		printf("\tid       = %d\n", p.m_id);
		printf("\tdate     = %s\n", gregorianToString(p.m_gregorian).c_str());
		printf("\tjulian   = %.3f\n", p.m_julian);
		printf("\tcoords   = (%.3f, %.3f) deg\n", interp.getDegRA(), interp.getDegDEC());
		printf("\tstartPos = (%.2f, %.2f) mm\n", start.first, start.second);
		printf("\tmidPos   = (%.2f, %.2f) mm\n", x, y);
		printf("\tfinalPos = (%.2f, %.2f) mm\n", end.first, end.second);
		printf("\tdrift    = %.2f mm\n", drift);
		printf("\tmagnitude= %.2f\n", mag);
		printf("\tgrade    = %c\n", p.m_grade);
		printf("\texposure = %.1f mins\n", p.m_exp*14400);
		printf("-----------------------------------------\n");
	}
	
	/*
		Takes the plate prefix/emulsions/filters from the plate record and calculates the 
		magnitude of the faintest possible object on that plate.

		Return values pulled from http://www.roe.ac.uk/ifa/wfau/ukstu/telescope.html#fe
	*/
	static double limitingMagnitude(const string& buffer) {
		string prefix = buffer.substr(0, 2);
		boost::algorithm::trim(prefix);
		if (prefix == "U" || prefix == "B" || prefix == "V") 
			return 21.0;
		else if (prefix == "J" || prefix == "BJ" || prefix == "OR")
			return 22.5;
		else if (prefix == "R" || prefix == "HA")
			return 21.5;
		else if (prefix == "I")
			return 19.0;
		else if (prefix == "OR") {
			string emulsion = buffer.substr(40, 6);
			boost::algorithm::trim(emulsion);
			if (emulsion == "IIIa-F")
				return 21.5;
			else
				return 22.5;
		}
		return 23.0;	// default value
	}

	/*
		Converting a gregorian date string from "yymmdd" to "dd/mm/yyyy"
	*/
	static string gregorianToString(const string& greg) {
		int year  = stoi(greg.substr(0,2));
		int month = stoi(greg.substr(2,2));
		int day   = stoi(greg.substr(4,2));
		year += (year < 17) ? 2000 : 1900;
		string date = "";
		if (day < 10) date += '0';
		date += to_string(day) + '/';
		if (month < 10) date += '0';
		date += to_string(month) + '/';
		date += to_string(year);
		return date;
	}

	/*
		Calculates the xi/eta coordinates of the object at the start and end of the exposure time
		This is to account for the potential dragging across the image for longer exposure times.
		Returns the values as an std::pair<double> object, with xi as first and eta in second
	*/
	static void exposureBoundaries(const Plate& p, const vector<Coords>& coords, const vector<double>& times, 
	                               const double degree, pair<double, double>& start, pair<double, double>& end) {
		double expTime   = p.exposure();
		double startTime = p.julian() - (expTime/2.0);
		double endTime   = p.julian() + (expTime/2.0);
		Coords startCoords = Coords::polyInterp(coords, times, startTime, degree);
		Coords endCoords   = Coords::polyInterp(coords, times, endTime,   degree);

		double xiStart, xiEnd, etaStart, etaEnd;
		int status1, status2;
		Coords::gnomonic(startCoords, p.m_coords, xiStart, etaStart, status1);
		Coords::gnomonic(endCoords,   p.m_coords, xiEnd,   etaEnd,   status2);

		double degreesToMillimetres = 3600.0 / 67.12;
		double radsToDegs = 180.0 / M_PI;
		double x1 = (xiStart  * radsToDegs * degreesToMillimetres) + (354.5/2.0);
		double y1 = (etaStart * radsToDegs * degreesToMillimetres) + (354.5/2.0);
		double x2 = (xiEnd    * radsToDegs * degreesToMillimetres) + (354.5/2.0);
		double y2 = (etaEnd   * radsToDegs * degreesToMillimetres) + (354.5/2.0);
		start = {x1, y1};
		end   = {x2, y2};
	}

};

#endif