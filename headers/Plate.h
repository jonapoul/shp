#ifndef PLATE_H
#define PLATE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Coords.h"
using namespace std;

class Plate {
private:
	int m_id;				// plate identification number
	Coords m_coords;		// combines RA and DEC to one object
	string m_gregorian;		// UT gregorian date string, formatted as yymmdd
	double m_julian;		// julian date at the start of exposure
	char m_grade;			// graded plate quality, A being best and C being worst

public:
	Plate(int num=0, const Coords& c={}, const string& g="", const double j=0.0, const char grade=' ')
		: m_id(num), m_coords(c), m_gregorian(g), m_julian(j), m_grade(grade) { }
	Plate(const Plate& p)
		: m_id(p.m_id), m_coords(p.m_coords), m_gregorian(p.m_gregorian), m_julian(p.m_julian), m_grade(p.m_grade) { }

	int    getID() 		  const { return m_id; };
	Coords getCoords() 	  const { return m_coords; };
	string getGregorian() const { return m_gregorian; }
	double getJulian()    const { return m_julian; };
	char   getGrade()	  const { return m_grade; };

	void setID    	 (const int id) 	{ m_id = id; };
	void setCoords	 (const Coords& c)  { m_coords = c; };
	void setGregorian(const string& g)  { m_gregorian = g; }
	void setJulian 	 (const double j) 	{ m_julian = j; };
	void setGrade 	 (const char grade) { m_grade = grade; };

	bool parsePlateString(const string& buffer);
	void printPlate() const;
	static double convertDate(const string& gregorian, const string& lst);
	static double gregorianToJulian(float d, int m, int y);
	static double LSTtoGST(const int hour, const int min, const float sec);
	static double GSTtoUT(const double gst, const double JD);
	static void readPlateCatalog(vector<Plate>& plates, const string& filename);
	static void printMatch(const Plate& p, const Coords& interp, const double distance, const double x, const double y);

};

/*
	Reads a line from the plate catalog file and fills in the relevant fields of the Plate object
	Input: 
		buffer = full string representing all of one plate record
	Output:
		boolean flag to indicate whether the plate is valid for our purposes
*/
bool Plate::parsePlateString(const string& buffer) {
	// plate suffix column
	// if this isn't blank it indicates shenanigans while recording the image, such as a tracking shot or multiple images
	// T = tracked shot, M = multiple shots, P = full-aperture prism
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

	if (!isdigit(lst[0]) || 
	    !isdigit(lst[1]) || 
	    !isdigit(lst[2]) || 
	    !isdigit(lst[3]))
		return false;
	m_julian = convertDate(m_gregorian, lst);
	double exposureTime = stod(buffer.substr(52,4))/14400.0;
	m_julian += exposureTime;
	// plate quality grade, from A to C
	m_grade = buffer[56];

	return true;
}

/*
	Prints out a plate object's info. Used for testing
*/
void Plate::printPlate() const {
	printf("ID = %6d ", m_id);
	printf("RA = %s ", m_coords.RAtoString().c_str());
	printf("DEC = %s ", m_coords.DECtoString().c_str());
	printf("GD = %s ", m_gregorian.c_str());
	printf("JD = %.2f ", m_julian);
	printf("Grade = %c\n", m_grade);
}

/* 
	Takes in a 6-char date string, in the format "yymmdd", outputs a Julian date as a float.
	This assumes that all dates are between 1st Jan 1917 and 31st Dec 2016, which is fine for
	these plate catalogs since they don't go beyond ~2003

	Shamelessly pilfered from "Practical Astronomy with your Calculator or Spreadsheet"
*/
double Plate::convertDate(const string& gregorianDate, const string& lst) {
	int year = stoi(gregorianDate.substr(0, 2));
	if (year < 100) {
		year += (year < 17) ? 2000 : 1900;
	}
	int month = stoi(gregorianDate.substr(2, 2));
	int day = stoi(gregorianDate.substr(4, 2));
	int hour = stoi(lst.substr(0, 2));
	int mins = stoi(lst.substr(2, 2));

	double julian = gregorianToJulian(day, month, year);
	double gst = LSTtoGST(hour, mins, 0.0);
	double ut = GSTtoUT(gst, julian);
	double frac = (ut - 12.0) / 24.0;
	julian += frac;

	return julian;
}

double Plate::gregorianToJulian(float d, int m, int y) {
	if (m < 3) { 
		y--; 
		m += 12; 
	}
	int A = y / 100;
	int B = 2 - A + (A/4);
	int C = int(365.25 * y);
	int D = int(30.6001 * (m+1));

	return B + C + D + d + 1720994.5;
}

/*
	Takes the LST at Siding Springs observatory and returns the Greenwich Sidereal Time
*/
double Plate::LSTtoGST(const int hour, const int min, const float sec) {
	double lstF = hour + (min/60.0) + (sec/3600.0);
	double longitude = 149.07/15; 	// in hours, the 149.07 taken from UKSTU website
	double gst = lstF - longitude;
	while (gst > 24) gst -= 24;
	while (gst < 0) gst += 24;
	return gst;
}

/*
	Converts decimal Greenwich Sidereal Time to decimal Universal Time
	Also uses the precalculated julian date
*/
double Plate::GSTtoUT(const double gst, const double JD) {
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
void Plate::readPlateCatalog(vector<Plate>& plates, const string& filename) {
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
	}
	else {
		cout << "Plate file \"" << filename << "\" is not valid\n";
	}
}

void Plate::printMatch(const Plate& p, const Coords& interp, const double distance, const double x, const double y) {
	printf("plateID = %5d\n", p.m_id);
	printf("\tplateRA  =%10.4f deg\n", p.m_coords.getDegRA());
	printf("\tplateDEC =%10.4f deg\n", p.m_coords.getDegDEC());
	printf("\tephemRA  =%10.4f deg\n", interp.getDegRA());
	printf("\tephemDEC =%10.4f deg\n", interp.getDegDEC());
	printf("\tdistance =%10.4f deg\n", distance);
	printf("\tUTdate   =%10s\n", p.m_gregorian.c_str());
	printf("\tXi       =%10.4f\n", x);
	printf("\tEta      =%10.4f\n", y);
}


#endif