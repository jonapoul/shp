#ifndef PLATE_H
#define PLATE_H

#include <vector>
#include <fstream>
#include <sstream>
#include "Coords.h"
using namespace std;

class Plate {
private:
	int m_id;				// plate identification number
	Coords m_coords;		// combines RA and DEC to one object
	double m_julianDate;		// julian date at the start of exposure
	char m_grade;			// graded plate quality, A being best and C being 

public:
	Plate(int num=0, const Coords& c={}, const double jd=0.f, const string& lst="", const char grade=' ')
		: m_id(num), m_coords(c), m_julianDate(jd), m_grade(grade) { }
	Plate(const Plate& p)
		: m_id(p.m_id), m_coords(p.m_coords), m_julianDate(p.m_julianDate), m_grade(p.m_grade) { }

	int    getID()     const { return m_id; };
	Coords getCoords() const { return m_coords; };
	RA     getRA()     const { return m_coords.getRA(); };
	DEC    getDEC()    const { return m_coords.getDEC(); };
	double  getDate()  const { return m_julianDate; };
	char   getGrade()  const { return m_grade; };

	void setID    (const int id) 	 { m_id = id; };
	void setCoords(const Coords& c)  { m_coords = c; };
	void setRA    (const RA& ra) 	 { m_coords.setRA(ra); };
	void setDEC   (const DEC& dec) 	 { m_coords.setDEC(dec); };
	void setDate  (const double d) 	 { m_julianDate = d; };
	void setGrade (const char grade) { m_grade = grade; };

	bool parsePlateString(const string& buffer);
	void printPlate() const;
	static double convertDate(const string& gregorian, const string& lst);
	static int dayOfTheWeek(const int d, const int m, const int y);
	static int hourShift(const int day, const int month, const int year);
	static double gregorianToJulian(float d, int m, int y);
	static float LSTtoGST(const int hour, const int min);
	static float GSTtoUT(const float gst, const float JD);
	static void readPlateCatalog(vector<Plate>& plates, const string& filename);
};

/*
	Reads a line from the plate catalog file and fills in the relevant fields of the Plate object
*/
bool Plate::parsePlateString(const string& buffer) {
	// plate suffix column
	// if this isn't blank it indicates shenanigans while recording the image, such as a tracking shot or multiple images
	// T = tracked shot, M = multiple shots, P = full-aperture prism
	if (buffer[7] == 'T' || 
	    buffer[7] == 'M' || 
	    buffer[7] == 'P') 
		return false;
	// checking for the word "TEST", so we can throw it out
	if (buffer.substr(16,4) == "TEST" || 
	    buffer.substr(15,4) == "TEST")
		return false;
	// plate number, for later reference
	m_id = stoi(buffer.substr(2, 5));
	// reading the RA/DEC coordinates in sexagesimal format
	m_coords.parseCoordsFromPlate(buffer);
	// julian date calculated from gregorian
	string gregorianDate = buffer.substr(30, 6);
	string lst = buffer.substr(36, 4);
	if (!isdigit(lst[0]) || 
	    !isdigit(lst[1]) || 
	    !isdigit(lst[2]) || 
	    !isdigit(lst[3]))
		return false;
	m_julianDate = convertDate(gregorianDate, lst);
	// plate quality grade, from A to C
	m_grade = buffer[56];

	return true;
}

/*
	Prints out a plate object's info
	Used for testing
*/
void Plate::printPlate() const {
	printf("ID = %6d RA = %s ", m_id, m_coords.getRA().toString().c_str());
	printf("DEC = %s Date = %.2f Quality = %c\n", m_coords.getDEC().toString().c_str(), m_julianDate, m_grade);
}

/*
	Calculates the day of the week as an integer fro 0-6, representing Sunday=0, Monday=1... Saturday=6
*/
int Plate::dayOfTheWeek(const int d, const int m, const int y) {
	float julian = gregorianToJulian(d, m, y);
	float A = (julian + 1.5) / 7.f;
	float f = A - int(A);
	int n = int(7*f + 0.5);
	return n;
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

	// calculating the time difference between GMT and local time, either +10 or +11h
	int timeDiff = hourShift(day, month, year);
	double julian = gregorianToJulian(day, month, year);
	float gst = LSTtoGST(hour, mins);
	float ut = GSTtoUT(gst, julian);
	if (ut < 12) ut += 24;
	float frac = (ut - 12) / 24;
	julian += frac;

	return julian;
}

int Plate::hourShift(const int day, const int month, const int year) {
	int firstSunday, firstOfTheMonth, hourShift = 10;
	if (month > 10 || month < 4) {
	 	return 11;
	} else if (month == 4) {
		if (day > 7) {
			return 10;
		} else {
			firstOfTheMonth = dayOfTheWeek(1, 4, year);
			firstSunday = (firstOfTheMonth == 0) ? 1 : 8-firstOfTheMonth;
			if (day < firstSunday)
				return 11;
		}
	} else if (month == 10) {
		if (day > 7) {
			return 11;
		} else {
			firstOfTheMonth = dayOfTheWeek(1, 10, year);
			firstSunday = (firstOfTheMonth == 0) ? 1 : 8-firstOfTheMonth;
			if (day >= firstSunday)
				return 11;
		}
	}
	return 10;
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
float Plate::LSTtoGST(const int hour, const int min) {
	float lstF = hour + (min / 60.f);
	float longitude = 149.0661 / 15; 	// in hours, the 149 is longitude of the observatory in degrees
	float gst = lstF - longitude;
	while (gst > 24) gst -= 24;
	while (gst < 0) gst += 24;
	return gst;
}

/*
	Converts decimal Greenwich Sidereal Time to decimal Universal Time
*/
float Plate::GSTtoUT(const float gst, const float JD) {
	float S = JD - 2451545;
	float T = S / 36525.f;
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

#endif