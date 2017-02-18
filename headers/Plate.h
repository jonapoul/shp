#ifndef PLATE_H
#define PLATE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <utility>
#include <boost/algorithm/string.hpp>
#include "Coords.h"
#include "Ephemeris.h"
using namespace std;

#define PLATE_SIZE 354.5

class Plate {
private:
	int m_id;				// plate identification number
	Coords m_coords;		// combines RA and DEC to one object
	string m_gregorian;		// UT gregorian date string, formatted as yymmdd
	double m_julian;		// julian date in the middle of exposure
	double m_exp;			// exposure time in secs
	char m_grade;			// graded plate quality, A being best and C being worst
	double m_countLimit;	// lowest number of photon counts on a plate below which an object isn't visible

public:
	Plate(const int num=0, 
	      const Coords& c={}, 
	      const string& g="", 
	      const double j=0.0, 
	      const double exp=0.0, 
	      const char grade=' ', 
	      const double lim=0.0)
		: m_id(num), m_coords(c), m_gregorian(g), m_julian(j), m_exp(exp), m_grade(grade), m_countLimit(lim) { }
	Plate(const Plate& p)
		: m_id(p.m_id), m_coords(p.m_coords), m_gregorian(p.m_gregorian), m_julian(p.m_julian), m_exp(p.m_exp), m_grade(p.m_grade), m_countLimit(p.m_countLimit) { }

	inline int id() const { return m_id; };
	inline Coords coords() const { return m_coords; };
	inline string gregorian() const { return m_gregorian; }
	inline double julian() const { return m_julian; };
	inline double exposure() const { return m_exp; }
	inline char grade() const { return m_grade; };
	inline double countLimit() const { return m_countLimit; }
	inline void setID(const int id) { m_id = id; };
	inline void setCoords(const Coords& c) { m_coords = c; };
	inline void setGregorian(const string& g) { m_gregorian = g; }
	inline void setJulian(const double j) { m_julian = j; };
	inline void setExposure(const double e) { m_exp = e; }
	inline void setGrade(const char grade) { m_grade = grade; };
	inline void setCountLimit(const double lim) { m_countLimit = lim; }

	/*
		Reads a line from the plate catalog file and fills in the relevant fields of the Plate object
		buffer = full string representing all of one plate record
		Outputs a boolean flag to indicate whether the plate is valid for our purposes
	*/
	bool parsePlateString(const string& buffer) {
		// plate suffix column
		// if this isn't blank, it indicates shenanigans while recording the image
		// T = tracked shot, M = multiple shots, P = full-aperture prism (distorted)
		if (buffer[7] == 'T' || buffer[7] == 'M' || buffer[7] == 'P') 
			return false;
		// checking for the word "TEST", so we can throw it out
		if (buffer.substr(16,4) == "TEST" || buffer.substr(15,4) == "TEST")
			return false;
		// plate number, for later reference
		try { m_id = stoi(buffer.substr(2, 5)); }
		catch (...) { return false; }
		// reading the RA/DEC coordinates in sexagesimal format
		try { m_coords.parseCoordsFromPlate(buffer); }
		catch (...) { return false; }
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
		try { m_exp = stod(buffer.substr(52,4)) * 6.0; }
		catch (...) { return false; }
		m_julian += (m_exp / 86400.0) / 2.0;

		// plate quality grade, from A to C
		m_grade = buffer[56];
		// calculating the faintest apparent magnitude that would be visible on the plate
		m_countLimit = limitingCounts(buffer);

		// if no problems were found at any point, return true
		return true;
	}

	/*
		Checks whether the plate is missing from the plate archive room
	*/
	bool isMissing(const vector<int>& blacklist) {
		for (int id : blacklist) {
			if (this->id() == id) return true;
		}
		return false;
	}

	/* 
		Takes in a 6-char date string, in the format "yymmdd", outputs a Julian date as a double.
		This assumes that all dates are between 1st Jan 1917 and 31st Dec 2016, which is fine for
		these plate catalogs since they don't go beyond ~2003

		From "Practical Astronomy with your Calculator or Spreadsheet"
	*/
	static double convertDate(const string& gregorianDate, 
	                          const string& lst,
	                          bool isDebug = false) {
		int year = stoi(gregorianDate.substr(0, 2));
		if (year < 100)
			year += (year < 17) ? 2000 : 1900;
		int month = stoi(gregorianDate.substr(2, 2));
		int day   = stoi(gregorianDate.substr(4, 2));
		int hour  = stoi(lst.substr(0, 2));
		int mins  = stoi(lst.substr(2, 2));

		double julian = gregorianToJulian(double(day), month, year);
		double gst = LSTtoGST(hour, mins, 0.0, 149.0644);//149.07);
		double ut = GSTtoUT(gst, julian);
		double frac = (ut) / 24.0;
		julian += frac;

		if (isDebug) {		// date conversion debugging
			printf("\nJD   = %.5f\n", julian-frac);
			printf("GST  = %.3f\n", gst);
			printf("UT   = %.3f\n", ut);
			printf("FRAC = %.3f\n", frac);
			printf("JD   = %.5f\n", julian);
		}
		return julian;
	}

	/*
		Converts a UT Gregorian date to a floating point Julian date
	*/
	static double gregorianToJulian(const double day, 
	                                int month, 
	                                int year) {
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
	static double LSTtoGST(const int hour, 
	                       const int min, 
	                       const double sec, 
	                       const double longitude) {
		double lstDecimal = hour + (min/60.0) + (sec/3600.0);
		double longitudeHours = longitude/15.0; 	// decimal hours
		double gst = lstDecimal - longitudeHours;
		while (gst > 24) gst -= 24;
		while (gst < 0) gst += 24;
		return gst;
	}

	/*
		Converts decimal Greenwich Sidereal Time to decimal Universal Time
		Also uses the precalculated julian date
	*/
	static double GSTtoUT(const double gst, 
	                      const double JD) {
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
		Extracts the UT time for a given julian date, as a string in the format "hh:mm:ss.s"
	*/
	static string julianToUT(const double JD) {
		double zeroHour = floor(JD - 0.5) + 0.5;
		double remainder = (JD - zeroHour) * 24;

		int h    = int(remainder);
		remainder = (remainder - h) * 60.0;
		int m    = int(remainder);
		double s = (remainder - m) * 60.0;
		stringstream ss;
		if (h < 10) ss << '0';
		ss << h << ':';
		if (m < 10) ss << '0';
		ss << m << ':';
		if (s < 10) ss << '0';
		char sec[20];
		sprintf(sec, "%.0f", s);
		ss << sec;
		return ss.str();
	}

	/*
		Takes a Plate array reference and a filename
		Opens the filename and reads all valid plate records into the array
	*/
	static void readPlateCatalog(vector<Plate>& plates, 
	                             const string& filename) {
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
		Goes through all matched plates and prints the relevant info about them all
		This is a LITTLE BIT OF A MESS but it works
	*/
	static void printMatches(const vector<Plate>& p, 
	                         const vector<Coords>& c, 
	                         const vector<int>& count, 
	                         const vector<double>& mag, 
	                         const vector<double>& countLim, 
	                         const vector<pair<double,double>>& start, 
	                         const vector<pair<double,double>>& mid, 
	                         const vector<pair<double,double>>& end) {

		vector<pair<double,double>> middle;
		for (auto m : mid) {
			double x = radsToMM(m.first);
			double y = radsToMM(m.second);
			middle.push_back({x, y});
		}
		const int SIZE = 42;	// defines width of each printed data box
		size_t length;

		for (int i = 0; i < int(p.size()); i += 2) {
			bool canPrint = !(i == p.size() - 1 && p.size() % 2 == 1);
			stringstream ss;

			// Match number and Plate ID
			char buffer1[50], buffer2[50], buffer3[50], buffer4[50];
			sprintf(buffer1, "%03d", count[i]);
			sprintf(buffer2, "\tplate ID       = %d", p[i].id());
			ss << buffer1 << buffer2;
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length+3, ' ');
				sprintf(buffer3, "%03d", count[i+1]);
				sprintf(buffer4, "\tplate ID       = %d", p[i+1].id());
				ss << "| " << buffer3 << buffer4;
			}
			cout << ss.str() << '\n';
			ss.str("");

			// UT Date
			ss << "\tUT date/time   = " << gregorianToString(p[i].gregorian());
			ss << ", " << julianToUT(p[i].julian());
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				ss << "| " << "\tUT date        = " << gregorianToString(p[i+1].gregorian());
				ss << ", " << julianToUT(p[i+1].julian());
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Julian Date
			char buffer5[50], buffer6[50];
			sprintf(buffer5, "\tJulian date    = %.3f", p[i].julian());
			ss << buffer5;
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				sprintf(buffer6, "\tJulian date    = %.3f", p[i+1].julian());
				ss << "| " << buffer6;
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Object's RA/DEC coordinates in degrees
			char buffer7[50], buffer8[50];
			sprintf(buffer7, "\tObject Coords  = (%.3f, %.3f) deg", c[i].getDegRA(), c[i].getDegDEC());
			ss << buffer7;
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				sprintf(buffer8, "\tObject Coords  = (%.3f, %.3f) deg", c[i+1].getDegRA(), c[i+1].getDegDEC());
				ss << "| " << buffer8;
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Object's positional coordinates on the plate, in millimetres from the bottom left corner
			char buffer9[50], buffer10[50];
			sprintf(buffer9, "\tPlate Position = (%.3f, %.3f) mm", middle[i].first, middle[i].second);
			ss << buffer9;
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				sprintf(buffer10, "\tPlate Position = (%.3f, %.3f) mm", middle[i+1].first, middle[i+1].second);
				ss << "| " << buffer10;
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Length of the object's exposure trail on the plate in millimetres, plus the directional vector
			char buffer11[50], buffer12[50], buffer13[50], buffer14[50];
			double dx = end[i].first  - start[i].first;
			double dy = end[i].second - start[i].second;
			double drift = sqrt(dx*dx + dy*dy);
			sprintf(buffer11, "\tDrift length   = %.2f mm", drift);
			sprintf(buffer13, "\tDrift vector   = (%.2f, %.2f)", dx, dy);
			ss << buffer11;
			if (canPrint) {
				dx = end[i+1].first  - start[i+1].first;
				dy = end[i+1].second - start[i+1].second;
				drift = sqrt(dx*dx + dy*dy);
				sprintf(buffer12, "\tDrift length   = %.2f mm", drift);
				sprintf(buffer14, "\tDrift vector   = (%.2f, %.2f)", dx, dy);
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				ss << "| " << buffer12;
			}
			cout << ss.str() << '\n';
			ss.str("");
			ss << buffer13;
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				ss << "| " << buffer14;
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Object's apparent magnitude on the plate
			char buffer15[50], buffer16[50];
			sprintf(buffer15, "\tMagnitude      = %.2f", mag[i]);
			ss << buffer15;
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				sprintf(buffer16, "\tMagnitude      = %.2f", mag[i+1]);
				ss << "| " << buffer16;
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Plate quality grade
			ss << "\tPlate Grade    = " << p[i].grade();
			if (canPrint) {
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				ss << "| \tPlate Grade    = " << p[i+1].grade();
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Plate exposure time in minutes
			char buffer17[50], buffer18[50];
			sprintf(buffer17, "\tExposure       = %.1f mins", p[i].exposure() / 60.0);
			ss << buffer17;
			if (canPrint) {
				sprintf(buffer18, "\tExposure       = %.1f mins", p[i+1].exposure() / 60.0);
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				ss << "| " << buffer18;
			}
			cout << ss.str() << '\n';
			ss.str("");

			// Sort-of signal to noise ratio
			char buffer19[50], buffer20[50];
			double snr1 = Ephemeris::counts(p[i].exposure(), mag[i]) / p[i].countLimit();
			double snr2 = Ephemeris::counts(p[i+1].exposure(), mag[i+1]) / p[i+1].countLimit();
			sprintf(buffer19, "\tSNR            = %.2f", snr1);
			ss << buffer19;
			if (canPrint) {
				sprintf(buffer20, "\tSNR            = %.2f", snr2);
				length = SIZE-ss.str().length();
				ss << string(length, ' ');
				ss << "| " << buffer20;
			}
			cout << ss.str() << '\n';
			ss.str("");

			cout << string(SIZE*2.4, '-') << '\n';
		}
		/*for (int i = 0; i < middle.size(); i++) {
			printf("%d %.3f (%.3f %.3f) (%.3f %.3f) (%.3f %.3f)\n", p[i].id(), p[i].julian(), start[i].first, start[i].second, middle[i].first, middle[i].second, end[i].first, end[i].second);
		}*/
	}
	
	/*
		Takes the plate prefix/emulsions/filters from the plate record and calculates the 
		number of photon counts of the faintest possible object on that plate.

		Return values pulled from http://www.roe.ac.uk/ifa/wfau/ukstu/telescope.html#fe
	*/
	static double limitingCounts(const string& buffer) {
		string prefix = buffer.substr(0, 2);
		string emulsion = buffer.substr(40, 6);
		boost::algorithm::trim(prefix);
		boost::algorithm::trim(emulsion);
		double magLim, expLim;
		if (prefix == "U") {
			expLim = 180;
			magLim = 21;
		} else if (prefix == "B") {
			expLim = 60;
			magLim = 21;
		} else if (prefix == "BJ" || prefix == "J") {
			expLim = 60;
			magLim = 22.5;
		} else if (prefix == "V") {
			expLim = 60;
			magLim = 21;
		} else if (prefix == "OR") {
			expLim = 60;
			if (emulsion == "IIIaF")
				magLim = 21.5;
			else
				magLim = 22.5;
		} else if (prefix == "R") {
			expLim = 90;
			magLim = 21.5;
		} else if (prefix == "HA") {
			expLim = 180;
			magLim = 21.5;
		} else if (prefix == "I") {
			expLim = 90;
			magLim = 19.5;
		} else {
			expLim = 100;		// these default values were pulled out of the ether
			magLim = 22;
		}
		expLim *= 60;	// convert minutes to seconds
		return expLim * pow(10.0, -magLim / 2.5);
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
	static void exposureBoundaries(const Plate& p, 
	                               const pair<Coords,Coords>& coords, 
	                               const pair<double,double>& times, 
	                               pair<double, double>& start, 
	                               pair<double, double>& end) {
		double expTime   = p.exposure() / 86400.0;	// in days
		double startTime = p.julian() - (expTime/2.0);
		double endTime   = p.julian() + (expTime/2.0);
		Coords startCoords = Coords::linInterp(coords.first, times.first, coords.second, times.second, startTime);
		Coords endCoords   = Coords::linInterp(coords.first, times.first, coords.second, times.second, endTime);

		double xiStart, xiEnd, etaStart, etaEnd;
		int status1, status2;
		Coords::gnomonic(startCoords, p.m_coords, xiStart, etaStart, status1);
		Coords::gnomonic(endCoords,   p.m_coords, xiEnd,   etaEnd,   status2);
		if (status1 != 0 || status2 != 0) printf("Plate %d failed at exposureBoundaries()!\n", p.id());

		double x1 = radsToMM(xiStart);
		double y1 = radsToMM(etaStart);
		double x2 = radsToMM(xiEnd);
		double y2 = radsToMM(etaEnd);
		start = {x1, y1};
		end   = {x2, y2};
	}

	/*
		Converts the output of a gnomonic transformation to mm
			1) convert radians to arcseconds
			2) convert arcsecs to millimetres, using the given rate of 67.12"/mm
			3) add a shift of half the plate size
		This gives the output as relative to the bottom left corner of the plate
	*/
	static double radsToMM(const double rads) {
		return (rads * 3600.0 * 180.0) / (67.12 * M_PI) + (PLATE_SIZE / 2.0);
	}

	/*
		Inverse of the above
	*/
	static double mmToRads(const double mm) {
		return (mm - (PLATE_SIZE / 2.0)) * (67.12 * M_PI) / (3600.0 * 180.0);
	}

	/*
		Prints a summary of which plates were matched, but are known to be unviewable
		This is either due to the plate being known as broken/missing, or the signal-to-noise
		ratio of the plate is too low to spot the asteroid
	*/
	static void printMissingAndFaint(const vector<unsigned>& missingPlate, 
	                                 const vector<unsigned>& tooFaint, 
	                                 const int matchCount, 
	                                 const double closest) {
		int missingPlateCount = missingPlate.size(); 
		int tooFaintCount     = tooFaint.size();
		if (missingPlateCount > 0) {
			cout << missingPlateCount << (matchCount>0?" other":"") << " plate" << (missingPlateCount==1?"":"s");
			cout << " matched, but " << (missingPlateCount==1 ? "isn't" : "aren't") << " in the plate room:\n\t";
			for (int j = 0; j < missingPlateCount; j++) {
				printf("%5d ", missingPlate[j]);
				if ((j+1) % 5 == 0 && (j+1) < missingPlateCount) cout << "\n\t";
			}
			cout << endl;
		}
		if (tooFaintCount > 0) {
			cout << tooFaintCount << (matchCount>0||tooFaintCount>0?" other":"") << " plate" << (tooFaintCount==1?"":"s");
			cout << " matched, but " << (tooFaintCount==1 ? "is" : "are") << " too faint to be seen on the plate:\n\t";
			for (int j = 0; j < tooFaintCount; j++) {
				printf("%5d ", tooFaint[j]);
				if ((j+1) % 5 == 0 && (j+1) < tooFaintCount) cout << "\n\t";
			}
			cout << endl;
		}
	}

	static void printSummary(const double firstEphDate, 
	                         const double lastEphDate, 
	                         const int matchCount, 
	                         const string& objectName) {
		string firstDate = Plate::julianToGregorian(firstEphDate);
		string lastDate  = Plate::julianToGregorian(lastEphDate);
		string s = objectName;
		for (int i = 0; i < s.length(); i++) s[i] = toupper(objectName[i]);
		if (matchCount > 0) {
			cout << matchCount << " matching plate" << (matchCount>1?"s":"") << " found for " << s;
			cout << " between " << firstDate << " and " << lastDate << '\n';
		} else {
			cout << "No matches found for " << s << "!\n";
		}
	}
};

#endif