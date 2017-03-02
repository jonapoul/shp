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
#include "Parameters.h"

using std::cout;

class Plate {
private:
	int id_;						// plate identification number
	Coords coords_;			// combines RA and DEC to one object
	std::string greg_;	// UT gregorian date string, formatted as yymmdd
	double julian_;			// julian date in the middle of exposure
	double exp_;				// exposure time in secs
	char grade_;				// graded plate quality, A being best and C being worst
	double countLim_;		// lowest number of photon counts on a plate below which an object isn't visible
	int filterWL_;			// peak wavelength in angstroms

public:
	Plate(const int num = 0, 
	      const Coords& c = {}, 
	      const std::string& g = "", 
	      const double j = 0.0, 
	      const double exp = 0.0, 
	      const char grade = ' ', 
	      const double lim = 0.0,
	      const int wl = 0)
			: id_(num), coords_(c), greg_(g), julian_(j), 
				exp_(exp), grade_(grade), countLim_(lim), filterWL_(wl) { }
	Plate(const Plate& p)
			: id_(p.id_), coords_(p.coords_), greg_(p.greg_), julian_(p.julian_), 
				exp_(p.exp_), grade_(p.grade_), countLim_(p.countLim_), filterWL_(p.filterWL_) { }

	inline int id() const { return id_; };
	inline Coords coords() const { return coords_; };
	inline std::string gregorian() const { return greg_; }
	inline double julian() const { return julian_; };
	inline double exposure() const { return exp_; }
	inline char grade() const { return grade_; };
	inline double countLimit() const { return countLim_; }
	inline int wavelength() const { return filterWL_; }

	/*
		Reads a line from the plate catalog file and fills in the relevant fields
		of the Plate object.

		buffer = full string representing all of one plate record
		Outputs a boolean flag to indicate whether the plate is valid for our purposes
	*/
	bool parsePlateString(const std::string& buffer) {
		// T = tracked shot, M = multiple shots, P = full-aperture prism (distorted)
		if (buffer[7] == 'T' || buffer[7] == 'M' || buffer[7] == 'P') 			return false;
		if (buffer.substr(16,4) == "TEST" || buffer.substr(15,4) == "TEST") return false;

		try { id_ = stoi(buffer.substr(2, 5)); }
		catch (...) { return false; }

		try { coords_.parseCoordsFromPlate(buffer); }
		catch (...) { return false; }

		greg_ = buffer.substr(30, 6);
		std::string lst  = buffer.substr(36, 4);
		for (auto& c : lst)
			if (!isdigit(c)) c = '0';

		bool isDebug = (PRINT_DATE_BREAKDOWN && id_ == PLATE_TO_PRINT);
		julian_ = convertDate(greg_, lst, isDebug);

		try { exp_ = stod(buffer.substr(52,4)) * 6.0; }
		catch (...) { return false; }
		julian_ += (exp_ / 86400.0) / 2.0;

		grade_ = buffer[56];
		countLim_ = limitingCounts(buffer);
		
		std::string prefix = buffer.substr(0, 2);
		std::string filter = buffer.substr(46, 6);
		boost::algorithm::trim(filter);
		filterWL_ = findWavelength(prefix, filter);

		return true;
	}

	/*
		Checks whether the plate is missing from the plate archive room
	*/
	bool isMissing(const std::vector<int>& missingList) {
		for (const int id : missingList) {
			if (this->id() == id) return true;
		}
		return false;
	}

	/* 
		Takes in a 6-char date string, in the format "yymmdd", outputs a Julian 
		date as a double. This assumes that all dates are between 1st Jan 1917 
		and 31st Dec 2016, which is fine for these plate catalogs since they 
		don't go beyond ~2003

		From "Practical Astronomy with your Calculator or Spreadsheet"
	*/
	static double convertDate(const std::string& gregorianDate, 
	                          const std::string& lst,
	                          const bool isDebug = false) {
		int year = stoi(gregorianDate.substr(0, 2));
		if (year < 100)
			year += (year < 17) ? 2000 : 1900;
		int month = stoi(gregorianDate.substr(2, 2));
		int day   = stoi(gregorianDate.substr(4, 2));
		int hour  = stoi(lst.substr(0, 2));
		int mins  = stoi(lst.substr(2, 2));

		double julian = gregorianToJulian(double(day), month, year);
		double gst = LSTtoGST(hour, mins, 0.0, SITE_LONGITUDE);
		double ut = GSTtoUT(gst, julian);
		double frac = ut / 24.0;
		julian += frac;

		if (isDebug) {		// date conversion debugging
			printf("Date = %s\n", gregorianDate.c_str());
			printf("LST  = %s\n", lst.c_str());
			printf("JD   = %.5f\n", julian-frac);
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
	static std::string julianToGregorian(double jd) {
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
		
		std::string output;
		if (d < 10) output += '0';
		output += std::to_string(d) + '/';
		if (m < 10) output += '0';
		output += std::to_string(m) + '/';
		output += std::to_string(y);
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
	static std::string julianToUT(const double JD) {
		double zeroHour = floor(JD - 0.5) + 0.5;
		double remainder = (JD - zeroHour) * 24;

		int h    = int(remainder);
		remainder = (remainder - h) * 60.0;
		int m    = int(remainder);
		double s = (remainder - m) * 60.0;
		std::stringstream ss;
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
	static void readPlateCatalog(std::vector<Plate>& plates, 
	                             const std::string& filename) {
		std::ifstream platesFile(filename);
		if (platesFile.is_open()) {
			while (!platesFile.eof()) {
				std::string buffer;
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
		Takes the plate prefix/emulsions/filters from the plate record and calculates the 
		number of photon counts of the faintest possible object on that plate.

		Return values pulled from http://www.roe.ac.uk/ifa/wfau/ukstu/telescope.html#fe
	*/
	static double limitingCounts(const std::string& buffer) {
		std::string prefix = buffer.substr(0, 2);
		std::string emulsion = buffer.substr(40, 6);
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
	static std::string gregorianToString(const std::string& greg) {
		int year  = stoi(greg.substr(0,2));
		int month = stoi(greg.substr(2,2));
		int day   = stoi(greg.substr(4,2));
		year += (year < 17) ? 2000 : 1900;
		std::string date = "";
		if (day < 10) date += '0';
		date += std::to_string(day) + '/';
		if (month < 10) date += '0';
		date += std::to_string(month) + '/';
		date += std::to_string(year);
		return date;
	}

	/*
		Calculates the xi/eta coordinates of the object at the start and end of the exposure time
		This is to account for the potential dragging across the image for longer exposure times.
		Returns the values as an std::pair<double> object, with xi as first and eta in second
	*/
	static void exposureBoundaries(const Plate& p, 
	                               const std::pair<Coords,Coords>& coords, 
	                               const std::pair<double,double>& times, 
	                               std::pair<double, double>& start, 
	                               std::pair<double, double>& end) {
		double expTime   = p.exposure() / 86400.0;	// in days
		double startTime = p.julian() - (expTime/2.0);
		double endTime   = p.julian() + (expTime/2.0);
		Coords startCoords = Coords::linInterp(coords.first, times.first, coords.second, times.second, startTime);
		Coords endCoords   = Coords::linInterp(coords.first, times.first, coords.second, times.second, endTime);

		double xiStart, xiEnd, etaStart, etaEnd;
		int status1, status2;
		Coords::gnomonic(startCoords, p.coords_, xiStart, etaStart, status1);
		Coords::gnomonic(endCoords,   p.coords_, xiEnd,   etaEnd,   status2);
		if (status1 != 0 || status2 != 0) printf("Plate %d failed at exposureBoundaries()!\n", p.id());

		double x1 = Coords::radsToMM(xiStart);
		double y1 = Coords::radsToMM(etaStart);
		double x2 = Coords::radsToMM(xiEnd);
		double y2 = Coords::radsToMM(etaEnd);
		start = {x1, y1};
		end   = {x2, y2};
	}

	/*
		Prints a summary of which plates were matched, but are known to be unviewable
		This is either due to the plate being known as broken/missing, or the signal-to-noise
		ratio of the plate is too low to spot the asteroid
	*/
	static void printMissingAndFaint(const std::vector<unsigned>& missingPlate, 
	                                 const std::vector<unsigned>& tooFaint, 
	                                 const int matchCount, 
	                                 const double closest, 
	                                 const bool filterSNR) {
		int missingPlateCount = missingPlate.size(); 
		int tooFaintCount     = tooFaint.size();
		if (missingPlateCount > 0) {
			cout << missingPlateCount << (matchCount>0?" other":"") << " plate" << (missingPlateCount==1?"":"s");
			cout << " matched, but " << (missingPlateCount==1 ? "isn't" : "aren't") << " in the plate room:\n\t";
			for (int j = 0; j < missingPlateCount; j++) {
				printf("%5d ", missingPlate[j]);
				if ((j+1) % 5 == 0 && (j+1) < missingPlateCount) cout << "\n\t";
			}
			cout << '\n';
		}
		if (tooFaintCount > 0) {
			cout << tooFaintCount << (matchCount>0||tooFaintCount>0?" other":"") << " plate" << (tooFaintCount==1?"":"s");
			cout << " matched, but " << (tooFaintCount==1 ? "has" : "have");
			cout << " a Signal-to-Noise Ratio below " << SNR_LIMIT << ":\n\t";
			for (int j = 0; j < tooFaintCount; j++) {
				printf("%5d ", tooFaint[j]);
				if ((j+1) % 7 == 0 && (j+1) < tooFaintCount) cout << "\n\t";
			}
			if (filterSNR) {
				cout << "\nAdd a -snr flag to the command to ignore SNR filtering\n";
			}
		} 
		if (!filterSNR) {
			cout << "SNR filtering has been ignored\n";
		}
		if (matchCount == 0 && tooFaintCount == 0 && missingPlateCount == 0) {
			cout << "Closest angular distance to a plate centre was " << closest << "°\n";
		}
	}

	static void printSummary(const double firstEphDate, 
	                         const double lastEphDate, 
	                         const int matchCount, 
	                         const std::string& objectName) {
		std::string firstDate = Plate::julianToGregorian(firstEphDate);
		std::string lastDate  = Plate::julianToGregorian(lastEphDate);
		std::string s = objectName;
		for (int i = 0; i < s.length(); i++) s[i] = toupper(objectName[i]);
		if (matchCount > 0) {
			cout << matchCount << " matching plate" << (matchCount>1?"s":"") << " found for " << s;
			cout << " between " << firstDate << " and " << lastDate << '\n';
			cout << "The shown dates/times represent the middle of the exposure, not the start\n";
			cout << "Uncertainties represent an ellipse up to 3 sigma, i.e. 99.7% certainty\n";
		} else {
			cout << "No matches found for " << s << "!\n";
		}
	}

	/*
		Returns the filter wavelength based on the plate prefix and filter code strings
		
		Values taken from:
		http://www.roe.ac.uk/ifa/wfau/ukstu/pltcat.html#filter
			and 
		http://www.roe.ac.uk/ifa/wfau/ukstu/pltcat.html#prefix
	*/
	static int findWavelength(const std::string& prefix, 
	                          const std::string& filter) {

		// special filters
		std::string filterPrefix = filter.substr(0, 3);
		if (filterPrefix == "AAO") {
			if (filter == "AAO372") return 3727;
			if (filter == "AAO460") return 4600;
			if (filter == "AAO468") return 4686;
			if (filter == "AAO486") return 4860;
			if (filter == "AAO500") return 5007;
			if (filter == "AAO540") return 5400;
			if (filter == "AAO562") return 5620;
			if (filter == "AAO587") return 5870;
			if (filter == "AAO630") return 6300;
			if (filter == "AAO643") return 6430;
			if (filter == "AAO672") return 6725;
			if (filter == "AAO656") return 6567;
		} else if (filterPrefix == "MB ") {
			if (filter == "MB 486") return 4861;
			if (filter == "MB 675") return 6751;
			if (filter == "MB 569") return 5696;
			if (filter == "MB 666") return 6668;
			if (filter == "MB 657") return 6567;
		} else if (filterPrefix == "HZ ") {
			if (filter == "HZ 580") return 5800;
			if (filter == "HZ 600") return 6000;
			if (filter == "HZ 620") return 6200;
			if (filter == "HZ 640") return 6400;
			if (filter == "HZ 660") return 6600;
		}
		if (filter == "RG630H") return 6567;

		// singular prefixes
		if (prefix == " U") return 3550;
		if (prefix == " B") return 3850;
		if (prefix == " J") return 3950;
		if (prefix == " V") return 4950;
		if (prefix == " R") return 6300;
		if (prefix == " I") return 7150;
		if (prefix == " Z") return 10000;

		// two letter prefixes
		if (prefix[0] == 'U') return 3050;
		if (prefix[0] == 'B') return 3850;
		if (prefix[0] == 'J') return 3950;
		if (prefix[0] == 'Y') return 4550;
		if (prefix[0] == 'V') return 4950;
		if (prefix[0] == 'O') return 5900;
		if (prefix[0] == 'R') return 6300;
		if (prefix[0] == 'I') return 7150;
		if (prefix[0] == 'W') return 8300;
		if (prefix[0] == 'Z') return 10000;
		
		// default, assumes visible wavelength
		return 5300;
	}
};

#endif