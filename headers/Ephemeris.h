#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "Coords.h"

using namespace std;

class Ephemeris {
private:
	double m_day;		// time in Julian days
	Coords m_coords;	// RA/DEC coordinates
	double m_lst;		// Apparent Local Sidereal Time from the observing point
	double m_mag;		// Apparent magnitude of the object in the sky
	double m_diam;		// Angular diameter of the object in the sky
	double m_dRA;		// 3sigma error in RA in arcseconds
	double m_dDEC;		// 3sigma error in DEC in arcseconds

public:
	Ephemeris(float day=0.0, Coords c={}, double lst=0.0, double mag=0.0, double diam=0.0, double dra=0.0, double ddec=0.0)
		: m_day(day), m_coords(c), m_lst(lst), m_mag(mag), m_diam(diam), m_dRA(dra), m_dDEC(ddec) { }
	Ephemeris(const Ephemeris& e)
		: m_day(e.m_day), m_coords(e.m_coords), m_lst(e.m_lst), m_mag(e.m_mag), m_diam(e.m_diam), m_dRA(e.m_dRA), m_dDEC(e.m_dDEC) { }

	inline double julian() const { return m_day; }
	inline Coords coords() const { return m_coords; }
	inline double lst()	const { return m_lst; }
	inline double mag() const { return m_mag; }
	inline double diam() const { return m_diam; }
	inline double dRA() const { return m_dRA; }
	inline double dDEC() const { return m_dDEC; }
	inline void setJulian(const double d) { m_day = d; }
	inline void setCoords(const Coords& c) { m_coords = c; }
	inline void setLST(const double lst) { m_lst = lst; }
	inline void setMag(const float mag) { m_mag = mag; }
	inline void setDiam(const double diam) { m_diam = diam; }
	inline void setdRA(const double dra) { m_dRA = dra; }
	inline void setdDEC(const double ddec) { m_dDEC = ddec; }

	/*
		Takes a full ephemeris record string and pulls the info from it, such as julian day, 
		RA/DEC coordinates (plus their errors), LST, apparent magnitude
	*/
	bool parseEphemerisString(const string& s, const bool isSurfBrt) {
		if (s.length() == 0) return false;
		string s2 = s.substr(0,18) + s.substr(21);
		stringstream ss(s2);
		string temp, dRA_str, dDEC_str, diamStr, lstStr, magStr;
		double ra, dec;
		ss >> m_day >> ra >> dec >> lstStr >> magStr; 
		if (isSurfBrt) ss >> temp;
		ss >> diamStr >> dRA_str >> dDEC_str;

		// converting the two double values (in degrees) to full RA/DEC objects
		m_coords = Coords(ra, dec);

		// checking whether some of the values are valid
		m_mag  = (magStr   == "n.a.") ? -100 : stod(magStr);
		m_lst  = (lstStr   == "n.a.") ?  0.0 : stod(lstStr);
		m_diam = (diamStr  == "n.a.") ? -1.0 : stod(diamStr);
		m_dRA  = (dRA_str  == "n.a.") ?  0.0 : stod(dRA_str);
		m_dDEC = (dDEC_str == "n.a.") ?  0.0 : stod(dDEC_str);
		return true;
	}

	/*
		Prints all the data stored in an Ephemeris object. Used for testing
	*/
	void printEphemeris() const {
		printf("JD = %.2f ", m_day);
		printf("RA = %8.4f ", m_coords.getDegRA());
		printf("DEC = %8.4f ", m_coords.getDegDEC());
		printf("LST = %7.4f ", m_lst);
		printf("ApMag = %5.2f ", m_mag);
		printf("Diam = %.5f ", m_diam);
		printf("dRA = %.4f ", m_dRA);
		printf("dDEC = %.4f\n", m_dDEC);
	}

	/*
		Goes through the ephemeris file to pick out all relevant data, then stores them all in 
		a vector of Ephemeris objects. The $$SOE and $$EOE tags signify the start and end of the 
		data lines.
	*/
	static void readEphemerisFile(vector<Ephemeris>& eph, const string& filename) {
		ifstream ephemerisFile(filename);
		if (ephemerisFile.is_open()) {
			bool canReadEntries = false;
			bool isSurfBrt = false;
			while (!ephemerisFile.eof()) {
				string buffer;
				getline(ephemerisFile, buffer);
				// if the column headers include "S-brt", pass that flag to the parsing function
				if (!canReadEntries && buffer.find("S-brt") != string::npos)
					isSurfBrt = true;
				if (buffer == "$$EOE") 
					break;
				if (canReadEntries) {
					Ephemeris e;
					if (e.parseEphemerisString(buffer, isSurfBrt)) {
						eph.push_back(e);
					}
				}
				if (buffer == "$$SOE") 
					canReadEntries = true;
			}
			ephemerisFile.close();
		}
		else {
			cout << "Ephemeris file \"" << filename << "\" is not valid\n";
			exit(1);
		}
	}

	/*
		Linearly interpolates a floating point number between two points
		Intended for apparent magnitude estimation, but can be used for anything else
	*/
	static double linInterp(const double mag0, const double t0, const double mag1, const double t1, const double t) {
		return mag0 + (t-t0)*(mag1-mag0)/(t1-t0);
	}

	/*
		Takes an array of ephemeris objects, a number of records to return and an array index to surround.
		Returns an array of numCoords x Coords objects surrounding index i, and the corresponding julian dates.

		Used to find x/y values for polynomial least-squares fitting.
	*/
	static void findNearbyEphs(const vector<Ephemeris>& e, const int numCoords, const int i,
								vector<Coords>& c, vector<double>& t) {
		size_t N = e.size();
		if (i > N) {
			cout << "Error in Ephemeris::findNearbyCoords(): index " << i << " > vector size " << N << ".\n";
			exit(1);
		}
		int shift = (numCoords % 2 == 1) ? numCoords/2 : (numCoords-1)/2;
		int diff = (i-shift < 0) ? shift-i : 0;
		for (int j = 0; j < numCoords; j++) {
			c[j] = e[i-shift+j+diff].coords();
			t[j] = e[i-shift+j+diff].julian();
		}
	}

};

#endif