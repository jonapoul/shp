#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include "Coords.h"

using namespace std;

class Ephemeris {
private:
	double m_day;		// time in Julian days
	Coords m_coords;	// RA/DEC coordinates
	double m_lst;		// Apparent Local Sidereal Time from the observing point
	float m_mag;		// Apparent magnitude of the object in the sky
	double m_dRA;		// 3sigma error in RA in arcseconds
	double m_dDEC;		// 3sigma error in DEC in arcseconds

public:
	Ephemeris(float day=0.f, Coords c={}, double lst=0.f, float mag=0.f, double dra=0.f, double ddec=0.f)
		: m_day(day), m_coords(c), m_lst(lst), m_mag(mag), m_dRA(dra), m_dDEC(ddec) { }
	Ephemeris(const Ephemeris& e)
		: m_day(e.m_day), m_coords(e.m_coords), m_lst(e.m_lst), m_mag(e.m_mag), m_dRA(e.m_dRA), m_dDEC(e.m_dDEC) { }

	double 	getJulian()   const { return m_day; };
	Coords 	getCoords()   const { return m_coords; };
	RA 		getRA() 	  const { return m_coords.getRA(); };
	DEC 	getDEC() 	  const { return m_coords.getDEC(); };
	double 	getLST()	  const { return m_lst; };
	float 	getApMag()	  const { return m_mag; };
	double 	getErrorRA()  const { return m_dRA; };
	double 	getErrorDEC() const { return m_dDEC; };

	void setJulian  (const double d) 	{ m_day = d; };
	void setCoords  (const Coords& c) 	{ m_coords = c; };
	void setRA      (const RA& r) 		{ m_coords.setRA(r); };
	void setDEC     (const DEC& d) 	 	{ m_coords.setDEC(d); };
	void setApLST   (const double lst)  { m_lst = lst; };
	void setApMag   (const float mag) 	{ m_mag = mag; };
	void setErrorRA (const double dra)  { m_dRA = dra; };
	void setErrorDEC(const double ddec) { m_dDEC = ddec; };

	bool parseEphemerisString(const string& s);
	void printEphemeris() const;
	static void readEphemerisFile(vector<Ephemeris>& eph, const string& filename);
};


bool Ephemeris::parseEphemerisString(const string& s) {
	if (s.length() == 0) return false;

	stringstream ss(s);
	string temp, dRA_str, dDEC_str;
	double ra, dec;
	ss >> m_day >> temp >> ra >> dec >> m_lst >> m_mag >> temp >> dRA_str >> dDEC_str;

	// converting the two double values (in degrees) to full RA/DEC objects
	m_coords = { ra, dec };

	// checking whether the errors in RA/DEC are valid
	m_dRA  = (dRA_str  == "n.a.") ? 0.f : stod(dRA_str);
	m_dDEC = (dDEC_str == "n.a.") ? 0.f : stod(dDEC_str);
	return true;
}

void Ephemeris::printEphemeris() const {
	printf("JD = %.2f ", m_day);
	printf("RA = %.3f ", m_coords.getRA().getDegrees());
	printf("DEC = %.3f ", m_coords.getDEC().getDegrees());
	printf("LST = %.4f ", m_lst);
	printf("ApMag = %.2f ", m_mag);
	printf("dRA = %.2f ", m_dRA);
	printf("dDEC = %.2f\n", m_dDEC);
}

void Ephemeris::readEphemerisFile(vector<Ephemeris>& ephemerides, const string& filename) {
	ifstream ephemerisFile("mars.txt");
	if (ephemerisFile.is_open()) {
		bool canReadEntries = false;
		while (!ephemerisFile.eof()) {
			string buffer;
			getline(ephemerisFile, buffer);
			if (buffer == "$$EOE") break;
			if (canReadEntries) {
				Ephemeris e;
				if (e.parseEphemerisString(buffer))
					ephemerides.push_back(e);
			}
			if (buffer == "$$SOE") canReadEntries = true;
		}
		ephemerisFile.close();
	}
	else {
		cout << "Ephemeris file \"" << filename << "\" is not valid\n";
	}
}

#endif