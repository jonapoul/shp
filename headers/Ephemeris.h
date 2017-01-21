#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include "Coords.h"

class Ephemeris {
private:
	float m_day;		// time in Julian days
	Coords m_coords;	// RA/DEC coordinates
	double m_lst;		// Apparent Local Sidereal Time from the observing point
	float m_mag;		// Apparent magnitude of the object in the sky
	double m_dRA;		// 3sigma error in RA in arcseconds
	double m_dDEC;		// 3sigma error in DEC in arcseconds

public:
	Ephemeris(float day=0.f, Coords c={}, double lst=0.f, float mag=0.f, double dra=0.f, double ddec=0.f);
	Ephemeris(const Ephemeris& e);
	Ephemeris(const std::string& buffer);

	float getJulianDay() const;
	Coords getCoords() const;
	RA getRA() const;
	DEC getDEC() const;
	double getLST() const;
	float getApMag() const;
	double getErrorRA() const;
	double getErrorDEC() const;

	void setJulianDay(const float d);
	void setCoords(const Coords& c);
	void setRA(const RA& r);
	void setDEC(const DEC& d);
	void setApLST(const double lst);
	void setApMag(const float mag);
	void setErrorRA(const double dra);
	void setErrorDEC(const double ddec);

	bool parseEphemerisString(const std::string& s);
	static void readEphemerisFile(std::vector<Ephemeris>& eph, const std::string& filename);
	void printEphemeris() const;
};

Ephemeris::Ephemeris(float day, Coords c, double lst, float mag, double dra, double ddec)
	: m_day(day), m_coords(c), m_lst(lst), m_mag(mag), m_dRA(dra), m_dDEC(ddec) { }
Ephemeris::Ephemeris(const Ephemeris& e)
	: m_day(e.m_day), m_coords(e.m_coords), m_lst(e.m_lst), m_mag(e.m_mag), m_dRA(e.m_dRA), m_dDEC(e.m_dDEC) { }
Ephemeris::Ephemeris(const std::string& buffer) {

}

float Ephemeris::getJulianDay() const { return m_day; }
Coords Ephemeris::getCoords() const { return m_coords; }
RA Ephemeris::getRA() const { return m_coords.getRA(); }
DEC Ephemeris::getDEC() const { return m_coords.getDEC(); }
double Ephemeris::getLST() const { return m_lst; }
float Ephemeris::getApMag() const { return m_mag; }
double Ephemeris::getErrorRA() const { return m_dRA; }
double Ephemeris::getErrorDEC() const { return m_dDEC; }

void Ephemeris::setJulianDay(const float d) { m_day = d; }
void Ephemeris::setCoords(const Coords &c) { m_coords = c; }
void Ephemeris::setRA(const RA &r) { m_coords.setRA(r); }
void Ephemeris::setDEC(const DEC &d) { m_coords.setDEC(d); }
void Ephemeris::setApLST(const double lst) { m_lst = lst; }
void Ephemeris::setApMag(const float mag) { m_mag = mag; }
void Ephemeris::setErrorRA(const double dra) { m_dRA = dra; }
void Ephemeris::setErrorDEC(const double ddec) { m_dDEC = ddec; }

bool Ephemeris::parseEphemerisString(const std::string& s) {
	if (s.length() == 0) return false;

	std::stringstream ss(s);
	std::string temp, dRA_str, dDEC_str;
	double ra, dec;

	// the two "temp" insertions represent irrelevant info
	// first = whether it's daytime or not
	// second = surface brightness
	ss >> m_day >> temp >> ra >> dec >> m_lst >> m_mag >> temp >> dRA_str >> dDEC_str;

	// converting the two double values (in degrees) to sexagesimal RA/DEC values
	m_coords = { {ra}, {dec} };

	// checking whether the errors in RA/DEC are valid
	m_dRA  = (dRA_str  == "n.a.") ? 0.f : std::stod(dRA_str);
	m_dDEC = (dDEC_str == "n.a.") ? 0.f : std::stod(dDEC_str);
	return true;
}

void Ephemeris::readEphemerisFile(std::vector<Ephemeris>& eph, const std::string& filename) {
	std::ifstream ephemerisFile("mars.txt");
	if (ephemerisFile.is_open()) {
		bool canReadEntries = false;
		while (!ephemerisFile.eof()) {
			std::string buffer;
			getline(ephemerisFile, buffer);
			if (buffer == "$$EOE") canReadEntries = false;
			if (canReadEntries) {
				Ephemeris e;
				if (e.parseEphemerisString(buffer))
					eph.push_back(e);
			}
			if (buffer == "$$SOE") canReadEntries = true;
		}
		ephemerisFile.close();
	}
}

void Ephemeris::printEphemeris() const {
	printf("Day = %.1f RA = %.3f DEC = %.3f ", m_day, m_coords.getRA().getDegrees(), m_coords.getDEC().getDegrees());
	printf("LST = %.2f ApMag = %.2f dRA = %.1f dDEC = %.1f\n", m_lst, m_mag, m_dRA, m_dDEC);
}

#endif