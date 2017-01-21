#ifndef EPHEMERIS_H
#define EPHEMERIS_H

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
	void setCoords(const RA& r, const DEC& d);
	void setRA(const RA& r);
	void setDEC(const DEC& d);
	void setApLST(const double lst);
	void setApMag(const float mag);
	void setErrorRA(const double dra);
	void setErrorDEC(const double ddec);
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
void Ephemeris::setCoords(const RA &r, const DEC &d) { m_coords = {r, d}; }
void Ephemeris::setRA(const RA &r) { m_coords.setRA(r); }

#endif