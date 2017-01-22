#ifndef COORDS_H
#define COORDS_H

#include <iostream>
#include <string>
#include <math.h>
#include "RA.h"
#include "DEC.h"
using namespace std;

class Coords {
private:
	RA  m_ra;
	DEC m_dec;

public:
	Coords(const RA& r = {}, const DEC& d = {});
	Coords(const Coords& c);

	RA getRA() const;
	DEC getDEC() const;

	void setRA(const RA& r);
	void setDEC(const DEC& d);

	void parseFromPlateRecord(const string& record);
	static float angularDistance(const Coords& a, const Coords& b, const bool returnValueInDegrees = true);

	friend ostream& operator<<(ostream& os, const Coords& c);
	friend Coords operator+ (const Coords& a, const Coords& b);
	friend Coords operator- (const Coords& a, const Coords& b);
};

Coords::Coords(const RA& r, const DEC& d) 
	: m_ra(r), m_dec(d) { }
Coords::Coords(const Coords& c) 
	: m_ra(c.m_ra), m_dec(c.m_dec) { }

RA Coords::getRA() const { return m_ra; }
DEC Coords::getDEC() const { return m_dec; }

void Coords::setRA(const RA& r) { m_ra = r; }
void Coords::setDEC(const DEC &d) { m_dec = d; }

void Coords::parseFromPlateRecord(const string& record) {
	if (record.length() < 31) {
		cerr << "string passed to Coords::parseFromPlateRecord() is too short\n";
		return;
	}
	m_ra = {
		stoi(record.substr(20,2)),
		stoi(record.substr(22,2)),
		stoi(record.substr(24,1)) * 6.f
	};
	m_dec = {
		!(record[25] == '-'),		// if the sign column is '-', set isPositive to true
		stoi(record.substr(26,2)),
		stoi(record.substr(28,2)),
		0.f
	};
}

float Coords::angularDistance(const Coords& a, const Coords& b, const bool returnValueInDegrees) {
	float a_ra  = a.m_ra.getRadians();
	float a_dec = a.m_dec.getRadians() * (a.m_dec.isPositive() ? 1 : -1);
	float b_ra  = b.m_ra.getRadians();
	float b_dec = b.m_dec.getRadians() * (b.m_dec.isPositive() ? 1 : -1);;

	float cosx   = sin(a_dec)*sin(b_dec) + cos(a_dec)*cos(b_dec)*cos(a_ra-b_ra);
	if (returnValueInDegrees)	return acos(cosx) * 180.f / M_PI;
	else						return acos(cosx);
}

ostream& operator<<(ostream& os, const Coords& c) {
	os << c.m_ra.toString() << ", " << c.m_dec.toString();
    return os;
}

Coords operator+(const Coords&a, const Coords& b) {
	return {a.m_ra+b.m_ra, a.m_dec+b.m_dec};
}

Coords operator-(const Coords& a, const Coords& b) {
	return {a.m_ra-b.m_ra, a.m_dec-b.m_dec};
}

#endif