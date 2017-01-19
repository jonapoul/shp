#ifndef COORDS_H
#define COORDS_H

#include <iostream>
#include <string>
#include <math.h>
#include "RA.h"
#include "DEC.h"
using std::cout;

class Coords {
public:
	RA  m_ra;
	DEC m_dec;

	Coords(const RA& r = {}, const DEC& d = {});
	Coords(const Coords& c);
	float toDecimal();

	static float angularDistance(Coords& a, Coords& b, bool returnValueInDegrees = true);

	friend std::ostream& operator<<(std::ostream& os, Coords& c);
	friend Coords operator+ (Coords& a, Coords& b);
	friend Coords operator- (Coords& a, Coords& b);
};

Coords::Coords(const RA& r, const DEC& d) : m_ra(r), m_dec(d) { }

Coords::Coords(const Coords& c) : m_ra(c.m_ra), m_dec(c.m_dec) { }

float Coords::toDecimal() {
	return sqrt( m_ra.toDecimal()*m_ra.toDecimal() + m_dec.toDecimal()*m_dec.toDecimal() );
}

float Coords::angularDistance(Coords& a, Coords& b, bool returnValueInDegrees) {
	float a_ra   = a.m_ra.toRadians();
	float a_dec  = a.m_dec.toRadians() * (a.m_dec.m_isPositive ? 1 : -1);
	float b_ra   = b.m_ra.toRadians();
	float b_dec  = b.m_dec.toRadians() * (b.m_dec.m_isPositive ? 1 : -1);;

	float cosx   = sin(a_dec)*sin(b_dec) + cos(a_dec)*cos(b_dec)*cos(a_ra-b_ra);
	if (returnValueInDegrees)	return acos(cosx) * 180.f / M_PI;
	else						return acos(cosx);
}

std::ostream& operator<<(std::ostream& os, Coords& c) {
	os << c.m_ra.toString() << ", " << c.m_dec.toString();
    return os;
}

Coords operator+(Coords&a, Coords& b) {
	RA ra;
	DEC dec;
	ra 	= a.m_ra  + b.m_ra;
	dec = a.m_dec + b.m_dec;
	return {ra, dec};
}

Coords operator-(Coords& a, Coords& b) {
	RA ra;
	DEC dec;
	ra 	= a.m_ra  - b.m_ra;
	dec = a.m_dec - b.m_dec;
	return {ra, dec};
}

#endif