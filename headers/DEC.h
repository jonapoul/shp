#ifndef DEC_H	
#define DEC_H

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <sstream>
using namespace std;
class RA;

class DEC {
private:
	bool m_isPositive;
	int m_deg;
	int m_min;
	float m_sec;
	double m_degrees;
	double m_radians;

public:
	DEC(const bool isPos = true, const int d = 0, const int m = 0, const float s = 0.f);
	DEC(const DEC& dec);
	DEC(double decimal);

	bool 	isPositive() const { return m_isPositive; };
	int  	getDeg() 	 const { return m_deg; };
	int  	getMins() 	 const { return m_min; };
	float 	getSecs() 	 const { return m_sec; };
	double 	getDegrees() const { return m_degrees; };
	double 	getRadians() const { return m_radians; };
	
	void setIsPositive(const bool isPos) { m_isPositive = isPos; };
	void setDeg 	  (const int d) 	 { m_deg = d; };
	void setMins	  (const int m) 	 { m_min = m; };
	void setSecs	  (const float s) 	 { m_sec = s; };
	void setDegrees   (const double d) 	 { m_degrees = d; m_radians = d*M_PI/180; };
	void setRadians   (const double r) 	 { m_radians = r; m_degrees = r*180/M_PI; };

	void   fixDEC();
	string toString()  const;
	float  toDecimal() const;
	float  toRadians() const;

	friend ostream& operator<<(ostream& os,  const DEC& ra);
	friend DEC 		operator+ (const DEC& a, const DEC& b);
	friend DEC 		operator- (const DEC& a, const DEC& b);
	friend bool 	operator==(const DEC& a, const DEC& b);
	friend bool 	operator!=(const DEC& a, const DEC& b);
	friend bool 	operator< (const DEC& a, const DEC& b);
	friend bool 	operator> (const DEC& a, const DEC& b);
	friend bool 	operator<=(const DEC& a, const DEC& b);
	friend bool 	operator>=(const DEC& a, const DEC& b);
};

DEC::DEC(const bool isPos, const int d, const int m, const float s) 
	: m_isPositive(isPos), m_deg(d), m_min(m), m_sec(s) {
	if (d < 0) {
		m_isPositive = false;
		m_deg *= -1;
	}
	m_degrees = toDecimal();
	m_radians = m_degrees * M_PI / 180.f;
}

DEC::DEC(const DEC& dec) 
	: m_isPositive(dec.m_isPositive), m_deg(dec.m_deg), m_min(dec.m_min), m_sec(dec.m_sec) {
	m_degrees = toDecimal();
	m_radians = m_degrees * M_PI / 180.f;
}

DEC::DEC(double decimal) {
	while (decimal < -360.f) decimal += 360.f;
	while (decimal > 360.f) decimal -= 360.f;

	m_degrees = decimal;
	m_radians = m_degrees * M_PI / 180.f;

	if (decimal < 0.f) {
		m_isPositive = false;
		decimal *= -1;
	}
	else 
		m_isPositive = true;

	m_deg = int(decimal);
	decimal = (decimal - m_deg) * 60.f;
	m_min = int(decimal);
	decimal = (decimal - m_min) * 60.f;
	m_sec = decimal;
}

void DEC::fixDEC() {
	while (m_sec < 0.f) 	{ m_sec += 60.f; 	m_min -= 1; }
	while (m_sec > 60.f)	{ m_sec -= 60.f; 	m_min += 1; }
	while (m_min < 0)		{ m_min += 60; 		m_deg -= 1; }
	while (m_min > 60)		{ m_min -= 60; 		m_deg += 1; }
	while (m_deg < -90) 	{ m_deg += 90; }
	while (m_deg > 90)		{ m_deg -= 90; }
}

string DEC::toString() const {
	string output;
	if (m_isPositive) output += '+';
	else			  output += '-';
	if (m_deg < 10) output += '0';
	output += to_string(m_deg) + '\370';
	if (m_min < 10) output += '0';
	output += to_string(m_min) + '\'';
	if (m_sec < 10) output += '0';
	output += to_string(int(m_sec)) + '\"';
	return output;
}

float DEC::toDecimal() const {
	return (m_isPositive ? 1 : -1) * (m_deg + m_min*(1.f/60.f) + m_sec*(1.f/3600.f));
}

float DEC::toRadians() const {
	float pi = 3.1415926535;
	return toDecimal() * (pi/180.f) * (m_isPositive ? 1 : -1);
}

ostream& operator<<(ostream& os, const DEC& dec) {
	os << dec.toString();
    return os;
}

bool operator==(const DEC& a, const DEC& b) { 
	return abs(a.m_degrees - b.m_degrees) < 1e-6;
}

bool operator!=(const DEC& a, const DEC& b) { 
	return !(a == b); 
}

bool operator<(const DEC& a, const DEC& b) { 
	return (a.m_degrees < b.m_degrees); 
}

bool operator>(const DEC& a, const DEC& b) { 
	return (abs(a.m_degrees - b.m_degrees) < 1e-6) ? false : !(a < b); 
}

bool operator<=(const DEC& a, const DEC& b) { 
	return (a.m_degrees <= b.m_degrees); 
}

bool operator>=(const DEC& a, const DEC& b) { 
	return (abs(a.m_degrees - b.m_degrees) < 1e-6) ? false : !(a < b); 
}

#endif