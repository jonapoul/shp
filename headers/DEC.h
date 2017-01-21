#ifndef DEC_H	
#define DEC_H

#include <string>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <sstream>
using std::cout;

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
	DEC(bool isPos = true, int d = 0, int m = 0, float s = 0.f);
	DEC(const DEC& dec);
	DEC(double decimal);

	bool isPositive() const;
	int getDeg() const;
	int getMins() const;
	float getSecs() const;
	double getDegrees() const;
	double getRadians() const;
	
	void setIsPositive(const bool isPos);
	void setDeg(const int d);
	void setMins(const int m);
	void setSecs(const float s);
	void setDegrees(const double d);
	void setRadians(const double r);

	std::string toString() const;
	void fixDEC();
	float toDecimal() const;
	float toRadians() const;

	friend std::ostream& operator<<(std::ostream& os, const RA& ra);
	friend DEC operator+ (const DEC& a, const DEC& b);
	friend DEC operator- (const DEC& a, const DEC& b);

};

DEC::DEC(bool isPos, int d, int m, float s) 
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

bool DEC::isPositive() const { return m_isPositive; }
int DEC::getDeg() const { return m_deg; }
int DEC::getMins() const { return m_min; }
float DEC::getSecs() const { return m_sec; }
double DEC::getDegrees() const { return m_degrees; }
double DEC::getRadians() const { return m_radians; }

void DEC::setIsPositive(const bool isPos) { m_isPositive = isPos; }
void DEC::setDeg(const int d) { m_deg = d; }
void DEC::setMins(const int m) { m_min = m; }
void DEC::setSecs(const float s) { m_sec = s; }
void DEC::setDegrees(const double d) { m_degrees = d; }
void DEC::setRadians(const double r) { m_radians = r; };

std::string DEC::toString() const {
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2);
	if (m_isPositive)
		ss << '+';
	else 
		ss << '-';
	if (m_deg < 10 && m_deg > -10) 	
		ss << '0';
	ss << abs(m_deg) << '\370';
	if (m_min < 10)
		ss << '0';
	ss << m_min << '\'';
	if (m_sec < 10)
		ss << '0';
	ss << int(m_sec) << '\"';
	return ss.str();
}

void DEC::fixDEC() {
	while (m_sec < 0.f) 	{ m_sec += 60.f; 	m_min -= 1; }
	while (m_sec > 60.f)	{ m_sec -= 60.f; 	m_min += 1; }
	while (m_min < 0)		{ m_min += 60; 		m_deg -= 1; }
	while (m_min > 60)		{ m_min -= 60; 		m_deg += 1; }
	while (m_deg < -90) 	{ m_deg += 90; }
	while (m_deg > 90)		{ m_deg -= 90; }
}

float DEC::toDecimal() const {
	return (m_isPositive ? 1 : -1) * (m_deg + m_min*(1.f/60.f) + m_sec*(1.f/3600.f));
}

float DEC::toRadians() const {
	float pi = 3.1415926535;
	return toDecimal() * (pi/180.f) * (m_isPositive ? 1 : -1);
}

std::ostream& operator<<(std::ostream& os, const DEC& dec) {
	os << dec.toString();
    return os;
}

DEC operator+(const DEC& a, const DEC& b) {
	DEC output(a.m_deg+b.m_deg, a.m_min+b.m_min, a.m_sec+b.m_sec); 
	output.fixDEC();
	return output;
}

DEC operator-(const DEC& a, const DEC& b) {
	DEC output(a.m_deg-b.m_deg, a.m_min-b.m_min, a.m_sec-b.m_sec); 
	output.fixDEC();
	return output;
}

#endif