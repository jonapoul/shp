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
public:
	bool m_isPositive;
	int m_deg;
	int m_min;
	float m_sec;

	DEC(bool isPos = true, int d = 0, int m = 0, float s = 0.f);
	DEC(const DEC& dec);
	DEC(float decimal);
	std::string toString();
	void fixDEC();
	float toDecimal();
	float toRadians();
	void parseFromString(std::string s);

	friend std::ostream& operator<<(std::ostream& os, RA& ra);
	friend DEC operator+ (DEC& a, DEC& b);
	friend DEC operator- (DEC& a, DEC& b);

};

DEC::DEC(bool isPos, int d, int m, float s) : m_isPositive(isPos), m_deg(d), m_min(m), m_sec(s) {
	if (d < 0) {
		m_isPositive = false;
		m_deg *= -1;
	}
	fixDEC();
}

DEC::DEC(const DEC& dec) : m_isPositive(dec.m_isPositive), m_deg(dec.m_deg), m_min(dec.m_min), m_sec(dec.m_sec) { }

DEC::DEC(float decimal) {
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

std::string DEC::toString() {
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
	ss << m_sec << '\"';
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

float DEC::toDecimal() {
	return (m_isPositive ? 1 : -1) * (m_deg + m_min*(1.f/60.f) + m_sec*(1.f/3600.f));
}

float DEC::toRadians() {
	float pi = 3.1415926535;
	return toDecimal() * (pi/180.f) * (m_isPositive ? 1 : -1);
}

void DEC::parseFromString(std::string s) {
	if (s.length() < 31) {
		cout << "string passed to DEC::parseFromString() is too short\n";
		return;
	}
	char sign = s[25];
	if (s[25] == '-')
		m_isPositive = false;
	else if (s[25] == '+' || s[25] == ' ') 
		m_isPositive = true;
	else {
		cout << "s[25] isn't a sign, backing out of DEC::parseFromString()\n";
		return;
	}
	m_deg = stoi(s.substr(26,2));
	m_min = stoi(s.substr(28,2));
	m_sec = 0.f;	
}

std::ostream& operator<<(std::ostream& os, DEC& dec) {
	os << dec.toString();
    return os;
}

DEC operator+(DEC& a, DEC& b) {
	DEC output(a.m_deg+b.m_deg, a.m_min+b.m_min, a.m_sec+b.m_sec); 
	output.fixDEC();
	return output;
}

DEC operator-(DEC& a, DEC& b) {
	DEC output(a.m_deg-b.m_deg, a.m_min-b.m_min, a.m_sec-b.m_sec); 
	output.fixDEC();
	return output;
}

#endif