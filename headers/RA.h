#ifndef RA_H
#define RA_H

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <sstream>
#include "DEC.h"
using namespace std;

class RA {
private:
	int m_hour;
	int m_mins;
	float m_secs;
	double m_degrees;
	double m_radians;

public:
	RA(const int h = 0, const int m = 0, const float s = 0.f);
	RA(const RA& ra);
	RA(const double decimal);

	int    getHours() 	const { return m_hour; };
	int    getMins() 	const { return m_mins; };
	float  getSecs() 	const { return m_secs; };
	double getDegrees() const { return m_degrees; };
	double getRadians() const { return m_radians; };

	void setHours  (const int h) 	{ m_hour = h; };
	void setMins   (const int m)	{ m_mins = m; };
	void setSecs   (const float s)	{ m_secs = s; };
	void setDegrees(const double d)	{ m_degrees = d; m_radians = d*M_PI/180; };
	void setRadians(const double r)	{ m_radians = r; m_degrees = r*180/M_PI; };

	void   fixRA();
	string toString()  const;
	float  toDecimal() const;
	float  toRadians() const;

	friend ostream& operator<<(ostream& os, const RA& ra);
	friend bool 	operator==(const RA& a, const RA& b);
	friend bool 	operator!=(const RA& a, const RA& b);
	friend bool 	operator< (const RA& a, const RA& b);
	friend bool 	operator> (const RA& a, const RA& b);
	friend bool 	operator<=(const RA& a, const RA& b);
	friend bool 	operator>=(const RA& a, const RA& b);

};

RA::RA(const int h, const int m, const float s) 
	: m_hour(h), m_mins(m), m_secs(s) { 
	m_degrees = toDecimal();
	m_radians = m_degrees * M_PI / 180.f;
}

RA::RA(const RA& ra) 
	: m_hour(ra.m_hour), m_mins(ra.m_mins), m_secs(ra.m_secs), m_degrees(ra.m_degrees), m_radians(ra.m_radians) {
	m_degrees = toDecimal();
	m_radians = m_degrees * M_PI / 180.f;
}

RA::RA(const double decimal) : m_degrees(decimal) {
	while (m_degrees < 0.f) 
		m_degrees += 360.f;

	m_radians = m_degrees * M_PI / 180.f;

	float hourDecimal = m_degrees * 24.f/360.f;
	m_hour = int(hourDecimal);
	hourDecimal = (hourDecimal - m_hour) * 60.f;
	m_mins = int(hourDecimal);
	hourDecimal = (hourDecimal - m_mins) * 60.f;
	m_secs = hourDecimal;
}

void RA::fixRA() {
	while (m_secs < 0.f)	{ m_secs += 60.f; 	m_mins -= 1; }
	while (m_secs > 60.f)	{ m_secs -= 60.f; 	m_mins += 1; }
	while (m_mins < 0)		{ m_mins += 60; 	m_hour -= 1; }
	while (m_mins > 60)		{ m_mins -= 60;	 	m_hour += 1; }
	while (m_hour < 0)		{ m_hour += 24; }
	while (m_hour > 24)		{ m_hour -= 24; }
}

string RA::toString() const {
	string output;
	if (m_hour < 10) output += '0';
	output += to_string(m_hour) + 'h';
	if (m_mins < 10) output += '0';
	output += to_string(m_mins) + 'm';
	if (m_secs < 10) output += '0';
	output += to_string(int(m_secs)) + 's';
	return output;
}

float RA::toDecimal() const {
	return m_hour*(360/24) + m_mins*0.25 + m_secs*(1.f/240.f);
}

float RA::toRadians() const {
	float pi = 3.1415926535;
	return toDecimal() * (pi/180.f);
}

ostream& operator<<(ostream& os, const RA& ra) {
	os << ra.toString();
    return os;
}

bool operator==(const RA& a, const RA& b) { 
	return abs(a.m_degrees - b.m_degrees) < 1e-6; 
}

bool operator!=(const RA& a, const RA& b) { 
	return !(a == b); 
}

bool operator<(const RA& a, const RA& b) { 
	return (a.m_degrees < b.m_degrees); 
}

bool operator>(const RA& a, const RA& b) { 
	return (abs(a.m_degrees - b.m_degrees) < 1e-6) ? false : !(a < b); 
}

bool operator<=(const RA& a, const RA& b) { 
	return (a.m_degrees <= b.m_degrees); 
}

bool operator>=(const RA& a, const RA& b) { 
	return (abs(a.m_degrees - b.m_degrees) < 1e-6) ? false : !(a < b); 
}

#endif