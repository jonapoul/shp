#ifndef RA_H
#define RA_H

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include "DEC.h"
using std::cout;

class DEC;

class RA {
private:
	int m_hour;
	int m_mins;
	float m_secs;
	double m_degrees;
	double m_radians;

public:
	RA(int h = 0, int m = 0, float s = 0.f);
	RA(const RA& ra);
	RA(double decimal);

	int getHours() const;
	int getMins() const;
	float getSecs() const;
	double getDegrees() const;
	double getRadians() const;

	void setHours(const int h);
	void setMins(const int m);
	void setSecs(const float s);
	void setDegrees(const double d);
	void setRadians(const double r);

	std::string toString() const;
	void fixRA();
	float toDecimal() const;
	float toRadians() const;

	friend std::ostream& operator<<(std::ostream& os, const RA& ra);
	friend RA operator+ (const RA& a, const RA& b);
	friend RA operator- (const RA& a, const RA& b);
};

RA::RA(int h, int m, float s) 
	: m_hour(h), m_mins(m), m_secs(s) { 
	m_degrees = toDecimal();
	m_radians = m_degrees * M_PI / 180.f;
}

RA::RA(const RA& ra) 
	: m_hour(ra.m_hour), m_mins(ra.m_mins), m_secs(ra.m_secs), m_degrees(ra.m_degrees), m_radians(ra.m_radians) {
	m_degrees = toDecimal();
	m_radians = m_degrees * M_PI / 180.f;
}

RA::RA(double decimal) {
	while (decimal < 0.f) 
		decimal += 360.f;

	m_degrees = decimal;
	m_radians = decimal * M_PI / 180.f;

	float hourDecimal = decimal * 24.f/360.f;
	m_hour = int(hourDecimal);
	hourDecimal = (hourDecimal - m_hour) * 60.f;
	m_mins = int(hourDecimal);
	hourDecimal = (hourDecimal - m_mins) * 60.f;
	m_secs = hourDecimal;
}

int RA::getHours() const { return m_hour; }
int RA::getMins()const { return m_mins; }
float RA::getSecs() const { return m_secs; }
double RA::getDegrees() const { return m_degrees; }
double RA::getRadians() const { return m_radians; }

void RA::setHours(const int h) { m_hour = h; }
void RA::setMins(const int m) { m_mins = m; }
void RA::setSecs(const float s) { m_secs = s; }
void RA::setDegrees(const double d) { m_degrees = d; }
void RA::setRadians(const double r) { m_radians = r; }

std::string RA::toString() const {
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2);
	if (m_hour < 10)	ss << '0';
	ss << m_hour << 'h';
	if (m_mins < 10)	ss << '0';
	ss << m_mins << 'm';
	if (m_secs < 10) 	ss << '0';
	ss << int(m_secs) << 's';
	return ss.str();
}

void RA::fixRA() {
	while (m_secs < 0.f)	{ m_secs += 60.f; 	m_mins -= 1; }
	while (m_secs > 60.f)	{ m_secs -= 60.f; 	m_mins += 1; }
	while (m_mins < 0)		{ m_mins += 60; 	m_hour -= 1; }
	while (m_mins > 60)		{ m_mins -= 60;	 	m_hour += 1; }
	while (m_hour < 0)		{ m_hour += 24; }
	while (m_hour > 24)		{ m_hour -= 24; }
}

float RA::toDecimal() const {
	return m_hour*(360/24) + m_mins*0.25 + m_secs*(1.f/240.f);
}

float RA::toRadians() const {
	float pi = 3.1415926535;
	return toDecimal() * (pi/180.f);
}

std::ostream& operator<<(std::ostream& os, const RA& ra) {
	os << ra.toString();
    return os;
}

RA operator+(const RA& a, const RA& b) {
	RA output(a.m_hour+b.m_hour, a.m_mins+b.m_mins, a.m_secs+b.m_secs); 
	output.fixRA();
	return output;
}

RA operator-(const RA& a, const RA& b) {
	RA output(a.m_hour-b.m_hour, a.m_mins-b.m_mins, a.m_secs-b.m_secs); 
	output.fixRA();
	return output;
}

#endif