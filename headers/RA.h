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
public:
	int m_hour;
	int m_mins;
	float m_secs;
	double m_degrees;
	double m_radians;

	RA(int h = 0, int m = 0, float s = 0.f, double deg = 0.f, double rad = 0.f);
	RA(const RA& ra);
	RA(float decimal);
	int getHour();
	int getMins();
	float getSecs();
	std::string toString();
	void fixRA();
	float toDecimal();
	float toRadians();
	void parseFromString(std::string s);

	friend std::ostream& operator<<(std::ostream& os, RA& ra);
	friend RA operator+ (RA& a, RA& b);
	friend RA operator- (RA& a, RA& b);
};

RA::RA(int h, int m, float s, double deg, double rad) 
	: m_hour(h), m_mins(m), m_secs(s), m_degrees(deg), m_radians(rad) { 
	if (abs(m_degrees - 0.f) < 1e-10) {
		m_degrees = toDecimal();
		m_radians = m_degrees * M_PI / 180.f;
	}
}

RA::RA(const RA& ra) 
	: m_hour(ra.m_hour), m_mins(ra.m_mins), m_secs(ra.m_secs), m_degrees(ra.m_degrees), m_radians(ra.m_radians) {
	if (abs(m_degrees - 0.f) < 1e-10) {
		m_degrees = toDecimal();
		m_radians = m_degrees * M_PI / 180.f;
	}
}

RA::RA(float decimal) {
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

std::string RA::toString() {
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2);
	if (m_hour < 10)	ss << '0';
	ss << m_hour << 'h';
	if (m_mins < 10)	ss << '0';
	ss << m_mins << 'm';
	if (m_secs < 10) 	ss << '0';
	ss << m_secs << 's';
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

float RA::toDecimal() {
	return m_hour*(360/24) + m_mins*0.25 + m_secs*(1.f/240.f);
}

float RA::toRadians() {
	float pi = 3.1415926535;
	return toDecimal() * (pi/180.f);
}

void RA::parseFromString(std::string s) {
	if (s.length() < 31) {
		cout << "string passed to RA::parseFromString() is too short\n";
		return;
	}

	m_hour = stoi(s.substr(20,2));
	m_mins = stoi(s.substr(22,2));
	m_secs = stoi(s.substr(24,1)) * 6.f;
}

std::ostream& operator<<(std::ostream& os, RA& ra) {
	os << ra.toString();
    return os;
}

RA operator+(RA& a, RA& b) {
	RA output(a.m_hour+b.m_hour, a.m_mins+b.m_mins, a.m_secs+b.m_secs); 
	output.fixRA();
	return output;
}

RA operator-(RA& a, RA& b) {
	RA output(a.m_hour-b.m_hour, a.m_mins-b.m_mins, a.m_secs-b.m_secs); 
	output.fixRA();
	return output;
}

#endif