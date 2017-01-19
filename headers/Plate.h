#ifndef PLATE_H
#define PLATE_H

#include "Coords.h"
using std::cout;

class Plate {
public:
	int m_number;		// plate identification number
	Coords m_coords;	// combines RA and DEC to one object
	std::string m_date;	// Date of the plate being recorded, format = yymmdd
	std::string m_lst;	// LST at the start of the exposure, format = hhmm

	Plate(int num = 0, const Coords& c = {});
	Plate(const Plate& p);
	Plate(int num, const RA& ra, const DEC& dec);
	void parsePlateData(std::string buffer);
};

Plate::Plate(int num, const Coords& c) : m_number(num), m_coords(c) { }

Plate::Plate(const Plate& p) : m_number(p.m_number), m_coords(p.m_coords) { }

Plate::Plate(int num, const RA& ra, const DEC& dec) : m_number(num), m_coords({ra},{dec}) { }

void Plate::parsePlateData(std::string buffer) {
	m_coords.m_ra.parseFromString(buffer);
	m_coords.m_dec.parseFromString(buffer);

	m_number = stoi(buffer.substr(2, 5));
	m_date = buffer.substr(30, 6);
	m_lst  = buffer.substr(36, 4);
	if (m_lst == "   0") m_lst = "0000";
	//cout << m_date << ' ' << m_lst << '\t';
}

#endif