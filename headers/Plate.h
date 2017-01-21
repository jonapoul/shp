#ifndef PLATE_H
#define PLATE_H

#include "Coords.h"
using std::cout;

#define LIMITING_MAGNITUDE 23.0

class Plate {
private:
	int m_id;			// plate identification number
	Coords m_coords;	// combines RA and DEC to one object
	std::string m_date;	// Date of the plate being recorded, format = yymmdd
	std::string m_lst;	// LST at the start of the exposure, format = hhmm
	char m_grade;		// graded plate quality, A being best and C being 

public:
	Plate(int num = 0, const Coords& c = {}, const std::string& d = "", const std::string& lst = "", char grade = ' ');
	Plate(const Plate& p);
	Plate(const std::string& buffer);

	int getID() const;
	Coords getCoords() const;
	std::string getDate() const;
	std::string getLST() const;
	char getGrade() const;

	void setID(const int id);
	void setCoords(const Coords& c);
	void setDate(const std::string& d);
	void setLST(const std::string& lst);
	void setGrade(const char grade);

	bool parsePlateString(const std::string& buffer);
	void printPlate() const;
};

Plate::Plate(const int num, const Coords& c, const std::string& d, const std::string& lst, char grade) 
	: m_id(num), m_coords(c), m_date(d), m_lst(lst), m_grade(grade) { }
Plate::Plate(const Plate& p) 
	: m_id(p.m_id), m_coords(p.m_coords), m_date(p.m_date), m_lst(p.m_lst), m_grade(p.m_grade) { }
Plate::Plate(const std::string& buffer) {
	
}

int Plate::getID() const { return m_id; }
Coords Plate::getCoords() const { return m_coords; }
std::string Plate::getDate() const { return m_date; }
std::string Plate::getLST() const { return m_lst; }
char Plate::getGrade() const { return m_grade; }

void Plate::setID(const int id) { m_id = id; }
void Plate::setCoords(const Coords& c) { m_coords = c; }
void Plate::setDate(const std::string& d) { m_date = d; }
void Plate::setLST(const std::string& lst) { m_lst = lst; }
void Plate::setGrade(const char grade) { m_grade = grade; }

bool Plate::parsePlateString(const std::string& buffer) {
	// plate suffix column
	// if this isn't blank it indicates shenanigans while recording the image, such as a tracking shot or multiple images
	if (buffer[7] != ' ') return false;
	// checking for the word "TEST", s we can throw it out
	if (buffer.substr(16,4) == "TEST" || buffer.substr(15,4) == "TEST")
		return false;
	// plate number, for later reference
	m_id = stoi(buffer.substr(2, 5));
	// reading the RA/DEC coordinates in sexagesimal format
	m_coords.parseFromPlateRecord(buffer);
	// gregorian date string recorded as "yymmdd"
	m_date = buffer.substr(30, 6);
	// Local sidereal time when the plate was recorded
	m_lst = buffer.substr(36, 4);
	// invalid recorded LST => useless plate
	if (m_lst == "   0") 
		return false;
	// plate quality grade, from A to C
	m_grade = buffer[56];

	return true;
}

void Plate::printPlate() const {
	printf("ID = %8d RA = %s DEC = %s ", m_id, m_coords.getRA().toString().c_str(), m_coords.getDEC().toString().c_str());
	printf("Date = %s LST = %s Quality = %c\n", m_date.c_str(), m_lst.c_str(), m_grade);
}

#endif