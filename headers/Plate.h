#ifndef PLATE_H
#define PLATE_H

#include <vector>
#include <fstream>
#include <sstream>
#include "Coords.h"
using namespace std;

class Plate {
private:
	int m_id;				// plate identification number
	Coords m_coords;		// combines RA and DEC to one object
	float m_julianDate;		// julian date at the start of exposure
	string m_lst;			// LST at the start of the exposure, format = hhmm
	char m_grade;			// graded plate quality, A being best and C being 

public:
	Plate(int num=0, const Coords& c={}, const float jd=0.f, const string& lst="", const char grade=' ');
	Plate(const Plate& p);

	int getID() const;
	Coords getCoords() const;
	RA getRA() const;
	DEC getDEC() const;
	float getDate() const;
	string getLST() const;
	char getGrade() const;

	void setID(const int id);
	void setCoords(const Coords& c);
	void setRA(const RA& ra);
	void setDEC(const DEC& dec);
	void setDate(const float d);
	void setLST(const string& lst);
	void setGrade(const char grade);

	bool parsePlateString(const string& buffer);
	void printPlate() const;
	static float gregorianToJulian(const string& dateStr);
	static string julianToGregorian(float julian);
	static void readPlateCatalog(vector<Plate>& plates, const string& filename);
	static int partition(vector<Plate>& plates, const int p, const int q);
	static void quickSort(vector<Plate>& plates, const int p, const int q);
};

Plate::Plate(const int num, const Coords& c, const float jd, const string& lst, const char grade) 
	: m_id(num), m_coords(c), m_julianDate(jd), 
	  m_lst(lst), m_grade(grade) { }
Plate::Plate(const Plate& p) 
	: m_id(p.m_id), m_coords(p.m_coords), m_julianDate(p.m_julianDate),
	  m_lst(p.m_lst), m_grade(p.m_grade) { }

int Plate::getID() const { return m_id; }
Coords Plate::getCoords() const { return m_coords; }
RA Plate::getRA() const { return m_coords.getRA(); }
DEC Plate::getDEC() const {return m_coords.getDEC(); }
float Plate::getDate() const { return m_julianDate; }
string Plate::getLST() const { return m_lst; }
char Plate::getGrade() const { return m_grade; }

void Plate::setID(const int id) { m_id = id; }
void Plate::setCoords(const Coords& c) { m_coords = c; }
void Plate::setRA(const RA& ra) { m_coords.setRA(ra); }
void Plate::setDEC(const DEC &dec) { m_coords.setDEC(dec); }
void Plate::setDate(const float d) { m_julianDate = d; }
void Plate::setLST(const string& lst) { m_lst = lst; }
void Plate::setGrade(const char grade) { m_grade = grade; }

bool Plate::parsePlateString(const string& buffer) {
	// plate suffix column
	// if this isn't blank it indicates shenanigans while recording the image, such as a tracking shot or multiple images
	// T = tracked shot, M = multiple shots, P = full-aperture prism
	if (buffer[7] == 'T' || buffer[7] == 'M' || buffer[7] == 'P') 
		return false;
	// checking for the word "TEST", s we can throw it out
	if (buffer.substr(16,4) == "TEST" || buffer.substr(15,4) == "TEST")
		return false;
	// checking for an invalid LST string (if any of the 4 digits are blank)
	m_lst = buffer.substr(36, 4);
	if (stoi(m_lst) < 1000) 
		return false;
	// plate number, for later reference
	m_id = stoi(buffer.substr(2, 5));
	// reading the RA/DEC coordinates in sexagesimal format
	m_coords.parseFromPlateRecord(buffer);
	// julian date calculated from gregorian
	m_julianDate = gregorianToJulian(buffer.substr(30, 6));
	// plate quality grade, from A to C
	m_grade = buffer[56];

	return true;
}

void Plate::printPlate() const {
	printf("ID = %6d RA = %s DEC = %s ", m_id, m_coords.getRA().toString().c_str(), m_coords.getDEC().toString().c_str());
	printf("Date = %.1f LST = %s Quality = %c\n", m_julianDate, m_lst.c_str(), m_grade);
}

float Plate::gregorianToJulian(const string& gregorian) {
	int date = stoi(gregorian);
	int y = stoi(gregorian.substr(0, 2));
	int m = stoi(gregorian.substr(2, 2));
	int d = stoi(gregorian.substr(4));

	y += (y < 17) ? 2000 : 1900;
	if (m < 3) { y--; m += 12; }

	int A = y / 100;
	int B = 2 - A + (A/4);
	int C = int(365.25 * y);
	int D = int(30.6001 * (m+1));

	return B + C + D + d + 1720994.5;
}

string Plate::julianToGregorian(float julian) {
	julian += 0.5;
	int I = int(julian);
	float F = julian - I;

	int B;
	if (I > 2299160) {
		int A = int( (I-1867216.25) / 36524.25 );
		B = I + A + 1 - (A/4);
	}
	else B = I;
	int C = B + 1524;
	int D = int( (C-122.1) / 365.25 );
	int E = int(365.25 * D);
	int G = int( (C-E) / 30.6001);

	float d = C - E + F - int(30.6001 * G);
	int m = (G < 13.5) ? G-1 : G-13;
	int y = (m < 2.5) ? D-4715 : D-4716;

	// DO SOMETHING TO GET LST OUT HERE

	stringstream ss;
	y -= ((y < 16) ? 1900 : 2000);
	ss << ((y < 10) ? "0" : "");
	ss << y;
	ss << ((m < 10) ? "0" : "");
	ss << m;
	ss << ((d < 10) ? "0" : "");
	ss << int(d);
	return ss.str();
}

void Plate::readPlateCatalog(vector<Plate>& plates, const string& filename) {
	ifstream platesFile(filename);
	if (platesFile.is_open()) {
		while (!platesFile.eof()) {
			string buffer;
			getline(platesFile, buffer);
			if (buffer.length() > 0) {
				Plate p;
				if (p.parsePlateString(buffer)) 
					plates.push_back(p);
			}
		}
		platesFile.close();
	}
	else {
		cout << "Plate file \"" << filename << "\" is not valid\n";
	}
}

int Plate::partition(vector<Plate>& plates, const int p, const int q) {
	Plate x = plates[p];
	int i = p;
	for (int j = p+1; j < q; j++) {
		if (plates[j].getDate() <= x.getDate()) {
			i++;
			swap(plates[i], plates[j]);
		}
	}
	swap(plates[i], plates[p]);
	return i;
}

void Plate::quickSort(vector<Plate>& plates, const int p, const int q) {
	int r;
	if (p < q) {
		r = partition(plates, p, q);
		quickSort(plates, p, r);  
		quickSort(plates, r+1, q);
	}
}

#endif