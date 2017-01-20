#include "headers/Plate.h"
#include <fstream>
#include <vector>
#include <cctype>
#include <stdio.h>
using std::cout;

int main(int argc, char* argv[]) {
	
	std::vector<Plate> plates;
	std::ifstream file("catlog_ukstu.txt");
	if (file.is_open()) {
		while (!file.eof()) {
			std::string buffer;
			getline(file, buffer);
			if (buffer.length() > 0) {
				Plate p;
				if (p.parsePlateData(buffer))
					plates.push_back(p);
			}
		}
		file.close();
	}

	std::ifstream ephemerisFile("mars.txt");
	if (ephemerisFile.is_open()) {
		while (!ephemerisFile.eof()) {
			std::string buffer;
			getline(ephemerisFile, buffer);
			if (buffer.length() > 5) {
				if ( isdigit(buffer[0]) && isdigit(buffer[1]) && isdigit(buffer[2]) ) {
					double julianDate = std::stod(buffer.substr(0,17));
					//cout << julianDate << '\n';
				}
			}
		}
		ephemerisFile.close();
	}

	// for (auto p : plates) {
	// 	///////////////////////////////////////////////////////////////
	// 	// Can't print date or lst for some reason?
	// 	/////////////////////////////////////////////////////////////
	// 	cout << p.m_number << ' ' << p.m_date << ' ' << p.m_lst << '\n';
	// }

	
}

// sort out projection for an image to coordinates
// linear interpolation between two coordinates of (ra,dec,time), maybe quadratic or cubic?
// convert date into julian calendar as a double value
// convert any ra/dec value into decimal degrees and radians to be stored in the object AS ITS READ IN, to save computation
// look up what LST/julian calendar are
// around 6x6 degrees shown on each plate
// ignore tracked/multiple plates from suffix option
// carry on with reading an ephemeris line from the mars file
// check for the word "TEST", so the plate can be discarded. It comes up a lot in the last catalog file