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
				p.parsePlateData(buffer);
				plates.push_back(p);
			}
		}
		file.close();
	}

	std::ifstream ephemerisFile("mars_ephemeris.txt");
	if (ephemerisFile.is_open()) {
		while (!ephemerisFile.eof()) {
			std::string buffer;
			getline(ephemerisFile, buffer);
			if (buffer.length() > 5) {
				std::string dateStr = buffer.substr(1,4);
				auto isAllDigits = [](std::string s) {
					for (auto c : s)
						if (!isdigit(c))	return false;
					return true;
				};
				if ( isAllDigits(dateStr) ) {
					//cout << dateStr << '\n';
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

// carry on with reading an ephemeris line from the mars file
// also find a good timespace between records, is 4h too much? Jumps by ~2arcmin per record for DEC 
// carry on with Particle simulation functions
// check for the word "TEST", so the plate can be discarded. It comes up a lot in the last catalog file
// test gnomonic projection, how fast/efficient is it?