#include "headers/Plate.h"
#include "headers/Ephemeris.h"
#include <fstream>
#include <vector>
#include <cctype>
#include <stdio.h>
using std::cout;

int main(int argc, char* argv[]) {
	
	std::vector<Plate> plates;
	std::ifstream platesFile("catlog_ukstu.txt");
	if (platesFile.is_open()) {
		while (!platesFile.eof()) {
			std::string buffer;
			getline(platesFile, buffer);
			if (buffer.length() > 0) {
				Plate p;
				if (p.parsePlateString(buffer)) 
					plates.push_back(p);
			}
		}
		platesFile.close();
	};

	std::vector<Ephemeris> ephemerides;
	std::ifstream ephemerisFile("mars.txt");
	if (ephemerisFile.is_open()) {
		bool canReadEntries = false;
		while (!ephemerisFile.eof()) {
			std::string buffer;
			getline(ephemerisFile, buffer);
			if 		(buffer == "$$SOE") canReadEntries = true;
			else if (buffer == "$$EOE") canReadEntries = false;
			if (canReadEntries) {
				Ephemeris e;
				if (e.parseEphemerisString(buffer))
					ephemerides.push_back(e);
			}
		}
		ephemerisFile.close();
	}
	
}

// sort out projection for an image to coordinates
// linear interpolation between two coordinates of (ra,dec,time), maybe quadratic or cubic?
// convert date into julian calendar as a double value
// look up what LST/julian calendar are
// around 6x6 degrees shown on each plate
// carry on with reading an ephemeris line from the mars file