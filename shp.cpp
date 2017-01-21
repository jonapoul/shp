#include "headers/Plate.h"
#include "headers/Ephemeris.h"
#include <vector>
#include <stdio.h>
using std::cout;

int main(int argc, char* argv[]) {
	
	std::vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catlog_ukstu.txt");

	std::vector<Ephemeris> ephemerides;
	Ephemeris::readEphemerisFile(ephemerides, "mars.txt");
	

	return 0;	
}

// sort out projection for an image to coordinates
// linear interpolation between two coordinates of (ra,dec,time), maybe quadratic or cubic?
// convert date into julian calendar as a double value
// look up what LST/julian calendar are
// around 6x6 degrees shown on each plate