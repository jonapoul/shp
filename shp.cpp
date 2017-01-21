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

// calculate julian day using both date plus LST (on plates)
	// bottom of page 9 in book
// make ephemeris/plate constructors based on a string buffer
	// have to check whether the buffer is valid first
// sort out projection for an image to coordinates
// linear interpolation between two coordinates of (ra,dec,time), maybe quadratic or cubic?
// look up what LST/julian calendar are
// around 6x6 degrees shown on each plate