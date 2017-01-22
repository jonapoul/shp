#include <iostream>
#include <chrono>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using std::cout;

int main(int argc, char* argv[]) {	
	std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();	


	std::vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catlog_ukstu.txt");
	std::vector<Ephemeris> ephemerides;
	Ephemeris::readEphemerisFile(ephemerides, "mars.txt");

	Plate::quickSort(plates, 0, (int)plates.size());

	for (auto p : plates)
		p.printPlate();

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}

// calculate julian day using both date plus LST (on plates)
	// bottom of page 9 in book
// make ephemeris/plate constructors based on a string buffer
	// have to check whether the buffer is valid first before sending it to the constructor
// sort out projection for an image to coordinates
// interpolation between two coordinates of (ra,dec,time), linear/quadratic/cubic?
// look up what LST/julian calendar are
// around 6x6 degrees along each plate axis
// actually work out the limiting magnitude for each filter/emulsion, defaulting to 23 or if no matches