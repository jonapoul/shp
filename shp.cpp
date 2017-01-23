#include <iostream>
#include <chrono>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using namespace std;

int main(int argc, char* argv[]) {	
	chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
	


	vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catlog_ukstu.txt");
	for (auto p : plates)
		p.printPlate();
	vector<Ephemeris> ephemerides;
	Ephemeris::readEphemerisFile(ephemerides, "mars.txt");


	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	cout << "\nElapsed time: " << elapsed_seconds.count() << "s\n\n";
}


// carry on with Plate::convertDate()


// make ephemeris/plate constructors based on a string buffer
	// have to check whether the buffer is valid first before sending it to the constructor
// sort out projection for an image to coordinates
// interpolation between two coordinates of (ra,dec,time), linear/quadratic/cubic?
// around 6x6 degrees along each plate axis
// actually work out the limiting magnitude for each filter/emulsion, defaulting to 23 or if no matches