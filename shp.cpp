#include <iostream>
#include <chrono>
#include <math.h>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using namespace std;


int main(int argc, char* argv[]) {	
	chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
	
	// reading all the relevant info from the two files
	vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catlog_ukstu.txt");
	vector<Ephemeris> ephemerides;
	Ephemeris::readEphemerisFile(ephemerides, "mars.txt");

	// defining the first date in the ephemeris file, for later reference
	double firstDate = ephemerides[0].getJulian();
	// 6.4/sqrt(pi) is the radius of the circle that overlaps the most with a 6.4x6.4 square
	// im using this as a stopgap for now, until I can figure out the proper square one
	double distanceThreshold = 6.4 / sqrt(M_PI);
	int count = 0;

	// through every ephemeris record
	for (int p = 0; p < (int)plates.size()-1; p++) {
		double plateDate  = plates[p].getJulian();

		// i = index of the ephemeris record immediately before the plate date
		int i 			  = int( plateDate-firstDate );
		Coords before 	  = ephemerides[i].getCoords();
		double beforeDate = ephemerides[i].getJulian();
		Coords after 	  = ephemerides[i+1].getCoords();
		double afterDate  = ephemerides[i+1].getJulian();

		// linearly interpolating the ephemeris coordinates between the two records immediately
		// before and after the plate date
		Coords interped	  = Coords::interpolate(before, beforeDate, after, afterDate, plateDate);
		double angularDistance = Coords::angularDistance(interped, plates[p].getCoords(), true);

		if (angularDistance < distanceThreshold) {
			printf("plateID = %5d\n", plates[p].getID());
			printf("\tplateRA  =%10.4f deg\n", plates[p].getCoords().getDegRA());
			printf("\tplateDEC =%10.4f deg\n", plates[p].getCoords().getDegDEC());
			printf("\tephemRA  =%10.4f deg\n", interped.getDegRA());
			printf("\tephemDEC =%10.4f deg\n", interped.getDegDEC());
			printf("\tdistance =%10.4f deg\n", angularDistance);
			printf("\tUTdate   =%10s\n", plates[p].getGregorian().c_str());
			count++;
		}
	}
	if (count > 0) cout << count << " matching plates found within " << distanceThreshold << " degrees\n";
	else 		   cout << "No matches found!\n";



	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	cout << "\nElapsed time: " << elapsed_seconds.count() << "s\n\n";
}

// take into account the latitude of the observatory??? (-31.27 degrees)
// work out the limiting magnitude for each filter/emulsion, defaulting to 23 or if no matches
// 6.4x6.4 degrees along each plate axis
// make an attempt at quadratic/cubic fitting instead of linear
	// too time-consuming probably?