#include <iostream>
#include <chrono>
#include <math.h>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using namespace std;


int main(int argc, char* argv[]) {	
	chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
	
	// reading all the relevant info from the plate catalog file and the ephemeris
	vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catalog.txt");
	vector<Ephemeris> ephemerides;
	string objectName;
	Ephemeris::readEphemerisFile(ephemerides, "ephemeris.txt", objectName);
	// defining the first date in the ephemeris file, for later reference
	double firstEphDate = ephemerides[0].getJulian();
	// This is the threshold used for the first check, to make sure the angular distance is reasonable
	// before i do the gnomonic and inverse transformations to get more specific
	double distanceThreshold = sqrt(2*3.2*3.2);
	int matchCount = 0, tooFaintCount = 0;

	// through every plate record
	for (int p = 0; p < (int)plates.size()-1; p++) {
		double plateDate  = plates[p].getJulian();

		// i = index of the ephemeris record immediately before the plate date
		int i 			  = int( plateDate-firstEphDate );
		Coords before 	  = ephemerides[i].getCoords();
		double beforeDate = ephemerides[i].getJulian();
		Coords after 	  = ephemerides[i+1].getCoords();
		double afterDate  = ephemerides[i+1].getJulian();

		// linearly interpolating the ephemeris coordinates between the two records immediately
		// before and after the plate date
		Coords interpedCoords  = Coords::interpolate(before, beforeDate, after, afterDate, plateDate);
		Coords plateCoords	   = plates[p].getCoords();
		double angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
		double interpedApMag   = Ephemeris::interpolate(ephemerides[i].getApMag(), beforeDate, 
		                                                ephemerides[i+1].getApMag(), afterDate, plateDate);

		if (angularDistance < distanceThreshold && interpedApMag <= plates[p].getMagLimit()) {
			double x, y;
			int status;
			// calculate the x/y coordinates of the interpolated point, with the plate centre used as
			// the tangent point for the projection
			Coords::gnomonic(interpedCoords, plateCoords, x, y, status);
			// if the transformation went ok...
			if (status == 0) {
				Coords fromXAxis, fromYAxis;
				// calculate RA/DEC coordinates of the points that lie on the x/y axes respectively
				Coords::inverseGnomonic(x, 0.0, interpedCoords, fromXAxis);
				Coords::inverseGnomonic(0.0, y, interpedCoords, fromYAxis);
				// calculate the angular distance between the points (in degrees)
				double distanceToXAxis = Coords::angularDistance(interpedCoords, fromXAxis);
				double distanceToYAxis = Coords::angularDistance(interpedCoords, fromYAxis);
				// if the distances are both less than 3.2 degrees, the point will be on the plate
				if (distanceToXAxis < 3.2 && distanceToYAxis < 3.2) {
					Plate::printMatch(plates[p], interpedCoords, angularDistance, x, y);
					printf("\tdtX=%f dtY=%f\n", distanceToXAxis, distanceToYAxis);
					matchCount++;
				}
			}
		} else if (interpedApMag > plates[p].getMagLimit()) {
			tooFaintCount++;
		}
	}
	if (matchCount > 0) {
		cout << matchCount << " matching plates found for " << objectName << "!\n";
	} else {
		cout << "No matches found for " << objectName << "!\n";
	}
	if (tooFaintCount > 0) {
		cout << tooFaintCount << (matchCount>0 ? " other" : "") << " plate" << (tooFaintCount==1 ? "" : "s");
	 	cout << " matched, but " << (tooFaintCount==1 ? "was" : "were") << " too faint to show up on the plate!\n";
	}



	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	printf("\nElapsed time: %.4fs\n\n", elapsed_seconds.count());
}


// do scaling to find position on the plate as millimetres from top and left sides

// try to incorporate the errors in RA/DEC from the ephemeris somewhere
// make an attempt at quadratic/cubic fitting instead of linear
	// too time-consuming probably?