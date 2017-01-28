#include <iostream>
#include <chrono>
#include <math.h>
#include <iomanip>
#include <vector>
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

		// interpolating the ephemeris coordinates using the ephemerides surrounding it
		int numInterp = 6;
		vector<Coords> nearbyCoords(numInterp);
		vector<double> nearbyTimes(numInterp);
		Ephemeris::findNearbyEphs(ephemerides, numInterp, i, nearbyCoords, nearbyTimes);
		printf("i=%6d\n", i);
		for (int j = 0; j < numInterp; j++) {
			printf("\tt%d=%11.1f ra=%8.2f dec=%8.2f\n",j,nearbyTimes[j], nearbyCoords[j].getDegRA(), nearbyCoords[j].getDegDEC());
		}
		/*Coords interpedCoords  = Coords::cartesianInterp(nearbyCoords, nearbyTimes, plateDate, 2);
		// UNCOMMENT THIS ONE IF THE ABOVE ONE DOESNT WORK
		//Coords interpedCoords  = Coords::cartesianInterp(before, beforeDate, after, afterDate, plateDate);
		Coords plateCoords	   = plates[p].getCoords();
		double angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
		double interpedApMag   = Ephemeris::linInterp(ephemerides[i].getApMag(), beforeDate, 
		                                              ephemerides[i+1].getApMag(), afterDate, 
		                                              plateDate);

		if (angularDistance < distanceThreshold && interpedApMag <= plates[p].getMagLimit()) {
			double xi, eta;
			int status;
			// calculate the x/y coordinates of the interpolated point, with the plate centre used as
			// the tangent point for the projection
			Coords::gnomonic(interpedCoords, plateCoords, xi, eta, status);
			// if the transformation went ok...
			if (status == 0) {
				Coords fromXAxis, fromYAxis;
				// calculate RA/DEC coordinates of the points that lie on the x/y axes respectively
				Coords::inverseGnomonic(xi,  0.0, interpedCoords, fromXAxis);
				Coords::inverseGnomonic(0.0, eta, interpedCoords, fromYAxis);
				// calculate the angular distance between the points (in degrees)
				double toXAxis = Coords::angularDistance(interpedCoords, fromXAxis);
				double toYAxis = Coords::angularDistance(interpedCoords, fromYAxis);
				// if the distances are both less than 3.2 degrees, the point will be on the plate
				if (toXAxis < 3.2 && toYAxis < 3.2) {
					Plate::printMatch(plates[p], interpedCoords, xi, eta, matchCount+1,
									  interpedApMag, toXAxis, toYAxis);
					matchCount++;
				}
			}
		} else if (interpedApMag > plates[p].getMagLimit()) {
			tooFaintCount++;
		}*/
	}

	// printing a summary of how many plates matched
	if (matchCount > 0) {
		cout << matchCount << " matching plates found for " << objectName << "!\n";
	} else {
		cout << "No matches found for " << objectName << "!\n";
	}

	// if any playes matched in terms of coordinates, but the object was too faint to show up,
	// tooFaintCount increments. If this has happened more than once, this bit is printed
	if (tooFaintCount > 0) {
		cout << tooFaintCount << (matchCount>0 ? " other" : "") << " plate" << (tooFaintCount==1 ? "" : "s");
		cout << " matched, but " << (tooFaintCount==1 ? "was" : "were") << " too faint to show up on the plate!\n";
	}






	// printing the total time taken when running the program
	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	printf("\nElapsed time: %.4fs\n\n", elapsed_seconds.count());
}


// CARRY ON CHECKING FINDNEARBYEPHS()
	// then check the new cartesianInterp() with test cases

// try to incorporate the errors in RA/DEC from the ephemeris somewhere
	// error propagation through the transformation?
	// even error circle around the point on the plate? 
	// probably distorted?