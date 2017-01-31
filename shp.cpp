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
	double secondEphDate = ephemerides[1].getJulian();
	double step = secondEphDate - firstEphDate;

	// sqrt(3.2*3.2 + 3.2*3.2) is the minimum circular radius for the object to possibly
	// be on the plate. This is used as a first check before performing a polynomial fit
	// to get more accurate, to save processing time
	double distanceThreshold = sqrt(2*3.2*3.2);
	int matchCount = 0, tooFaintCount = 0;

	// through every plate record
	for (int p = 0; p < (int)plates.size()-1; p++) {
		double plateDate  = plates[p].getJulian();

		// i = index of the ephemeris record immediately before the plate date
		int i 			  = int( (plateDate-firstEphDate)/step );
		Coords before 	  = ephemerides[i].getCoords();
		double beforeDate = ephemerides[i].getJulian();
		Coords after 	  = ephemerides[i+1].getCoords();
		double afterDate  = ephemerides[i+1].getJulian();

		// initial interpolation
		Coords interpedCoords  = Coords::linInterp(before, beforeDate, after, afterDate, plateDate);
		
		// calculating the angular distance between the linearly interpolated coordinates and
		// the plate centre
		Coords plateCoords	   = plates[p].getCoords();
		double angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
		double interpedApMag   = Ephemeris::linInterp(ephemerides[i].getApMag(), beforeDate, ephemerides[i+1].getApMag(), afterDate, plateDate);

		// if the object is within a reasonable angular distance of the plate, do more accurate tests
		if (angularDistance < 2*distanceThreshold && interpedApMag <= plates[p].getMagLimit()) {
			// interpolating the ephemeris coordinates using the ephemeris records surrounding it
			int numInterp = 10;
			vector<Coords> nearbyCoords(numInterp);
			vector<double> nearbyTimes(numInterp);
			Ephemeris::findNearbyEphs(ephemerides, numInterp, i, nearbyCoords, nearbyTimes);

			// performing the polynomial least-squares fit using the nearby records
			interpedCoords = Coords::polyInterp(nearbyCoords, nearbyTimes, plateDate, 2);

			// recalculating the angular distance to make sure the object is definitely
			// somewhere on the plate
			angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
			if (angularDistance < distanceThreshold) {
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
						Plate::printMatch(plates[p], interpedCoords, xi, eta, matchCount+1, interpedApMag, toXAxis, toYAxis);
						matchCount++;
					}
				}
			}
		} else if (interpedApMag > plates[p].getMagLimit()) {
			tooFaintCount++;
		}
	}

	// printing a summary of how many plates matched
	if (matchCount > 0)
		cout << matchCount << " matching plates found for " << objectName << "!\n";
	else
		cout << "No matches found for " << objectName << "!\n";

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


// check the new cartesianInterp() with test cases
	// lu_substitute in polyfit() with id=6590, all others are fine

// try to incorporate the errors in RA/DEC from the ephemeris somewhere
	// error propagation through the transformation?
	// even error circle around the point on the plate? 
	// probably distorted?