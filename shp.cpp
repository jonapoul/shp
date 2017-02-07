#include <iostream>
#include <chrono>
#include <math.h>
#include <vector>
#include <utility>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using namespace std;

int main(int argc, char* argv[]) {	
	chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();

	// determining the program parameters from command line arguments
	int degree, numInterp;
	string objectName;
	Ephemeris::determineParameters(argc, argv, objectName, numInterp, degree);
	// an array of plate IDs that I couldn't find in the archive room
	vector<int> blacklist = {4910, 14587};
	// reading all valid plate records
	vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catalog.txt");
	// reading all ephemeris records
	vector<Ephemeris> eph;
	Ephemeris::readEphemerisFile(eph, "ephemeris/"+objectName+".txt");
	objectName[0] = toupper(objectName[0]);

	// defining the first date in the ephemeris file, for later reference
	double firstEphDate = eph[0].julian();
	double secondEphDate = eph[1].julian();
	double step = secondEphDate - firstEphDate;

	// sqrt(3.2*3.2 + 3.2*3.2) is the minimum circular radius for the object to possibly be on the plate. This is used as a first check before performing a polynomial fit, to save processing time
	double distanceThreshold = sqrt(2 * 3.2*3.2);
	int matchCount = 0, tooFaintCount = 0, missingPlateCount = 0;
	int loopLimit = plates.size() - 1;

	vector<Plate> matchPlates;
	vector<Coords> matchCoords;
	vector<int> matchCounts;
	vector<double> matchMags;
	vector<pair<double,double>> matchStart, matchMid, matchEnd;

	// through every plate record
	for (int p = 0; p < loopLimit; p++) {
		double plateDate  = plates[p].julian();

		// i = index of the ephemeris record immediately before the plate date
		int i = int( (plateDate-firstEphDate)/step );
		if (i < 0) 	// this happens if the first ephemeris record is after the first plate
			continue;
		Coords before 	  = eph[i].coords();
		double beforeDate = eph[i].julian();
		Coords after 	  = eph[i+1].coords();
		double afterDate  = eph[i+1].julian();

		// approximate linear interpolation
		Coords interpedCoords  = Coords::linInterp(before, beforeDate, after, afterDate, plateDate);

		// calculating the angular distance between the interpolated coordinates and the plate centre
		Coords plateCoords	   = plates[p].coords();
		double angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
		double magnitude 	   = Ephemeris::linInterp(eph[i].mag(), beforeDate, eph[i+1].mag(), afterDate, plateDate);

		// if the object is within a reasonable angular distance of the plate, do more accurate tests
		if (angularDistance < 1.5*distanceThreshold && magnitude <= plates[p].magLimit()) {

			// interpolating the ephemeris coordinates using the ephemeris records surrounding it
			vector<Coords> nearbyCoords(numInterp);
			vector<double> nearbyTimes(numInterp);
			Ephemeris::findNearbyEphs(eph, numInterp, i, nearbyCoords, nearbyTimes);

			// performing the polynomial least-squares fit using the nearby ephemeris records
			interpedCoords = Coords::polyInterp(nearbyCoords, nearbyTimes, plateDate, degree);

			// recalculating the angular distance to make sure the object is definitely somewhere on the plate
			angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
			if (angularDistance < distanceThreshold) {
				double xi, eta;
				int status;
				// calculate the x/y coordinates of the interpolated point, with the plate centre as the tangent point
				Coords::gnomonic(interpedCoords, plateCoords, xi, eta, status);
				// if the transformation went ok
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
						// checks whether the matched plate is known to be missing
						if ( Plate::plateIsMissing(plates[p], blacklist) ) {
							missingPlateCount++;
							continue;
						}
						// calculating the xi/eta coords of the points at the start and end of the exposure
						pair<double,double> start, end, mid(xi, eta);
						Plate::exposureBoundaries(plates[p], nearbyCoords, nearbyTimes, degree, start, end);
						// adding the relevant values to arrays to be printed at the end
						matchCounts.push_back(++matchCount);
						matchPlates.push_back(plates[p]);
						matchCoords.push_back(interpedCoords);
						matchMags.push_back(magnitude);
						matchStart.push_back(start);
						matchMid.push_back(mid);
						matchEnd.push_back(end);
						double x = (xi  * (3600.0 / 67.12) * (180.0/M_PI)) + (355.0/2.0);
						double y = (eta * (3600.0 / 67.12) * (180.0/M_PI)) + (355.0/2.0);
						printf("%.3f %.3f\n", x, y);
					}
				}
			}
		} else if (magnitude > plates[p].magLimit()) {
			tooFaintCount++;
		}
	}

	// printing a summary of how many plates matched
	//Plate::printMatches(matchPlates, matchCoords, matchCounts, matchMags, matchStart, matchMid, matchEnd);
	string firstDate = Plate::julianToGregorian(firstEphDate);
	string lastDate  = Plate::julianToGregorian(eph[eph.size()-1].julian());
	if (matchCount > 0) {
		printf("%d matching plates found for %s between %s and %s\n",
		       	matchCount,	objectName.c_str(), firstDate.c_str(), lastDate.c_str() );
	}
	else
		printf("No matches found for %s!\n", objectName.c_str());

	// if any playes matched in terms of coordinates, but the object was too faint to show up, tooFaintCount increments. If this has happened more than once, this bit is printed
	if (tooFaintCount > 0) {
		cout << tooFaintCount << (matchCount>0 ? " other" : "") << " plate" << (tooFaintCount==1 ? "" : "s");
		cout << " matched, but " << (tooFaintCount==1 ? "was" : "were") << " too faint to show up on the plate\n";
	}
	if (missingPlateCount > 0) {
		cout << missingPlateCount << (matchCount>0 ? " other" : "") << " plate" << (missingPlateCount==1 ? "" : "s");
		cout << " matched, but " << (missingPlateCount==1 ? "isn't" : "aren't") << " in the plate room\n";
	}
	string sign = "th";
	if 		(degree % 10 == 1) sign = "st";
	else if (degree % 10 == 2) sign = "nd";
	else if (degree % 10 == 3) sign = "rd";
	printf("Interpolated using %d surrounding ephemerides in a %d%s order polynomial\n", numInterp, degree, sign.c_str());

	// printing the total time taken when running the program	
	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	printf("Elapsed time: %.4fs\n", elapsed_seconds.count());
}


// TO DO
	// proper lininterp instead of polyfit
	// test 		double frac = (ut) / 24.0; instead of ut-12.0 on convertDate()
	// regularise time points??
	// use comparison shot from another sky survey to tell if the object im looking at is actually what i think it is
	// scale magnitude limits based on exposure times (logarithmically?)
	// try to incorporate the errors in RA/DEC from the ephemeris somewhere
		// error propagation through the transformation?
		// even error circle around the point on the plate? 
		// probably distorted?
		// comes from extrapolating objects paths when they were recently discovered, eg Eris/Makemake