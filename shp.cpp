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
	string objectName;
	bool couldntRead;
	Ephemeris::determineParameters(argc, argv, objectName, couldntRead);
	// an array of plate IDs that I couldn't find in the archive room
	vector<int> blacklist = {4910, 14587};
	// reading all valid plate records
	vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catalog.txt");
	// reading all ephemeris records
	vector<Ephemeris> eph;
	Ephemeris::readEphemerisFile(eph, objectName);
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

	double closest = 180.0;

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
		// this is only used if zero plates match, to give the user an idea of how close any of them were
		closest = (angularDistance < closest) ? angularDistance : closest;

		// if the object is within a reasonable angular distance of the plate, do more accurate tests
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
					// checks whether the plate's magnitude limit is sufficient to spot the object
					double magnitude = Ephemeris::linInterp(eph[i].mag(), beforeDate, eph[i+1].mag(), afterDate, plateDate);
					if (magnitude > plates[p].magLimit()) {
						tooFaintCount++;
						continue;
					}
					// xi/eta coords of the points at the start, middle and end of the exposure
					pair<double,double> start, mid(xi, eta), end;
					Plate::exposureBoundaries(plates[p], {before,after}, {beforeDate,afterDate}, start, end);
					// adding the relevant values to arrays to be printed at the end
					matchCounts.push_back(++matchCount);
					matchPlates.push_back(plates[p]);
					matchCoords.push_back(interpedCoords);
					matchMags.push_back(magnitude);
					matchStart.push_back(start);
					matchMid.push_back(mid);
					matchEnd.push_back(end);
				}
			}
		}
	}

	// printing a summary of how many plates matched
	Plate::printMatches(matchPlates, matchCoords, matchCounts, matchMags, matchStart, matchMid, matchEnd);
	string firstDate = Plate::julianToGregorian(firstEphDate);
	string lastDate  = Plate::julianToGregorian(eph[eph.size()-1].julian());
	if (argc < 2) {
		cout << "Defaulted to CERES\n";
	} else if (couldntRead) {
		cout << "Filename \"" << argv[1] << ".txt\" doesn't exist in /ephemeris, defaulted to CERES\n";
	}
	for (auto& c : objectName) c = toupper(c);
	if (matchCount > 0) {
		cout << matchCount << " matching plate" << (matchCount>1?"s":"") << " found for " << objectName;
		cout << " between " << firstDate << " and " << lastDate << '\n';
	} else {
		cout << "No matches found for " << objectName << "!\n";
	}

	// if any playes matched in terms of coordinates, but the object was too faint to show up, tooFaintCount increments. If this has happened more than once, this bit is printed
	if (tooFaintCount > 0) {
		cout << tooFaintCount << (matchCount>0 ? " other" : "") << " plate" << (tooFaintCount==1 ? "" : "s");
		cout << " matched, but " << (tooFaintCount==1 ? "was" : "were") << " too faint to show up on the plate\n";
	} if (missingPlateCount > 0) {
		cout << missingPlateCount << (matchCount>0 ? " other" : "") << " plate" << (missingPlateCount==1 ? "" : "s");
		cout << " matched, but " << (missingPlateCount==1 ? "isn't" : "aren't") << " in the plate room\n";
	} if (tooFaintCount == 0 && matchCount == 0) {
		printf("Closest approach was %.3f degrees from plate centre\n", closest);
	}
	// printing the total time taken when running the program	
	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	printf("Elapsed time: %.4fs\n", elapsed_seconds.count());
}


// TO DO
	// add option on spreadsheet to add plate size field, plus adjusted position

	// scale magnitude limits based on exposure times (logarithmically?)
	// try to incorporate the errors in RA/DEC from the ephemeris somewhere
		// error propagation through the transformation?
		// even error circle around the point on the plate? 
		// probably distorted?
		// comes from extrapolating objects paths backwards if they were recently discovered, eg Eris/Makemake