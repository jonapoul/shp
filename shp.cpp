#include <iostream>
#include <chrono>
#include <math.h>
#include <vector>
#include <utility>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using namespace std;

int main(int argc, char* argv[]) {	
	// determining the program parameters from command line arguments
	string objectName;
	Ephemeris::determineParameters(argc, argv, objectName);

	chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();

	// an array of plate IDs that I couldn't find in the archive room
	vector<int> blacklist = {4910, 14587, 15320, 17266};
	// reading all valid plate records
	vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catalog.txt");
	// reading all ephemeris records
	vector<Ephemeris> eph;
	Ephemeris::readEphemerisFile(eph, objectName);

	// defining the first date in the ephemeris file, for later reference
	double firstEphDate  = eph[0].julian();
	double secondEphDate = eph[1].julian();
	double lastEphDate   = eph[eph.size()-1].julian();
	double step = secondEphDate - firstEphDate;

	// sqrt(3.2*3.2 + 3.2*3.2) is the minimum circular radius for the object to possibly be on the plate. This is used as a first check before performing a polynomial fit, to save processing time
	double distanceThreshold = sqrt(2 * 3.2*3.2);
	int matchCount = 0;

	// arrays for holding data to be printed at the end
	vector<unsigned> tooFaint, missingPlate, tooLowSNR;
	vector<Plate> matchPlates;
	vector<Coords> matchCoords;
	vector<int> matchCounts;
	vector<double> matchMags, matchLimMags, matchSNR;
	vector<pair<double,double>> matchStart, matchMid, matchEnd;

	// used to determine the closest angular distance between plate centre and object coordinates
	// only used in case there are no matches, eg mercury/venus/makemake
	double closest = 180.0;

	// through every plate record
	for (int pl = 0; pl < plates.size()-1; pl++) {
		Plate p = plates[pl];
		double plateDate = p.julian();

		// i = index of the ephemeris record immediately before the plate date
		int i = int( (plateDate - firstEphDate) / step );
		// this happens if the first ephemeris record is after the first plate
		if (i < 0) continue;

		Coords before 	  = eph[i].coords();
		double beforeDate = eph[i].julian();
		Coords after 	  = eph[i+1].coords();
		double afterDate  = eph[i+1].julian();

		// approximate linear interpolation
		Coords interpedCoords  = Coords::linInterp(before, beforeDate, after, afterDate, plateDate);

		// calculating the angular distance between the interpolated coordinates and the plate centre
		Coords plateCoords	   = p.coords();
		double angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
		// this is only used if zero plates match, to give the user an idea of how close any of them were
		closest = (angularDistance < closest) ? angularDistance : closest;

		// if the object is within a reasonable angular distance of the plate, do more accurate tests
		if (angularDistance < distanceThreshold) {	
			double xi, eta;
			int status;
			// calculate the x/y coordinates of the interpolated point, with the plate centre as the tangent point
			Coords::gnomonic(interpedCoords, plateCoords, xi, eta, status);
			// if the transformation was valid
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
					if ( p.isMissing(blacklist) ) {
						missingPlate.push_back(p.id());;
						continue;
					}
					// checks whether the plate's magnitude limit is sufficient to spot the object
					double magnitude = Ephemeris::linInterp(eph[i].mag(), beforeDate, eph[i+1].mag(), afterDate, plateDate);
					if (magnitude > p.magLimit()) {
						tooFaint.push_back(p.id());
						continue;
					}
					// calculating the relative signal-to-noise ratio of the plate
					// uses a baseline photon countrate 'n0' and calculates a relative countrate 'n' for this magnitude
					// then S/N = sqrt(n*t)
					double snr = Plate::signalToNoiseRatio(magnitude, p.exposure());
					if (snr < 100) {
						tooLowSNR.push_back(p.id());
						continue;
					}
					// xi/eta coords of the points at the start, middle and end of the exposure
					pair<double,double> start, mid(xi, eta), end;
					Plate::exposureBoundaries(p, {before,after}, {beforeDate,afterDate}, start, end);
					// adding the relevant values to arrays to be printed at the end
					matchCounts.push_back(++matchCount);
					matchPlates.push_back(p);
					matchCoords.push_back(interpedCoords);
					matchMags.push_back(magnitude);
					matchLimMags.push_back(p.magLimit());
					matchStart.push_back(start);
					matchMid.push_back(mid);
					matchEnd.push_back(end);
					matchSNR.push_back(snr);
				}
			}
		}
	}

	// printing a summary of how many plates matched, and how many were too faint/missing
	Plate::printMatches(matchPlates, matchCoords, matchCounts, matchMags, matchLimMags, matchStart, matchMid, matchEnd, matchSNR);
	Plate::printSummary(firstEphDate, lastEphDate, matchCount, objectName);
	Plate::printMissingAndFaint(missingPlate, tooFaint, tooLowSNR, matchCount, closest);

	// printing the total time taken when running the program	
	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	printf("Elapsed time: %.4fs\n", elapsed_seconds.count());
}


// TO DO
	// TEST NEREUS, IT HAS ONE AT 17.75 MAG
		// downloaded a fits file of plate 7233, it's backwards in the x direction but otherwise ok
		// need to load it into gaia to check ra/dec accuracy
	// magnitude scaling using signal to noise over time
	// look at data reduction manual for what i'll be doing after
		// precision 
	// supercosmos catalogue download for images
	// scale magnitude limits based on exposure times (logarithmically?)
	// try to incorporate the errors in RA/DEC from the ephemeris somewhere
		// rough estimate for error region on the plate
		// use this to pipe back more accurate measurements to minor planet centre maybe?
		// error comes from extrapolating object paths backwards if they were recently discovered, eg Eris/Makemake