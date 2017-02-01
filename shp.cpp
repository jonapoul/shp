#include <iostream>
#include <chrono>
#include <math.h>
#include <vector>
//#include <boost/filesystem.hpp>
#include <experimental/filesystem>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using namespace std;

void determineParameters(int argc, char* argv[], string& name, int& num, int& power) {
	name = "";
	num = power = 0;
	for (int i = 1; i < argc; i++) {
		string param = string(argv[i]);
		if (param == "help") {
			cout << "\nOPTIONS:\n";
			cout << "\t-p  :\tflags the string name of the object you want to load the ephemerides for\n";
			cout << "\t-t :\tflags the integer number of surrounding ephemerides to use in the interpolation process (>=2)\n";
			cout << "\t-d :\tflags the integer number of polynomial coefficients to generate for the interpolation (>=0)\n";
			cout << "\ne.g. : " << argv[0] << " -n ceres -t 10 -p 3\n";
			cout << "This gives the matching plates for Ceres, using 10 surrounding ephemeris records, using a cubic fit\n\n";
			exit(1);
		} else if (param == "-n") {
			name = argv[i+1];
			if (!experimental::filesystem::exists("ephemeris/"+name+".txt")) {
				cout << "ephemeris/" << name << ".txt doesn't exist. Exiting...\n";
				exit(1);
			}
		} else if (param.substr(0,2) == "-t") {
			num = stoi(argv[i+1]);
		} else if (param.substr(0,2) == "-p") {
			power = stoi(argv[i+1]);
		}
	}

	// default values
	if (name == "") name = "mars";
	if (num == 0)   num = 6;
	if (power == 0) power = 2;
}


int main(int argc, char* argv[]) {	
	chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();

	// reading all valid plate records
	vector<Plate> plates;
	Plate::readPlateCatalog(plates, "catalog.txt");
	vector<Ephemeris> eph;
	
	// determining the program parameters from command line arguments
	int degree, numInterp;
	string objectName;
	determineParameters(argc, argv, objectName, numInterp, degree);

	// reading all ephemeris records
	Ephemeris::readEphemerisFile(eph, "ephemeris/"+objectName+".txt");
	objectName[0] = toupper(objectName[0]);
	

	// defining the first date in the ephemeris file, for later reference
	double firstEphDate = eph[0].julian();
	double secondEphDate = eph[1].julian();
	double step = secondEphDate - firstEphDate;

	// sqrt(3.2*3.2 + 3.2*3.2) is the minimum circular radius for the object to possibly
	// be on the plate. This is used as a first check before performing a polynomial fit
	// to get more accurate, to save processing time
	double distanceThreshold = sqrt(2*3.2*3.2);
	int matchCount = 0, tooFaintCount = 0;

	// through every plate record
	for (int p = 0; p < (int)plates.size()-1; p++) {
		double plateDate  = plates[p].julian();

		// i = index of the ephemeris record immediately before the plate date
		int i 			  = int( (plateDate-firstEphDate)/step );
		Coords before 	  = eph[i].coords();
		double beforeDate = eph[i].julian();
		Coords after 	  = eph[i+1].coords();
		double afterDate  = eph[i+1].julian();

		// approximate linear interpolation
		Coords interpedCoords  = Coords::linInterp(before, beforeDate, after, afterDate, plateDate);

		// calculating the angular distance between the linearly interpolated coordinates and
		// the plate centre
		Coords plateCoords	   = plates[p].coords();
		double angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
		double magnitude 	   = Ephemeris::linInterp(eph[i].mag(), beforeDate, eph[i+1].mag(), afterDate, plateDate);

		// if the object is within a reasonable angular distance of the plate, do more accurate tests
		if (angularDistance < 2*distanceThreshold && magnitude <= plates[p].magLimit()) {
			// interpolating the ephemeris coordinates using the ephemeris records surrounding it
			vector<Coords> nearbyCoords(numInterp);
			vector<double> nearbyTimes(numInterp);
			Ephemeris::findNearbyEphs(eph, numInterp, i, nearbyCoords, nearbyTimes);

			// performing the polynomial least-squares fit using the nearby ephemeris records
			interpedCoords = Coords::polyInterp(nearbyCoords, nearbyTimes, plateDate, degree);

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
						double diameter = Ephemeris::linInterp(eph[i].diam(), beforeDate, eph[i+1].diam(), afterDate, plateDate);
						Plate::printMatch(plates[p], interpedCoords, xi, eta, matchCount+1, magnitude, toXAxis, toYAxis, diameter);
						matchCount++;
					}
				}
			}
		} else if (magnitude > plates[p].magLimit()) {
			tooFaintCount++;
		}
	}

	// printing a summary of how many plates matched
	string firstDate = Plate::julianToGregorian(firstEphDate);
	string lastDate  = Plate::julianToGregorian(eph[eph.size()-1].julian());
	if (matchCount > 0)
		printf("%d matching plates found for %s between %s and %s\n",matchCount,objectName.c_str(),firstDate.c_str(),lastDate.c_str());
	else
		printf("No matches found for %s!\n", objectName.c_str());

	// if any playes matched in terms of coordinates, but the object was too faint to show up,
	// tooFaintCount increments. If this has happened more than once, this bit is printed
	if (tooFaintCount > 0) {
		cout << tooFaintCount << (matchCount>0 ? " other" : "") << " plate" << (tooFaintCount==1 ? "" : "s");
		cout << " matched, but " << (tooFaintCount==1 ? "was" : "were") << " too faint to show up on the plate!\n";
	}
	string sign = "th";
	if 		(degree % 10 == 1) sign = "st";
	else if (degree % 10 == 2) sign = "nd";
	else if (degree % 10 == 3) sign = "rd";
	printf("Used %d surrounding records in a %d%s order polynomial\n", numInterp, degree, sign.c_str());

	// printing the total time taken when running the program	
	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	printf("Elapsed time: %.4fs\n", elapsed_seconds.count());
}


// maybe add a filter for objects with too smallof a diameter? 
	// ask what would be a good threshold
	// add a message at the end to say how many of those were filtered out

// try to incorporate the errors in RA/DEC from the ephemeris somewhere
	// error propagation through the transformation?
	// even error circle around the point on the plate? 
	// probably distorted?