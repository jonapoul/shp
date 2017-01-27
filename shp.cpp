#include <iostream>
#include <chrono>
#include <math.h>
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
using namespace std;

void buildCorners(vector<Coords>& c, double ra, double dec, double dra, double ddec) {
	Coords centre = Coords(ra, dec), buf;								// centre
	c.push_back(centre);
	buf = centre + Coords( 0.0, ddec);									// top right
	buf = buf + Coords( dra/cos(buf.getDEC().getRadians()) , 0.0 );
	c.push_back(buf);
	buf = centre;														// middle right
	buf = buf + Coords( dra/cos(buf.getDEC().getRadians()) , 0.0 );
	c.push_back(buf);
	buf = centre + Coords( 0.0, -ddec );								// bottom right
	buf = buf + Coords( dra/cos(buf.getDEC().getRadians()) , 0.0 );
	c.push_back(buf);
	buf = centre + Coords( 0.0, -ddec );								// bottom middle
	c.push_back(buf);
	buf = centre + Coords( 0.0, -ddec );								// bottom left
	buf = buf + Coords(-dra/cos(buf.getDEC().getRadians()) , 0.0 );
	c.push_back(buf);
	buf = centre;														// middle left
	buf = buf + Coords(-dra/cos(buf.getDEC().getRadians()) , 0.0 );
	c.push_back(buf);
	buf = centre + Coords( 0.0, ddec );									// top left
	buf = buf + Coords(-dra/cos(buf.getDEC().getRadians()) , 0.0 );
	c.push_back(buf);
	buf = centre + Coords( 0.0, ddec );									// top middle
	c.push_back(buf);
}

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
		
		double angularDistance = Coords::angularDistance(interped, plates[p].getCoords());
		if (angularDistance < distanceThreshold) {
			printf("plateID = %5d\n", plates[p].getID());
			printf("\tplateRA  =%10.4f deg\n", plates[p].getRA().getDegrees());
			printf("\tplateDEC =%10.4f deg\n", plates[p].getDEC().getDegrees());
			printf("\tephemRA  =%10.4f deg\n", interped.getRA().getDegrees());
			printf("\tephemDEC =%10.4f deg\n", interped.getDEC().getDegrees());
			printf("\tdistance =%10.4f deg\n", angularDistance);
			printf("\tUTdate   =%10s\n", plates[p].getGregorian().c_str());
			count++;
		}
	}
	if (count > 0)
		cout << count << " matching plates found within " << distanceThreshold << " degrees\n";
	else 
		cout << "No matches found!\n";

	/*double r = 90, d = 60;
	double dAngle = 3.2;
	vector<Coords> coords;
	buildCorners(coords, r, d, dAngle, dAngle);
	vector<string> labels = {"mm", "tr", "mr", "br", "bm", "bl", "ml", "tl", "tm"};
	double x, y;
	int status;
	
	printf("%5s %10s %10s %10s %10s %10s %10s\n", "label", "orig-ra", "orig-dec", "x-coord", "y-coord", "invert-ra", "invert-dec");
	for (int i = 0; i < (int)coords.size(); i++) {
		printf("%5s %10.4f %10.4f ", labels[i].c_str(), coords[i].getRA().getDegrees(), coords[i].getDEC().getDegrees());
		Coords::gnomonic(coords[i], coords[0], x, y, status);
		printf("%10.4f %10.4f ", x*100, y*100);
		Coords buf;
		Coords::inverseGnomonic(x, y, coords[0], buf);
		printf("%10.4f %10.4f\n", buf.getRA().getDegrees(), buf.getDEC().getDegrees());
	}*/

	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - start;
	cout << "\nElapsed time: " << elapsed_seconds.count() << "s\n\n";
}

// take into account the latitude of the observatory??? (-31.27 degrees)
// work out the limiting magnitude for each filter/emulsion, defaulting to 23 or if no matches
// 6.4x6.4 degrees along each plate axis
// make an attempt at quadratic/cubic fitting instead of linear
	// too time-consuming probably?