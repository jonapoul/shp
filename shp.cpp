#include <math.h>
#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <utility>
#include "headers/Coords.h"
#include "headers/Ephemeris.h"
#include "headers/Plate.h"
#include "headers/Match.h"
#include "headers/Definitions.h"
namespace chr = std::chrono;

int main(int argc, char* argv[]) {  
  /* determining the asteroid name from command line arguments */
  std::string objectName;
  bool filterSNR = true;
  Ephemeris::determineParameters(argc, argv, objectName, filterSNR);

  chr::time_point<chr::system_clock> start = chr::system_clock::now();

  /* an array of plate IDs that I couldn't find in the archive room, just so i 
     can disregard these plates */
  std::vector<int> missingList = {
    2377, 2380, 4910, 13772, 13944, 14587, 15320, 17266, 18890
  };
  /* reading all valid plate and ephemeris records */
  std::vector<Plate> plates;
  Plate::readPlateCatalog(plates, "catalog.txt");
  std::vector<Ephemeris> eph;
  Ephemeris::readEphemerisFile(eph, objectName);

  /* defining the first date in the ephemeris file, for later reference */
  double firstEphDate  = eph[0].julian();
  double secondEphDate = eph[1].julian();
  double lastEphDate   = eph[eph.size()-1].julian();
  double step = secondEphDate - firstEphDate;

  /* arrays and counters to keep track of how many matches we have */
  unsigned matchCount = 0;
  std::vector<unsigned> tooFaint, missingPlate;
  std::vector<Match> matches;

  /* used to determine the closest angular distance between plate centre and 
     object coordinates, only used in case there are no matches
     eg mercury/venus/makemake */
  double closest = 180.0;

  /* through every plate record */
  for (int pl = 0; pl < plates.size()-1; pl++) {
    Plate p = plates[pl];
    double plateDate = p.julian();

    // i = index of the ephemeris record immediately before the plate date
    int i = int( (plateDate - firstEphDate) / step );
    // this happens if the first ephemeris record is after the first plate
    if (i < 0) continue;

    Coords before     = eph[i].coords();
    double beforeDate = eph[i].julian();
    Coords after      = eph[i+1].coords();
    double afterDate  = eph[i+1].julian();
    Coords interpedCoords = Coords::linInterp(before, beforeDate, after, afterDate, plateDate);
    /* calculating the angular distance between the interpolated coordinates 
       and the plate centre */
    Coords plateCoords     = p.coords();
    double angularDistance = Coords::angularDistance(interpedCoords, plateCoords);
    /* this is only used if zero plates match, to give the user an idea of how
       close any of them were */
    if (angularDistance < closest) 
      closest = angularDistance;

    /* if the object is within a reasonable angular distance of the plate, do 
       more accurate tests */
    double distanceThreshold = sqrt(2 * 3.2*3.2);
    if (angularDistance < distanceThreshold) {  
      double xi, eta;
      /* calculate the x/y coordinates of the interpolated point, with the plate
         centre as the tangent point */
      int status = Coords::gnomonic(interpedCoords, plateCoords, xi, eta);
      /* if the transformation was valid */
      if (status == 0) {
        Coords fromXAxis, fromYAxis;
        /* calculate RA/DEC coordinates of the points that lie on the x/y axes 
           respectively */
        Coords::inverseGnomonic(xi,  0.0, interpedCoords, fromXAxis);
        Coords::inverseGnomonic(0.0, eta, interpedCoords, fromYAxis);
        // calculate the angular distance between the points (in degrees)
        double toXAxis = Coords::angularDistance(interpedCoords, fromXAxis);
        double toYAxis = Coords::angularDistance(interpedCoords, fromYAxis);
        /* if the distances are both less than 3.2 degrees, the point will be 
           on the plate */
        if (toXAxis < 3.2 && toYAxis < 3.2) {
          // checks whether the matched plate is known to be missing
          if ( p.isMissing(missingList) ) {
            missingPlate.push_back(p.id());;
            continue;
          }
          /* checks whether the plate's signal to noise ratio is sufficient to
             spot the object */
          double mag = Ephemeris::linInterp(eph[i].mag(), beforeDate, 
                                            eph[i+1].mag(), afterDate, plateDate);
          double counts = Ephemeris::counts(p.exposure(), mag);
          double snr = counts / p.countLimit();
          if (filterSNR && snr < SNR_LIMIT) {
            tooFaint.push_back(p.id());
            continue;
          }
          std::pair<double,double> start, end;
          std::pair<double,double> mid(Coords::radsToMM(xi), Coords::radsToMM(eta));
          std::pair<Coords,Coords> boundCoords(before,after);
          std::pair<double,double> boundDates(beforeDate, afterDate);
          // getting the start, mid and end plate positions of the asteroid
          Plate::exposureBoundaries(p, boundCoords, boundDates, start, end);
          std::string uncer = Ephemeris::uncertainties(eph[i], eph[i+1], plateDate);
          Match m(matchCount++, p, interpedCoords, mag, start, mid, end, uncer);
          // adding this match to the array, to be printed at the end
          matches.push_back(m);
        }
      }
    }
  }

  // printing a summary of how many plates matched, and how many were too faint/missing
  Match::printMatches(matches);
  Match::printSummary(firstEphDate, lastEphDate, matchCount, objectName);
  Match::printMissingAndFaint(missingPlate, tooFaint, matchCount, closest, filterSNR);
  // printing the total time taken when running the program
  chr::duration<double> elapsed_seconds = chr::system_clock::now() - start;
  printf("Elapsed time: %.4fs\n", elapsed_seconds.count());
}