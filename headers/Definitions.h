#ifndef PARAMETERS_H
#define PARAMETERS_H

// lowest signal to noise ratio of an asteroid before a match is not printed as output
#define SNR_LIMIT 12

// scaling value of the plates, could be 67.14 according to the SuperCOSMOS Fits headers?
#define ARCSECS_PER_MM 67.12

// Millimetre size of the plates, varies by around +-0.5mm, probably not that relevant
#define PLATE_SIZE 354.5

// conversions between radians and degrees, taken from wolfram alpha
#define DEG_TO_RAD 0.017453292519943295769236907684886127134428718885417254560971914401710091
#define RAD_TO_DEG 57.29577951308232087679815481410517033240547246656432154916024386120284714

// value to assign to an object for which no apparent magnitude is given in the ephemeris
#define UNKNOWN_MAGNITUDE -100.0

// Longitude of Siding Springs observatory
#define SITE_LONGITUDE 149.070222

// decides whether to print out the final data in a spreadsheet copypastable format at the end
#define PRINT_FOR_SPREADSHEET false

// decides whether to print a list of the steps in the date conversion process
#define PRINT_DATE_BREAKDOWN false
// if the above is true, prints this plate's date conversion process
#define PLATE_TO_PRINT 14591

// used to specify which kind of angle you're working with when getting/setting RA/DEC
enum AngleUnit { 
  DEG, 
  RAD 
};

// used to pass around the x/y coords of start/end points of asteroid
struct Cartesian {
  double x;
  double y;
};

#endif