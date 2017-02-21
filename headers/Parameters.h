#ifndef PARAMETERS_H
#define PARAMETERS_H


// lowest signal to noise ratio of an asteroid before a match is not printed as output
#define SNR_LIMIT 12

// scaling value of the plates, could be 67.14 according to the SuperCOSMOS Fits headers?
#define ARCSECS_PER_MM 67.12

// Millimetre size of the plates, varies between this and 355, probably not that relevant
#define PLATE_SIZE 354.5

// conversions between radians and degrees, taken from wolfram alpha
#define DEG_TO_RAD 0.017453292519943295769236907684886127134428718885417254560971914401710091
#define RAD_TO_DEG 57.29577951308232087679815481410517033240547246656432154916024386120284714


#endif