#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <experimental/filesystem>
#include <boost/algorithm/string.hpp>
#include "Coords.h"
#include "Parameters.h"

using namespace std;
namespace fs = experimental::filesystem;


class Ephemeris {
private:
	double m_day;		// time in Julian days
	Coords m_coords;	// RA/DEC coordinates
	double m_mag;		// Apparent magnitude of the object in the sky
	double m_dRA;		// 3sigma error in RA in arcseconds
	double m_dDEC;		// 3sigma error in DEC in arcseconds

public:
	Ephemeris(const double day=0.0, 
	          const Coords c={}, 
	          const double mag=0.0, 
	          const double dra=0.0, 
	          const double ddec=0.0)
		: m_day(day), m_coords(c), m_mag(mag), m_dRA(dra), m_dDEC(ddec) { }
	Ephemeris(const Ephemeris& e)
		: m_day(e.m_day), m_coords(e.m_coords), m_mag(e.m_mag), m_dRA(e.m_dRA), m_dDEC(e.m_dDEC) { }

	inline double julian() const { return m_day; }
	inline Coords coords() const { return m_coords; }
	inline double mag() const { return m_mag; }
	inline double dRA() const { return m_dRA; }
	inline double dDEC() const { return m_dDEC; }
	inline void setJulian(const double d) { m_day = d; }
	inline void setCoords(const Coords& c) { m_coords = c; }
	inline void setMag(const double mag) { m_mag = mag; }
	inline void setdRA(const double dra) { m_dRA = dra; }
	inline void setdDEC(const double ddec) { m_dDEC = ddec; }

	/*
		Takes a full ephemeris record string and pulls the info from it, such as julian day, 
		RA/DEC coordinates (plus their errors), LST, apparent magnitude
	*/
	bool parseEphemerisString(const string& s, 
	                          const bool isSurfBrt) {
		if (s.length() == 0) 
			return false;
		string s2 = s.substr(0,18) + s.substr(21);
		stringstream ss(s2);
		string temp, dRA_str, dDEC_str, magStr;
		double ra, dec;
		ss >> m_day >> ra >> dec >> temp >> magStr; 
		if (isSurfBrt) ss >> temp;
		ss >> dRA_str >> dDEC_str;

		// converting the two double values (in degrees) to full RA/DEC objects
		m_coords = Coords(ra, dec);

		// checking whether some of the values are valid
		m_mag  = (magStr   == "n.a.") ? UNKNOWN_MAGNITUDE : stod(magStr);
		m_dRA  = (dRA_str  == "n.a.") ? 0.0 : stod(dRA_str);
		m_dDEC = (dDEC_str == "n.a.") ? 0.0 : stod(dDEC_str);
		return true;
	}

	/*
		Goes through the ephemeris file to pick out all relevant data, then stores them all in 
		a vector of Ephemeris objects. The $$SOE and $$EOE tags signify the start and end of the 
		data lines.
	*/
	static void readEphemerisFile(vector<Ephemeris>& eph, 
	                              const string& filename) {
		fs::path path;
		bool fileDoesntExist = true;
		for (auto& itr : fs::recursive_directory_iterator("./ephemeris")) {
			fs::path dir = itr.path();
			fs::path file = dir.string() +  '/' + filename + ".txt";
			if (is_directory(dir) && fs::exists(file)) {
				path = file;
				fileDoesntExist = false;
				break;
			}
		}
		if (fileDoesntExist) {
			cout << filename << ".txt couldn't be found in /ephemeris/, exiting...\n";
			exit(1);
		}
		ifstream ephemerisFile(path);
		if (ephemerisFile.is_open()) {
			bool canReadEntries = false;
			bool isSurfBrt = false;
			while (!ephemerisFile.eof()) {
				string buffer;
				getline(ephemerisFile, buffer);
				// if the column headers include "S-brt", pass that flag to the parsing function
				// so that it knows to ignore the surface brightness value when reading in
				if (!canReadEntries && buffer.find("S-brt") != string::npos)
					isSurfBrt = true;
				// End Of Entries, back out of the loop since all data is done
				if (buffer == "$$EOE")
					break;
				if (canReadEntries) {
					Ephemeris e;
					if (e.parseEphemerisString(buffer, isSurfBrt)) {
						eph.push_back(e);
					}
				}
				// Start Of Entries, enable parsing of all lines after this one
				if (buffer == "$$SOE") 
					canReadEntries = true;
			}
			ephemerisFile.close();
		} else {
			cout << "Ephemeris file \"" << filename << "\" is not valid\n";
			exit(1);
		}
	}

	/*
		Linearly interpolates a floating point number between two points
		Intended for apparent magnitude
	*/
	static double linInterp(const double mag0, 
	                        const double t0, 
	                        const double mag1, 
	                        const double t1, 
	                        const double t) {
		return mag0 + (t-t0)*(mag1-mag0)/(t1-t0);
	}

	/*
		Takes an array of ephemeris objects, a number of records to return and an array index to surround.
		Returns an array of numCoords x Coords objects surrounding index i, and the corresponding julian dates.

		Used to find x/y values for polynomial least-squares fitting.
	*/
	static void findNearbyEphs(const vector<Ephemeris>& e, 
	                           const int numCoords, 
	                           const int i, 
	                           vector<Coords>& c, vector<double>& t) {
		size_t N = e.size();
		if (i > N) {
			cout << "Error in Ephemeris::findNearbyCoords(): index " << i << " > vector size " << N << ".\n";
			exit(1);
		}
		int shift = (numCoords % 2 == 1) ? numCoords/2 : (numCoords-1)/2;
		int diff = (i-shift < 0) ? shift-i : 0;
		for (int j = 0; j < numCoords; j++) {
			c[j] = e[i-shift+j+diff].coords();
			t[j] = e[i-shift+j+diff].julian();
		}
	}

	/*
		Takes the argument passed to the program and finds what ephemeris file to load by recursively searching 
		through the ./ephemeris/ folder and testing if the argument exists as a .txt filename
	*/
	static void determineParameters(const int argc, 
	                                char* argv[], 
	                                string& name,
	                                bool& filterSNR) {
		name = "";
		// if not enough 
		if (argc < 2) {
			printFiles(name, filterSNR);
			return;
		}
		// picking up whether the user wants to ignore signal-to-noise ratio filtering
		for (int i = 1; i < argc; i++) {
			if (string(argv[i]) == "-snr")
				filterSNR = false; 
		}

		string param = string(argv[1]);
		for (auto& itr : fs::recursive_directory_iterator("./ephemeris")) {
			fs::path file = itr.path().string() +  '/' + param + ".txt";
			if (fs::exists(file)) {
				name = param;
				return;
			}
		}
		// if the argument doesnt exist as a filename, default to Ceres and pass a flag back to the main program
		printFiles(name, filterSNR);
	}

	static void printFiles(string& name, bool& filterSNR) {
		vector< vector<string> > folders;
		vector<string> folderNames;
		size_t maxLength = 0;

		for (auto& itr : fs::directory_iterator("./ephemeris")) {
			if (is_directory(itr.path())) {
				string name = itr.path().stem().string();
				
				// go through the recently-found folder and add all filename to an array
				fs::path path = itr.path();				
				vector<string> files;
				for (auto& itr : fs::recursive_directory_iterator(path)) {
						if (is_regular_file(itr.path())) 
						files.push_back(itr.path().stem().string());
				}
				// sorting the files in alphabetical order
				sort(files.begin(), files.end(), less<string>());

				// finding the longest filename length
				for (auto f : files) maxLength = (f.length() > maxLength) ? f.length() : maxLength;

				// Switch the first letter to a capital to it reads a bit better
				name[0] = toupper(name[0]);
				folderNames.push_back(name);
				folders.push_back(files);
			}
		}
		// sorting the folders into alphabetical order
		for (int i = 0; i < folders.size()-1; i++) {
			for (int j = 0; j < folders.size()-i-1; j++) {
				if (folderNames[j] > folderNames[j+1]) {
					string tempStr = folderNames[j];
					vector<string> tempVec = folders[j];

					folderNames[j] = folderNames[j+1];
					folders[j] = folders[j+1];

					folderNames[j+1] = tempStr;
					folders[j+1] = tempVec;
				}
			}
		}

		// printing out an ordered list of all files in each folder under ./ephemeris/
		// yes I know it's a bit messy but it comes out nice
		int FILES_PER_LINE = 7;
		for (size_t i = 0; i < folders.size(); i++) {
			cout << "   " << folderNames[i] << ":\n";
			for (size_t j = 0; j < folders[i].size(); j += FILES_PER_LINE) {
				cout << '\t';
				for (int k = 0; j+k < folders[i].size() && k < FILES_PER_LINE; k++) {
					cout << folders[i][j+k] << string(maxLength+3-folders[i][j+k].length(), ' ');
				}
				cout << '\n';
			}
			cout << '\n';
		}

		// taking input from the user
		cout << "Enter option: ";
		string temp, output;		
		while (true) {
			getline(cin, temp);
			// convert input to all lowercase, for convenience
			for (auto& c : output) c = tolower(c);
			size_t space = temp.find(" ");
			// if the input is two or more words
			if (space != string::npos) {
				output = temp.substr(0, space);
				string secondWord = temp.substr(space+1);
				// if the second word is -snr, set filtering flag to false
				filterSNR = !(secondWord == "-snr");
			} else {
				filterSNR = true;
				output = temp;
			}
			// go through every filename and compare the input word
			for (const auto& folder : folders) {
				for (const auto& file : folder) {
					if (file == output) {
						// if it matches, return it
						name = output;
						return;
					}
				}
			}
			// if no matches, ask for another and try again
			cout << "That filename doesn't exist, try again: ";
		}
	}

	/*
		Calculates how many photon counts would be recieved for a given exposure time (in seconds)
	*/
	static double counts(const double exposure,
	                     const double magnitude) {
		return exposure * pow(10.0, -magnitude / 2.5);
	}

	/*
		Function that outputs the 3 sigma uncertainties in RA and DEC directions, in units of mm
		Interpolates the dRA and dDEC values from two Ephemeris objects and converts them from 
		arcseconds to mm, then formats it approtiately
	*/
	static string uncertainties(const Ephemeris& before, 
	                            const Ephemeris& after, 
	                            const double t) {
		// these two in arcsecs
		double dra  = linInterp(before.dRA(),  before.julian(), after.dRA(),  after.julian(), t);
		double ddec = linInterp(before.dDEC(), before.julian(), after.dDEC(), after.julian(), t);
		// these two in mm
		double dx = dra  / ARCSECS_PER_MM;
		double dy = ddec / ARCSECS_PER_MM;
		// formatting
		char buf[50];
		if (dx > 1e5 || dy > 1e5) {
			sprintf(buf, "(±%.2e, ±%.2e)", dx, dy);
		} else {
			sprintf(buf, "(±%.3f, ±%.3f)", dx, dy);
		}
		return string(buf);
	}
};

#endif