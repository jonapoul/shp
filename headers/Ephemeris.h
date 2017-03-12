#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <experimental/filesystem>
#include <boost/algorithm/string.hpp>
#include "Coords.h"
#include "Definitions.h"

namespace fs = std::experimental::filesystem;
using std::cout;


class Ephemeris {
private:
	double jd_;			// time in Julian days
	Coords coords_;	// RA/DEC coordinates
	double mag_;		// Apparent magnitude of the object in the sky
	double dra_;		// 3sigma error in RA in arcseconds
	double ddec_;		// 3sigma error in DEC in arcseconds

public:
	Ephemeris(const double day = 0.0, 
	          const Coords c = {}, 
	          const double mag = 0.0, 
	          const double dra = 0.0, 
	          const double ddec = 0.0)
			: jd_(day), coords_(c), mag_(mag), 
				dra_(dra), ddec_(ddec) { }
	Ephemeris(const Ephemeris& e)
			: jd_(e.jd_), coords_(e.coords_), 
				mag_(e.mag_), dra_(e.dra_), ddec_(e.ddec_) { }

	inline double julian() const { return jd_; }
	inline Coords coords() const { return coords_; }
	inline double mag() const { return mag_; }
	inline double dRA() const { return dra_; }
	inline double dDEC() const { return ddec_; }

	/*
		Takes a full ephemeris record string and pulls the info from it, such as julian day, 
		RA/DEC coordinates (plus their errors), LST, apparent magnitude
	*/
	bool parseEphemerisString(const std::string& s, 
	                          const bool isSurfBrt) {
		if (s.length() == 0) 
			return false;
		std::string s2 = s.substr(0,18) + s.substr(21);
		std::stringstream ss(s2);
		std::string temp, dRA_str, dDEC_str, magStr;
		double ra, dec;
		ss >> jd_ >> ra >> dec >> temp >> magStr; 
		if (isSurfBrt) ss >> temp;
		ss >> dRA_str >> dDEC_str;

		// converting the two double values (in degrees) to full RA/DEC objects
		coords_ = Coords(ra, dec);

		// checking whether some of the values are valid
		mag_  = (magStr   == "n.a.") ? UNKNOWN_MAGNITUDE : stod(magStr);
		dra_  = (dRA_str  == "n.a.") ? 0.0 : stod(dRA_str);
		ddec_ = (dDEC_str == "n.a.") ? 0.0 : stod(dDEC_str);
		return true;
	}

	/*
		Goes through the ephemeris file to pick out all relevant data, then stores them all in 
		a vector of Ephemeris objects. The $$SOE and $$EOE tags signify the start and end of the 
		data lines.
	*/
	static void readEphemerisFile(std::vector<Ephemeris>& eph, 
	                              const std::string& filename) {
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
		std::ifstream ephemerisFile(path);
		if (ephemerisFile.is_open()) {
			bool canReadEntries = false;
			bool isSurfBrt = false;
			while (!ephemerisFile.eof()) {
				std::string buffer;
				getline(ephemerisFile, buffer);
				// if the column headers include "S-brt", pass that flag to the parsing function
				// so that it knows to ignore the surface brightness value when reading in
				if (!canReadEntries && buffer.find("S-brt") != std::string::npos)
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
	template<typename T>
	static T linInterp(const T y0, 
										 const T x0, 
	                   const T y1, 
	                   const T x1, 
	                   const T x) {
		return y0 + (x-x0)*(y1-y0)/(x1-x0);
	}

	/*
		Takes an array of ephemeris objects, a number of records to return and an array index to surround.
		Returns an array of numCoords x Coords objects surrounding index i, and the corresponding julian dates.

		Used to find x/y values for polynomial least-squares fitting.
	*/
	static void findNearbyEphs(const std::vector<Ephemeris>& e, 
	                           const int numCoords, 
	                           const int i, 
	                           std::vector<Coords>& c, 
	                           std::vector<double>& t) {
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
	                                std::string& name,
	                                bool& filterSNR) {
		name = "";
		// if not enough 
		if (argc < 2) {
			printFiles(name, filterSNR);
			return;
		}
		// picking up whether the user wants to ignore signal-to-noise ratio filtering
		for (int i = 1; i < argc; i++) {
			if (std::string(argv[i]) == "-snr")
				filterSNR = false; 
		}

		std::string param = std::string(argv[1]);
		for (auto& itr : fs::recursive_directory_iterator("./ephemeris")) {
			fs::path file = itr.path().string() +  '/' + param + ".txt";
			if (fs::exists(file)) {
				name = param;
				return;
			}
		}
		// if the argument doesnt exist as a filename, print out a list and ask which the user wants
		name = std::string(argv[1]);
		printFiles(name, filterSNR);
	}

	static void printFiles(std::string& name, 
	                       bool& filterSNR) {
		// finding the most recently updated ephemeris
		if (name == "new") {
			std::chrono::system_clock::time_point latestTime = std::chrono::system_clock::time_point::min();
			fs::path latestPath;
			std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
			for (auto& itr : fs::recursive_directory_iterator("./ephemeris")) {
				if (is_regular_file(itr.path())) {
					std::chrono::system_clock::time_point time = fs::last_write_time( itr.path() );
					if (time > latestTime) {
						latestTime = time;
						latestPath = itr.path();
					}
				}
			}
			name = latestPath.stem().string();
			return;
		}

		std::vector< std::vector<std::string> > folders;
		std::vector<std::string> folderNames;
		size_t maxLength = 0;
		for (auto& itr : fs::directory_iterator("./ephemeris")) {
			if (is_directory(itr.path())) {
				std::string filename = itr.path().stem().string();
				
				// go through the recently-found folder and add all filename to an array
				fs::path path = itr.path();				
				std::vector<std::string> files;
				for (auto& itr : fs::recursive_directory_iterator(path)) {
					if (is_regular_file(itr.path())) 
						files.push_back(itr.path().stem().string());
				}
				// sorting the files in alphabetical order
				std::sort(files.begin(), files.end(), std::less<std::string>());

				// finding the longest filename length
				for (auto f : files) maxLength = (f.length() > maxLength) ? f.length() : maxLength;

				// Switch the first letter to a capital to it reads a bit better
				filename[0] = toupper(filename[0]);
				folderNames.push_back(filename);
				folders.push_back(files);
			}
		}
		// sorting the folders into alphabetical order
		for (int i = 0; i < folders.size()-1; i++) {
			for (int j = 0; j < folders.size()-i-1; j++) {
				if (folderNames[j] > folderNames[j+1]) {
					std::string tempStr = folderNames[j];
					std::vector<std::string> tempVec = folders[j];

					folderNames[j] = folderNames[j+1];
					folders[j] = folders[j+1];

					folderNames[j+1] = tempStr;
					folders[j+1] = tempVec;
				}
			}
		}

		// printing out an ordered list of all files in each folder under ./ephemeris/
		// yes I know it's a bit messy but it comes out nice
		for (size_t i = 0; name.length() == 0 && i < folders.size(); i++) {
			cout << "   " << folderNames[i] << ":\n";
			for (size_t j = 0; j < folders[i].size(); j += FILES_PER_LINE) {
				cout << '\t';
				for (int k = 0; j+k < folders[i].size() && k < FILES_PER_LINE; k++) {
					cout << folders[i][j+k] << std::string(maxLength+1-folders[i][j+k].length(), ' ');
				}
				cout << '\n';
			}
			cout << '\n';
		}

		// taking input from the user
		if (name.length() == 0) 
			cout << "Enter option: ";
		std::string temp, output;
		while (true) {
			if (name.length() == 0) 
				getline(std::cin, temp);
			else 
				temp = name;
			// getting the year of the asteroid
			std::string year = (temp.length() >= 4) ? temp.substr(0,4) : temp;
			std::vector<std::string> possibilities;
			bool hasNumber = (!isalpha(temp[0]));
			// all to lower case
			for (auto& c : temp) c = tolower(c);
			char firstLetter = temp[0];
			size_t space = temp.find(" ");
			// if the input is two or more words
			if (space != std::string::npos) {
				output = temp.substr(0, space);
				std::string secondWord = temp.substr(space+1);
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
					} else if (hasNumber) {
						std::string fileYear = (temp.length() >= 4) ? file.substr(0,4) : file.substr(0, temp.length());
						if (year == fileYear)
							possibilities.push_back(file);
					} else {
						if (file[0] == firstLetter)
							possibilities.push_back(file);
					}
				}
			}
			// if no matches, ask for another and try again
			if (possibilities.size() > 0) {
				cout << "Did you mean:\n\t";
				size_t longest = 0;
				for (int i = 1; i <= possibilities.size(); i++) {
					cout << possibilities[i-1] << std::string(maxLength-possibilities[i-1].length()+1, ' ');
					if (i % FILES_PER_LINE == 0) cout << "\n\t";
				}
				cout << '\n';
			}
			cout << "Try again: ";
			name = "";
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
	static std::string uncertainties(const Ephemeris& before, 
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
		return std::string(buf);
	}
};

#endif