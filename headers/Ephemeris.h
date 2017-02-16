#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <experimental/filesystem>
#include <boost/algorithm/string.hpp>
#include "Coords.h"

using namespace std;
namespace fs = experimental::filesystem;


class Ephemeris {
private:
	double m_day;		// time in Julian days
	Coords m_coords;	// RA/DEC coordinates
	double m_lst;		// Apparent Local Sidereal Time from the observing point
	double m_mag;		// Apparent magnitude of the object in the sky
	double m_dRA;		// 3sigma error in RA in arcseconds
	double m_dDEC;		// 3sigma error in DEC in arcseconds

public:
	Ephemeris(const double day=0.0, 
	          const Coords c={}, 
	          const double lst=0.0, 
	          const double mag=0.0, 
	          const double dra=0.0, 
	          const double ddec=0.0)
		: m_day(day), m_coords(c), m_lst(lst), m_mag(mag), m_dRA(dra), m_dDEC(ddec) { }
	Ephemeris(const Ephemeris& e)
		: m_day(e.m_day), m_coords(e.m_coords), m_lst(e.m_lst), m_mag(e.m_mag), m_dRA(e.m_dRA), m_dDEC(e.m_dDEC) { }

	inline double julian() const { return m_day; }
	inline Coords coords() const { return m_coords; }
	inline double lst()	const { return m_lst; }
	inline double mag() const { return m_mag; }
	inline double dRA() const { return m_dRA; }
	inline double dDEC() const { return m_dDEC; }
	inline void setJulian(const double d) { m_day = d; }
	inline void setCoords(const Coords& c) { m_coords = c; }
	inline void setLST(const double lst) { m_lst = lst; }
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
		string temp, dRA_str, dDEC_str, lstStr, magStr;
		double ra, dec;
		ss >> m_day >> ra >> dec >> lstStr >> magStr; 
		if (isSurfBrt) ss >> temp;
		ss >> dRA_str >> dDEC_str;

		// converting the two double values (in degrees) to full RA/DEC objects
		m_coords = Coords(ra, dec);

		// checking whether some of the values are valid
		m_mag  = (magStr   == "n.a.") ? 0.0 : stod(magStr);
		m_lst  = (lstStr   == "n.a.") ? 0.0 : stod(lstStr);
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
	                              string& filename) {
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
			cout << filename << ".txt couldn't be found anywhere in /ephemeris/, exiting...\n";
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
				// Start Of Entries, enable parsing of the buffers after this
				if (buffer == "$$SOE") 
					canReadEntries = true;
			}
			ephemerisFile.close();
		} else {
			cout << "Ephemeris file \"" << filename << "\" is not valid\n";
			exit(1);
		}
		filename[0] = toupper(filename[0]);
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
		Takes the argument passed to the program and finds what ephemeris file to load by recursively searching through the ./ephemeris/ folder and testing if the argument exists as a .txt filename
	*/
	static void determineParameters(const int argc, 
	                                char* argv[], 
	                                string& name) {
		name = "";
		if (argc < 2) {
			name = printFiles();
			return;
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
		name = printFiles();
	}

	static string printFiles() {
		vector<string> comets, planets;
		// picking up all filenames from tboth folders
		for (auto& itr : fs::recursive_directory_iterator("./ephemeris/comets")) {
			if (is_regular_file(itr.path())) comets.push_back(itr.path().stem().string());
		}
		for (auto& itr : fs::recursive_directory_iterator("./ephemeris/planets")) {
			if (is_regular_file(itr.path())) planets.push_back(itr.path().stem().string());
		}
		sort(planets.begin(), planets.end(), less<string>());
		sort(comets.begin(),  comets.end(),  less<string>());
		

		// finding the length of the longest filename
		size_t maxLength = 0;
		for (auto p : planets) maxLength = (p.length() > maxLength) ? p.length() : maxLength;
		for (auto c : comets)  maxLength = (c.length() > maxLength) ? c.length() : maxLength;

		// printing all the file options
		cout << "   Planets:\n";
		for (size_t j = 0; j < planets.size(); j += 3) {
			cout << '\t' << planets[j] << string(maxLength+2-planets[j].length(), ' ');
			if (j+1 < planets.size()) 
				cout << '\t' << planets[j+1] << string(maxLength+2-planets[j+1].length(), ' ');
			if (j+2 < planets.size())
				cout << '\t' << planets[j+2] << string(maxLength+2-planets[j+2].length(), ' ');
			cout << '\n';
		}
		cout << "\n   Comets:\n";
		for (size_t j = 0; j < comets.size(); j += 3) {
			cout << '\t' << comets[j] << string(maxLength+2-comets[j].length(), ' ');
			if (j+1 < comets.size()) 
				cout << '\t' << comets[j+1] << string(maxLength+2-comets[j+1].length(), ' ');
			if (j+2 < comets.size())
				cout << '\t' << comets[j+2] << string(maxLength+2-comets[j+2].length(), ' ');
			cout << '\n';
		}

		// taking input from the user
		cout << "\nEnter option: ";
		string output;
		cin >> output;
		for (auto& c : output) c = tolower(c);
		bool choiceIsValid = false;
		while (!choiceIsValid) {
			for (auto p : planets) if (p == output) choiceIsValid = true;
			for (auto c : comets)  if (c == output) choiceIsValid = true;
			if (!choiceIsValid) {
				cout << "That file doesn't exist, try again: ";
				cin >> output;
				for (auto& c : output) c = tolower(c);
			}
		}
		return output;
	}

};

#endif