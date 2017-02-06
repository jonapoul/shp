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
	Ephemeris(double day=0.0, Coords c={}, double lst=0.0, double mag=0.0, double dra=0.0, double ddec=0.0)
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
	bool parseEphemerisString(const string& s, const bool isSurfBrt) {
		if (s.length() == 0) return false;
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
		Prints all the data stored in an Ephemeris object. Used for testing
	*/
	void printEphemeris() const {
		printf("JD = %.2f ", m_day);
		printf("RA = %8.4f ", m_coords.getDegRA());
		printf("DEC = %8.4f ", m_coords.getDegDEC());
		printf("LST = %7.4f ", m_lst);
		printf("ApMag = %5.2f ", m_mag);
		printf("dRA = %.4f ", m_dRA);
		printf("dDEC = %.4f\n", m_dDEC);
	}

	/*
		Goes through the ephemeris file to pick out all relevant data, then stores them all in 
		a vector of Ephemeris objects. The $$SOE and $$EOE tags signify the start and end of the 
		data lines.
	*/
	static void readEphemerisFile(vector<Ephemeris>& eph, const string& filename) {
		ifstream ephemerisFile(filename);
		if (ephemerisFile.is_open()) {
			bool canReadEntries = false;
			bool isSurfBrt = false;
			while (!ephemerisFile.eof()) {
				string buffer;
				getline(ephemerisFile, buffer);
				// if the column headers include "S-brt", pass that flag to the parsing function
				// so that it knows to pass over the surface brightness value when reading in
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
		}
		else {
			cout << "Ephemeris file \"" << filename << "\" is not valid\n";
			exit(1);
		}
	}

	/*
		Linearly interpolates a floating point number between two points
		Intended for apparent magnitude estimation, but can be used for anything else
	*/
	static double linInterp(const double mag0, const double t0, const double mag1, const double t1, const double t) {
		return mag0 + (t-t0)*(mag1-mag0)/(t1-t0);
	}

	/*
		Takes an array of ephemeris objects, a number of records to return and an array index to surround.
		Returns an array of numCoords x Coords objects surrounding index i, and the corresponding julian dates.

		Used to find x/y values for polynomial least-squares fitting.
	*/
	static void findNearbyEphs(const vector<Ephemeris>& e, const int numCoords, const int i, vector<Coords>& c, vector<double>& t) {
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

	static void determineParameters(int argc, char* argv[], string& name, int& num, int& power) {
		name = "";
		num = power = 0;
		for (int i = 1; i < argc; i++) {
			string param = string(argv[i]);
			if (param == "help") {
				cout << "\n   OPTIONS:\n";
				cout << "\t-o :\tthe string name of the object you want to load the ephemerides for\n";
				cout << "\t-n :\tthe integer number of surrounding ephemerides to use in the interpolation process (>=2)\n";
				cout << "\t-p :\tthe integer power of polynomial coefficients to generate for the interpolation (>=0)\n";
				cout << "\n   e.g. " << argv[0] << " -o ceres -n 10 -p 3\n";
				cout << "   This gives the matching plates for Ceres, using 10 surrounding ephemeris records, using a cubic fit\n\n";
				exit(1);
			} else if (param == "files") {
				printFiles();
				exit(1);
			} else if (param == "-o") {
				if (i+1 >= argc) {
					cout << "You need another argument after -o. Exiting...\n";
					exit(1);
				}
				name = argv[i+1];
				if (!fs::exists("ephemeris/"+name+".txt")) {
					cout << "ephemeris/" << name << ".txt doesn't exist. Exiting...\n";
					exit(1);
				}
			} else if (param == "-n") {
				if (i+1 >= argc) {
					cout << "You need another argument after -n. Exiting...\n";
					exit(1);
				}
				string numStr = argv[i+1];
				if (numStr.find_last_of("0123456789") == string::npos) {
					cout << "-n must be followed by an integer number. Exiting...\n";
					exit(1);
				}
				num = stoi(numStr);
				if (num < 2) {
					cout << "-n must be larger than 1. Exiting...\n";
					exit(1);
				}
			} else if (param == "-p") {
				if (i+1 >= argc) {
					cout << "You need another argument after -p. Exiting...\n";
					exit(1);
				}
				string powStr = argv[i+1];
				if (powStr.find_last_of("0123456789") == string::npos) {
					cout << "-p must be followed by an integer number. Exiting...\n";
					exit(1);
				}
				power = stoi(powStr);
				if (power < 1) {
					cout << "-p must be larger than 0. Exiting...\n";
					exit(1);
				}
			}
		}
		// default values
		if (name  == "") name  = "ceres";
		if (num   == 0)  num   = 30;
		if (power == 0)  power = 20;

		if (power > num) num = power;
	}

	static void printFiles() {
		fs::path ephPath = "./ephemeris/";
		fs::directory_iterator end;
		vector<string> files;
		for (fs::directory_iterator itr(ephPath); itr != end; itr++) {
			if (is_regular_file(itr->path()))
				files.push_back(itr->path().stem().string());
		}
		cout << "\n   FILES:\n";
		cout << "   The available ephemerides for the -o option are:\n";

		size_t maxLength = 0, size = files.size();
		for (auto f : files)
			maxLength = (f.length() > maxLength) ? f.length() : maxLength;
		for (size_t j = 0; j < files.size(); j += 3) {
			cout << '\t' << files[j] << string(maxLength+2-files[j].length(), ' ');
			if (j+1 < size) 
				cout << '\t' << files[j+1] << string(maxLength+2-files[j+1].length(), ' ');
			if (j+2 < size)
				cout << '\t' << files[j+2] << string(maxLength+2-files[j+2].length(), ' ');
			cout << '\n';
		}
		cout << "\n";
	}

};

#endif