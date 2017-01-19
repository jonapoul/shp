#include "headers/Plate.h"
#include "headers/Particle.h"
#include <fstream>
#include <vector>
using std::cout;

int main(int argc, char* argv[]) {
	
	std::vector<Plate> plates;
	std::ifstream file("catlog_ukstu.txt");
	if (file.is_open()) {
		while (!file.eof()) {
			std::string buffer;
			getline(file, buffer);
			if (buffer.length() > 0) {
				Plate p;
				p.parsePlateData(buffer);
				plates.push_back(p);
			}
		}
		file.close();
	}

	Particle p{"Earth", {}, {}, 1.f, 6371e3};
	cout << p.toString() << '\n';

	// for (auto p : plates) {
	// 	///////////////////////////////////////////////////////////////
	// 	// Can't print date or lst for some reason?
	// 	/////////////////////////////////////////////////////////////
	// 	cout << p.m_number << ' ' << p.m_date << ' ' << p.m_lst << '\n';
	// }

	
}

// carry on with reading an ephemeris line from the mars file
// also find a good timespace between records, is 4h too much? Jumps by ~2arcmin per record for DEC 
// carry on with Vector and Particle
// check for the word "TEST", so the plate can be discarded. It comes up a lot in the last catalog file
// test gnomonic projection, how fast/efficient is it?