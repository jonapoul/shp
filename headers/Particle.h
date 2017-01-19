#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <string>
#include <sstream>
#include "Vector.h"
#include "Coords.h"
using std::cout;

class Particle {
private:
	std::string m_name;
	Vector m_pos;	// units of AU			 = 1.496e11m
	Vector m_vel;	// units of AU/year		 = 4740.57m/s
	double m_mass;	// units of Earth masses = 5.97e24kg
	float m_radius;	// units of m
	Coords m_coord;	// RA/DEC

public:
	Particle(std::string n="", Vector pos={}, Vector vel={}, double m=0.f, float r=0.f);
	Particle(const Particle& p);
	std::string getName();
	Vector getR();
	Vector getV();
	double getMass();
	float getRadius();
	std::string toString();

};

Particle::Particle(std::string n, Vector pos, Vector vel, double m, float r)
	: m_name(n), m_pos(pos), m_vel(vel), m_mass(m), m_radius(r) { }
Particle::Particle(const Particle& p) 
	: m_name(p.m_name), m_pos(p.m_pos), m_vel(p.m_vel), m_mass(p.m_mass), m_radius(p.m_radius) { }
std::string Particle::getName() { return m_name; }
Vector Particle::getR() { return m_pos; }
Vector Particle::getV() { return m_vel; }
double Particle::getMass() { return m_mass; }
float Particle::getRadius() { return m_radius; }

std::string Particle::toString() {
	std::stringstream ss;
	ss << '[' << m_name << ",POS=" << m_pos.toString() << ",VEL=" << m_vel.toString() 
	   << ",RA=(" << m_coord.m_ra << "),DEC=(" << m_coord.m_dec << ")]";
	return ss.str();
}

#endif