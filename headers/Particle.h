#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <string>
#include "Vector.h"

class Particle {
private:
	std::string name;
	Vector m_pos;	// units of AU			 = 1.496e11m
	Vector m_vel;	// units of AU/year		 = 4740.57m/s
	double mass;		// units of Earth masses = 5.97e24kg
	float radius;		// units of m

public:
	Particle(std::string n="", Vector pos={}, Vector vel={}, double m=0.f, float r=0.f);
	Particle(const Particle& p);
	std::string getName();
	Vector getR();
	Vector getV();
	double getMass();
	float getRadius();

};

Particle::Particle(std::string n, Vector pos, Vector vel, double m, float r)
	: name(n), m_pos(pos), m_vel(vel), mass(m), radius(r) { }
Particle::Particle(const Particle& p) 
	: name(p.name), m_pos(p.m_pos), m_vel(p.m_vel), mass(p.mass), radius(p.radius) { }
std::string Particle::getName() { return name; }
Vector Particle::getR() { return m_pos; }
Vector Particle::getV() { return m_vel; }
double Particle::getMass() { return mass; }
float Particle::getRadius() { return radius; }

#endif