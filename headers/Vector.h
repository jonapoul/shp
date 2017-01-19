#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <string>
#include <sstream>

class Vector {
private:
	double m_x;
	double m_y;
	double m_z;	// all in units of AU = 1.496e11m

public:
	Vector(double XX = 0.f, double YY = 0.f, double ZZ = 0.f);
	Vector(const Vector& v);
	double getX();
	double getY();
	double getZ();
	void setX(double XX);
	void setY(double YY);
	void setZ(double ZZ);
	double magnitude();
	double magSq();
	std::string toString();
	static Vector cross(Vector& a, Vector& b);
	static double dot(Vector &a, Vector& b);

	friend Vector operator+(Vector& a, Vector& b);
	friend Vector operator-(Vector& a, Vector& b);
	template <typename T>
	friend Vector operator*(Vector& a, T b);
	template <typename T>
	friend Vector operator/(Vector& a, T b);
	friend bool operator==(Vector& a, Vector& b);
};

Vector::Vector(double XX, double YY, double ZZ) : m_x(XX), m_y(YY), m_z(ZZ) { }
Vector::Vector(const Vector& v) : m_x(v.m_x), m_y(v.m_y), m_z(v.m_z) { }
double Vector::getX() { return m_x; }
double Vector::getY() { return m_y; }
double Vector::getZ() { return m_z; }
void Vector::setX(double XX) { m_x = XX; }
void Vector::setY(double YY) { m_y = YY; }
void Vector::setZ(double ZZ) { m_z = ZZ; }

double Vector::magnitude() { 
	return sqrt(m_x*m_x + m_y*m_y + m_z*m_z); 
}

double Vector::magSq() { 
	return m_x*m_x + m_y*m_y + m_z*m_z; 
}

std::string Vector::toString() {
	std::stringstream ss;
	ss << '(' << m_x << ", " << m_y << ", " << m_z << ')';
	return ss.str();
}

Vector Vector::cross(Vector& a, Vector& b) {
	double xx = a.m_y*b.m_z - a.m_z*b.m_y;
	double yy = a.m_z*b.m_x - a.m_x*b.m_z;
	double zz = a.m_x*b.m_y - a.m_y*b.m_x;
	return {xx, yy, zz};
}

double Vector::dot(Vector& a, Vector& b) {
	return a.m_x*b.m_x + a.m_y*b.m_y + a.m_z*b.m_z;
}

Vector operator+(Vector& a, Vector& b) { return { a.m_x+b.m_x, a.m_y+b.m_y, a.m_z+b.m_z }; }
Vector operator-(Vector& a, Vector& b) { return { a.m_x-b.m_x, a.m_y-b.m_y, a.m_z-b.m_z }; }
template <typename T>
Vector operator*(Vector& a, T b) { return { a.m_x*b, a.m_y*b, a.m_z*b }; }
template <typename T>
Vector operator/(Vector& a, T b) { return { a.m_x/b, a.m_y/b, a.m_z/b }; }
bool operator==(Vector& a, Vector& b) {
	double threshold = 1e-10;
	double XX = abs(a.m_x - b.m_x);
	double YY = abs(a.m_y - b.m_y);
	double ZZ = abs(a.m_z - b.m_z);
	return (XX < threshold && YY < threshold && ZZ < threshold);
}

#endif