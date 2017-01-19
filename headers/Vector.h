#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <string>

class Vector {
private:
	double x;
	double y;
	double z;	// all in units of AU = 1.496e11m

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

Vector::Vector(double XX, double YY, double ZZ) : x(XX), y(YY), z(ZZ) { }
Vector::Vector(const Vector& v) : x(v.x), y(v.y), z(v.z) { }
double Vector::getX() { return x; }
double Vector::getY() { return y; }
double Vector::getZ() { return z; }
void Vector::setX(double XX) { x = XX; }
void Vector::setY(double YY) { y = YY; }
void Vector::setZ(double ZZ) { z = ZZ; }

double Vector::magnitude() { 
	return sqrt(x*x + y*y + z*z); 
}

double Vector::magSq() { 
	return x*x + y*y + z*z; 
}

Vector Vector::cross(Vector& a, Vector& b) {
	double xx = a.y*b.z-a.z*b.y;
	double yy = a.z*b.x-a.x*b.z;
	double zz = a.x*b.y-a.y*b.x;
	return {xx, yy, zz};
}

double Vector::dot(Vector& a, Vector& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector operator+(Vector& a, Vector& b) { return { a.x+b.x, a.y+b.y, a.z+b.z }; }
Vector operator-(Vector& a, Vector& b) { return { a.x-b.x, a.y-b.y, a.z-b.z }; }
template <typename T>
Vector operator*(Vector& a, T b) { return { a.x*b, a.y*b, a.z*b }; }
template <typename T>
Vector operator/(Vector& a, T b) { return { a.x/b, a.y/b, a.z/b }; }
bool operator==(Vector& a, Vector& b) {
	double threshold = 1e-10;
	double XX = abs(a.x-b.x);
	double YY = abs(a.y-b.y);
	double ZZ = abs(a.z-b.z);
	return (XX < threshold && YY < threshold && ZZ < threshold);
}

#endif