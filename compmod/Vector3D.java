/**
 * A class for vectors, complete with constructors, setters and getters, instance methods to
 * calculate magnitudes, addition and subtraction, cross products, dot products, scalar
 * multiplication, division of vectors, and a method used to compare two vectors. There is also a 
 * method to sort a vector array by their X values.
 *
 * @author J. Poulton
 * @author A. Hussnain
 * @version "10/2015"
 *
 */

import java.text.DecimalFormat;

public class Vector3D {

    /*
     * Private doubles to be used in the Vector3D class.
     */
    private double X, Y, Z;

    /*
     * Setter methods
     */
    /** Sets the x value of a vector.
     * @param xx a double to set the x value.
     */
    public void setX(double xx) { X = xx; }

    /** Sets the y value of a vector.
     * @param yy a double to set the y value.
     */   
    public void setY(double yy) { Y = yy; }

     /** Sets the z value of a vector.
     * @param zz a double to set the z value.
     */
    public void setZ(double zz) { Z = zz; }

    /*
     * Getter methods
     */
    /** Gets the x value of a vector.
     * @return a double instance representing the x value.
     */
    public double getX() { return X; }

    /** Gets the y value of a vector.
     * @return a double instance representing the y value.
     */
    public double getY() { return Y; }

    /** Gets the z value of a vector.
     * @return a double instance representing the z value.
     */
    public double getZ() { return Z; }

    /*
     *  Constructors
     */
    /** A default constructor that sets all elements in a vector to zero.
     */
    public Vector3D () {
	X = 0.0;
	Y = 0.0;
	Z = 0.0;
    }

    /** A constructor that creates a vector with three specified elements.
     * @param xx a double giving an x value to the vector.
     * @param yy a double giving an y value to the vector.
     * @param zz a double giving an z value to the vector.
     */
    public Vector3D (double xx, double yy, double zz) {
	X = xx;
	Y = yy;
	Z = zz;
    }

    /** Implementing a constructor that duplicates a vector.
     * @param original a vector giving a template to copy over to the new vector.
     */
    public Vector3D (Vector3D original) {
	X = original.X;
	Y = original.Y;
	Z = original.Z;
    }

    /*
     * Instance methods.
     */
    /** Instance method that computes the magnitude-squared of a vector.
     * @return a double representing the square of the vector's magnitude.    
     */
    public double magSquared() {
	return (X*X) + (Y*Y) + (Z*Z);
    }

    /** Instance method that returns the magnitude of a vector.
     * @return a double representing the magnitude of the vector.
     */
    public double magnitude() {
	return Math.sqrt(this.magSquared());
    }

    /** Instance method that prints out the vector in the appropriate format as a string.
     * @return a string showing each of the vector components in a (x,y,z) format.
     */
    public String toString() {

	// Shortening each coordinate to a 4-decimal-place number to save on file space.
	DecimalFormat df = new DecimalFormat("0.0000");
	String xx = df.format(X);
        String yy = df.format(Y);
        String zz = df.format(Z);
	return xx + " " + yy + " " + zz;
    }

    /** Instance method for multiplying a vector with a specified scalar.
     * @param a a double to be multiplied by the vector. 
     * @return the original vector multiplied in each component by the scalar a.
     */
    public Vector3D multiply(double a) {
	return new Vector3D(a*X, a*Y, a*Z);
    }

    /** Instance method for dividing a vector with a scalar.
     * @param a a double to divide the vector by.
     * @return the original vector divided by the scalar a in each component.
     */
    public Vector3D divide(double a) {
	return new Vector3D(X/a, Y/a, Z/a);
    }


    /*
     * Static Methods
     */
    /** Static method computing the addition of two vectors.
     * @param a the first vector to be added. 
     * @param b the second vector to be added.
     * @return vector a + vector b.
     */
    public static Vector3D add(Vector3D a, Vector3D b) {
	return new Vector3D(a.X+b.X, a.Y+b.Y, a.Z+b.Z);
    }

    /** Static method computing the addition of three vectors.
     * @param a the first vector to be added. 
     * @param b the second vector to be added.
     * @param c the third vector to be added.
     * @return vector a + vector b + vector c.
     */
    public static Vector3D add(Vector3D a, Vector3D b, Vector3D c) {
	double xx = a.X + b.X + c.X;
	double yy = a.Y + b.Y + c.Y;
	double zz = a.Z + b.Z + c.Z;
	return new Vector3D(xx, yy, zz);
    }

    /** Static method for the subtraction of two vectors, returns a - b.
     * @param a a vector to be subtracted by the latter.
     * @param b a vector to be subtracted from the former.
     * @return the vector a minus the vector b.
     */
    public static Vector3D subtract(Vector3D a, Vector3D b) {
	return new Vector3D(a.X-b.X, a.Y-b.Y, a.Z-b.Z);
   }

    /** Static method for the cross product of two vectors.
     * @param a the first factor.
     * @param b the second factor.
     * @return the vector product of (a x b).
     */
    public static Vector3D cross(Vector3D a, Vector3D b) {
	double xx = a.Y*b.Z-a.Z*b.Y;
	double yy = a.Z*b.X-a.X*b.Z;
	double zz = a.X*b.Y-a.Y*b.X;
	return new Vector3D(xx, yy, zz);
    }

    /** Static method for the scalar product of two vectors.
     * @param a the first vector of the dot product.
     * @param b the second vector in the dot product.
     * @return a double scalar product of a and b.
     */
    public static double dot(Vector3D a, Vector3D b) {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z;
    }

    /** Static method to compare two vectors and determine
     * whether they are equivalent.
     * @param a the first vector to be compared.
     * @param b the second vector to be compared.
     * @return True if the vectors are equivalent, False otherwise.
     */
    public static boolean compare(Vector3D a, Vector3D b) {
	if (a.X==b.X && a.Y==b.Y && a.Z==b.Z) {
	    return true;
	}
	else {
	    return false;
	}
    }

    /** A static method for sorting a Vector3D array by their X values in descending order.
     * @param f a Vector3D array to be sorted.
     */
    public static void bubbleSort(Vector3D[] f) {
	boolean flag = true;  
	Vector3D temp;   
	while (flag) {
	    flag = false;
            for(int i = 0; i < f.length-1; i++) {
		if (f[i].getX() < f[i+1].getX()) {
		    temp = new Vector3D(f[i]);
		    f[i] = new Vector3D(f[i+1]);
		    f[i+1] = new Vector3D(temp);
		    flag = true;
		} 
            } 
	}
    }


}
