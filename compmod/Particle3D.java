/**
 * Computer Modelling: Exercise 3 
 *
 * A class to define several operations to be done with a new constructor named Particle3D, 
 * for the purpose of simulating a particle under a given potential. This includes setters 
 * and getters for mass, position and name, methods to calculate the particle's kinetic
 * energy, methods to calculate the time evolution of a particle's force (and therefore
 * trajectory), a method to calculate the straight line distance between two particles,
 * and finally a method to output key information about a particle as a String.
 *
 * @author J. Poulton
 * @author A. Hussnain
 * @version "2/2016"
 *
 */

// Imports the Scanner we need for readParticle().
import java.util.Scanner;
import java.io.*;

public class Particle3D {

    /* Basic private variables (position and velocity vectors, mass, label) for the Particle3D 
     * class. */
    private Vector3D r;
    private Vector3D v;
    private double mass;
    private String label;


    /* Basic setter and getter methods to be used on particle properties. */

    /** Sets the mass of a particle.
     * @param m a double to set the mass.
     */
    public void setMass(double m) {mass = m;}

    /** Sets the name given to a particle.
     * @param label a string to set the name of a particle.
     */
    public void setLabel(String name) {label = name;}

    /** Sets the position vector of a particle.
     * @param pos the position vector to be a applied to the particle.
     */
    public void setR(Vector3D pos) {r = new Vector3D(pos);}

    /** Sets the velocity vector of a particle.
     * @param vel the velocity vector to be a applied to the particle.
     */
    public void setV(Vector3D vel) {v = new Vector3D(vel);}

    /** Gets the mass of a particle.
     * @return a double instance representing the mass.
     */
    public double getMass() {return mass;}

    /** Gets the label of a particle.
     * @return a string instance representing the particle label.
     */
    public String getLabel() {return label;}
    
    /** Gets the position vector of a particle.
     * @return a Vector3D representing the position vector of a particle.
     */
    public Vector3D getR() {return r;}

    /** Gets the velocity vector of a particle.
     * @return a Vector3D representing the velocity vector of a particle.
     */
    public Vector3D getV() {return v;}

    /*
     * Constructors
     */

    /** A default constructor that sets all values to 0.
     */
    public Particle3D () {
	r = new Vector3D();
	v = new Vector3D();
	mass = 0; 
	label = "";
    }

    /** A Constructor that sets particle values to the given vectors, double and string.
     * @param pos a vector to represent position at that time.
     * @param vel a vector to represent velocity at that time.
     * @param m a double to represent the particle mass.
     * @param name a string to represent the particle's name.
     */
    public Particle3D (Vector3D pos, Vector3D vel, double m, String name) {
	r = new Vector3D(pos);
	v = new Vector3D(vel);
	mass = m; 
	label = name;
    }

    /** A Constructor that copies an existing particle's values into a new one.
     * @param p the particle to be copied over.
     */
    public Particle3D (Particle3D p) {
	r = new Vector3D(p.r);
	v = new Vector3D(p.v);
        mass = p.mass;
	label = p.label;
    }

    /*
     * Instance Methods
     */
    /** An instance method to convert the particle position and name into a string format.
     * @return a string containing the name and particle position.
     */
    public String toString() {
	return label + " " + r.toString();
    }

    /** An instance method to read in initial particle information from a given input file, then 
     * returns all this as a Particle3D object.
     * @param input a Scanner file from which to read.
     */
    public void readParticle(Scanner input) {
	r.setX( input.nextDouble() );
	r.setY( input.nextDouble() );
	r.setZ( input.nextDouble() );
	v.setX( input.nextDouble() );
	v.setY( input.nextDouble() );
	v.setZ( input.nextDouble() );
	mass = input.nextDouble();
	label = input.next();
    }

    /** An instance method to evolve the velocity of a particle with time according to a force.
     * v = v + (f/m)*dt
     * @param dt a double to define the step to be used in the time evolution process.
     * @param force a Vector3D object to define the force acting on the particle at a particular time.
     */
    public void leapVelocity(double dt, Vector3D force) {
	this.v = Vector3D.add(v, force.multiply(dt/mass));
    }

    /** An instance method to evolve the position of a particle according to a velocity.
     * r = r + v*dt
     * @param dt a double to define the step to be used in the time evolution process.
     */
    public void leapPosition(double dt) {
	this.r = Vector3D.add(r, v.multiply(dt));
    }

    /** An instance method to evolve the position of a particle according to a velocity as well as 
     * an external force.
     * r = r + v*dt + (f/2m)*dt*dt
     * @param dt a double to define the step to be used in the time evolution process.
     * @param force a Vector3D object to define the external force on the particle.
     */
    public void leapPosition(double dt, Vector3D force) {
	Vector3D a = v.multiply(dt);
	Vector3D b = force.multiply(dt*dt*0.5/mass);
	r = Vector3D.add(r, a, b);
    }

    /** An instance method to compute and return the kinetic energy of a particle based on its mass 
     * and velocity vector.
     * E_k = (1/2)*m*v^2
     * @return a double that represents the kinetic energy of the particle.
     */
    public double kineticEnergy() {
	return 0.5 * mass * Vector3D.dot(v, v);
    }

    /** An instance method to compute the gravitational potential energy of a particle with 
     * respect to another particle.
     * E_p = (-G*M*m)/(r)
     * @param p a Particle3D object to represent the first particle.
     * @param G the graviational constant.
     * @return a double to represent the potential energy between the two particles.
     */
    public double potentialEnergy(Particle3D p, double G) {
	double r = Particle3D.separation(this, p).magnitude();
	return -1.0 * G * mass * p.mass / r;
    }

    /** An instance method to compute a gravitational force vector between the given particle 
     * and one other.
     * @param p a Particle3D object to define the central (other) particle.
     * @param G the gravitaional constant.
     * @return a Vector3D object to show the direction and magnitude of the force between the two.
     */
    public Vector3D gravForce(Particle3D p, double G) {

	// Position vector between the two particles
	Vector3D vectorDiff = Vector3D.subtract(this.r, p.r);

	// Unit vector between the two particles
	Vector3D unitDiff = vectorDiff.divide( vectorDiff.magnitude() );

	// Magnitude of the force between the two particles
	double f = -1.0 * G * this.mass * p.mass / vectorDiff.magSquared();
	return unitDiff.multiply(f);
    }
    
    /** An instance method to compute the total energy of a single particle, summing the kinetic
     * energy plus all the potential energies from the other bodies in the simulation.
     * @param particles a Particle3D array of all particles in the system.
     * @param G the gravitational constant.
     * @return the total energy of the particle.
     */
    public double particleEnergy(Particle3D[] particles, double G) {
	double kinetic = 0.0, potential = 0.0;
	for (int i = 0; i < particles.length; i++) {
	    // Filtering out the possibility of calculating an infinite potential
	    if ( !compare(this, particles[i]) ) {
		kinetic += this.kineticEnergy();
		potential += this.potentialEnergy(particles[i], G);
	    }
	}
	return kinetic + potential;	
    }


    /*
     * Static Methods
     */

    /** A static method to compute the relative straight line distance between two Particle3D
     * objects.
     * @param a a Particle3D object to represent the first particle.
     * @param b a Particle3D object to represent the second particle.
     * @return a Vector3D to represent the vector between the two particles.
     */
    public static Vector3D separation(Particle3D a, Particle3D b) {
	return Vector3D.subtract(a.r, b.r);
    }

    /** Static method to compare two particles and determine
     * whether they are equivalent.
     * @param a the first particle to be compared.
     * @param b the second particle to be compared.
     * @return True if the particles are equivalent, False otherwise.
     */
    public static boolean compare(Particle3D a, Particle3D b) {
	boolean pos = Vector3D.compare(a.r,b.r);
	boolean vel = Vector3D.compare(a.v,b.v);

	// If positions, velocities, masses and labels are equal:
	if ( pos && vel && (a.mass == b.mass) && (a.label.equals(b.label)) ) {
	    return true;
	}
	else {
	    return false;
	}
    }

    /** A static method to compute the total energy of a system of particles, summing the kinetic
     * energies plus all the potential energies from the other bodies in the simulation, making
     * sure not to double count the potential energies of particles already taken account of.
     * @param particles a Particle3D array of all particles in the system.
     * @param G the gravitational constant.
     * @return the total energy of the system.
     */
    public static double totalEnergy(Particle3D[] particles, double G) {
	double kinetic = 0.0, potential = 0.0;
	for (int i = 0; i < particles.length; i++) {
	    for (int j = i+1; j < particles.length; j++) {
		potential += particles[i].potentialEnergy(particles[j], G);
	    }
	    kinetic += particles[i].kineticEnergy();
	}
	return kinetic + potential;
    }

    /** A static method to evolve the velocities of an array of particles with time according to a
     * force.
     * v(t+dt) = v(t) + (f/m)*dt
     * @param particles a Particle3D array representing the particles to be updated.
     * @param dt a double to define the step to be used in the time evolution process.
     * @param force a Vector3D array to define the forces acting on the particle at a particular 
     * time. These forces are the averages of f(t) and f(t+dt).
     */
    public static void leapVelocityArray(Particle3D[] particles, double dt, Vector3D[] force) {
	for (int i = 0; i < particles.length; i++) {
	    particles[i].leapVelocity( dt, force[i] );
    	}
    }

    /** A static method to evolve the positions of an array of particles with time according 
     * to a force.
     * r(t+dt) = r(t) + v(t)*dt + (f(t)*dt*dt)/(2*m)
     * @param particles a Particle3D array representing the particles to be updated.
     * @param dt a double to define the step to be used in the time evolution process.
     * @param force a Vector3D array to define the forces acting on the particle at a particular time. These forces are the averages of f(t) and f(t+dt).
     */
    public static void leapPositionArray(Particle3D[] particles, double dt, Vector3D[] force) {
	for (int i = 0; i < particles.length; i++) {
	    particles[i].leapPosition( dt, force[i] );
    	}
    }

    /** A static method to sum the total force vectors acting on each element of an array of
     * particles. 
     * @param particles an array of Particle3D objects to be acted on.
     * @param G the gravitational constant.
     * @param force the array of force vectors to be set to zero, then given new values.
     */
    public static void gravForceArray(Particle3D[] particles, double G, Vector3D[] force) {
	for (int i = 0; i < force.length; i++) {
	    force[i] = new Vector3D();
	    for (int j = 0; j < particles.length; j++) {
		// To avoid adding an infinite force between the particle and itself:
		if (i != j) {
		    force[i] = Vector3D.add(force[i], particles[i].gravForce(particles[j], G)); }
	    }
	}
    }

    /** A static method to calculate the mean force vectors of the pre-timestep and
     * post-timestep force vectors.
     * F_mean = 0.5 * ( F_old + F_new )
     * @param old the force vector array from before the timestep update.
     * @param updated the force vector array from after the timestep update.
     * @param mean the mean vector array of the two previous arrays.
     */
    public static void meanForceArray(Vector3D[] old, Vector3D[] updated, Vector3D[] mean){
	for (int i = 0; i < old.length; i++) {
	    mean[i] = Vector3D.add(old[i], updated[i]).multiply(0.5);
	}
    }

    /** A static method for determining which particle the simulation will count orbits around. For 
     * example, we could count the moon's orbits around the earth or the sun. This is done by
     * calculating the force magnitudes between particle p and every other particle in the system,
     * and if none of them are greater than 1/100th of the largest (probably the sun), the method
     * outputs the array index of the particle with the largest force. If there is one greater than
     * 1/100th, the method checks whether one of the two options orbits around the other, and if so,
     * outputs the index of the former. For example, for the Moon, this method finds that the forces
     * from the Earth and Sun are similar, and finds that Earth orbits around the Sun, so outputs 
     * the index of Earth to the main program since the data from that is more relevant to the user.
     * @param particles the array of all particles in the system.
     * @param p the particle for which we want to find an orbit centre.
     * @param G the gravitational constant.
     * @return the integer index of p's central particle in the array.
     */
    public static int orbitCentre(Particle3D[] particles, Particle3D p, double G) {
        int N = particles.length;
	boolean check = false;

	/* Here we are treating the combination of the force magnitude from one particle and the
	 * array index of that particle as a Vector3D object, since it helpfully combines multiple
	 * doubles to one object, and we can just ignore the Z component. This is so that we can sort
	 * the array by the force magnitudes whilst also keeping track of which particle is causing
	 * each particular force. */
	Vector3D[] f = new Vector3D[N];
	for (int i = 0; i < N; i++) {
	    // Adds the particle index integer to the Y component of the "vector"
	    f[i] = new Vector3D(0, (double)i, 0);

	    /* Returns a negative index number to tell the program that particle p is too massive to
	     * orbit meaningfully around any other particles. This happens when p is more than 1000
	     * times more massive than any other particle in the system. */

	    // If p and particles[i] are not the same particle:
	    if ( !compare(p, particles[i]) ) { 

		// If particle p's mass is less than 1000 times particles[i]'s mass:
		if ( p.getMass() < 1000 * particles[i].getMass() ) {

		    // Tell the system that p is not the sun.
		    check = true;
		}
		// Setting the X component as the i'th particle's force magnitude on p
		f[i].setX( p.gravForce(particles[i], G).magnitude() );
	    }
	}
	
	/* Returns an index outside the array bounds if particle p is too massive. This is so that
	 * the main program knows that the particle doesn't have a proper orbital centre, for
	 * example the Sun doesn't really orbit around any other object.
	 */
	if ( !check ) return -1;

	// Sorting each of the "vectors" by their forces in descending order.
	Vector3D.bubbleSort(f);

	/* Checking whether the force of any particle is >0.01 times the largest force, then counting
	 * how many there are and putting their array indices into the orbitCentres[] array to keep
	 * track of which particles are orbital centre possibilities.
	 */
	int length = 1;
	for (int i = 1; i < N; i++) {
	    if ( f[i].getX() > 0.01*f[0].getX() ) { 
		length++; 
	    }
	}
	int[] orbitCentres = new int[length];
	int j = 0;
	for (int i = 0; i < N; i++) {
	    if ( f[i].getX() > 0.01*f[0].getX() ) {
		orbitCentres[j] = (int)f[i].getY();
		j++;
	    }
	}

	/* If there is more than one possible candidate, returns the index of the particle that has
	 * the smallest mass */
	if (length > 1) {
	    int minIndex = 0;
	    double minMass = Double.MAX_VALUE;
	    for (int i = 0; i < length; i++) {
		if (particles[orbitCentres[i]].getMass() < minMass) {
		    minMass = particles[orbitCentres[i]].getMass();
		    minIndex = orbitCentres[i];
		}
	    }
	    return minIndex;
	}
	// If there is only one candidate, simply returns the array index of that one.
	return orbitCentres[0];
    }

    /** A static method to change the velocities of each particle in the array so that each one 
     * travels in a perfectly circular orbit around its central particle.
     * @param p the array of particles in the system.
     * @param c the array of integers representing the particles that each particle orbits around.
     * @param G the gravitational constant.
     */
    public static void circularVelocity(Particle3D[] p, int[] c, double G) {
	int N = p.length;
	Vector3D rVec = new Vector3D();
	Vector3D vVec = new Vector3D();
	Vector3D planeVector = new Vector3D();
	Vector3D vDirection = new Vector3D();
	Vector3D vUnit = new Vector3D();
	double r, mu, vMag;
	for (int i = 0; i < N; i++) {

	    // If p[i] orbits around anything:
	    if (c[i] < N) {	
		// Initial position vector with respect to the orbit centre
		rVec = Vector3D.subtract(p[i].getR(), p[c[i]].getR());

		// Initial velocity vector with respect to the orbit centre
		vVec = Vector3D.subtract(p[i].getV(), p[c[i]].getV());

		// Vector that defines the plane of the initial (elliptical) orbit
		planeVector = Vector3D.cross(rVec,vVec);

		// Direction vector that the circular orbit will initially follow
		vDirection = Vector3D.cross(planeVector, rVec);

		// Unit vector in that circular direction
		vUnit = vDirection.divide(vDirection.magnitude());

		// Magnitude of distance to orbit centre
		r = rVec.magnitude();

		// Summed mass of the two bodies
		mu = p[i].getMass() + p[c[i]].getMass();

		// Magnitude of tangential velocity
		vMag = Math.sqrt(G*mu/r);

		p[i].setV( Vector3D.add(vUnit.multiply(vMag), p[c[i]].getV()) );
	    }
	}
    }

    /** A static method to adjust the velocities of an array of Particle3D objects based on the total
     * initial momentum vector of every particle in the array.
     * Change in velocity = -(total momentum)/(total mass)
     * @param p an array of all particles in the system.
     */
    public static void momentumCorrection(Particle3D[] p) {
	Vector3D momentum = new Vector3D();
	Vector3D vCOM = new Vector3D();
	Vector3D v = new Vector3D();
	double m = 0.0;
	double massTotal = 0.0;
	int N = p.length;
	for (int i = 0; i < N; i++) {
	    v = p[i].getV();
	    m = p[i].getMass();
	    // Summing up all the momentum vectors
	    momentum = Vector3D.add( momentum, v.multiply(m) );

	    // Summing up all the mass
	    massTotal += m;
	}
	// Dividing the momentum vector by the negative of the summed mass
	vCOM = momentum.divide(massTotal);

	// Applying this adjustment to each velocity vector 
	for (int i=0; i<N; i++) { 
	    p[i].setV(Vector3D.subtract(p[i].getV(), vCOM)); 
	}
    }

    /** A static method to count the number of lines in a given text file. This is so that the
     * simulation can allocate the correct number of spaces in the Particle3D array. This requires
     * that the particle details input file is arranged with each particles details being on a
     * separate line to all the others.
     * @param filename the input file to be counted.
     * @return the number of lines in the given file as an integer.
     */
    public static int countLines(String filename) throws IOException {
	LineNumberReader reader  = new LineNumberReader(new FileReader(filename));
	while (reader.readLine() != null) {}
	reader.close();
	return reader.getLineNumber();
    }

}

