/**
 * This class is used to simulate the orbital movement paths of an arbitrary number of
 * point particles with time, taking into account all gravitational forces between each
 * particle. It uses the Verlet algorithm to update each particles velocity and position
 * after each user-specified timestep.
 * It also tracks the total energy of the entire system, as well as outputting useful orbital
 * information such as apoapsis/periapsis distances, number of completed and partial orbits,
 * eccentricity and orbital period. Period is calculated separately via Kepler's 2nd Law and
 * 3rd Law, for comparison between the methods.
 * This program assumes the particle input data to be in units of earth sidereal years
 * (for time), earth masses (for mass), Astronomical Units (AU) for distance and AU/year 
 * (for velocity) and outputs in the same units.
 * The program is run by typing:
 * "java ParticleManyBody parameters.input particles.input traj.xyz"
 * into the console, with an optional flag of "circle" on the end if you want to run circular
 * orbits.
 *
 * @author J. Poulton
 * @author A. Hussnain
 * @version "3/2016"
 */

import java.io.*;
import java.util.*;

public class ParticleManyBody  { 

    public static void main (String[] argv) throws IOException {

	// Used for comparison in the countdown timer at the end of each loop.
	double timeInitial = System.currentTimeMillis(), timeCurrent, timeEnd, timeTaken;

	// Checking whether the user added the correct arguments to the program
	if (argv.length < 3) {
	    System.err.println("There should be at least three arguments to the program:\n (1) the system parameters\n (2) the initial particle information\n (3) the output trajectory file\n (4) the optional \"circle\" signifier to use only circular orbits");
	    System.exit(0);
	}

	/* Opening scanners to read input data from two files. The first command line argument 
	 * is for the three simulation parameters: "timestep dt" "number of steps n" 
	 * "gravitational constant G". The second argument is for the initial particle states: 
	 * position vector, velocity vector, mass and label, in that order. The third argument
	 * is for the output trajectory file. */
	Scanner parameters = new Scanner(new BufferedReader(new FileReader(argv[0])));
	Scanner particleInfo = new Scanner(new BufferedReader(new FileReader(argv[1])));
	PrintWriter trajectory = new PrintWriter(new FileWriter(argv[2]));
	PrintWriter energy = new PrintWriter(new FileWriter("energy.dat"));

	/* Using a method to count the number of lines in the initial particle details so the 
	 * system can create the correct number of particles in an array. It is important to make 
	 * sure that there aren't any blank lines underneath all the particle information in the
	 * file, since this will lead to an out of bounds exception later on. */
	int N = Particle3D.countLines(argv[1]);
	
	// Null array of particles to be simulated.
	Particle3D[] particles = new Particle3D[N];

	// oldForce[] holds the force vectors from the previous timestep.
	Vector3D[] oldForce = new Vector3D[N];

	// newForce[] holds the new force vectors on each particle after updating.
	Vector3D[] newForce = new Vector3D[N];

	// meanForce[] holds the mean value of the two previous arrays.
	Vector3D[] meanForce = new Vector3D[N];

	for (int i = 0; i < N; i++) {
	    //Initialising the various arrays
	    particles[i] = new Particle3D();
	    oldForce[i] = new Vector3D();
	    meanForce[i] = new Vector3D();
	    newForce[i] = new Vector3D();
	    
	    // Attempts to read particle information from the input file. 
	    try { particles[i].readParticle(particleInfo); }
	    /* If the particle input file was formatted incorrectly, this section will show the 
	     * user what they did wrong and how to correct it. */
	    catch (NoSuchElementException e) {
	    	System.err.println("\nERROR: There is a problem with your "+argv[1]+" file, check if there is a blank or partially-filled line at the bottom of the file. This file must be in the following format for each line (with single word labels for each particle):\n    <positionX1> <positionY1> <positionZ1> <velocityX1> <velocityY1> <velocityZ1> <mass1> <label1>\n    <positionX2> <positionY2> <positionZ2> <velocityX2> <velocityY2> <velocityZ2> <mass2> <label2>\n");
	    	System.exit(0);
	    }
	}

	// Assigning the timestep size dt, the number of time steps n, and gravity constant G. 
	double dt = parameters.nextDouble();
	int n = parameters.nextInt();
	double G = parameters.nextDouble();

	// Closing both input feeds, since we've finished reading in information.
	parameters.close();
	particleInfo.close();

	int[] centres = new int[N];
	for (int i = 0; i < N; i++) {
	    /* Calculating which particle each other particle orbits around, eg. Mercury orbits 
	     * around the Sun instead of Jupiter.
	     * centres[] is an array of integers, where centres[i] is the INDEX (in the main 
	     * particles[] array) of the particle around which we will count the orbits of
	     * particles[i].
	     * In other words, particles[i] orbits around particles[centres[i]]. */
	    centres[i] = Particle3D.orbitCentre(particles, particles[i], G);
	}
	
	// Checking whether any two particles are in the same position.
	for (int i = 0; i < N; i++) {
	    for (int j = 0; j < N; j++) {
		// If the two position vectors are equal, print an error and exit the program.
		if ( i != j && Vector3D.compare( particles[i].getR(), particles[j].getR() ) ) {
		    System.err.println("ERROR: \""+particles[i].getLabel()+"\" and \""+particles[j].getLabel()+"\" are in the same position, please fix this before continuing.");
		    System.exit(0);
	    	}
	    }
	}

	// Defines the initial system time as 0.
	double t = 0;

	// Calculating the initial forces on each particle.
	Particle3D.gravForceArray(particles, G, oldForce);        
	
	// In case the user wants to use circular orbits for each particle
	if (argv.length == 4) {
	    /* The optional 4th argument determines the shape of the orbit, if it's "circle" it 
	     * changes all velocities to circular tangential velocities. Otherwise, the program 
	     * treats the orbits as ellipses and leaves all velocites as they were read in from the
	     * input file. */
	    String shape = argv[3].toLowerCase().trim();
	    if (shape.equals("circle")) { 
		Particle3D.circularVelocity(particles, centres, G); 
	    }

	    // In case the user mispells the 4th argument or types something else:
	    if (!shape.equals("circle")) System.err.println("The fourth argument you gave to the program will only do anything if it is \"circle\". Continuing with simulating elliptical orbits...");
	}

	// Booleans for checking whether the apo/periapses have been found yet.
	boolean[] apoCheck = new boolean[N];  
	boolean[] periCheck = new boolean[N];

 	int count = 0;                         // for printing the correct frame number to file
	String[] apsisTest = new String[N];    // the particle is moving towards apoapsis or periapsis
	double rCurrent;                       // current distance to orbit centre
	double rApo = 0.0;                     // current maximum distance over the orbital cycle
	double rPeri = 0.0;                    // current minimum distance over the orbital cycle
	double[] r0 = new double[N];           // initial orbital radii before the simulation starts
	double[] h = new double[N];	       // |r x v| with respect to the orbit centre
	double[] areaTotal = new double[N];    // total area of the ellipse
	double[] areaSwept = new double[N];    // elliptical area swept out in a timestep dt 
	double[] numTimesteps = new double[N]; // number of timesteps to complete an orbit
	double[] period2nd = new double[N];    // period calculated via Kepler's 2nd law
	double[] period3rd = new double[N];    // period calculated via Kepler's 3rd law
	double[] a = new double[N]; 	       // semi-major axes
	double[] b = new double[N];	       // semi-minor axes
	double[] e = new double[N];            // eccentricities
	Vector3D[] rCentre = new Vector3D[N];  // position of particle wrt its orbit centre
	Vector3D[] vCentre = new Vector3D[N];  // velocity of particle wrt its orbit centre
	Particle3D[] apoapsis = new Particle3D[N];   // the saved particle at maximum distance
	Particle3D[] periapsis = new Particle3D[N];  // the saved particle at minimum distance

	// Initialising the various arrays ready for use
	for (int i = 0; i < N; i++) {
	    if (centres[i] >= 0) {
		r0[i] = Particle3D.separation(particles[i], particles[centres[i]]).magnitude();
		apsisTest[i] = "";
		apoapsis[i] = new Particle3D(particles[i]);
		periapsis[i] = new Particle3D(particles[i]);

		// Set to min and max possible values so they can only be overwritten.
		apoapsis[i].setR( new Vector3D(0.0, 0.0, 0.0) );
		periapsis[i].setR( new Vector3D(Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE) );
		rCentre[i] = new Vector3D();
		vCentre[i] = new Vector3D();
	    }
	}

	// Adjusting the velocities depending on the initial total momentum vector of the system.
	Particle3D.momentumCorrection(particles);
	

	/*****************************************************************************
	 ******************** Beginning of main simulation loop **********************
	 *****************************************************************************/
	for (int i = 0; i < n; i++){

	    for (int j = 0; j < N; j++) {
	    	// Copying the new force vectors to the old force array.
	    	oldForce[j] = new Vector3D(newForce[j]);

		// Reinitialising the newForce array ready to be recalculated.
	    	newForce[j] = new Vector3D();
	    }

	    // Prints every 10 frames to the trajectory file to save on space.
	    if (i % 5 == 0) {
		trajectory.printf("%d\nPoint = %d\n", N, count);
		for (int j = 0; j < N; j++) {
		    trajectory.println(particles[j].toString());
		}
		count++;
	    }

	    // Calculating and printing the total system energy against time to file.
	    if (i > 0) {
		energy.printf("%f %.20f\n", t, Particle3D.totalEnergy(particles, G));
	    }

	    // Updating the positions using current velocity.
	    Particle3D.leapPositionArray(particles, dt, oldForce);

	    // Updating the forces using the new positions.
	    Particle3D.gravForceArray(particles, G, newForce);

	    // Calculating the mean force vectors from the previous and current forces.
	    Particle3D.meanForceArray(oldForce, newForce, meanForce);

            // Updating the velocities using the mean forces.
            Particle3D.leapVelocityArray(particles, dt, meanForce);

	    // Apoapsis/periapsis calculation section
	    for (int j = 0; j < N; j++) {

	    	/* If the particle actually orbits around something (see Particle3D.orbitCentre 
		 * for more explanation) */
	    	if (centres[j] >= 0) {

		    // The positions and velocities relative to the orbit centre
		    rCentre[j] = new Vector3D(Particle3D.separation(particles[centres[j]], particles[j]));
		    vCentre[j] = new Vector3D(Vector3D.subtract(particles[centres[j]].getV(), particles[j].getV()));

		    // Current distance from the particle to its orbit centre.
		    rCurrent = rCentre[j].magnitude();

		    // Maximum saved radial distance so far.
		    rApo = apoapsis[j].getR().magnitude();

		    // Minimum saved radial distance so far.
		    rPeri = periapsis[j].getR().magnitude();

		    /* If the first time-evolved radius is greater than the saved one, the system 
		     * is initially heading towards apoapsis. */
		    if (i == 0 && rCurrent > r0[j]) { apsisTest[j] = "apoapsis"; }

		    /* If the first time-evolved radius is less than the saved one, the system is 
		     * initially heading towards periapsis. */
		    if (i == 0 && rCurrent < r0[j]) { apsisTest[j] = "periapsis"; }

		    // If the particle is known to be heading towards apoapsis:
		    if (i > 0 && apsisTest[j].equals("apoapsis") && !apoCheck[j]) {

			// If the current distance is greater than the previous maximum:
			if (rCurrent > rApo) {

			    // Save the current particle state as the new maximum.
			    apoapsis[j] = new Particle3D(particles[j]);
			    apoapsis[j].setR( rCentre[j] );
			}

			// If the current distance is less than the previous maximum:
			if (rCurrent < rApo) {

			    /* Maximum radial distance has been found, set apoCheck to true to 
			     * signify that apoapsis has been found. */
			    apoCheck[j] = true;
			}
		    }
		    // If particle is hading towards periapsis:
		    if (i > 0 && apsisTest[j].equals("periapsis") && !periCheck[j]) {

			// If the current radial distance is less than the saved minimum:
			if (rCurrent < rPeri) {

			    // Save the current particle state as the new minimum.
			    periapsis[j] = new Particle3D(particles[j]);
			    periapsis[j].setR( rCentre[j] );
			}

			// If the current distance is greater than the saved minimum:
			if (rCurrent > rPeri) {

			    /* Minimum radial distance has been found, set periCheck to true to 
			     * signify that periapsis has been found. */
			    periCheck[j] = true;
			}
		    }

		    /* If either apoapsis or periapsis (but not both) has been found, set the 
		     * apsisTest to the other option so the program knows what to look for next. */
		    if (apoCheck[j] && !periCheck[j]) { apsisTest[j] = "periapsis"; }
		    if (!apoCheck[j] && periCheck[j]) { apsisTest[j] = "apoapsis"; }

		    /* These two if() statements account for any local maxima/minima by checking 
		     * whether the current distance is further away/closer than the saved
		     * apo/periapsis respectively and overwriting with this new value. This is 
		     * especially necessary for Neptune since it has several local maxima and
		     * minima in its orbital distance over time. */
		    if (rCurrent > 1.001*rApo) {
			apoCheck[j] = false;
			apsisTest[j] = "apoapsis";
			apoapsis[j] = new Particle3D(particles[j]);
			apoapsis[j].setR( rCentre[j] );
		    }
		    if (rCurrent < rPeri/1.001) {
			periCheck[j] = false;
			apsisTest[j] = "periapsis";
			periapsis[j] = new Particle3D(particles[j]);
			periapsis[j].setR( rCentre[j] );
		    }
		}
	    }

	    // Incrementing the time
	    t += dt;

	    // Prints a countdown of the approximate time to the end of the main simulation loop
	    if ( i != 0 ) {
		timeCurrent = System.currentTimeMillis();
		timeTaken = timeCurrent - timeInitial;

		// Calculates the average amount of time taken per loop and extrapolates to the end.
		double timePerLoop = timeTaken/(1000*(double)i);
		double stepsLeft = (double)(n-i);
		System.out.printf( "Approx simulation time remaining: %6.1fs   \r", stepsLeft*timePerLoop );
	    }
	}
	/*****************************************************************************
	 *********************** End of main simulation loop *************************
	 *****************************************************************************/


	/**
	 * Section to calculate the orbital periods by two methods, for comparison
	 */
	for (int i = 0; i < N; i++) {

	    // If both the apoapsis and periapsis have been found
	    if (apoCheck[i] && periCheck[i]) {

		// Calculating the semimajor axis a: R_apo + R_peri = a(1+e) + a(1-e) = 2a
		a[i] = ( periapsis[i].getR().magnitude() + apoapsis[i].getR().magnitude() )/2.0;
		
		// Summed mass of particle[j] and the particle it orbits around
		double m = particles[i].getMass() + particles[centres[i]].getMass();
		
		/* 
		 * Calculating the orbital period of each particle via Kepler's 3rd Law:
		 * T^2 = (4 * pi^2 * a^3) / ( G * M ) 
		 */
		double tSquared = Math.pow(a[i], 3) * 4 * Math.PI * Math.PI / (G*m);
		period3rd[i] = Math.sqrt( tSquared );
		
		/* 
		 * Calculating the orbital period via Kepler's 2nd Law 
		 */
		// Cross product of position and velocity vectors wrt the orbit centre
		h[i] = Vector3D.cross(rCentre[i], vCentre[i]).magnitude();

		// Semiminor axis b
		b[i] = Math.sqrt(apoapsis[i].getR().magnitude()*periapsis[i].getR().magnitude());

		// Eccentricity e
		e[i] = Math.sqrt(1 - (b[i]*b[i])/(a[i]*a[i]));

		// Total area of the orbital ellipse
		areaTotal[i] = Math.PI * a[i] * b[i];

		// Swept area during the timestep dt
		areaSwept[i] = h[i] * dt / 2;

		// Number of timesteps to complete an orbit
		numTimesteps[i] = areaTotal[i]/areaSwept[i];

		// Orbital period from Kepler's 2nd Law
		period2nd[i] = numTimesteps[i] * dt;
	    }
	}

	// Outputting out the calculated orbital data of each orbiting particle in the system.
	System.out.printf("Printing data to file...                    ");
	PrintWriter data = new PrintWriter(new FileWriter("orbitaldata.dat"));
	for (int i=0; i<N; i++) {
	    if (centres[i] >= 0) {
		data.printf("%s (orbiting around %s):\n",particles[i].getLabel(), particles[centres[i]].getLabel());
		if (a[i] != 0 && b[i] != 0) {
		    data.printf("Apoapsis at r      = %f AU\n", 				apoapsis[i].getR().magnitude());
		    data.printf("Periapsis at r     = %f AU\n", 				periapsis[i].getR().magnitude());
		    data.printf("Eccentricity e     = %f\n", 					e[i]);
		    data.printf("Orbital period T   = %f years (2nd Law)\n",	period2nd[i]);
		    data.printf("                   = %f years (3rd Law)\n",	period3rd[i]);
		    data.printf("Number of orbits N = %f (2nd Law)\n",			((double)n*dt)/period2nd[i]);
		    data.printf("                   = %f (3rd Law)\n\n",		((double)n*dt)/period3rd[i]);
		}
		else {
		    data.printf("Not enough time was simulated.\n\n");
		}
	    }
	}

        // Closing the various output file streams.
        trajectory.close();
        energy.close();
        data.close();

        timeEnd = System.currentTimeMillis();
	timeTaken = (timeEnd - timeInitial)/1000;

	// Telling the user what has been simulated and which outputs have been created
	System.out.println("\rDone!                                             ");
	System.out.printf(" - "+ N +" particles simulated for %.2f years\n", n*dt);
	System.out.println(" - Trajectory data output to "+argv[2]);
	System.out.println(" - Total energy of the system output to energy.dat");
	System.out.println(" - Orbital counts and other data output to orbitaldata.dat");
	System.out.printf(" - %.6f seconds of total simulation time\n", timeTaken);
    }
}
