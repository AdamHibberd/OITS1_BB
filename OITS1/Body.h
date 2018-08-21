#pragma once
#define PI_OITS  (4 * atan(1))
#include <string>
#include <cstring>
//#include <cstdio>
//#include <cstdlib>
#include <iostream>
#include "Ephemeris.h"
#include "Orbit.h"



using namespace std;

#include <..\SPICE\cspice\include\SpiceUsr.h>


class Body
{
public:
	char  name[50];		// String identifier for Celestial Body, eg 'Earth'
	char* ID;			// String defining Celestial Body SPICE ID
	double	time;		// Current time(secs)
	double	mu;			// Gravitational Mass of Body(m3 / s2)
	double	GM;         // Gravitational Mass of Attracting Centre(Sun) (m3 / s2)
	double	radius;     // Equatorial radius of Body(m)
	Orbit	orbit;		// Orbit of Celestial Body(Orbit object)
	Ephemeris	ephem0; // Ephemeris of Celestial Body at Epoch(normally t = 0)
	Ephemeris	ephemt;	// Ephemeris of Celestial Body at t = time
	double	state[6];   // State of Body as furnished by SPICE - state(1:3) = position (km)
						//										 state(4:6) = velocity(km / s)
	double	SpoI;       // Sphere of Influence as suggested by Laplace
	int Fixed_Point = 0; //  = 0 Means Body is in orbit, =1 Means Intermediate Point, Otherwise a Fixed Point
	Body();
	int Orbit_To_Ephem(double);
	double Cz(double);
	double Sz(double);
	double dCdz(double);
	double dSdz(double);
	int Compute_Ephem_At_t(double, int, double);
	int Calculate_Orbit_From_Ephem(double);
	int Sphere_Of_Influence();
	int Calculate_True_Anomaly();
	int Sign(double);
	~Body();
};

