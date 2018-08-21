#pragma once
#define PI_OITS  (4 * atan(1))
#include "Body.h"

class Transfer_Orbit
{
public:
		int ni[2] = {2 * 0};	//  2 element array corresponding to number of iterations for convergence for both ways
		int nmax;               // Maximum number of iterations
		Body bodyd;             // Departure Planet / Body
		Body bodya;             // Arrival Planet / Body
		Body transfer_body[2];  // Transfer Body - One for short way and the other for long way
		Ephemeris ephemd[2];     // Departure ephemeris - Long way and short way
		Ephemeris ephema[2];     // Arrival ephemeris - Long way and Short way
		double true_anom_dep[2];// True anomaly at Departure - Long way and Short Way
		double true_anom_arr[2];// True anomaly at Arrival - Long way and Short Way
		double perihelion[2];   // Perihelia of Transfer ARC - Long Way and Short way
		double td = 0;          // Departure Time
		double tar;             // Arrival Time
		Vector3d dVd[2] = { 2 * Vector3d::Zero() };  // Delta - V at Departure - for long way and short way
		Vector3d dVa[2] = { 2 * Vector3d::Zero() };  // Delta - V at Arrival - for long way and short way
	Transfer_Orbit();
	int Calculate_Transfer(double, double, double, int, int); // Calculate Transfer Orbits
	double fzero(double, double, double, int);
	int Calculate_Perihelion(); // Calculate closest approach of Transfer body to Sun i.e. Perihelia
	~Transfer_Orbit();
};

