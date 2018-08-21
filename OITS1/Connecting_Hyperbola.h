#pragma once
#define PI_OITS  (4 * atan(1))
#include <iomanip>
#include "Body.h"
class Connecting_Hyperbola
{
public:
		Vector3d VA = Vector3d::Zero();	// Arrival Velocity
		Vector3d VD = Vector3d::Zero();	// Required Departure Velocity
		double alpha;	// Angle between VD & VA
		double alpha_thresh = 10.0 * PI_OITS / 180.0; // threshold on alpha
		Body Planet;	// Central Body
		Body Probe;		// Body travelling on hyperbola wrt Central Body(Arriving)
		Body Probe2;	// Body travelling on hyperbola wrt Central Body(Departing)
		double Per;		// Periapsis Distance
		double beta;	// Angle in B - Plane, Beta
		double DV;		// DeltaV required at Periapsis
		Matrix3d Trans1;	// Transformation Matrices Based on Arrival Velocity
		Matrix3d Trans2;	// Transformation Matrices Based on Arrival Velocity
		Matrix3d BTrans;	// Overall Transformation Matrix
		int LOW_PER_ERROR; // FLAG for Low Periapsis
		int NO_ENCOUNTER; // FLAG for no encounter dynamics at SSO
		int NIT;		 //Number of Iterations Of fzero 

	Connecting_Hyperbola();
	int Transform_Matrix();
	int Compute_Hyperbola();
	double fzero(double, double, int, int);
	int Orbits_From_Hyperbolas();
	Vector3d  Calculate_Departure_Velocity();
	~Connecting_Hyperbola();
};

