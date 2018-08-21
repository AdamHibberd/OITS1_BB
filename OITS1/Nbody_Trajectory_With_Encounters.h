#pragma once
#include "Connecting_Hyperbola.h"
#include "Nbody_Trajectory.h"
#include <algorithm>

class Nbody_Trajectory_With_Encounters :
	public Nbody_Trajectory
{
public:
		Connecting_Hyperbola **HYPERB;         //# Array of Hyperbolas
		Connecting_Hyperbola **Hyperbola;      //# Connecting Hyperbolas
		int DIVERGING;      //# Interplanetary Trajectory Has unacceptably low Periapsis here
		int NO_ENCOUNTER;   //# Interplanetary Trajectory has no encounter dynamics for this SSO
	
	Nbody_Trajectory_With_Encounters();	
		
	Nbody_Trajectory_With_Encounters(int, Body*);
	int Compute_Total_Deltav(double*, int *, double, int, int, int);
	~Nbody_Trajectory_With_Encounters();
};

