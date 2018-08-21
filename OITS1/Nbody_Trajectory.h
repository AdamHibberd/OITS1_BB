#pragma once
#include "Body.h"
#include "Transfer_orbit.h"
#include <algorithm>

class Nbody_Trajectory
{
public:
	string Name;			// Name of Nbody_Trajectory Object
	int	Nbody;			    // Number of Celestial Bodies to visit including first and last body
	int	Ntrans;				// Number of transfers
	int	NP;					// Possible Permutations
	int	**perm;				// Array of All Permutations
	Body *Body_Set;			// List of Bodies to Visit
	Transfer_Orbit	*Trans_Set;      // List of Intermediate Transfers
	double	*DeltaV;        // Total Connecting DeltaVs
	double BestDeltaV;		// Best Total DeltaV
	int	Best;				// Best Permutation
	Vector3d **deltaV;		// Array of deltav vectors
	double **dV;			// Array of Magnitudes of deltaV
	Vector3d **VD;			// Array of Departure Velocities
	Vector3d **VA;			// Array of Arrival Velocities


	Nbody_Trajectory();
	Nbody_Trajectory(int, Body*);
	int Initialise_Transfers();
	int Calculate_Nbody_Ephem(int *, double);
	int Calculate_Nbody_Transfers(double, int, int);
	int Compute_Total_Deltav(double *, int *, double , int , int , int);
	~Nbody_Trajectory();
};

