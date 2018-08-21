#pragma once
#include "Nbody_Trajectory_With_Encounters.h"
class Mission
{
public:
	
	
		Nbody_Trajectory_With_Encounters Trajectory;    // Trajectory Of Mission
		double *Mission_Times;      // Array of relative Encounter Times
		double *Absolute_Times;      // Array of Absolute Encounter Times
		double *Spice_Min_Times;    // Array of Minimum Times for SPICE Calculation
		double *Spice_Max_Times;    // Array of Maximum Times for SPICE Calculation
		int Out_Of_Spice_Bounds = 0; // Flag to indicate SPICE was unable to be used for a body i
		double TotaldV;            // Total DeltaV of Mission
		int FlybyRendez = 0;      // Flyby or Rendezvous Flag = 0 (FLYBY) = 1 (RENDEZVOUS)
		int wayflag = 1;		// Default to Prograde only
		int *mode;				// Array of NB modes = 2 (Calculated by Spice) = 1 (Otherwise)
		

		Mission();

		Mission(int, Body *, double *, double *, double *, int, int);
		
		double Compute_DeltaV(double *);
		int Set_Absolute_Times(double *);
		int Set_Mission_Times(double *);

		~Mission();
};

