#include "stdafx.h"
#include "Nbody_Trajectory.h"


Nbody_Trajectory::Nbody_Trajectory()
{
}
////# Nbody_Trajectory is the class constructor method
Nbody_Trajectory::Nbody_Trajectory(int NB, Body *bodies)
{
	int i,j;
	// # Nbody_Trajectory is the class constructor
	//#
	//# INPUT :
	//#
	//# bodies : List of Celestial Bodies to be connected(each of Body class)
	// #

	Nbody = NB;
	Ntrans = Nbody - 1;
	Body_Set = new Body[Nbody];

	// Initialise set of Encounter Bodies
	for (i = 0; i < Nbody; i++) Body_Set[i] = bodies[i];

	// Initialise set of Transfer Orbits
	Trans_Set = new Transfer_Orbit[Ntrans];
	Initialise_Transfers();

	// Initialise DeltaVs
	deltaV = new Vector3d *[Ntrans];
	for (i = 0; i < Ntrans; i++)
	{
		deltaV[i] = new Vector3d[2];
		for (j = 0; j < 2; j++)deltaV[i][j]= Vector3d::Zero();
	}
	dV = new double*[NP];
	for (i = 0; i < NP; i++)
	{
		dV[i] = new double[Nbody];
		for (j = 0; j < Nbody; j++) dV[i][j] = 0.0;
	}

	DeltaV = new double[NP];
	
	BestDeltaV = 0.0;

	VD = new Vector3d*[Nbody];
	VA = new Vector3d*[Nbody];
	for (i = 0; i < Nbody; i++)
	{
		VD[i] = new Vector3d[2];
		VA[i] = new Vector3d[2];
		for (j = 0; j < 2; j++)
		{
			VD[i][j] = Vector3d::Zero();
			VA[i][j] = Vector3d::Zero();
		}
	}

	return;

}

////# Initialise_Transfers initialises permutation information and transfer orbits
int Nbody_Trajectory::Initialise_Transfers()
{
// # Initialise_Transfers initialises permutation information and transfer orbits
//#
	int i, j, bytes;	// index counters

	// Initialise Number of Permutations
	NP = 1 <<  Ntrans;
	perm = new int*[NP];
	for (i=0 ; i<NP ; i++) perm[i] = new int[Ntrans];
	
	for (bytes = 0; bytes< NP; bytes++)
	{
		for (i = 0; i < Ntrans; i++) perm[bytes][i] = (bytes & ( 1 << i )) / ( 1 << i );
	}


	// Initialise set of Transfer Orbits

	for (i = 0; i < Ntrans; i++)
	{
		Trans_Set[i].bodyd = Body_Set[i];
		Trans_Set[i].bodya = Body_Set[i + 1];
		for (j = 0; j < 2; j++)
		{
			Trans_Set[i].transfer_body[j].ephem0 = Body_Set[i].ephemt;
			Trans_Set[i].transfer_body[j].ephemt = Body_Set[i + 1].ephemt;
		}
	}
	return 0;

} //# Initialise_Transfers


////# Calculate_Nbody_ephem computes list of Ephemeris for each body of the Body_Set
int Nbody_Trajectory::Calculate_Nbody_Ephem(int * mode, double thresh)
// # Calculate_Nbody_ephem computes list of Ephemeris for each body of the Body_Set
//#
//# INPUT:
//#
//# mode : Array of Modes of operation for calculating Ephemeris(see Body class method compute_ephem_at_t)
// # thresh : Threshold for calculating Ephemeris(see Body class method compute_ephem_at_t)
// #
//# OUTPUT :
//#
//# obj : The Nbody_Trajectory in question, with all the Ephemeris of the Bodies calculated
//#
{
	int i;
	for (i = 0; i < Nbody; i++)
		Body_Set[i].Compute_Ephem_At_t(Body_Set[i].ephemt.t, mode[i], thresh);
	return 0;
}

//# Calculate_Nbody_ephem

////# Calculate_Nbody_transfers computes list of Transfer Orbits Connecting each of the Bodies successively
int Nbody_Trajectory::Calculate_Nbody_Transfers(double thresh, int itmax, int wayflag)
// # Calculate_Nbody_transfers computes list of Transfer Orbits Connecting each of the Bodies successively
//#
//# INPUT:
//#
//# thresh : Threshold for calculating Transfer Orbits(see Transfer_orbit class method Calculate_transfer)
// # itmax : Maximum number of iterations for Calculate_transfer
//#
//# OUTPUT :
//#
//#  The Nbody_Trajectory in question, with all the Transfer_orbits calculated
//#
{
	int i;
	double td, tar; // Departure and Arrival Times
	// Compute list of Transfer Orbits

	for (i = 0; i < Ntrans; i++)
	{
		td = Trans_Set[i].bodyd.ephemt.t;
		tar = Trans_Set[i].bodya.ephemt.t;
		Trans_Set[i].Calculate_Transfer(td, tar, thresh, itmax, wayflag);
		
	}
	return 0;
} //# Calculate_Nbody_transfers

////# Compute_Total_Deltav computes Overall DeltaV for different times for a multi - planet Trajectory
int Nbody_Trajectory::Compute_Total_Deltav(double *t, int *mode, double thresh, int maxit, int wayflag, int flybyrendez)
// # Compute_Total_Deltav computes Overall DeltaV for different times for a multi - planet Trajectory
//#
//# INPUT:
//#
//# obj : Current Nbody_Trajectory in question
//# t : Array of times t(1) = Launch Time but t(2)etc are traj durations
//# mode : Array of Mode of operation for calculating Ephemeris(see Body class method compute_ephem_at_t)
// # thresh : Threshold for calculating Ephemeris(see Body class method compute_ephem_at_t)
// # maxit : Maximum number of iterations for Calculate_transfer
//# flybyrendez : Flyby or rendezvous with destination planet
//#
//# OUTPUT :
//#
//#  The Nbody_Trajectory in question, with all the Delta V variables calculated
//#
{
	int i, j, k;
	Vector3d Vdiff;

	for (i = 0; i < Nbody; i++)
		Body_Set[i].ephemt.t = t[i];
	
	
	Calculate_Nbody_Ephem(mode, thresh);

	Initialise_Transfers();

	Calculate_Nbody_Transfers(thresh, maxit, wayflag);

	for (j = 0; j < Nbody; j++)
	{
		for (k = 0; k < 2; k++)
		{
			if (j == 0)
			{
				// Start with Departure Planet
				VA[j][k] = Vector3d::Zero();
				VD[j][k] = Trans_Set[j].transfer_body[k].ephem0.v - Body_Set[j].ephemt.v;
				deltaV[j][k] = VD[j][k];
				continue;
			}

			if (j == Nbody - 1)
			{
				// End With Final Planet
				VD[j][k] = Vector3d::Zero();
				VA[j][k] = Trans_Set[j - 1].transfer_body[k].ephemt.v - Body_Set[j].ephemt.v;
				continue;
			}

			VD[j][k] = Trans_Set[j].transfer_body[k].ephem0.v - Body_Set[j].ephemt.v;
			VA[j][k] = Trans_Set[j - 1].transfer_body[k].ephemt.v - Body_Set[j].ephemt.v;
			deltaV[j][k] = VD[j][k] - VA[j][k];
		}
	}


	for (i = 0; i < NP; i++)
	{
		DeltaV[i] = 0.0;

		for (j = 0; j < Nbody; j++)
		{
			if (j == 0)
				dV[i][j] = VD[j][perm[i][j]].norm();
			else if (j == Nbody - 1)
			{
				if (flybyrendez > 0)
					dV[i][j] = VA[j][perm[i][j - 1]].norm();
				else
					dV[i][j] = 0;

			}
			else
			{
				Vdiff = VD[j][perm[i][j]] - VA[j][perm[i][j - 1]];
				dV[i][j] = Vdiff.norm();
			}
			if (isnan(dV[i][j]))dV[i][j] = 1e50; // ignore if calculation wrong
			DeltaV[i] = DeltaV[i] + dV[i][j];

		}


	}
	/* Find minimum of all DeltaV[i] */

	BestDeltaV = DeltaV[0];
	Best = 0;

	for (i = 1; i < NP; i++)
	{
		if (DeltaV[i] < BestDeltaV)
		{
			BestDeltaV = DeltaV[i];
			Best = i;
		}
	}
//	cout << Best << endl;
	// Best = min_element(DeltaV, DeltaV+NP)- DeltaV;
	// BestDeltaV = DeltaV[Best];

	return 0;
}

Nbody_Trajectory::~Nbody_Trajectory()
{
}
