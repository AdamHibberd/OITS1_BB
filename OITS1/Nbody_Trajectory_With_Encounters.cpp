#include "stdafx.h"
#include "Nbody_Trajectory_With_Encounters.h"


Nbody_Trajectory_With_Encounters::Nbody_Trajectory_With_Encounters()
{
}

Nbody_Trajectory_With_Encounters::Nbody_Trajectory_With_Encounters(int NB, Body *bodies) : Nbody_Trajectory( NB, bodies)
{
////# Nbody_Trajectory_With_Encounters is class constructor
//#
//# INPUT :
//#
//# bodies : List of Celestial Bodies to be connected(each of Body class)
// #
//# OUTPUT :
//#
//# The Nbody_Trajectory in question, Initialised
//#
	int i;

//	Nbody_Trajectory(NB, bodies);
//	cout << Ntrans << endl;
	Hyperbola = new Connecting_Hyperbola*[NP];
	for (i = 0; i < NP; i++) Hyperbola[i] = new Connecting_Hyperbola[Ntrans];
	HYPERB = new Connecting_Hyperbola*[4];
	for (i = 0; i < 4; i++) HYPERB[i] = new Connecting_Hyperbola[Ntrans];

}

////# Compute_Total_Deltav computes Overall DeltaV for different times for a
////# multi - planet Trajectory with Planetary Encounters
int Nbody_Trajectory_With_Encounters::Compute_Total_Deltav(double *t, int *mode, double thresh, int maxit, int wayflag, int flybyrendez)
// # Compute_Total_Deltav computes Overall DeltaV for different times for a multi - planet Trajectory with Planetary Encounters
//# INPUT:
//#
//# t : Array of times t(1) = Launch Time but t(2)etc are traj durations
//# array of modes : Mode of operation for calculating Ephemeris(see Body class method compute_ephem_at_t)
// # thresh : Threshold for calculating Ephemeris(see Body class method compute_ephem_at_t)
// # maxit : Maximum number of iterations for Calculate_transfer
//# flybyrendez : Flyby or rendezvous with destination planet
//#
//# OUTPUT :
//#
//#  The Nbody_Trajectory_With_Encounters in question, with all the Delta V variables calculated
//#    
{
	int i, j, k;
	double VelPer;

	Nbody_Trajectory::Compute_Total_Deltav(t, mode, thresh, maxit, wayflag, flybyrendez);

	for (i = 0; i < NP; i++) DeltaV[i] = 0.0;

	if (wayflag == 0)
	{
		for (j = 1; j < Ntrans; j++)
		{
			HYPERB[0][j].VA = VA[j][0];
			HYPERB[0][j].VD = VD[j][0];
			HYPERB[0][j].Planet = Body_Set[j];
			//     HYPERB[0][j].Transform();
			HYPERB[0][j].Compute_Hyperbola();


			HYPERB[1][j].VA = VA[j][1];
			HYPERB[1][j].VD = VD[j][0];
			HYPERB[1][j].Planet = Body_Set[j];
			//      HYPERB[1][j].Transform();
			HYPERB[1][j].Compute_Hyperbola();


			HYPERB[2][j].VA = VA[j][0];
			HYPERB[2][j].VD = VD[j][1];
			HYPERB[2][j].Planet = Body_Set[j];
			//       HYPERB[2][j].Transform();
			HYPERB[2][j].Compute_Hyperbola();


			HYPERB[3][j].VA = VA[j][1];
			HYPERB[3][j].VD = VD[j][1];
			HYPERB[3][j].Planet = Body_Set[j];
			//       HYPERB[3][j].Transform();
			HYPERB[3][j].Compute_Hyperbola();
		}

		for (i = 0 ; i < NP ; i++)
		{
			for (j = 0 ; j < Nbody; j++)
			{
				if (j == 0)
				{
					VelPer = VD[j][perm[i][j]].norm();
					dV[i][j] = VelPer;
				}
				else if (j == Nbody-1)
				{
					if (flybyrendez > 0)
						dV[i][j] = VA[j][perm[i][ j - 1]].norm();
					else
						dV[i][ j] = 0;
				}
				else
				{
					if ((perm[i][j] == 0) && (perm[i][j - 1] == 0))
						Hyperbola[i][j] = HYPERB[0][j];
					else if ((perm[i][j] == 0) && (perm[i][j - 1] == 1))
						Hyperbola[i][j] = HYPERB[1][j];
					else if ((perm[i][j] == 1) && (perm[i][j - 1] == 0))
						Hyperbola[i][j] = HYPERB[2][j];
					else if ((perm[i][j] == 1) && (perm[i][j - 1] == 1))
						Hyperbola[i][j] = HYPERB[3][j];
					
					dV[i][j] = abs(Hyperbola[i][j].DV);
										
				}
				if (isnan(dV[i][j]))dV[i][j] = 1e7; // ignore if calculation wrong
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

		//Best = min_element(DeltaV, DeltaV + NP) - DeltaV;
		//BestDeltaV = DeltaV[Best];
		
	}
	else
	{
		i = Best;

		for (j = 0; j<Nbody; j++)
		{
			if (j == 0)
			{
				VelPer = VD[j][ perm[i][j]].norm() ;
				dV[i][j] = VelPer;
			}
			else if (j == Nbody-1)
			{
				if (flybyrendez > 0)
					dV[i][j] = VA[j][ perm[i][j - 1]].norm();
				else
					dV[i][j] = 0;
			}
			else
			{
				Hyperbola[i][j].VA = VA[j][ perm[i][j - 1]];
				Hyperbola[i][j].VD = VD[j][ perm[i][j]];
				Hyperbola[i][j].Planet = Body_Set[j];
		//     Hyperbola[i][j].Transform();
				Hyperbola[i][j].Compute_Hyperbola();
				dV[i][j] = abs(Hyperbola[i][j].DV);
				
			}
			if (isnan(dV[i][j]))dV[i][j] = 1e7; // ignore if calculation wrong
			DeltaV[i] = DeltaV[i] + dV[i][j];
//			cout << DeltaV[i] << endl;
		}

		BestDeltaV = DeltaV[i];
	}

	DIVERGING = 0;
	NO_ENCOUNTER = 0;

	if (Nbody > 2)
	{
		for (j = 1; j< Ntrans; j++)
		{
			if (Hyperbola[Best][j].LOW_PER_ERROR == 1)
//				cout << Body_Set[j - 1].ephemt.r << "\n";
//			cout << 180 / PI_OITS * acos(Body_Set[j - 1].ephemt.r.dot(Body_Set[j].ephemt.r) / Body_Set[j].ephemt.R / Body_Set[j - 1].ephemt.R);
//				cout << "\n\n";
//				   
//				   cout << Body_Set[j].ephemt.r << "\n";
//				   cout << 180 / PI_OITS * acos(Body_Set[j].ephemt.r.dot(Body_Set[j + 1].ephemt.r) / Body_Set[j].ephemt.R / Body_Set[j + 1].ephemt.R);
//					   cout << "\n\n";
//
//				   cout << Body_Set[j+1].ephemt.r << "\n";

				DIVERGING = DIVERGING + (1 << j);

			if (Hyperbola[Best][j].NO_ENCOUNTER == 1)
				NO_ENCOUNTER = NO_ENCOUNTER + (1 << j);
		}
	}
	return 0;
}

	 //# Compute_Total_Deltav

Nbody_Trajectory_With_Encounters::~Nbody_Trajectory_With_Encounters()
{
}
