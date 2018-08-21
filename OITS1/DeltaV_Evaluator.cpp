#include "stdafx.h"
#include "DeltaV_Evaluator.h"


/* DeltaV_Evaluator::DeltaV_Evaluator()
{
} */
DeltaV_Evaluator::DeltaV_Evaluator( const NOMAD::Parameters &t ) : NOMAD::Evaluator(t)
{
}
bool DeltaV_Evaluator::eval_x(NOMAD::Eval_Point &x, const NOMAD::Double &h_max, bool & count_eval) const
{

	int i;
	int pointer;	// Array pointer for Maximum Perihelia
	double *time;  // Array of Times for Calculation of DeltaV
	double DeltaV; // DeltaV
	int Best;		// Best Permutation of Transfers
	double *Min_Per_Error = NULL; // Array of Minimum Perisapsis errors
	double *Max_Per_Error = NULL; // Array of Maximum Perisapsis errors
	double *Min_Perihelia_Error = NULL; // Array of Minimum Perihelia Errors
	double Max_Duration_Error;	// Error in MAximum Duration if necessary
	double AA;		// Array of Pointers to Angles to be optimized for each Fixed Point
	Mission *MISS;	// Mission to be optimized


	MISS = &(PROJ->Solution);
	time = new double[MISS->Trajectory.Nbody];

	for (i = 0; i < MISS->Trajectory.Nbody; i++)
	{
		time[i] = x[i].value();
		AA = MISS->Trajectory.Nbody;
		if (MISS->Trajectory.Body_Set[i].Fixed_Point > 0)
		{
			MISS->Trajectory.Body_Set[i].ephem0.r(0) = PROJ->Min_Per[i] * cos(x[AA].value())*cos(x[AA + 1].value());
			MISS->Trajectory.Body_Set[i].ephem0.r(1) = PROJ->Min_Per[i] * sin(x[AA].value())*cos(x[AA + 1].value());
			MISS->Trajectory.Body_Set[i].ephem0.r(2) = PROJ->Min_Per[i] * sin(x[AA + 1].value());
			AA += 2;
		}
	}


	// Compute Objective Function First
	DeltaV = MISS->Compute_DeltaV(time);

	// Best way to go
	Best = MISS->Trajectory.Best;
	//Now Constraints
	if (PROJ->Nconstraints > 0)
	{

		Min_Per_Error = new double[MISS->Trajectory.Ntrans - 1];

		// First Work out array of Minimum Periapsis Constraint Errors
		for (i = 1; i < MISS->Trajectory.Ntrans; i++)
		{
			if ((MISS->Trajectory.NO_ENCOUNTER)&(1 << i))
				Min_Per_Error[i - 1] = 0.0;
			else
				Min_Per_Error[i - 1] = MISS->Trajectory.Hyperbola[Best][i].Planet.radius
				+ PROJ->Min_Per[i]
				- MISS->Trajectory.Hyperbola[Best][i].Per;
			//			cout << MISS->Trajectory.Hyperbola[Best][i].Planet.radius << " " << Min_Per_Error[i - 1] << endl;
		}

		Max_Per_Error = new double[MISS->Trajectory.Ntrans - 1];

		// Secondly the array of Maximum Periapsis Errors
		for (i = 1; i < MISS->Trajectory.Ntrans; i++)
		{
			if ((MISS->Trajectory.NO_ENCOUNTER)&(1 << i))
				Max_Per_Error[i - 1] = 0.0;
			else
			{
				MISS->Trajectory.Body_Set[i].Sphere_Of_Influence();
				Max_Per_Error[i - 1] = MISS->Trajectory.Hyperbola[Best][i].Per - MISS->Trajectory.Body_Set[i].SpoI;
			}
			//		cout << MISS->Trajectory.Hyperbola[Best][i].Per << " " << MISS->Trajectory.Body_Set[i].SpoI << endl;
		}
	}

	// Thirdly Work out Array of Minimum Perihelia Constraints if necessary
	if (PROJ->NPerihelia > 0)
	{
		Min_Perihelia_Error = new double[PROJ->NPerihelia];
		for (i = 0; i < MISS->Trajectory.Ntrans; i++)
		{
			MISS->Trajectory.Trans_Set[i].Calculate_Perihelion();
		}
		for (i = 0; i < PROJ->NPerihelia; i++)
		{
			Min_Perihelia_Error[i] = PROJ->Perihelia[i] - MISS->Trajectory.Trans_Set[PROJ->Per_Pointer[i]].perihelion[MISS->Trajectory.perm[Best][PROJ->Per_Pointer[i]]];
		}
	}

	// Finally work out Maximum Duration Constraint if Necessary
	if (PROJ->Max_Duration < 1e50)
	{
		Max_Duration_Error = MISS->Absolute_Times[MISS->Trajectory.Nbody - 1] - MISS->Absolute_Times[0] - PROJ->Max_Duration;
	}

	// Set up array of Objective and Constraints

	x.set_bb_output(0, NOMAD::Double(DeltaV));	// Objective is DeltaV

	if (PROJ->Nconstraints > 0)
	{
		for (i = 1; i < MISS->Trajectory.Ntrans; i++)
		{
			x.set_bb_output(i, NOMAD::Double(Min_Per_Error[i - 1]));
			x.set_bb_output(i + MISS->Trajectory.Ntrans - 1, NOMAD::Double(Max_Per_Error[i - 1]));
		}
	}
	if (PROJ->NPerihelia > 0)
	{
		for (i = 2 * MISS->Trajectory.Ntrans - 1; i < 2 * MISS->Trajectory.Ntrans + PROJ->NPerihelia - 1; i++)
			x.set_bb_output(i, NOMAD::Double(Min_Perihelia_Error[i + 1 - 2 * MISS->Trajectory.Ntrans]));
	}

	if (PROJ->Max_Duration < 1e50)
		x.set_bb_output(2 * MISS->Trajectory.Ntrans + PROJ->NPerihelia - 1, NOMAD::Double(Max_Duration_Error));
//	cout << DeltaV << "\n";

	count_eval = true;

	return true;

}

DeltaV_Evaluator::~DeltaV_Evaluator()
{
}
