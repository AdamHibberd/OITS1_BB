#include "stdafx.h"
#include "Mission.h"

Mission::Mission()
{
}

Mission::Mission(int NB, Body *bodies, double *times, double *MinSpice, double *MaxSpice, int flag1, int flag2)
{
	int i;
	double DeltaV;

	FlybyRendez = flag1;
	wayflag = flag2;

//	cout << times[0] << times[1] << endl;
	//char x;
	//cin >> x;
	// Initialise SPICE Min & MAx Times
	Spice_Min_Times = MinSpice;
	Spice_Max_Times = MaxSpice;

	mode = new int[NB];
	
	// Initialise Positions and Osculating Orbits of Various Bodies
	
	for (i = 0; i < NB; i++)
	{
		if (times[i] < Spice_Min_Times[i])
		{
			mode[i] = 1;
			bodies[i].Compute_Ephem_At_t(Spice_Min_Times[i], 2, 1e-4);
			bodies[i].Calculate_Orbit_From_Ephem(Spice_Min_Times[i]);
			bodies[i].ephem0 = bodies[i].ephemt;
		}
		else if (times[i] > Spice_Max_Times[i])
		{
			mode[i] = 1;
			bodies[i].Compute_Ephem_At_t(Spice_Max_Times[i], 2, 1e-4);
			bodies[i].Calculate_Orbit_From_Ephem(Spice_Max_Times[i]);
			bodies[i].ephem0 = bodies[i].ephemt;
		}
		else
		{
			mode[i] = 2;
			bodies[i].Compute_Ephem_At_t(times[i], mode[i], 1e-4);
			bodies[i].Calculate_Orbit_From_Ephem(times[i]);
		}
	}
	
	// Create an object of Type Nbody_Trajectory_With_Encounters
	
	Trajectory = Nbody_Trajectory_With_Encounters(NB, bodies);

	
	Mission_Times = new double[NB] { 0.0 };
	Absolute_Times = new double[NB] { 0.0 };

	Set_Mission_Times(times);

	DeltaV = Compute_DeltaV(Mission_Times);
	TotaldV = DeltaV;

/*	Min_time = new double[NB] {0.0};
	Max_time = new double[NB] {0.0};
	Min_Per = new double[NB] {0.0};
	Perihelia = new double[NB] {0.0};

	for (i = 0; i < Trajectory.Nbody; i++)
	{
		Min_time[i] = Mission_Times[i] - 60 * 60 * 24 * 28;
		Max_time[i] = Mission_Times[i] + 60 * 60 * 24 * 28;
		Min_Per[i] = 0.0;
	} */
}

int Mission::Set_Mission_Times(double *times)
{

	// Function to Initialise Mission_Times
	// Input times is array of absolute times
	// Output has first item as launch time and
	// remaining times as cruise durations
	int i;

	Mission_Times[0] = times[0];

	for (i = 1; i < Trajectory.Nbody; i++)
		Mission_Times[i] = times[i] - times[i - 1];

	return 0;
}

int Mission::Set_Absolute_Times(double *times)
{
	// Reverse of Set_Mission_Times
	int i;
	Absolute_Times[0]= times[0];

	for (i = 1; i < Trajectory.Nbody; i++)
		Absolute_Times[i] = times[i] + Absolute_Times[i - 1];

	return 0;
}

double Mission::Compute_DeltaV(double *times)
{
	int i;
	double DeltaV;

	Set_Absolute_Times(times);
	
	for (i = 0; i < Trajectory.Nbody; i++)
	{
		if (Absolute_Times[i] < Spice_Min_Times[i])
		{
			mode[i] = 1;
			Out_Of_Spice_Bounds = i + 1;
			Trajectory.Body_Set[i].Compute_Ephem_At_t(Spice_Min_Times[i], mode[i], 1e-4);
			Trajectory.Body_Set[i].Calculate_Orbit_From_Ephem(Spice_Min_Times[i]);
			Trajectory.Body_Set[i].ephem0 = Trajectory.Body_Set[i].ephemt;
		}
		else if (Absolute_Times[i] > Spice_Max_Times[i])
		{
			mode[i] = 1;
			Out_Of_Spice_Bounds = i + 1;
			Trajectory.Body_Set[i].Compute_Ephem_At_t(Spice_Max_Times[i], mode[i], 1e-4);
			Trajectory.Body_Set[i].Calculate_Orbit_From_Ephem(Spice_Max_Times[i]);
			Trajectory.Body_Set[i].ephem0 = Trajectory.Body_Set[i].ephemt;
		}
		else
			mode[i] = 2;
	}

	Trajectory.Compute_Total_Deltav(Absolute_Times, mode, 1e-4, 1000, wayflag, FlybyRendez);
	
	DeltaV = Trajectory.BestDeltaV;
	TotaldV = DeltaV;
	return DeltaV;
}


Mission::~Mission()
{
}
