// OITS1.cpp : Defines the entry point for the console application.
//

#define PI_OITS  3.141592654
#define ASTRO_UNIT 149597870700    // Astronomical Unit in Metres
#define NBODIES 4
#include "stdafx.h"

#include <iostream>
// #include <iomanip>
#include <fstream>
// #include <cmath>
// #include <cstdlib>
// #include <ctime>
// #include <sstream>

#include "Ephemeris.h"
#include "Orbit.h"
#include "Transfer_Orbit.h"
#include "Nbody_Trajectory.h"
#include "Connecting_Hyperbola.h"
#include "Project.h"

using namespace std;


int main(int argc,char **argv)
{
	int i, j, N;
	char y;
	double *times, *abstimes;
	double x[100];

	ifstream in (argv[1]);
	
	i = 0;
	while ((!in.eof())&&i<100)
	{
		in >> x[i];
		i++;
		//		cout << x[i] << " ";
	}
	N = i;
	in.close();

	Project Test;
	Test.Read_Data("inputfile.txt");


//	x = new double[N];

//	for (i = 0; i < N; i++) {
//		in >> x[i];
//		cout << x[i] << " ";
//	}
	in.close();

	Test.Initialize_SPICE();
	
	Test.Get_SPICE_List();

	for (i=0;i<Test.Body_Number;i++)
			if (strcmp(Test.Body_Chosen[i].ID, "IP")==0)Test.Add_Intermediate_Point(i);	
	
	Test.Merge_Data();



	times = new double[Test.Body_Number];
	abstimes = new double[Test.Body_Number];

	for (i = 0; i < Test.Body_Number; i++)
	{
		times[i] = x[i];
		abstimes[i] = times[i];
		if (i > 0)
			abstimes[i] = abstimes[i - 1] + times[i];
	}
	Test.Nconstraints = Test.Body_Number-2;

//	cout << Test.Nconstraints;

//	cout << "HEre";
//	cin >> x;
//	Test.Current_Mission = Mission::Mission(Test.Body_Number, Test.Body_Chosen, abstimes, Test.Min_Spice_Select, Test.Max_Spice_Select, Test.FlybyRendez, Test.wayflag);

//	Test.Initialize_Mission(times, Test.FlybyRendez, Test.wayflag);

//	cout << "HEre2";
//	cin >> x;


//	Test.Solution = Test.Current_Mission;
//	Test.Optimize_Mission();
//	Test.Current_Mission = Test.Solution;

	//	Test.Current_Mission = NewMiss;

	int pointer;	// Array pointer for Maximum Perihelia
	double *time;  // Array of Times for Calculation of DeltaV
	double DeltaV; // DeltaV
	int Best;		// Best Permutation of Transfers
	double *Min_Per_Error = NULL; // Array of Minimum Perisapsis errors
	double *Max_Per_Error = NULL; // Array of Maximum Perisapsis errors
	double *Min_Perihelia_Error = NULL; // Array of Minimum Perihelia Errors
	double Max_Duration_Error;	// Error in MAximum Duration if necessary
	int AA;		// Array of Pointers to Angles to be optimized for each Fixed Point
	Mission *MISS;	// Mission to be optimized

//	Test.Current_Mission = Mission::Mission(Test.Body_Number, Test.Body_Chosen, abstimes, Test.Min_Spice_Select, Test.Max_Spice_Select, Test.FlybyRendez, Test.wayflag);
//	MISS = &(Test.Current_Mission);
	
	for (i = 0; i < Test.Body_Number; i++)
	{
		AA = Test.Body_Number;
		if (Test.Body_Chosen[i].Fixed_Point > 0)
		{
			Test.Body_Chosen[i].ephem0.r(0) = Test.Min_Per[i] / 1000 * Test.AU * cos(x[AA])*cos(x[AA + 1]);
			Test.Body_Chosen[i].ephem0.r(1) = Test.Min_Per[i] / 1000 * Test.AU * sin(x[AA])*cos(x[AA + 1]);
			Test.Body_Chosen[i].ephem0.r(2) = Test.Min_Per[i] / 1000 * Test.AU * sin(x[AA + 1]);
			AA += 2;
		}
	}
	Test.Current_Mission = Mission::Mission(Test.Body_Number, Test.Body_Chosen, abstimes, Test.Min_Spice_Select, Test.Max_Spice_Select, Test.FlybyRendez, Test.wayflag);
	MISS = &(Test.Current_Mission);
	// Compute Objective Function First
//	DeltaV = MISS->Compute_DeltaV(times);
	DeltaV = MISS->TotaldV;
	// Best way to go
	Best = MISS->Trajectory.Best;
	//Now Constraints
	if (Test.Nconstraints > 0)
	{

		Min_Per_Error = new double[MISS->Trajectory.Ntrans - 1];

		// First Work out array of Minimum Periapsis Constraint Errors
		for (i = 1; i < MISS->Trajectory.Ntrans; i++)
		{
			if ((MISS->Trajectory.NO_ENCOUNTER)&(1 << i))
				Min_Per_Error[i - 1] = 0.0;
			else
				Min_Per_Error[i - 1] = MISS->Trajectory.Hyperbola[Best][i].Planet.radius
				+ Test.Min_Per[i]
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
	if (Test.NPerihelia > 0)
	{
		Min_Perihelia_Error = new double[Test.NPerihelia];
		for (i = 0; i < MISS->Trajectory.Ntrans; i++)
		{
			MISS->Trajectory.Trans_Set[i].Calculate_Perihelion();
		}
		for (i = 0; i < Test.NPerihelia; i++)
		{
			Min_Perihelia_Error[i] = Test.Perihelia[i] - MISS->Trajectory.Trans_Set[Test.Per_Pointer[i]].perihelion[MISS->Trajectory.perm[Best][Test.Per_Pointer[i]]];
		}
	}

	// Finally work out Maximum Duration Constraint if Necessary
	if (Test.Max_Duration < 1e50)
	{
		Max_Duration_Error = MISS->Absolute_Times[MISS->Trajectory.Nbody - 1] - MISS->Absolute_Times[0] - Test.Max_Duration;
	}

	// Set up array of Objective and Constraints

	// x.set_bb_output(0, NOMAD::Double(DeltaV));	// Objective is DeltaV

	int Duration_Flag = (Test.Max_Duration < 1e50);
	
	int Eval_Number = 1 + Test.Nconstraints * 2 + Test.NPerihelia + Duration_Flag;
	double *g;
	g = new double[Eval_Number];

	g[0] = DeltaV;

	if (Test.Nconstraints > 0)
	{
		for (i = 1; i < MISS->Trajectory.Ntrans; i++)
		{
			g[i]=Min_Per_Error[i - 1];
			g[i + MISS->Trajectory.Ntrans - 1]= Max_Per_Error[i - 1];
		}
	}
	if (Test.NPerihelia > 0)
	{
		for (i = 2 * MISS->Trajectory.Ntrans - 1; i < 2 * MISS->Trajectory.Ntrans + Test.NPerihelia - 1; i++)
			g[i] = Min_Perihelia_Error[i + 1 - 2 * MISS->Trajectory.Ntrans];
	}

	if (Test.Max_Duration < 1e50)
		g[2 * MISS->Trajectory.Ntrans + Test.NPerihelia - 1] = Max_Duration_Error;
	
	cout << std::fixed;
	cout << setprecision(16);
	
	for (i = 0; i < Eval_Number; i++)
	{
		cout << g[i] << " ";
	}


	return 0;


/*
	// JUICE
	Test.Body_Number = 7;
	double times[7] = { 7.073136691849103e+08 ,7.386660751349104e+08 ,7.512548871749103e+08 ,7.784238261749104e+08, 7.924849334749102e+08, 8.488985106749103e+08,9.420301846749103e+08 };
	Test.Body_List[1].mu = 4.869e24*(6.67259e-11);		//Venus
	Test.Body_List[2].mu = 3.98602e14;					//Earth
	Test.Body_List[3].mu = 0.6419e24*(6.67259e-11);		//Mars
	Test.Body_List[4].mu = 1898.6e24*(6.67259e-11);     //Jupiter
	Test.Body_List[5].mu = 568.46e24*(6.67259e-11);		//Saturn
	Test.Body_List[1].radius = 6052000;
	Test.Body_List[2].radius = 6378135;
	Test.Body_List[3].radius = 3397000;
	Test.Body_List[4].radius = 71492000;
	Test.Body_List[5].radius = 60268000;

	Test.Body_Chosen = new Body[Test.Body_Number];
	Test.Body_Chosen[0] = Test.Body_List[2];
	cout << Test.Body_Chosen[0].name<< endl;
	
	Test.Body_Chosen[1] = Test.Body_List[2];
	cout << Test.Body_Chosen[1].name<< endl;

	Test.Body_Chosen[2] = Test.Body_List[1];
	cout << Test.Body_Chosen[2].name << endl;

	Test.Body_Chosen[3] = Test.Body_List[2];
	cout << Test.Body_Chosen[3].name << endl;

	Test.Body_Chosen[4] = Test.Body_List[3];
	cout << Test.Body_Chosen[4].name << endl;

	Test.Body_Chosen[5] = Test.Body_List[2];
	cout << Test.Body_Chosen[5].name << endl;

	Test.Body_Chosen[6] = Test.Body_List[4];
	cout << Test.Body_Chosen[6].name << endl;
	cin >> x;
//	Test.Body_Chosen[3] = Test.Body_List[5];
//	cout << Test.Body_Chosen[3].name << endl;

	Test.Min_Spice_Select = new double[Test.Body_Number];
	Test.Max_Spice_Select = new double[Test.Body_Number];
	Test.Min_time = new double[Test.Body_Number];
	Test.Max_time = new double[Test.Body_Number];
	Test.Min_Per = new double[Test.Body_Number];

	Test.Min_time[0] = 7.073136691849103e+08;
	Test.Min_time[1] = 0.086400000000000e+08;
	Test.Min_time[2] = 0.043200000000000e+08;
	Test.Min_time[3] = 0.086400000000000e+08;
	Test.Min_time[4] = 0.043200000000000e+08;
	Test.Min_time[5] = 0.345600000000000e+08;
	Test.Min_time[6] = 0.864000000000000e+08;


	Test.Max_time[0] = 7.073136691849103e+08;
	Test.Max_time[1] = 0.432000000000000e+08;
	Test.Max_time[2] = 0.259200000000000e+08;
	Test.Max_time[3] = 0.432000000000000e+08;
	Test.Max_time[4] = 0.259200000000000e+08;
	Test.Max_time[5] = 0.777600000000000e+08;
	Test.Max_time[6] = 1.555200000000000e+08;


	Test.Min_Spice_Select[0] = Test.Min_Spice_Time[2];
	Test.Min_Spice_Select[1] = Test.Min_Spice_Time[2];
	Test.Min_Spice_Select[2] = Test.Min_Spice_Time[1];
	Test.Min_Spice_Select[3] = Test.Min_Spice_Time[2];
	Test.Min_Spice_Select[4] = Test.Min_Spice_Time[3];
	Test.Min_Spice_Select[5] = Test.Min_Spice_Time[2];
	Test.Min_Spice_Select[6] = Test.Min_Spice_Time[4];


	Test.Max_Spice_Select[0] = Test.Max_Spice_Time[2];
	Test.Max_Spice_Select[1] = Test.Max_Spice_Time[2];
	Test.Max_Spice_Select[2] = Test.Max_Spice_Time[1];
	Test.Max_Spice_Select[3] = Test.Max_Spice_Time[2];
	Test.Max_Spice_Select[4] = Test.Max_Spice_Time[3];
	Test.Max_Spice_Select[5] = Test.Max_Spice_Time[2];
	Test.Max_Spice_Select[6] = Test.Max_Spice_Time[4];

	

	Test.NPerihelia = 0;
	Test.Nconstraints = 5;
	cout << "HEre";
	cin >> x;
	Test.Initialize_Mission(times, 0, 1);
	cout << "HEre2";
	cin >> x;
	for (i = 0; i < Test.Body_Number; i++) Test.Min_Per[i] = 200000;

	Test.Solution = Test.Current_Mission;
	Test.Optimize_Mission();
	Test.Current_Mission = Test.Solution;

//	Test.Current_Mission = NewMiss;
// New_Mission = new Mission(2, Test.Body_Select, times, Test.Min_Spice_Select, Test.Max_Spice_Select);
	
	cout << Test.Current_Mission.Trajectory.Trans_Set[0].transfer_body[0].ephem0.V << endl;
	cout << Test.Current_Mission.Trajectory.Trans_Set[0].transfer_body[1].ephem0.V << endl;

	cout << Test.Current_Mission.Trajectory.Body_Set[0].ephemt.V << endl;
	cout << Test.Current_Mission.Trajectory.Body_Set[1].ephemt.V << endl;
	cout << Test.Current_Mission.TotaldV << endl;
	cin >> x; 
	

	// OUMUAMUA
Test.Body_Number = 4;
double times[4] = { 6.731424691852796e+08 ,1.080000000000000e+08 ,1.308960000000000e+08 ,1.944000000000000e+08 };

times[1] = times[0] + times[1];
times[2] = times[1] + times[2];
times[3] = times[2] + times[3];
Test.Body_List[1].mu = 4.869e24*(6.67259e-11);		//Venus
Test.Body_List[2].mu = 3.98602e14;					//Earth
Test.Body_List[3].mu = 0.6419e24*(6.67259e-11);		//Mars
Test.Body_List[4].mu = 1898.6e24*(6.67259e-11);     //Jupiter
Test.Body_List[5].mu = 568.46e24*(6.67259e-11);		//Saturn
Test.Body_List[1].radius = 6052000;
Test.Body_List[2].radius = 6378135;
Test.Body_List[3].radius = 3397000;
Test.Body_List[4].radius = 71492000;
Test.Body_List[5].radius = 60268000;

Test.Body_Chosen = new Body[Test.Body_Number];
Test.Body_Chosen[0] = Test.Body_List[2];
cout << Test.Body_Chosen[0].name << endl;

Test.Body_Chosen[1] = Test.Body_List[4];
cout << Test.Body_Chosen[1].name << endl;

Test.Body_Chosen[2] = Test.Body_List[Test.NBody_List-1];
cout << Test.Body_Chosen[2].name << endl;

Test.Body_Chosen[3] = Test.Body_List[Test.NBody_List-2];
cout << Test.Body_Chosen[3].name << endl;



Test.Min_Spice_Select = new double[Test.Body_Number];
Test.Max_Spice_Select = new double[Test.Body_Number];
Test.Min_time = new double[Test.Body_Number];
Test.Max_time = new double[Test.Body_Number];
Test.Min_Per = new double[Test.Body_Number];

Test.Min_time[0] = 6.705072691856551e+08;
Test.Min_time[1] = 0.086400000000000e+08;
Test.Min_time[2] = 0.025920000000000e+08;
Test.Min_time[3] = 0.432000000000000e+08;



Test.Max_time[0] = 6.757776691849042e+08;
Test.Max_time[1] = 2.073600000000000e+08;
Test.Max_time[2] = 2.592000000000000e+08;
Test.Max_time[3] = 3.456000000000000e+08;
/* for (i = 1; i < Body_Number; i++)
{
	Test.Min_time[i] += Test.Min_time[i - 1];
	Test.Max_time[i] += Test.Max_time[i - 1];
}



Test.Min_Spice_Select[0] = Test.Min_Spice_Time[2];
Test.Min_Spice_Select[1] = Test.Min_Spice_Time[4];
Test.Min_Spice_Select[2] = Test.Min_Spice_Time[Test.NBody_List-1];
Test.Min_Spice_Select[3] = Test.Min_Spice_Time[Test.NBody_List-2];



Test.Max_Spice_Select[0] = Test.Max_Spice_Time[2];
Test.Max_Spice_Select[1] = Test.Max_Spice_Time[4];
Test.Max_Spice_Select[2] = Test.Max_Spice_Time[Test.NBody_List-1];
Test.Max_Spice_Select[3] = Test.Max_Spice_Time[Test.NBody_List-2];




Test.NPerihelia = 0;
Test.Nconstraints = 2;
cout << "HEre";
cin >> x;
Test.Initialize_Mission(times, 0, 0);
cout << "HEre2";
cin >> x;
for (i = 0; i < Test.Body_Number; i++) Test.Min_Per[i] = 000000;
Test.Min_Per[2]=2.0944e+09;
Test.NPerihelia = 2;
Test.Per_Pointer = new int[Test.NPerihelia];
Test.Perihelia = new double[Test.NPerihelia];
Test.Per_Pointer[0] = 1;
Test.Per_Pointer[1] = 2;
Test.Perihelia[0] = 2.0000e+09;
Test.Perihelia[1] = 2.0000e+09;

Test.Solution = Test.Current_Mission;
Test.Optimize_Mission();
Test.Current_Mission = Test.Solution;

//	Test.Current_Mission = NewMiss;
// New_Mission = new Mission(2, Test.Body_Select, times, Test.Min_Spice_Select, Test.Max_Spice_Select);

cout << Test.Current_Mission.Trajectory.Trans_Set[0].transfer_body[0].ephem0.V << endl;
cout << Test.Current_Mission.Trajectory.Trans_Set[0].transfer_body[1].ephem0.V << endl;

cout << Test.Current_Mission.Trajectory.Body_Set[0].ephemt.V << endl;
cout << Test.Current_Mission.Trajectory.Body_Set[1].ephemt.V << endl;
cout << Test.Current_Mission.TotaldV << endl;
cin >> x;

return 0;
*/
}

