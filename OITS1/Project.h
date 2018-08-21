#pragma once
#define MAX_FILES 100
#define MAX_NBODY 20
#include <string>
#include <fstream>
#include <sstream>
#include <C:\Users\adamh\Documents\SPICE\cspice\include\SpiceUsr.h>
#include "Mission.h"

class Project
{
public:
		string name;       // Name of Project
		string naif_file;	// Filename for Leap Second File
		string BSP[MAX_FILES];        // Filename for Binary SPK file from JPL
		string MetaKernel;			// Leap Second File + List of all SPK files 
		int Number_BSP = 0;			// Total Number of Binary SPK Files Specified
		int Number_IP = 0;			// Total Number of intermediate points
		Body *Body_List;   // List of Bodies Available to Choose from
		int NBody_List=0;  // Number of Bodies Available to Choose from
		int Max_NBody = MAX_NBODY; // Maximum Number of Bodies allowed for Optimizer
		double *Min_time;	// Array of Minimum Times for Optimization
		double *Max_time;	// Array of Maximum Times for Optimization
		double *Min_Spice_Time; // Lower Limit of Spice Kernel Range
		double *Max_Spice_Time; // Upper Limit of Spice Kernel Range
		double *Min_Spice_Select; // Selected Values of Min_Spice_Time
		double *Max_Spice_Select;  // Selected Values of Max_Spice_Time
		Body *Body_Select;    // List of Selected Bodies
		Body *Body_Chosen;    // List of Chosen Bodies for Optimization
		int Body_Number;    // Number of Selected Body
		int FlybyRendez = 0;      // Flyby or Rendezvous Flag = 0 (FLYBY) = 1 (RENDEZVOUS)
		int wayflag = 1;		// Default to Prograde only
		Mission Current_Mission;
		Mission Global_Solution;    // Result of Running Global Optimizer
		Mission Local_Solution;    // Result of Running Local Optimzer
		Mission Solution;           // Result of Last Run of Optimizer
		int Run_Time = 20;       // Maximum Optimizer Run Time in minutes
		double Perihelia[MAX_NBODY];		// Array of Perihelia for Trajectory
		int Nconstraints=0;       // Number of Periapsis Constraints
		int NPerihelia=0;         // Number of Perihelia Constraints
		int Per_Pointer[MAX_NBODY-1];        // Pointer to Array of Perihelia Values
		double *Constr_Tol;         // Array of Constraint Tolerances
		double *Min_Per;            // Array of Minimum Periapsis for each Body
		double MAX_DURATION = 1e50;  // Upper Limit for Max_Duration
		double Max_Duration = 1e50;        // Maximum Duration Of Entire Mission(optional)
		double factor = 0.5;             // Factor used for Local Optimizer
		double AU = 149597870700;			// Astronical Unit (m)
	
	Project();
	int Read_Data(string);
	int Get_Line(std::ifstream& , string&);
	int Initialize_SPICE();
	int Get_SPICE_List();
	int Add_Intermediate_Point(int);
	int Merge_Data();
	int Initialize_Mission(double *, int, int);
	~Project();
};

