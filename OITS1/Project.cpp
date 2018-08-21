#include "stdafx.h"
#include "Project.h"


Project::Project()
{
}

int Project::Read_Data(string infile)
{
	string buffer;

	int i;
	int nlines = 0;
	int nfam;
	int nplan;
	char y;
	string family[5] = { "SPICE","TRAJ","BODY","TRANS", "END" };

	char x;
	int alldataflag = 0;
	// Open Input File
	ifstream input_file(infile);
	if (!(input_file.is_open()))
	{
		cout << " Cannot open input file :" << infile << endl;
		cin >> x;
		return 1;
	}

	while (alldataflag<15)
	{
		nfam = -1;
		if(Get_Line(input_file, buffer)==0)break;

	// Check for start of a family of data
		for (i = 0; i < 5; i++)
		{
			if (buffer.find(family[i]) != string::npos)
			{
//				cout << "Family found = " << family[i] << endl;
				nfam = i;			// Which Family has been found!
				break;
			}
		}
		switch (nfam)
		{
		case 0:						// SPICE Family
//			cout << "SPICE Family Found" << endl;
			if (Get_Line(input_file, buffer) == 0)break;
			naif_file = buffer;			// Get Leap Second File
			while (Get_Line(input_file, buffer)>0)
			{
				if ((buffer.find("END") ==0) || (Number_BSP >= MAX_FILES))break;
				BSP[Number_BSP] = buffer; // Get Binary SPK file
				Number_BSP++;
//				cout << " BSP File " << Number_BSP << " = " << BSP[Number_BSP-1] << endl;
			}
			if (Number_BSP >= MAX_FILES)cout << " Too Many Binary SPK Files, Limit = " << MAX_FILES << endl;
			MetaKernel = naif_file;
//			for (i = 0; i < Number_BSP; i++)MetaKernel = MetaKernel + "\n" + BSP[i];
//			cout << MetaKernel;
//			cin >> x;
			alldataflag += 1;
			break;
		case 1:						// TRAJ Family
//			cout << "TRAJ Family Found" << endl;
			if (Get_Line(input_file, buffer) == 0)break;
			name = buffer;			// Get name of Current Mission
			if (Get_Line(input_file, buffer) == 0)break;
			

					/*  INITIALISE TRAJECTORY DATA */

			Body_Number = atoi(buffer.c_str());	// Number of Bodies
			Body_Chosen = new Body[Body_Number];//  Allocate Space for Bodies
			Min_Per = new double[Body_Number];	// Allocate Space for Minimum Periapsis Altitude


			if (Get_Line(input_file, buffer) == 0)break;
			if ((buffer.find("Y") != string::npos) || (buffer.find("y") != string::npos))FlybyRendez = 1; // Set Mission to Rendezvous
			if (Get_Line(input_file, buffer) == 0)break;
			if ((buffer.find("N") != string::npos) || (buffer.find("n") != string::npos)) wayflag = 0; // Set Mission to Prograde and Retrograde
			if (Get_Line(input_file, buffer) == 0)break;
			if ((buffer.find("M") == string::npos) && (buffer.find("m") == string::npos)&&(buffer.find("END")==string::npos)&&(buffer.find("end")==string::npos))
				Max_Duration = 365.25*24*60*60*atof(buffer.c_str()); // Set Maximum Duration
//			cout << "Body_Number = " << Body_Number << " BR= " << FlybyRendez << " WF= " << wayflag << " MD = " << Max_Duration << endl;
			alldataflag += 2;
			break;
			

		case 2:						// BODY Family
//			cout << "BODY Family Found" << endl;
			for (i = 0; i < Body_Number; i++) 
			{
				if (Get_Line(input_file, buffer) == 0)break;
				nplan = atoi(buffer.c_str());	//Current Number of Body
				if (Get_Line(input_file, buffer) == 0)break;
				Body_Chosen[nplan-1].ID = new char[buffer.length()];
				Body_Chosen[nplan - 1].ID[buffer.length() ] = '\0';
				memcpy(Body_Chosen[nplan-1].ID, buffer.c_str(), buffer.length());// Current ID of Body
				if (strcmp(Body_Chosen[nplan - 1].ID, "IP") == 0) Number_IP++;
				if (Get_Line(input_file, buffer) == 0)break;
				Body_Chosen[nplan - 1].mu = 0.0;
				Min_Per[nplan - 1] = atof(buffer.c_str())*1000.0; // Minimum Perapsis Altitude
//				cout << "nplan = " << nplan << " ID = " << Body_Chosen[nplan - 1].ID << " Periapsis = " << Min_Per[nplan - 1] << endl;
			}
			alldataflag += 4;
			break;

		case 3:						// TRANS Family
//			cout << "TRANS Family Found" << endl;
						
			while ((Get_Line(input_file, buffer) > 0)&&(buffer.find("END") == string::npos) && (buffer.find("end") == string::npos)&& (NPerihelia<MAX_NBODY-1))
			{
				Per_Pointer[NPerihelia] = atoi(buffer.c_str());
				if (Get_Line(input_file, buffer) == 0)break;
				if (buffer.find("p") != string::npos)
					buffer = buffer.substr(0, buffer.find("p"));
				else if (buffer.find("P") != string::npos)
						buffer=buffer.substr(0, buffer.find("P"));
				
				Perihelia[NPerihelia] = AU*atof(buffer.c_str());
//				cout << "NPerihelia = " << (NPerihelia + 1) << " Transfer = " << Per_Pointer[NPerihelia] << " Perihelia = " << Perihelia[NPerihelia] << endl;
				NPerihelia++;
				
			}
			alldataflag += 8;
			break;

		case 4:						// END Found
			break;
			
		default:					// Unknown Family name
			cout << "Unknown Family Name : " << buffer << endl;
		}
		
	//	Remove_Comment(buffer, buffer.length()); 
	}
	if (alldataflag < 7) {
		cout << "All Data Not Present" << endl;
		cin >> x;
		return 1;
	}
	input_file.close();
	return 0;
	// First Read Name of Project

}
int Project::Get_Line(ifstream& input_file, string& first_bit)
{

	string buffer;
	string no_comment_bit;
	size_t pos1, pos2, pos3=0,pos4=0;
	int nlines = 0;
	
	while (pos3 == 0||pos4==0)
	{
		getline(input_file, buffer);
		if (input_file.eof())break;

		pos1 = buffer.find("%");			// Find Position of Comment in Line
		pos2 = buffer.length();				// Find Length of buffer string
		pos3 = pos1;
		if (pos2 < pos3)pos3 = pos2;		// Find Amount of useful data on line
		if (pos3 == 0) continue;			// Line Contains no data: Discard
		no_comment_bit = buffer.substr(0, pos3); // Remove Comments
		pos4 = pos3;
		while ((pos4 > 0)&&((no_comment_bit[pos4-1] == ' ') || (no_comment_bit[pos4-1] == '\t')))pos4--; // Remove Trailing blanks
		first_bit = no_comment_bit.substr(0, pos4);
		if (pos4 > 0)nlines++;				// Line Contains Data
		else continue;

	//	cout << first_bit << "E" << endl << " Nlines = " << nlines << " pos1 = " << pos1 << " pos2 = " << pos2 << " pos3 = " << pos3 << endl;
	}
	
	
	return nlines;
	
}

int Project::Initialize_SPICE()
{
	char x[1000];
	int k;
	char *Kern;

	tkvrsn_c("toolkit");
	

	/* Initialise Leap Second File*/
	furnsh_c((char *)naif_file.c_str());

// Initialise SPICE Kernels
	for (k = 0; k < Number_BSP; k++) 
	{
		furnsh_c((char * )BSP[k].c_str());
	}

	return 0;

}
int Project::Get_SPICE_List()
{
	int MAXIV = 1000;
	int WINSIZ = 2 * MAXIV;
	char x;

	SPICEDOUBLE_CELL(cover, 1000);
	SPICEINT_CELL(ids, 100);

	SpiceChar               timstr1[51];
	SpiceChar               timstr2[51];

	SpiceDouble             b;
	SpiceDouble             e;

	SpiceInt                i;
	SpiceInt                j;
	SpiceInt                k;
	SpiceInt                niv;
	SpiceInt                obj;
	SpiceInt				idcode;
	SpiceInt				count;

	SpiceBoolean  found;

	Body * BODYLIST, *BODYLISTN;
	char * Kern;
	int NBody;
	char buffer[50] = { " " };
	double *SPICE_TIME1, *SPICE_TIME2;

	// Initialise Minimum and Maximum SPICE Times;

	Min_Spice_Select = new double[Body_Number];
	Max_Spice_Select = new double[Body_Number];

	for (i = 0; i < Body_Number; i++)
	{
		if (strcmp(Body_Chosen[i].ID, "IP")==0)continue;
		prsint_c(Body_Chosen[i].ID, &idcode);
		
		for (j = 0; j < Number_BSP; j++)
		{
			Kern = (char *)(BSP[j].c_str());
			spkcov_c(Kern, idcode, &cover);
		}
		wnfetd_c(&cover, 0, &b, &e);
		Min_Spice_Select[i] = b;
		Max_Spice_Select[i] = e;

			/* .
				Convert the endpoints to TDB calendar
				format time strings and display them.
				. */
/*		timout_c(b,
					"YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
					51,
					timstr1);

		printf("\n"
				"Interval:  %d\n"
				"Start:     %s\n",
				1,
				timstr1);

		timout_c(e,
				"YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
				51,
				timstr2);
			printf("Stop:      %s\n", timstr2);
*/
	}
return 0;
}


/*	for (k = 0; k < Number_BSP; k++)
	{
		spkobj_c((char *)BSP[k].c_str(), &ids);



		/*
			Find the set of objects in the SPK file.
			




			// Number of separate bodies in Kernel

		NBody = card_c(&ids) - NBody_List;
		for (i = 0; i < NBody; i++)
		{
			cout << SPICE_CELL_ELEM_I(&ids, i);
		}

		BODYLISTN = new Body[NBody];
		SPICE_TIME1 = new double[NBody + NBody_List];
		SPICE_TIME2 = new double[NBody + NBody_List];

		if (NBody_List > 0)
		{
			memcpy(SPICE_TIME1, Min_Spice_Time, NBody_List * sizeof(double));
			memcpy(SPICE_TIME2, Max_Spice_Time, NBody_List * sizeof(double));
		}
		Min_Spice_Time = SPICE_TIME1;
		Max_Spice_Time = SPICE_TIME2;



		/*
			We want to display the coverage for each object.Loop over
			the contents of the ID code set, find the coverage for
			each item in the set, and display the coverage.
			
		for (i = 0; i < NBody; i++)
		{
			/*
				Find the coverage window for the current object.
				Empty the coverage window each time so we don't
				include data for the previous object.
				
			obj = SPICE_CELL_ELEM_I(&ids, i + NBody_List);

			scard_c(0, &cover);
			spkcov_c((char *)BSP[k].c_str(), obj, &cover);

			 Display a simple banner.
				

			printf("%s\n", "========================================");

			printf("Coverage for object %d\n", (int)obj);

			
				Convert the coverage interval start and stop times to TDB
				calendar strings.
				*/
				/*
					Get the endpoints of the jth interval.
					
			wnfetd_c(&cover, 0, &b, &e);

			/*
				Convert the endpoints to TDB calendar
				format time strings and display them.
				
			timout_c(b,
				"YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
				51,
				timstr1);

			printf("\n"
				"Interval:  %d\n"
				"Start:     %s\n",
				1,
				timstr1);

			timout_c(e,
				"YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
				51,
				timstr2);
							printf("Stop:      %s\n", timstr2);
			sprintf_s(buffer, "%d", SPICE_CELL_ELEM_I(&ids, i + NBody_List));
			//			printf("\n%s\n", buffer);
			BODYLISTN[i].ID = new char[strlen(buffer) + 1];
			memcpy(BODYLISTN[i].ID, buffer, strlen(buffer) + 1);
			bodc2n_c(atoi(BODYLISTN[i].ID), 50, BODYLISTN[i].name, &found);
			BODYLISTN[i].radius = 0;
			str2et_c(timstr1, &Min_Spice_Time[i + NBody_List]);
			str2et_c(timstr2, &Max_Spice_Time[i + NBody_List]);
		}

		if (NBody_List == 0)
		{
			Body_List = BODYLISTN;
			NBody_List = NBody;
		}
		else
		{
			BODYLIST = new Body[NBody_List + NBody];
			memcpy(BODYLIST, Body_List, NBody_List * sizeof(Body));
			memcpy(BODYLIST + NBody_List, BODYLISTN, NBody * sizeof(Body));
			NBody_List = NBody_List + NBody;
			Body_List = BODYLIST;
		}
	}*/


/*
	string doscommand;
	char   psBuffer[128];
	string temp;
	int i,j,l;
	int n;
	size_t timeposfirst, timeposlast;
	typedef struct spice_line
	{
		struct spice_line *next;
		string line;
	} SPICEL;
	SPICEL *listend, *head, *current, *oldcurrent, *startdata, *enddata;
	int NBody;
	FILE *pPipe;
	Body * BODYLIST, *BODYLISTN;
	double *Old_Min_Times, *Old_Max_Times;

	char *Kern = (char*)Kernel.c_str();

	// Set up DOS Command to Get Run Down of SPICE Data File

	doscommand = "C:/Users/adamh/Documents/SPICE/mice/mice/exe/brief -t " +  Kernel;
	char *doscommandc = (char*)doscommand.c_str();

	// Initialise Kernel
	furnsh_c(Kern);

	// Get Output of SPICE brief and pit in output
	if ((pPipe = _popen(doscommandc, "rt")) == NULL)
		return -1;

	/*	 Read pipe until end of file, or an error occurs. */
/*	head = new SPICEL;
	current = head;
	while (fgets(psBuffer, 128, pPipe))
	{
		current->line = psBuffer;
		current->next = new SPICEL;
		oldcurrent = current;
		current = current->next;
	}
	listend = current;
	current = head;
	while (current != listend)
	{
		// cout << current->line << endl;
		current = current->next;

	}
//	/* Close pipe and print return value of pPipe. */
/*	if (feof(pPipe))
	{
		printf("\nProcess returned %d\n", _pclose(pPipe));
	}
	else
	{
		printf("Error: Failed to read the pipe to the end.\n");
	}
	
	// Find the line where all the Body ID's are listed

	current = head;

	while (current != listend)
	{
		if ((current->line.find("------- ")) != std::string::npos)break;
		current = current->next;
	}
		// Cater for absence of Data in File	
	if (current == listend)
	{

		cout << "No Bodies Found" << endl;
		return 0;
	}

	NBody = 0;

	startdata = current->next;
	// cout << current->line << endl;
	timeposfirst = current->line.find("-----------------------------",0);
	timeposlast =  current->line.find("-----------------------------",timeposfirst+1);
	
	while (current != listend)
	{
		current = current->next;
		enddata = current;
		// Extract List of Bodies from Relevant Line
		if (current->line.length() == 0)break;
		NBody++;
	}

	if (NBody > 0)
		cout << "Data Found: Nbodies = " << NBody << endl;
	else
	{
		cout << "No Bodies Found" << endl;
		return 0;
	}
	cout << sizeof(Body) << endl;
	// Initialize Body_List
	current = startdata;
	j = 0;
		
	if (NBody_List > 0)
	{
		BODYLIST = new Body[NBody_List];
		memcpy(BODYLIST, Body_List, NBody_List * sizeof(Body));
		BODYLISTN = new Body[NBody];
//		Old_Min_Times = Min_Spice_Time;
//		Old_Max_Times = Max_Spice_Time;
		Old_Min_Times = new double[NBody_List];
		Old_Max_Times = new double[NBody_List];
		memcpy(Old_Min_Times, Min_Spice_Time, sizeof(double)*NBody_List);
		memcpy(Old_Max_Times, Max_Spice_Time, sizeof(double)*NBody_List);
		Min_Spice_Time = new double[NBody + NBody_List];
		Max_Spice_Time = new double[NBody + NBody_List];
		for (l = 0; l < NBody_List; l++)
		{
			/*	Body_List[l] = BODYLIST[l];
			Min_Spice_Time[l] = Old_Min_Times[l];
			Max_Spice_Time[l] = Old_Max_Times[l]; */
/*			cout << Body_List[l].ID << endl;
		}
		while (current != enddata)
		{
			if (current == startdata)
			{
				string time1str = current->line.substr(timeposfirst, 20);
				string time2str = current->line.substr(timeposlast, 20);
				// cout << time1str << endl;
				// cout << time2str << endl;
				char *time1ch = (char *)time1str.c_str();
				char *time2ch = (char *)time2str.c_str();
				str2et_c(time1ch, &Min_Spice_Time[j+NBody_List]);
				str2et_c(time2ch, &Max_Spice_Time[j+NBody_List]);
			}
			else
			{
				Min_Spice_Time[j+NBody_List] = Min_Spice_Time[j + NBody_List-1];
				Max_Spice_Time[j+NBody_List] = Max_Spice_Time[j + NBody_List-1];
			}	
			n = current->line.find(" ", 0);

			temp = current->line.substr(0, n);
			BODYLISTN[j].ID = (char *)temp.c_str();
			BODYLISTN[j].ID = new char[temp.length() + 1];
			BODYLISTN[j].ID[temp.length()] = '\0';
			memcpy(BODYLISTN[j].ID, (char *)temp.c_str(), temp.length() * sizeof(char));

			bodc2n_c(atoi(BODYLISTN[j].ID), 50, BODYLISTN[j].name, &n);
			BODYLISTN[j].radius = 0;

			j++;
			current = current->next;
		}
		Body_List = new Body[NBody_List + NBody];
		memcpy(Body_List, BODYLIST, NBody_List * sizeof(Body));
		for (l = 0; l < NBody_List; l++)
		{
			Body_List[l] = BODYLIST[l];
			Min_Spice_Time[l] = Old_Min_Times[l];
			Max_Spice_Time[l] = Old_Max_Times[l]; 
			cout << Body_List[l].ID << endl;
		}
		for (l = NBody_List; l < NBody_List + NBody; l++)
		{
			Body_List[l] = BODYLISTN[l - NBody_List];
			cout << Body_List[l].ID << endl;
		}
		NBody_List = NBody_List + NBody;
	}
	else
	{
		NBody_List = NBody;
		Body_List = new Body[NBody_List];
		Min_Spice_Time = new double[NBody];
		Max_Spice_Time = new double[NBody];
		while (current != enddata)
		{
			if (current == startdata)
			{
				string time1str = current->line.substr(timeposfirst, 20);
				string time2str = current->line.substr(timeposlast, 20);
				// cout << time1str << endl;
				// cout << time2str << endl;
				char *time1ch = (char *)time1str.c_str();
				char *time2ch = (char *)time2str.c_str();
				str2et_c(time1ch, &Min_Spice_Time[j]);
				str2et_c(time2ch, &Max_Spice_Time[j]);
			}
			else
			{
				Min_Spice_Time[j] = Min_Spice_Time[j - 1];
				Max_Spice_Time[j] = Max_Spice_Time[j - 1];
			}
			n = current->line.find(" ", 0);

			temp = current->line.substr(0, n);
			Body_List[j].ID = new char[temp.length() + 1];
			Body_List[j].ID[temp.length()] = '\0';
			memcpy(Body_List[j].ID, (char *)temp.c_str(), temp.length()* sizeof(char));
			bodc2n_c(atoi(Body_List[j].ID), 50, Body_List[j].name, &n);
			Body_List[j].radius = 0;

			j++;
			current = current->next;
		}

				for (l = 0; l < NBody_List; l++)
		{
			/*	Body_List[l] = BODYLIST[l];
			Min_Spice_Time[l] = Old_Min_Times[l];
			Max_Spice_Time[l] = Old_Max_Times[l]; */
/*			cout << Body_List[l].ID << endl;
		}


	}
	return 0;
}
*/

int Project::Add_Intermediate_Point(int num)

// Add_Intermediate_Point   :   Adds an INTERMEDIATE POINT to be Optimized to the Possible Bodies to Select From
{
	
	Body_Chosen[num].ID = "INTERMEDIATE POINT";
	memcpy(Body_Chosen[num].name, Body_Chosen[num].ID, strlen(Body_Chosen[num].ID) + 1);
	Body_Chosen[num].radius = 0;
	Body_Chosen[num].Fixed_Point = 1;
	Body_Chosen[num].ephem0.r = { cos(num/Body_Number*2*PI_OITS)*Min_Per[num]/1000*AU, sin(num / Body_Number * 2 * PI_OITS)*Min_Per[num] / 1000 * AU, 0 };
	Body_Chosen[num].ephem0.v = { 0, 0, 0 };
	Body_Chosen[num].ephem0.t = 0;
	Body_Chosen[num].mu = 0.0;
	double spice_min = -1e50;
	double spice_max = 1e50;


	Min_Spice_Select[num] = spice_min;
	Max_Spice_Select[num] = spice_max;
	return 0 ;
}
int Project::Merge_Data()
{

	//Merges Data on the Planets with Data extracted from Get_SPICE_List

	char *PLANET_NAMES[] = { "MERCURY", "VENUS", "EARTH", "MARS", "JUPITER", "SATURN", "URANUS", "NEPTUNE", "PLUTO" };
	double PLANET_MU[9] = { 6.67259e-11*0.3302e24, 6.67259e-11*4.869e24,
							6.67259e-11*5.9736e24, 6.67259e-11*0.6419e24,
							6.67259e-11*1898.6e24, 6.67259e-11*568.46e24,
							6.67259e-11*86.83e24, 6.67259e-11*102.43e24, 6.67259e-11*0.0125e24 };
	double PLANET_RADIUS[9] = { 2440e3, 6052e3, 6378135, 3397e3, 71492e3, 60268e3, 25559e3, 24766e3, 1137e3, };
	SpiceBoolean   found;
	int i, j;

	for (i = 0; i < Body_Number; i++)
	{

		bodc2n_c(atoi(Body_Chosen[i].ID), 50, Body_Chosen[i].name, &found);

		for (j = 0; j < 9; j++)
			if (strstr(Body_Chosen[i].name, PLANET_NAMES[j]) !=NULL )break;
		if (j < 9)
		{
			Body_Chosen[i].radius = PLANET_RADIUS[j];
			Body_Chosen[i].mu = PLANET_MU[j];
//			cout << Body_Chosen[i].name << endl;
		}
	}
	return 0;
}




int Project::Initialize_Mission(double *times, int flag , int flag2)
{
	Current_Mission = Mission::Mission(Body_Number, Body_Chosen, times, Min_Spice_Select, Max_Spice_Select, flag, flag2);
//	Current_Mission.Per_Pointer = Per_Pointer;
//	Current_Mission.NPerihelia = NPerihelia;
//	Current_Mission.Nconstraints = Nconstraints;
//	Current_Mission.Max_duration = Max_Duration;
//	Current_Mission.Perihelia = Perihelia;
//	Current_Mission.Min_time = Min_time;
//	Current_Mission.Max_time = Max_time;
	Current_Mission.wayflag = wayflag;
	return 0;

}


Project::~Project()
{
}
