#include "stdafx.h"
#include "Connecting_Hyperbola.h"

Connecting_Hyperbola::Connecting_Hyperbola()
{
}

////# Transform calculated Transformation Matrices from Ecliptic to B - Plane base on Arrival Velocity
int Connecting_Hyperbola::Transform_Matrix()
// # Transform calculated Transformation Matrices
//#
//# INPUT:
//#
//# Object in question with uninitialsied Transformation Matrices
//#
//# OUTPUT :
//#
//#  Object in question with initialsied Transformation Matrices
//#

{

	double Vij;	// Projection of VA onto ecliptic
	double Vmag; // Magnitude of VA(speed)
	double cosdel, sindel; // Cosine and Sine of angle 'delta'
	double cosgam, singam; // Cosine and Sine of angle 'gamma'

	Vij = sqrt(pow(VA(1), 2) + pow(VA(0), 2));

	// Magnitude of VA(speed)
	Vmag = VA.norm();

	// Compute angles
	cosdel = VA(1) / Vij;
	sindel = VA(0) / Vij;
	cosgam = VA(2) / Vmag;
	singam = Vij / Vmag;

	// Transformation Matrices
	Trans1 << cosdel, -sindel, 0, 
			  sindel, cosdel, 0, 
					0,     0, 1;

	Trans2 << 1,	  0,	  0,
			  0, cosgam, -singam, 
			  0,  singam, cosgam;
	return 0;

} //# Transform

////# Compute_Hyperbola connecting Arrival Velocity VA with Departure Velocity VD
////# and find Perapsis Distance, DeltaV required at Periapsis and Angle beta
int Connecting_Hyperbola::Compute_Hyperbola()
// # Compute_Hyperbola connecting Arrival Velocity VA with Departure Velocity VD
//#
//# INPUT:
//#
//#  Object in question with given VA, VD
//#
//# OUTPUT :
//#
//# Object in question with calculated alpha, beta, Per and DV
//#

//       options = optimset('Display', 'iter');
{
	double cosalpha;	// Cosine aplha angle
	double sinbeta;		// Cosine beta angle
	double VAV;
	Vector3d Vdiff;
	double X0;			// Initial Guess of True anomaly at Arrival
	double X;			// Solution True Anomaly at Arrival 
	double TolX = 1e-50;

	// Determine angle alpha between Departure Velocity and Arrival
	// Velocity

	cosalpha = VA.dot(VD) / VA.norm() / VD.norm();
	alpha = acos(cosalpha);
	VAV = sqrt(VA.norm()*VD.norm());

	if (isnan(alpha))
	{
		DV = 1e7;
		Per = 0;
	}
	if (Planet.mu == 0)
	{
		Per = 0;
		Vdiff = VD - VA;
		if (!isnan(alpha)) DV = Vdiff.norm();
		NO_ENCOUNTER = 1;
	}
	else if ((!isnan(alpha))&&(abs(alpha) > PI_OITS - alpha_thresh))
	{
		Per = Planet.mu*pow((PI_OITS - abs(alpha)), 2) / 8 / pow(VAV, 2);
		DV = abs(1.0 - Per / Planet.radius)*1e7;
	
		
		LOW_PER_ERROR = 1;
		NO_ENCOUNTER = 0;
		
//		cout << DV << endl;
	}
	else if (!isnan(alpha))
	{
		LOW_PER_ERROR = 0;
		NO_ENCOUNTER = 0;


		// Determine angle beta

		sinbeta = (VA(1)*VD(0) - VA(0)*VD(1))/VD.norm() / sqrt(pow(VA(0), 2) + pow(VA(1), 2)) / sin(alpha);
		beta = asin(sinbeta);

		// Guess Angle between periapsis and Arrival Asymptote

		X0 = (PI_OITS + alpha) / 2;

		X = fzero(X0, TolX, 0, 1000);
//		if (isnan(X))
//		{
//			LOW_PER_ERROR = 1;
//			cout << " ISNAN \n" << "X0= " << X0 << "VA=" << VA << " VD= " << VD << "\n";
//		}
//		double Per1, Per2,DV1,DV2;

//		Per1 = -Planet.mu / pow(VA.norm(), 2) * (1 + 1 / cos(X));
//		DV1 = sqrt(2 * Planet.mu / Per1 + pow(VD.norm(), 2)) - sqrt(2 * Planet.mu / Per1 + pow(VA.norm(), 2));
//		cout << "X1 = " << setprecision(20) << X << endl;
//		cout << " Perigee 1 : " << Per1 << endl;
//		cout << " DV 1 : " << DV1 << endl;


		Per = -Planet.mu / pow(VA.norm(), 2) * (1 + 1 / cos(X));
		DV = sqrt(2 * Planet.mu / Per + pow(VD.norm(), 2)) - sqrt(2 * Planet.mu / Per + pow(VA.norm(), 2));

//		X = fzero(X0, TolX, 1, 1000);


//		Per2 = -Planet.mu / pow(VD.norm(), 2) * (1 + 1 / cos(X));
//		DV2 = sqrt(2 * Planet.mu / Per2 + pow(VD.norm(), 2)) - sqrt(2 * Planet.mu / Per2 + pow(VA.norm(), 2));
//		cout << "X2 = " << X << endl;
//		cout << " Perigee 2 : " << Per2 << endl;
//		cout << " DV 2 : " << DV2 << endl;
		

	}
	return 0;
}


// Employs Newton Iteration to converge to solution for angle theta (output x)
double Connecting_Hyperbola::fzero(double x0, double xtol, int mode, int Maxit)
{
	double func,dfuncbydx;
	double VD2, VA2;
	double cal, sal;
	double x,xold,dx=0;
	double reduce = 1.0;
	double dxold = 0;
	double funcold = 0;
	int i;

	VD2 = pow(VD.norm(), 2);
	VA2 = pow(VA.norm(), 2);
	cal = cos(alpha);
	sal = sin(alpha);
	
	x = x0;
	for (i = 0; i < Maxit; i++)
	{
		if (mode == 0)
		{
			func = (VA2 + VD2*cal)*cos(x) + (VD2 - VA2)*sal*sin(x)*cos(x) + (VD2 - VA2)*pow(cos(x), 2) * cal + VD2*sal*sin(x);
			func = func / Planet.mu;
		//	func = 1 / Planet.mu * (((VD2 - VA2)*cos(x) + VD2)*cos(alpha - x) + VA2*cos(x));
			dfuncbydx = 1/Planet.mu*((VD2 - VA2)*sin(alpha - 2 * x) + VD2 * sin(alpha - x) - VA2 *sin(x));
		}
		else
		{	
			func = ((VA2 - VD2)*cos(alpha - x) + VD2) *cos(x) + VA2*cos(alpha-x);
			dfuncbydx = (VA2 - VD2)*sin(alpha - 2 * x) + VA2 * sin(alpha - x) - VD2 *sin(x);
		}

		if ((funcold * func) < 0.0)
		{
			reduce = reduce*0.1;
		}
		else
		{
			reduce = 1.0;
		}
		
		dxold = dx;
		funcold = func;
		xold = x;

		dx = reduce*func / dfuncbydx;
		x = x - dx;		
/*		cout << "x= " << setprecision(50) << x << endl;
		cout << "dx= " << dx << endl;
		cout << "DELTA= " << x - xold << endl;
		cout << "reduce = " << reduce << endl;
		cout << "fx = " << func << endl << endl; */
		if (abs(x-xold) < xtol) break;
	}
	NIT = i;
//	cout << NIT << endl;
	if (x >  2*PI_OITS) x = x -  2*PI_OITS*int(x / 2 / PI_OITS);
	else if (x < -2 * PI_OITS)
	{
		x = abs(x);
		x = PI_OITS - x + PI_OITS*int(x / PI_OITS);
	}
	return x;
}

// Calculates the orbital parameters from the Hyperbolas, assuming a solution has been found
int Connecting_Hyperbola::Orbits_From_Hyperbolas()
{
	double VAS, VDS;	// Arrival and Departure Speed
	double ENA, END;	// Arrival and Departure Energy
	double HA, HD;		// Arrival and Departure Angular Momentum Magnitudes
	Vector3d uH;		// Unit Angular Momentum Vector (Arrival and Departure are co-linear)
	Vector3d HAV, HDV;	// Angular Momentum Vector (Arrival and Departure) 
	double L;			// Longitude of Ascending Node (Same both ways)
	double I;			// Inclination (Same both ways)
	Vector3d EvA, EvD;	// Laplace - Range - Lenz vector
	double EA, ED;		// Magnitudes of EvA,EvD
	double sinwA, coswA;// sine and cos of Argument of Perigee
	double sinwD, coswD;//                  ""

	// Initialise Gravitational Masses
	Probe.orbit.GM = Planet.mu;
	Probe2.orbit.GM = Planet.mu;

	// Arrival Speed
	VAS = VA.norm();
	// Departure Speed
	VDS = VD.norm();

	// Arrival Energy
	ENA = VAS*VAS / 2;
	// Departure Energy
	END = VDS*VDS / 2;

	// Arrival Semi - major Axis
	Probe.orbit.a = -Probe.orbit.GM / 2 / ENA;
	Probe.orbit.arec = 1 / Probe.orbit.a;
	// Departure Semi - major Axis
	Probe2.orbit.a = -Probe2.orbit.GM / 2 / END;
	Probe2.orbit.arec = 1 / Probe2.orbit.a;

	// Arrival Eccentricity
	Probe.orbit.e = Per*VAS*VAS / Probe.orbit.GM + 1;
	// Departure Eccentricity
	Probe2.orbit.e = Per*VDS*VDS / Probe2.orbit.GM + 1;

	// Arrival Semi - latus Rectum or Parameter
	Probe.orbit.p = Probe.orbit.a*(1 - pow(Probe.orbit.e, 2));
	// Departure Semi - latus Rectum or Parameter
	Probe2.orbit.p = Probe2.orbit.a*(1 - pow(Probe2.orbit.e, 2));

	// Arrival Angular Momentum
	HA = sqrt(Probe.orbit.p*Probe.orbit.GM);
	// Departure Angular Momentum
	HD = sqrt(Probe2.orbit.p*Probe2.orbit.GM);

	// Unit Angular Momentum Vector
	uH = VA.cross(VD) / VAS / VDS / sin(alpha);

	// Arrival Angular Momenum Vector
	HAV = HA * uH;
	// Departure Angular Momentum Vector
	HDV = HD * uH;

	// Longitude of Ascending Node
	L = atan2(uH(0), -uH(1));
	Probe.orbit.loan = L;
	Probe2.orbit.loan = L;

	// Inclination
	I = atan2(sqrt(uH(0)*uH(0) + uH(1)*uH(1)), uH(2));
	Probe.orbit.I = I;
	Probe2.orbit.I = I;

	// Laplace - Range - Lenz vector for Arrival EvA
	EvA = VA.cross(HAV) + Probe.orbit.GM / VAS * VA;
	EA = EvA.norm();
	// Laplace - Range - Lenz vector for Departure EvD
	EvD = VD.cross(HDV) - Probe2.orbit.GM / VDS * VD;
	ED = EvD.norm();

	// Arrival  Argument of Perigee wA needs to be calculated
	sinwA = EvA(2) / EA / sin(I);
	coswA = (EvA(0) + EA * sinwA * cos(I) * sin(L)) / EA / cos(L);
	Probe.orbit.aop = atan2(sinwA, coswA);

	// Departure Argument of Perigee wD needs to be calculated
	sinwD = EvD(2) / ED / sin(I);
	coswD = (EvD(0) + ED * sinwD * cos(I) * sin(L)) / ED / cos(L);
	Probe2.orbit.aop = atan2(sinwD, coswD);


	// Initialise True Anomalies to equal zero
	Probe.orbit.ta0 = 0;
	Probe2.orbit.ta0 = 0;

	// Inialise Epochs to Zero
	Probe.orbit.epoch = 0;
	Probe2.orbit.epoch = 0;
	return 0;
}

////# Calculate_Departure_Velocity based on solution calculated by Compute_Hyperbola(should equal VD)
Vector3d  Connecting_Hyperbola::Calculate_Departure_Velocity()
// # Calculate_Departure_Velocity based on solution calculated by Compute_Hyperbola(should equal VD)
// #
//# INPUT:
//#
//# Connecting_Hyperbola Object with Computed Hyperbolas
//#
//# OUTPUT :
//#
//# Vdep : Departure Velocity should equal VD
//#
{
	double Vmag;
	double Periapsis;
	double Beta;
	double VPer;
	double e;
	double costheta, sintheta;
	double edash;
	double costhetadash, sinthetadash;
	double Vdep2;
	double Vdepmag;
	Vector3d Vdep;
	Vector3d Vtemp;

	Matrix3d T;

	// Compute Arrival Speed

	Vmag = VA.norm();

		
	// First Compute Velocity at Periapsis

	Periapsis = Per;
	Beta = beta;

	VPer = sqrt(2 * Planet.mu / Periapsis + Vmag * Vmag);

	// Calculate Arrival Eccentricity

	e = Periapsis*Vmag*Vmag / Planet.mu + 1;

	// Calculate Angle between Arrival Velocity and Periapsis

	costheta = -1 / e;

	sintheta = sqrt(1 - costheta *costheta);

	// Calculate Departing Eccentricity

	edash = Periapsis * pow((VPer + DV), 2) / Planet.mu - 1;

	// Calculate Angle Between Departure Velocity and Periapsis

	costhetadash = -1 / edash;

	sinthetadash = sqrt(1 - costhetadash *costhetadash);

	// Work out the Departure Velocity Try in the B - Plane First

	Vdep2 = pow((VPer + DV) , 2) - 2 * Planet.mu / Periapsis;

	Vdepmag = sqrt(pow((VPer + DV) , 2) - 2 * Planet.mu / Periapsis);

	Vdep << 0,
			-Vdepmag*(sintheta*costhetadash + costheta*sinthetadash), 
			-Vdepmag*(sintheta*sinthetadash - costhetadash*costheta);

	// Now Rotate to the Ecliptic

	T << cos(Beta), sin(Beta), 0,
		-sin(Beta), cos(Beta), 0,
		 0,			0,		   1;

	Vtemp = -Trans1.transpose() * Trans2.transpose() * T.transpose() * Vdep;

	Vdep = Vtemp;

	return Vdep;

} //#Calculate_Departure_Velocity

Connecting_Hyperbola::~Connecting_Hyperbola()
{
}
