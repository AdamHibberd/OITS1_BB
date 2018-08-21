#include "stdafx.h"
#include "Body.h"

Body::Body()
{
}

//# Calculate Ephemeris from Orbit
int Body::Orbit_To_Ephem(double t)
{
//# orbittoephem Calculates ephemeris ephem0 using Orbit parameters based upon the true anomaly
//#
//# INPUT:
//# t : Epoch for ephemeris ephem0.t
//#
//# OUTPUT :
//#
//# ephem0 : Ephemeris at time t
	double horizv;		// Horizontal Velocity
	double radialv;		// Radial Velocity

	// Set Transformation Matrix
	orbit.Set_Trans_PtoE();

	// Time for ephemeris data
	ephem0.t = t;

	// Radial distance
	ephem0.R = orbit.p / (1 + orbit.e * cos(orbit.ta0));

	// Position components in perifocal system
	ephem0.r(0) = ephem0.R * cos(orbit.ta0);
	ephem0.r(1) = ephem0.R * sin(orbit.ta0);
	ephem0.r(2) = 0;

	// Position components in ecliptic system
	ephem0.r = orbit.Trans_PtoE * ephem0.r;

	// Radial velocity
	radialv = sqrt(orbit.GM / orbit.p) * orbit.e * sin(orbit.ta0);

	// Horizontal velocity in perifocal system
	horizv = sqrt(orbit.GM / orbit.p) * (1 + orbit.e * cos(orbit.ta0));

	// Velocity components in perifocal system
	ephem0.v(0) = -sqrt(orbit.GM / orbit.p)*sin(orbit.ta0);
	ephem0.v(1) = sqrt(orbit.GM / orbit.p)*(orbit.e + cos(orbit.ta0));
	ephem0.v(2) = 0;
	ephem0.V = sqrt(pow(radialv,2) + pow(horizv,2));

// Velocity components in ecliptic system
	ephem0.v = orbit.Trans_PtoE * ephem0.v;

	return 0; //# orbittoephem
}


////# Function Cz using universal variable formulation
double Body::Cz( double z)
{
// # Cz Special Function of Auxiliary Variable z
//#
//# INPUT:
//#
//# z : Auxiliary variable z
//#
//# OUTPUT :
//#
//# val : value of Cz function at z
//#        

	double val; // Return Value
	
	if (abs(z) < 0.01)
		val = 1.0 / 2 - z / 4 / 3 / 2 + pow(z,2) / 6 / 5 / 4 / 3 / 2 - pow(z,3) / 8 / 7 / 6 / 5 / 4 / 3 / 2;
	else if (z <= -0.01)
		val = (1.0 - cosh(sqrt(-z))) / z;
	else
		val = (1.0 - cos(sqrt(z))) / z;
	return val;

} //# Cz

////# Function Sz using universal variable formulation
double Body::Sz(double z)
{
// # Sz Special Function of Auxiliary Variable z
//#
//# INPUT:
//#
//# z : Auxiliary variable z
//#
//# OUTPUT :
//#
//# val : value of Sz function at z
//#  

	double val; // Return Value

	if (abs(z) < 0.01)
		val = 1.0 / 3 / 2 - z / 5 / 4 / 3 / 2 + pow(z,2) / 7 / 6 / 5 / 4 / 3 / 2 - pow(z,3)/ 9 / 8 / 7 / 6 / 5 / 4 / 3 / 2 ;
	else if (z <= -0.01)
		val = (sinh(sqrt(-z)) - sqrt(-z)) / pow(sqrt(-z),3);
	else
		val = (sqrt(z) - sin(sqrt(z))) / pow(sqrt(z),3);
	return val;

} //# Sz

////# Function dCdz using universal variable formulation
double Body::dCdz( double z)
{
// # dCdz Derivative of Cz wrt z
//#
//# INPUT:
//#
//# z : Auxiliary variable z
//#
//# OUTPUT :
//#
//# val : value of Cz gradient function at z
//#

	double val; // Return Value

	if (abs(z) < 0.01)
		val = -1.0 / 4 / 3 / 2 + 2 * z / 6 / 5 / 4 / 3 / 2 - 3 * pow(z,2) / 8 / 7 / 6 / 5 / 4 / 3 / 2 + 4 * pow(z,3)/ 10 / 9 / 8 / 7 / 6 / 5 / 4 / 3 / 2 ;
	else
		val = 1.0 / 2 / z*(1.0 - z*Sz(z) - 2 * Cz(z));
	return val;

} //# dCdz

////# Function dSdz using universal variable formulation
double Body::dSdz( double z)
{
// # dSdz Derivative of Sz wrt z
//#
//# INPUT:
//#
//# z : Auxiliary variable z
//#
//# OUTPUT :
//#
//# val : value of Sz gradient function at z
//#

	double val; // Return Value

	if (abs(z) < 0.01)
		val = -1.0 / 5 / 4 / 3 / 2 + 2 * z / 7 / 6 / 5 / 4 / 3 / 2 - 3 * pow(z, 2) / 9 / 8 / 7 / 6 / 5 / 4 / 3 / 2 + 4 * pow(z, 3) / 11 / 10 / 9 / 8 / 7 / 6 / 5 / 4 / 3 / 2 ;
	else
		val = 1.0 / 2 / z*(Cz(z) - 3 * Sz(z));
	return val;

} //# dSdz

////# Compute ephem at time t
int Body::Compute_Ephem_At_t(double t, int mode, double thresh)
{
// # Compute_Ephem_At_t Calculates ephemeris ephemt at time t using Orbit parameters : method used depends upon mode
//#
//# INPUT :
//#
//# t : Time at which ephemeris ephemt.t will be calculated
//# mode : mode = 1 use Universal Formulation from Fixed Orbital Parameters
//#               : mode = 2 Use SPICE
//# thresh : time tolerance for numerical iteration required for mode = 1
// #
//# OUTPUT :
//#
//#            

	double time = t;	// Time at which ephemeris ephemt.t will be calculated
	double Dt;			// Overall change in time
	double rdotv0;		// r0.v0 = initial dot product of r0 and v0
	double xn;			// Universal Variable x
	double zn;			// Universal Variable z
	double tn;			// Time of Flight (at iteration step n)
	double dtdxn;		// Gradient of tn by dx
	double f, g;		// f and g values
	double fdot, gdot;	// Gradients of f & g
	int i;				// Iteration Step Counter
	SpiceDouble    state[6];
	SpiceDouble    et;
	SpiceDouble	lt;

	ephemt.t = t;
	et = t;      
	
	ConstSpiceChar  * versn;
	if (Fixed_Point > 0)
	{
		ephemt = ephem0;
		ephemt.t = t;
		ephemt.R = ephemt.r.norm();
		ephemt.V = ephemt.v.norm();
		return 0;
	}

	if (mode == 2)
	{
		
		spkezr_c(ID, et, "ECLIPJ2000", "NONE", "SUN", state, &lt);
		
		
		for (i = 0 ; i < 3; i++)
		{
			ephemt.r(i) = state[i]*1e+3;
			ephemt.v(i) = state[i+3]*1e+3;
		}
		ephemt.R = ephemt.r.norm();
		ephemt.V = ephemt.v.norm();
		return 0;
	}
	else

		// Overall change in time
		Dt = time - ephem0.t;
		if (abs(Dt) <= thresh)
		{
			ephemt = ephem0;
			return 0;
		}
		// Correct for Orbital Time Period if necessary

		// Dt = Dt - orbit.TP*fix(Dt / orbit.TP);

		// r0.v0 = initial dot product of r0 and v0
		rdotv0 = ephem0.r.dot(ephem0.v);

		// Initial Guess for x(universal variable)
		if (orbit.e <= 1)
			xn = sqrt(orbit.GM)*orbit.arec*Dt;
		else
			xn = Sign(Dt)*sqrt(-1 / orbit.arec)*log(-2 * orbit.arec*orbit.GM*Dt / (rdotv0 + Sign(Dt)*sqrt(-orbit.GM / orbit.arec)*(1 - orbit.arec*ephem0.R)));
		// Initial value of z (universal variable)
		zn = pow(xn, 2) * orbit.arec;
		// Initial guess for t, time of flight
		tn = rdotv0*pow(xn, 2) * Cz(zn) / orbit.GM + (1 - orbit.arec*ephem0.R)*pow(xn, 3) * Sz(zn) / sqrt(orbit.GM) + ephem0.R*xn / sqrt(orbit.GM);
		// Initial guess for gradient dt by dx
		dtdxn = rdotv0*xn*(1 - zn*Sz(zn)) / orbit.GM + pow(xn, 2) * Cz(zn) / sqrt(orbit.GM) + ephem0.R * (1 - zn*Cz(zn)) / sqrt(orbit.GM);

		// Do Newton iteration
		i = 0;
		while (abs(tn - Dt) > thresh)
		{
			i++;
			if (i > 1000)	break;
			xn = xn + (Dt - tn) / dtdxn;
			zn = pow(xn, 2) * orbit.arec;
			tn = rdotv0*pow(xn, 2) * Cz(zn) / orbit.GM + (1 - orbit.arec*ephem0.R)*pow(xn, 3) * Sz(zn) / sqrt(orbit.GM) + ephem0.R*xn / sqrt(orbit.GM);
			dtdxn = rdotv0*xn*(1 - zn*Sz(zn)) / orbit.GM + pow(xn, 2) * Cz(zn) / sqrt(orbit.GM) + ephem0.R * (1 - zn*Cz(zn)) / sqrt(orbit.GM);
		}

		// Compute f & g for position vector at time Dt

		f = 1 - pow(xn,2) * Cz(zn) / ephem0.R;
		g = Dt - pow(xn,3) * Sz(zn) / sqrt(orbit.GM);

		ephemt.r = f * ephem0.r + g * ephem0.v;
		ephemt.R = ephemt.r.norm();

		// Compute fdot & gdot for velocity vector at time Dt
	
		gdot = 1 - pow(xn, 2) * Cz(zn) / ephemt.R;
		fdot = sqrt(orbit.GM) / ephemt.R / ephem0.R*xn*(zn*Sz(zn) - 1);

		ephemt.v = fdot * ephem0.r + gdot * ephem0.v;
		ephemt.V = ephemt.v.norm();

		// Compute true anomaly

		orbit.ta = atan2(ephemt.r.dot(ephemt.v) / ephemt.R*sqrt(orbit.p / orbit.GM), orbit.p / ephemt.R - 1);

		// Correct for Orbital Time Period for Elliptical Orbits
		return 0;

} //# compute_ephem_at_t

////# calculate_orbit_from_ephem : Calculates orbit from ephemeris ephemt

int Body::Calculate_Orbit_From_Ephem(double time)
{
// # Calculates osculating orbital parameters from cartesians
//#
//# INPUT:
//#
//# t : Epoch for ephemeris ephemt.t
//#
//# OUTPUT :
//#
//#    

	double mu;		// Gravitational Mass
	Vector3d x;		// Position Vector
	Vector3d v;		// Velocity Vecor
	Vector3d h;		// Angular Momentum Vector
	Vector3d Ev;	// Laplace - Range - Lenz vector
	double R, V;	// Respectively Radial Distance and Speed
	double H;		// Magnitude of Angular Momentum Vector
	double E;		// Magnitude of Ev
	double L;		// Longitude of Ascending Node
	double I;		// Inclination of Orbit (rads)
	double sinw;	// Sine of Argument of Perigee
	double cosw;	// Cosine of Argument of Perigee
	double w;		// Argument of Perigee (rads)
	double Energy;	// Energy per unit Mass
	double e;		// Eccentricity of Orbit
	double rxd, ryd, rydd; // Used to Calculate True Anomaly
	double costa;	// Cosine True Anomaly
	double sinta;	// Sine True Anomaly
	double true_anom; // True Anomaly (rads)
	double X, Y;	// Used to calculate True Anomaly
	mu = orbit.GM;

	// Epoch for orbit

	orbit.epoch = time;

	// Radial Distance

	x = ephemt.r;
	v = ephemt.v;
	R = x.norm();

	// Speed

	V = v.norm();

	// Angular Momentum Vector h

	h = x.cross(v);

	// Magnitude of angular momentum H

	H = h.norm();

	// Longitude of Ascending Node

	L = atan2(h(0), -h(1));
	orbit.loan = L;

	// Inclination

	I = atan2(sqrt(pow(h(0), 2) + pow(h(1), 2)), h(2));
	orbit.I = I;

	// Laplace - Range - Lenz vector Ev

	Ev = v.cross(h) - mu / R * x;

	E = Ev.norm();

	// Argument of Perigee w needs to be calculated

//	sinw = Ev(2) / E / sin(I);

//	cosw = (Ev(0) + E * sinw * cos(I) * sin(L)) / E / cos(L);

//	w = atan2(sinw, cosw);
	w = atan2(Ev(2) * cos(L), ( sin(I) * Ev(0) + Ev(2) * cos(I) * sin(L) ));
	orbit.aop = w;

	// Energy per unit mass gives semi - major axis a

	Energy = V * V / 2 - mu / R;

	if (abs(Energy) > 0)
		orbit.a = -mu / 2 / Energy;

	orbit.arec = -Energy * 2 / mu;

	// Eccentricity, e

	e = sqrt(1 + 2 * Energy * pow((H / mu), 2));
	orbit.e = e;

	// Parameter, p(semi = latus rectum)

	orbit.p = pow(H, 2) / mu;

	// Time Period
	if (e < 1)
		orbit.TP = 2 * PI_OITS*sqrt(pow(orbit.a, 3) / mu);

	// True Anomaly

	rxd = x(0) * cos(L) + x(1) * sin(L);
	ryd = x(1) * cos(L) - x(0) * sin(L);
	rydd = ryd * cos(I) + x(2) * sin(I);

//	sinta = (rydd * cosw - rxd * sinw) / R;
//	costa = (rydd * sinw + rxd * cosw) / R;

//	true_anom = atan2(sinta, costa);

	Y = rydd * Ev(0)* sin(I) + rydd * Ev(2) * cos(I)*sin(L) - rxd * Ev(2)*cos(L);
	X = rydd * Ev(2)* cos(L) - rxd * Ev(2) * cos(I)*sin(L) - rxd * Ev(0)*sin(I);

	true_anom = atan2(Y, X);
	
	orbit.ta = true_anom;
	Calculate_True_Anomaly();
	return(0);

}//# calculate_orbit_from_ephem


 // Calculate Sphere of Influence

int Body::Sphere_Of_Influence()
{

	SpoI = orbit.a*pow(mu / orbit.GM,0.4);

	return 0;

}
int Body::Calculate_True_Anomaly()
{
	double MU;		// Centre of Attraction
	Vector3d x;		// Position Vector
	Vector3d v;		// Velocity Vecor
	Vector3d h;		// Angular Momentum Vector
	double R, H, p;	
	double Radial_Velocity;
	double true_anom;

	// Gravitational Mass
		MU = orbit.GM;

	// Radial Distance
		x = ephemt.r;
		v = ephemt.v;
		R = x.norm();

	// Angular Momentum Vector h

		h = x.cross(v);

	// Magnitude of angular momentum H

		H = h.norm();

	// Parameter, p(semi = latus rectum)

		p = H * H / MU;

		Radial_Velocity = x.dot(v) / R;


			true_anom = atan2(R*Radial_Velocity*sqrt(p / MU), (p - R));
			orbit.ta = true_anom;
		
		return 0;
}


int Body::Sign(double input)
{
	int output;
	if (input < 0.0)
		output = -1;
	else
		output = 1;
	return output;
}
Body::~Body()
{
}
