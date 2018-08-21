#include "stdafx.h"
#include "Transfer_Orbit.h"


Transfer_Orbit::Transfer_Orbit()
{
}

//# Calculate Transfer Orbit from Departure Body to Body 2
//# Uses Universal Variables Formulation as described in Fundamentals
//# of Astrodynamics, Bate, Mueller, White
int Transfer_Orbit::Calculate_Transfer(double tdin, double tarin, double thresh, int itmax, int wayflag)
// # Calculate_transfer Calculates 2 transfer orbits(long way & short way) between 2 Bodies
//#
//# INPUT:
//#
//# obj : Current Transfer_orbit object in question
//# td : Departure Time from bodyd(Departure Body)
//# tar : Arrival Time at bodya(Arrival Body)
//# thresh : Threshold for iteration process
//# itmax : Maximum number of iteration cycles before giving up
//# wayflag : = 0 return short way and longway
//#               : = 1 return prograde only
//#
//# OUTPUT 
//#
//# obj : Transfer_orbit object with departure and arrival ephemeris and orbit calculated
//#      
{
double rd;			// Departure radial distance
double ra;			// Arrival radial distance
double GM;			// Sun's Gravitational Mass
double cosdta;		// Cosine of Change in true anomaly
double dta;			// Change in True Anomaly
double A;			// Constant
double z0, zn;	// Values of the Universal Variable z
double yn, xn;		// Values of Univeral Variables y & x
double Sn, Cn;		// Values of Special Functions of z
double dSdzn, dCdzn; // Values of Gradient of Special Functions of z
double tn;			// Calculation of Time of Flight on the nth iteration
double dtdzn;		// Gradient of tn w.r.t. zn
double dtdznold;	// Value of dtdzn on last cycle of iteration
double dzn;			// Value of change of z for nth iteration step
double znmin=0;		// Minimum value of zn
double deltatnew;	// Difference between Required ToF and Calculated ToF
double deltatold;	// Difference between Required ToF and Calculated ToF
double deltatime;	// Change of delta t
double factor;		// Reduction factor in dzn
double f, g;		// Scalar Coefficients for Arrival Velocity Vector
double gdot;		// Rate of change of  g
int errflag, newtry, exitflag; // Flags for Newton Raphson Iteration 
int i, j;			// Iteration Step Counters

	nmax = itmax;

	// Set Departure and Arrival times amd Positions
	td = tdin;
	tar = tarin;
	for (i = 0; i < 2; i++)
	{
		ephemd[i] = bodyd.ephemt;
		ephema[i] = bodya.ephemt;
	}

	rd = bodyd.ephemt.R;
	ra = bodya.ephemt.R;

	// Sun's Gravitational Mass
	GM = bodya.orbit.GM;

	// Start to iterate to Solution based on j = 1 (one way) then
	//  try j = 2 (opposite way)

	for (j = 0; j < 2; j++)
	{
		// Change in true anomaly
		cosdta = bodya.ephemt.r.dot(bodyd.ephemt.r) / ra / rd;

		dta = acos(cosdta);

		if (j == 0) dta = 2 * PI_OITS - dta;

		// Calculate constant A
		A = sqrt(rd*ra)*sin(dta) / sqrt(1 - cos(dta));

		if (A > 0 && dta < PI_OITS)
		{
			z0 = 0;
			znmin = fzero(z0, 1e-13, A, 1000);
			znmin = znmin + 1e-11;
			zn = znmin;
		}
		else
		{
			// Guess inital value of z
			zn = pow((3 / 2 * PI_OITS), 2);
		}

		// Calculate Special Functions of zn and their gradient wrt zn
		Sn = bodya.Sz(zn);
		Cn = bodya.Cz(zn);
		dSdzn = bodya.dSdz(zn);
		dCdzn = bodya.dCdz(zn);

		// Calculate special universal variables for initial guess
		yn = ra + rd - A * (1 - zn*Sn) / sqrt(Cn);
		xn = sqrt(yn / Cn);

		// Calculate corresponding Time - of - Flight for initial guess
		tn = (pow(xn , 3) * Sn + A * sqrt(yn)) / sqrt(GM);

		// Calculate Gradient of tn wrt zn
		dtdzn = (pow(xn , 3) * (dSdzn - 3 * Sn*dCdzn / 2 / Cn) + A / 8 * (3 * Sn*sqrt(yn) / Cn + A / xn)) / sqrt(GM);

		// Calculate initial step dzn
		dzn = (tar - td - tn) / dtdzn;

		// Inialise Error in ToF wrt Required ToF
		deltatold = tar - td - tn;

		// Iniialise Variables for Newton interation
		i = 0;
		factor = 1;
		errflag = 0;
		newtry = 0;
		exitflag = 0;

		// Do Newton Iteration on zn
		while (exitflag == 0)
		{

			exitflag = (abs(tn - tar + td) < thresh);

			i++;

			// Break out if not converging
			if (i > itmax)	break;

			// Don't allow factor to get to small - try a different
			// guess for dzn
			if (abs(factor*dzn) < thresh / abs(dtdzn))
			{
				errflag = 1;
				factor = 1;
			}

			// Normal change in zn
			if (errflag == 0)
			{
				dzn = factor*deltatold / dtdzn;
			}
			// Else if Not Converging try a totally different guess
			else
			{
				zn = 0;
				dzn = pow((rand() * 2 * PI_OITS) , 2);
				errflag = 0;
				newtry = 1;
			}

			// Iterate to solution
			zn = zn + dzn;

			// Don't allow zn to get to large or small
			if (zn > pow((2 * PI_OITS) , 2))	zn = pow((2 * PI_OITS) , 2) - 0.01; // LOOK HERE !!!!!!


			// Calculate Special Functions of zn and their gradient wrt
			// zn
			Sn = bodya.Sz(zn);
			Cn = bodya.Cz(zn);
			dSdzn = bodya.dSdz(zn);
			dCdzn = bodya.dCdz(zn);


			// Calculate Universal Variable yn based on zn
			yn = ra + rd - A * (1 - zn*Sn) / sqrt(Cn);

			// Don't allow yn to become -ve
			if (yn < 0)
			{
				zn = zn - dzn;
				factor = factor*0.1;
				continue;
			}

			// Calculate Universal Variable xn based on yn and zn
			xn = sqrt(yn / Cn);

			// Calculate ToF
			tn = (pow(xn, 3) * Sn + A * sqrt(yn)) / sqrt(GM);

			// Calculate Error in ToF wrt required

			deltatnew = tar - td - tn;

			// Determine Change in ToF wrt previous iteration
			deltatime = deltatnew - deltatold;

			// If ToF has reached a local minimum try a new guess
			// for dzn
			if ((abs(deltatime) < thresh / 10) & (newtry == 0))
			{
				errflag = 1;
				continue;
			}
			else
			{
				errflag = 0;
			}

			// Don't allow the new guess to be further away from the
			// solution than the old guess
			//  if (abs(deltatold) - abs(deltatnew))<0 & newtry == 0
				// '3'
				//      rateflag = 1;
			//      zn = zn - dzn;
			//     factor = factor*0.1;
			//     continue;
			//  end

			//			rateflag = 0;

			// Don't go backwards in time
			if (tn < 0)
			{
				zn = zn - dzn;
				factor = factor*0.1;
				continue;
			}

			deltatold = deltatnew;
			dtdznold = dtdzn;

			// Calculate Gradient of ToF wrt zn
			dtdzn = (pow(xn, 3) * (dSdzn - 3 * Sn*dCdzn / 2 / Cn) + A / 8 * (3 * Sn*sqrt(yn) / Cn + A / xn)) / sqrt(GM);

			// Don't allow Gradient of ToF to become -ve
			//       if (dtdzn<0 & newtry == 0)
				// '5'
				//           zn = zn - dzn;
			//           factor = factor*0.1;
			//            dtdzn = dtdznold;
			//       end

			// Normal Operation of iteration
			newtry = 0;

			}


			// Update Number of iterations
			ni[j] = i;
//			cout << ni[j] << endl;
			// Compute f & g & gdot for velocity vectors time td and ta

			f = 1 - pow(xn, 2) * Cn / rd;
			f = 1 - yn / rd;
			g = tar - td - pow(xn, 3) * Sn / sqrt(GM);
			g = A * sqrt(yn / GM);
			gdot = 1 - pow(xn, 2) * Cn / ra;
			gdot = 1 - yn / ra;

			// Velocity vector at departure

			ephemd[j].v = (ephema[j].r - f * ephemd[j].r) / g;
			ephemd[j].V = ephemd[j].v.norm();

			// Velocity vector at arrival

			ephema[j].v = (-ephemd[j].r + gdot * ephema[j].r) / g;
			ephema[j].V = ephema[j].v.norm();

			// Calculate Transfer Orbit

			transfer_body[j].ephem0 = ephemd[j];

			transfer_body[j].ephemt = ephemd[j];
			transfer_body[j].Calculate_Orbit_From_Ephem(td);
			true_anom_dep[j] = transfer_body[j].orbit.ta;

			transfer_body[j].ephemt = ephema[j];
			transfer_body[j].Calculate_Orbit_From_Ephem(tar);
			true_anom_arr[j] = transfer_body[j].orbit.ta;

	}

	// Take action if Prograde to be considered only
	if (wayflag == 1)
	{
		// Firstly determine whether angle between first Transfer Orbit velocity and Body's velocity is > PI/2 (i.e. retrograde)
		if (ephemd[0].v.topRows(2).dot(bodyd.ephemt.v.topRows(2)) < 0)
		{
			// If so then determine whether second Transfer Velocity is prograde
			
			if (ephemd[1].v.topRows(2).dot(bodyd.ephemt.v.topRows(2)) > 0)
			{
				// If so then set first way to second way
				ephemd[0] = ephemd[1];
				ephema[0] = ephema[1];
				transfer_body[0] = transfer_body[1];
				true_anom_dep[0] = true_anom_dep[1];
				true_anom_arr[0] = true_anom_arr[1];
			}
			// If both ways are retrograde then do nothing
		}
		// Now if first way is prograde determine if second way is retograde

		else if (ephemd[1].v.topRows(2).dot(bodyd.ephemt.v.topRows(2)) < 0)
		{
			// If so then set second way equal to first way
			ephemd[1] = ephemd[0];
			ephema[1] = ephema[0];
			transfer_body[1] = transfer_body[0];
			true_anom_dep[1] = true_anom_dep[0];
			true_anom_arr[1] = true_anom_arr[0];
		}
	}
	return 0;
}
////# YN Calculates auxiliary variable y from z and other data
double Transfer_Orbit::fzero(double z0, double ztol, double A, int Maxit)
{
// # YN Calculates auxiliary variable y from z and other data
//#
//# INPUT:
//#
//# zn : Current value of auxiliary variable zn
//#
//# OUTPUT :
	//#
	//# fn : Value of y
	//#      
	int i;
	double S1,S2;
	double C1,C2;
	double rd = bodyd.ephemt.R;
	double ra = bodya.ephemt.R;
	double dz;
	double dzold = 0.0;
	double func1, func2, dfunc, dfuncdz,grad;
	double zn = z0;
	double dzg = ztol / 10;
	double reduce = 1.0;
	
	for (i = 0; i < Maxit; i++)
	{
		S1 = bodya.Sz(zn); 
		S2 = bodya.Sz(zn + dzg);
		C1 = bodya.Cz(zn);
		C2 = bodya.Cz(zn + dzg);
//		cout << " A = " << A << " S1 = " << S1 << " C1 = " << C1 << endl;
		func1 = rd + ra - A*(1.0 - zn*S1) / sqrt(C1);
		func2 = rd + ra - A*(1.0 - (zn + dzg)*S2) / sqrt(C2);
		dfunc = func2 - func1;
		grad = A / 4.0 * sqrt(C1) ;
		dfuncdz = dfunc / dzg;

		dz =  reduce*func1 / grad;
		
		zn = zn - dz;

		if ((dzold * dz) < 0.0)
			reduce = reduce*0.1;
		else
		{
			reduce = 1.0;

			//		cout << "func = " << func1 << " zn = " << zn << " grad " << grad << " dfuncdz = " << dfuncdz << endl;

			if (abs(dz) < ztol)break;
		}
		dzold = dz;
	}
//	cout << "i= " << i << endl;
	return zn;
}


//# Function to calculate the closes approach of the transfer orbit to the Sun during its passage
int Transfer_Orbit::Calculate_Perihelion()
{
	int j;		// Index

	// Treat each transfer orbit in turn
	for (j = 0; j < 2; j++)
	{
		// Firstly Determine if there has been a change in sign of true anomaly
		if (true_anom_dep[j] * true_anom_arr[j] < 0)
		{
			//  .... and the body is approaching the sun
			if (ephemd[j].r.dot(ephemd[j].v) < 0)
			{
				// Set closest approach as Perhielion of Orbit
				perihelion[j] = transfer_body[j].orbit.a * (1 - transfer_body[j].orbit.e);
			}
			//  .... and the body is going away from the sun
			else
			{
				// Set closest approach to minimum of radii of end points
				perihelion[j] = min(ephemd[j].R, ephema[j].R);
			}
		}

		// IF there has not been a change in sign of true anomaly 

		else
		{
			// .... If the orbit is elliptical 
			if (transfer_body[j].orbit.e < 1)
			{
				// ... If the body has spent more than half the time-period in the elliptical orbit then
				if (abs(tar - td) > transfer_body[j].orbit.TP / 2)
				{
					// Set closest approach as Perhielion of Orbit
					perihelion[j] = transfer_body[j].orbit.a * (1 - transfer_body[j].orbit.e);
				}
				else
				{
					// Set the Closest approach to be the minimum of the radial distances at Arrival & Departure
					perihelion[j] = min(ephemd[j].R, ephema[j].R);
				}
			}
			// .... If the Orbit is Hyperbolic		
			else
			{
				// Set the Closest approach to be the minimum of the radial distances at Arrival & Departure
				perihelion[j] = min(ephemd[j].R, ephema[j].R);
			}
		}
	}
	return 0;
}

Transfer_Orbit::~Transfer_Orbit()
{
}
