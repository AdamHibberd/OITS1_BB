#pragma once
#include <Eigen/Dense>
#ifndef _GREAT_BIG_CENTRE_
#define _GREAT_BIG_CENTRE_ 1.32712440018e20 //Assumes the Sun
#endif

using namespace Eigen;

class Orbit
{
public:

	double epoch;			// epoch = time of orbital parameters
	double	ta0;            // ta0 = true anomaly at epoch
	double	ta;             // ta = true anomaly
	double	p;              // p = parameter or semi - latus rectum
	double	a;              // a = semi - major axis
	double	arec;           // arec = 1 / a
	double	e;              // e = eccentricity
	double	aop;            // aop = argument of perigee
	double	loan;           // loan = longitude of ascending node
	double	I;              // I = inclination
	double	GM = _GREAT_BIG_CENTRE_; // GM = gravitational mass of attracting centre
	double	TP;             // TP = Orbital Time period for elliptical orbits
	Matrix3d Trans_PtoE;	// Trans_PtoE = Perifocal to Ecliptic Transformation Matrix
	Orbit();
	~Orbit();
	int Set_Trans_PtoE();
};

