#pragma once
#include <Eigen/Dense>

using namespace Eigen;

class Ephemeris
{
public:
	double t;				// Time
	Vector3d r;				// Cartesian Position Vector
	double R;				// Radial Distance i.e. magnitude of r
	Vector3d v;				// Cartesian Velocity Vector
	double V;				// Speed i.e. magnitude of v
	Ephemeris();
	~Ephemeris();
};

