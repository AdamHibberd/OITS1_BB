#include "stdafx.h"
#include "Orbit.h"


Orbit::Orbit()
{
}


Orbit::~Orbit()
{
}

int Orbit::Set_Trans_PtoE()
{
	Matrix3d  Trans_OtoE;
	Matrix3d  Trans_PtoO;
	
	Trans_OtoE(0, 0) = cos(loan);
	Trans_OtoE(0, 1) = -sin(loan)*cos(I);
	Trans_OtoE(0, 2) = sin(I)*sin(loan);
	Trans_OtoE(1, 0) = sin(loan);
	Trans_OtoE(1, 1) = cos(loan)*cos(I);
	Trans_OtoE(1, 2) = -sin(I)*cos(loan);
	Trans_OtoE(2, 0) = 0;
	Trans_OtoE(2, 1) = sin(I);
	Trans_OtoE(2, 2) = cos(I);
	Trans_PtoO(0, 0) = cos(aop);
	Trans_PtoO(0, 1) = -sin(aop);
	Trans_PtoO(0, 2) = 0;
	Trans_PtoO(1, 0) = sin(aop);
	Trans_PtoO(1, 1) = cos(aop);
	Trans_PtoO(1, 2) = 0;
	Trans_PtoO(2, 0) = 0;
	Trans_PtoO(2, 1) = 0;
	Trans_PtoO(2, 2) = 1;

	Trans_PtoE = Trans_OtoE * Trans_PtoO;
	return 0;
}
