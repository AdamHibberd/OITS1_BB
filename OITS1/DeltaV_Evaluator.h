#pragma once
#include "Project.h"
#include "c:\Users\adamh\Documents\nomad.3.8.1\builds\VisualStudio\src\Evaluator.hpp"


class DeltaV_Evaluator :
	public NOMAD::Evaluator
{

public:

	class Project * PROJ;
	
//	DeltaV_Evaluator();
	
	DeltaV_Evaluator(const NOMAD::Parameters &t );
	bool eval_x(NOMAD::Eval_Point &x, const NOMAD::Double &h_max, bool & count_eval) const;
	~DeltaV_Evaluator();
};

