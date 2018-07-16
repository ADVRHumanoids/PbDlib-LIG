/**
Copyright (C) 2014, Danilo Bruno, Sylvain Calinon

This file is part of PbDLib (Programming-by-demonstration C++ Library).

    PbDLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PbDLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PbDLib.  If not, see <http://www.gnu.org/licenses/>.
*/

/*! \file lqr.h
\brief lwr class
	This class is the implementation of a minimal intervention controller with 
	linear quadratic regulation, as described in:

	@inproceedings{Calinon14ICRA,
		author="Calinon, S. and Bruno, D. and Caldwell, D. G.",
		title="A task-parameterized probabilistic model with minimal intervention control",
		booktitle="Proc. {IEEE} Intl Conf. on Robotics and Automation ({ICRA})",
		year="2014",
		month="May-June",
		address="Hong Kong, China",
		pages="3339--3344"
	}

\author Danilo Bruno, Sylvain Calinon
\bug No known bugs.
*/

#ifndef LQR_H
#define LQR_H

//#include "pbdlib/gmm.h"
//#include "pbdlib/gmr.h"
//#include "datapoints.h"
#include "armadillo"
#include <vector>

using namespace arma;

namespace pbdlib
{

class LQR
{
	private:
		mat A;				//Matrix defining linear system
		mat B;				//Matrix defining linear system
		double dt;			//time step
		bool prob_set; 		//flag to 1 if problem matrices are set

		std::vector<mat> Q;	//Matrix defining variability
		mat R;				//Matrix defining command weights
		mat Target;			//Target trajectory

		int nbData;			//number of Data Points to evaluate
		int nbVar;			//Number of variables
		int nbCtrl;			//Number of controller variables;
		std::vector<mat> S;	//Riccati solution
		mat d;				//Linear DE solution

		std::vector<mat> L;	//LQR Gains
		mat M;				//Feedforward term

	public:
		LQR(mat,mat,double);
		//LQR(GMM_Model* model);
		//mat reproduction_finiteHorizon(mat,mat,colvec);
		bool evaluate_gains_finiteHorizon(mat,colvec);
		bool evaluate_gains_infiniteHorizon();
		//Datapoints* reproduction_finiteHorizon(Datapoints*, Datapoints*, Datapoint);
		//Datapoints* reproduction_infiniteHorizon(Datapoints*, Datapoints*, Datapoint);

		void setQ(std::vector<mat> _Q){this->Q = _Q;};
		void setR(mat _R){this->R = R;};
		void setTarget(mat _Target){this->Target = _Target;};

		std::vector<mat> getGains(){return this->L;};
		mat getFF(){return this->M;};

		bool setProblem(mat,std::vector<mat>,mat);

		mat diff(mat);
		mat solveAlgebraicRiccati(mat,mat,mat,mat);
};

} //end pbdlib namespace

#endif



