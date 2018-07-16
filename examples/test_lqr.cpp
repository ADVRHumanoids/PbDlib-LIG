
/**
Copyright (C) 2014, Danilo Bruno

This file is part of PbDLib.

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

/*! \file test_lqr
\brief Testing lqr

\author Danilo Bruno
\bug No known bugs.
*/

#include <iostream>
#include <sstream>
#include "pbdlib/lqr.h"
#include "armadillo"
#include <vector>


using namespace pbdlib;
using namespace arma;

int main(int argc, char **argv)
{
	mat A,B;
	float dt,rFactor;
	int nbData = 400;
	mat R,Target;
	std::vector<float> DataIn;
	std::vector<mat> Q;

	A << 0 << 1 << endr << 0 << 0 << endr; 	//matrix for second order system
	B << 0 << endr << 1 << endr;				//matrix for force controller

	dt = 0.01;
	rFactor = 0.1;
	R = rFactor;

	Target = zeros(2,nbData);

	Target.row(0) = 100*ones(1,nbData);

	for (int t=0;t<nbData;t++){
		DataIn.push_back(t*dt);
		mat QTmp = zeros(2,2);
		QTmp(0,0) = 1;
		Q.push_back(QTmp);
	}

	LQR test_lqr(A,B,dt);

	test_lqr.setProblem(R,Q,Target);

	//Test infinite horizon
	test_lqr.evaluate_gains_infiniteHorizon();

	//Test finite horizon

	//mat S;
	//colvec d;

	//S = zeros(2,2);
	//d = zeros(2,1);

	//test_lqr.evaluate_gains_finiteHorizon(S,d);

	//Output data
	for (int t=0;t<nbData;t++){
		std::cout << t << " , " << test_lqr.getGains().at(t) << std::endl;
	}
	std::cout << test_lqr.getFF() << std::endl;

	return 0;
}
