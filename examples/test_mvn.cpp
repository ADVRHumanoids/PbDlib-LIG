/**
Copyright (C) 2014, Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon

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

/*! \file learn_gmm.cpp
\brief Learning GMM model
Learning a GMM model from a demonstration saved in the file data_txyz.txt

\author Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon
\bug No known bugs.
*/

#include "pbdlib/mvn.h"
#include <sstream>

#define nbVar 4

using namespace pbdlib;

int main(int argc, char **argv)
{
	// Make new distribution
	//
	colvec Mu;
	mat SIGMA;

	SIGMA << 0.0453 << 0.0021 << 0.0682 << endr
	   	  << 0.0021 << 0.4463 << 0.5093 << endr
		  << 0.0682 << 0.5093 << 0.8964; 
	SIGMA = SIGMA;
	Mu << 0.4843 << 0.3222 << 0.5827;	
	GaussianDistribution G1(Mu,SIGMA);

	cout << "Mu: " << endl << G1.getMU() << endl;
	cout << "Sigma: " << endl << G1.getSIGMA() << endl;
	cout << "Lambda: " << endl << G1.getLAMBDA() << endl;

	cout << "getPDF(Mu): " << endl <<
		G1.getPDFValue(Mu);

	// Stochastic sampling
	mat samps = G1.stochasticSampling(20);
	cout << "Samples: " << endl;
	for (int i = 0;i<20;i++)
		cout << samps.col(i).t() << endl;

	cout << "Probs of samples" << endl;
	for (int i = 0;i<20;i++)
		cout << G1.getPDFValue(samps.col(i)) << endl;
}

