/**
Copyright (C) 2014, Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso, Sylvain Calinon

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

/*! \file gmr.h
\brief gmr class
This class is the implementation of the Gaussian mixture regression (GMR), as described in:

@article{Calinon07SMC,
	title="On Learning, Representing and Generalizing a Task in a Humanoid Robot",
	author="S. Calinon and F. Guenter and A. Billard",
	journal="IEEE Transactions on Systems, Man and Cybernetics, Part B. 
	Special issue on robot learning by observation, demonstration and imitation",
	year="2007",
	volume="37",
	number="2",
	pages="286--298"
}

\author Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso, Sylvain Calinon
\bug No known bugs.
*/

#ifndef GMR_H
#define GMR_H

#include "pbdlib/gmm.h"
#include "armadillo"

using namespace arma;

namespace pbdlib
{

class GMR
{
	private:
		GMM_Model *gmm;
		Datapoints *data_in;
		/*!
		COMPONENTS: A vector of Gaussian components, each component has PRIORS, MU and SIGMA
		*/
		std::vector<GaussianDistribution> COMPONENTS;

		// Regression Variables:
		uint i,t;
		colvec Mu, MuOut, MuOutTmp;
		mat Sigma, SigmaOut, SigmaOutTmp, beta;
		mat InvSigmaInIn; 
		mat Pxi;	
		GaussianDistribution* Gtmp;


	public:
		GMR(){}
		// The GMM_Model should be already learnt. (SIGMA, MU, PRIORS)
		GMR(GMM_Model* model);
		/*!
		setGMMModel(): Set the GMM model that will be used for regression
		*/
		void setGMMModel(GMM_Model* gmmmodel);
		/*!
		regression(): Gaussian mixture regression
		*/
		void regression(GMM_Model*, mat Data_in, urowvec Var_In, urowvec Var_Out);
		GMM_Model* regression(Datapoints* data_in);
		GMM_Model* regression(mat Data_in, urowvec Var_In, urowvec Var_Out);

};

} //end pbdlib namespace

#endif


