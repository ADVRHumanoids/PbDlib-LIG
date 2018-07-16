/**
Copyright (C) 2014, Milad S. Malekzadeh, Sylvain Calinon

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

/*! \file test_rewardWeightedRefinement
\brief Testing Expectation Maximization (EM) based reward-weighted refinement

The "reward" function can be defined based on problem.

\author Milad S. Malekzadeh, Sylvain Calinon
\bug No known bugs.
*/

#include "armadillo"
#include "pbdlib/mvn.h"

using namespace pbdlib;
using namespace arma;

// Reward function
rowvec reward( const mat& sample )
{
	vec mu = 8 * ones<vec>(2,1);
	mat sigma = 5 * eye<mat>(2,2);
	GaussianDistribution *gd;
	gd = new GaussianDistribution(mu, sigma);

	return trans(gd->getPDFValue( sample )) / as_scalar(gd->getPDFValue(mu));
}

int main(int argc, char **argv)
{
	// Refinement parameters
	double nEpisods = 150; // number of refinement iterations
	double minSigma = 1e-1; // minimum exploration noise
	uint nbImportanceSampling = 5; // number of importance sampling (choose 0 if there is no need)
	uint nbS = 10; // number of samples used in stochastic sampling
	uint nVars = 2; // number of refinement variables


	// Set Initial policy and noise covariance
	vec mu = ones<vec>(nVars,1);
	mat sigma = eye<mat>(nVars,nVars);

	GaussianDistribution *GD;
	GD = new GaussianDistribution(mu, sigma);

	mat minNoise, p, pNoisy;
	rowvec r;

	minNoise = zeros(nVars,nVars);
	minNoise.diag() = minSigma * ones(nVars, 1);
	//cout << endl << "minSigma:\n" << minNoise << endl;

	// EM-based reward weighted refinement
	for (uint i=0; i<=nEpisods; i++)
	{
		// Stochatic policy distribution sampling
		pNoisy = GD->stochasticSampling( nbS );
		//cout << endl << "pNoisy: " << endl << pNoisy; // DEBUG

		// Reward evaluation of the policies
		rowvec rNoisy = reward( pNoisy );
		//cout << endl << "rr:\n" << rNoisy << endl; // DEBUG

		// Add the new policy and reward to dataset
		p = join_rows(p, pNoisy);
		r = join_rows(r, rNoisy);

		// after the first iteration nbS should be 1
		nbS = 1;

		// Reward-weighted update of policy distribution
		//GD->setParamsFromData(p);
		//GD->setParamsFromData(p, r);
		GD->setParamsFromData(p, r, nbImportanceSampling);
		GD->setSIGMA( GD->getSIGMA() + minNoise );
	}

	cout << "//////////////////////////////////////////////////" << endl;
	cout << "/////////////// Before Refinement ////////////////" << endl;
	cout << "//////////////////////////////////////////////////" << endl;

	cout << endl << "Initial Mu: \n" << mu << endl;
	cout << endl << "Initial Sigma: \n" << sigma << endl;

	cout << "//////////////////////////////////////////////////" << endl;
	cout << "/////////////// After Refinement /////////////////" << endl;
	cout << "//////////////////////////////////////////////////" << endl;

	cout << endl << "Refined Mu: \n" << GD->getMU() << endl;
	cout << endl << "Refined Sigma: \n" << GD->getSIGMA() << endl;
	cout << endl << "Returns: \n" << r << endl;
	cout << endl << "Explored policies: \n" << p << endl;

	return 0;
}

