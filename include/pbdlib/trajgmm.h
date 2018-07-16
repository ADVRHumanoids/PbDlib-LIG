/**
Copyright (C) 2015, Leonel Rozo, Davide De Tommaso, Milad Malekzadeh

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

/*! \file trajgmm.h
\brief TrajGMM class
The class TrajGMM allows to use a trajectory GMM and to learn the parameters from task demonstrations

\author Leonel Rozo
\bug No known bugs.
\todo Function init_TrajGMM_timeBased() should be moved to GMM_Model class
*/

#ifndef TRAJGMM_H
#define TRAJGMM_H

#include "pbdlib/gmm.h"
#include "armadillo"
#include <math.h>

using namespace arma;

namespace pbdlib
{

class TrajGMM : public GMM_Model
{
	private:
		uint 		nVARPOS;	// Dimension of position data
		uint 		nDERIV;		// Number of static and dynamic features
		double	dt;				// Time step to be used in matrix PHI

		sp_mat PHI;				// Window matrix for all the demonstrations
		sp_mat PHI1;			// Window matrix for only one demonstrations

		void init_TrajGMM_timeBased(double regularization = 1E-5);	// Initialization with equal bins splitting
		template<class armaType>
			armaType circshift(const armaType &_matrix, int _shX, int _shY);	// Circular shift for matrices/vectors
		colvec getMuQ(vec _q);		// Function to get the super vector Mu containing the all the means given a sequence of states
		mat 	 getSigmaQ(vec _q);	// Function to get the super (sparse) matrix Sigma containing the all the covariances given a sequence of states

	public:
		// Constructors
		TrajGMM(std::vector<Demonstration> &demos, uint _nSTATES, uint _nVARPOS, uint _nDERIV, double _dt);
		TrajGMM(uint _nSTATES, uint _nVARPOS, uint _nDERIV, double _dt);
		TrajGMM(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path,
				const std::string &vars_path, const std::string &modelInfo_path);
		TrajGMM(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path, const std::string &modelInfo_path);
		// Destructor
		~TrajGMM(){}

		uint EM_learn(double regularization = 1E-5);		// Override EM learning that uses equal bins splitting
    void constructPHI(uint _nData, uint _nSamples);	// Function to construct the window matrix PHI
    mat  leastSquaresOptimalData(colvec _q);				// Function to carry out the WLS-like optimization for a given sequence of states
		std::vector<GaussianDistribution> leastSquaresOptimalProb(colvec _q);	// Function to carry out the WLS-like optimization for a given sequence of states

		uint 	 getNumVARSPOS(){return nVARPOS;};// Returns dimension of position data
		uint 	 getNumDERIV(){return nDERIV;};		// Returns number of static and dynamic features
		double getDt(){return dt;};							// Returns the time step used in the window matrix PHI
		sp_mat getPHI(){return PHI;};						// Returns window matrix for all the demonstrations
		sp_mat getPHI1(){return PHI1;};					// Returns window matrix for only one demonstration
	};

} //end of pbdlib namespace

#endif


