/**


Copyright (C) 2015, Leonel Rozo

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

/*! \file
\brief ADHSMM class
The class ADHSMM implements an adaptive duration Hidden semi-Markov model where both the model states and duration probabilities are represented
by a Gaussian mixture model. The training phase is carried out by a HMM-based approach.

\author Leonel Rozo
\bugs		No bugs known
\todo		Handling duration GMM with different number of components for every ADHSMM state
\todo 	Determine if a composition pattern may prevent inheritance for some methods of HSMM (e.g., getDurationCOMPONENTS())
*/

#ifndef ADHSMM_H
#define ADHSMM_H

#include "pbdlib/datapoints.h"
#include "pbdlib/hsmm.h"
#include "pbdlib/demonstration.h"
#include "pbdlib/mvn.h"
#include "pbdlib/gmm.h"
#include "pbdlib/gmr.h"
#include "armadillo"
#include <math.h>

using namespace arma;

namespace pbdlib
{

class ADHSMM : public HSMM
{
	public:
		// Constructors
		ADHSMM(const std::string &priors_path,
					 const std::string &mu_path,
					 const std::string &sigma_path,
					 const std::string &transition_path,
					 const std::string &durPriors_path,
					 const std::string &durMu_path,
					 const std::string &durSigma_path,
					 uint T);

		// Destructor
		~ADHSMM(){}

		// Set and get functions for GMM encoding duration probabilities
		void setNumSTATESDUR(uint _nSTATESDUR);
		uint getNumSTATESDUR(){return nSTATESDUR;};
		void setNumVARSDUR(uint _nVARSDUR);
		uint getNumVARSDUR(){return nVARSDUR;};
		void setDurationGMMs(std::vector<GMM_Model> durGMMs);
		std::vector<GMM_Model>& getDurationGMMs();
		GMM_Model& getDurationGMMs(uint idGMM);

		// -> Forward variables for recursive (readable) computation:
		void stepRecursiveForwardVariable(colvec _u, urowvec _varIn, urowvec _varOut, uint _tn);	// Forward variable computation without observation
		void stepRecursiveForwardVariable(colvec _u, urowvec _varIn, urowvec _varOut, uint _tn, mat& _obs);	// Forward variable computation with full observation
		void stepRecursiveForwardVariable(colvec _u, urowvec _varIn, urowvec _varOut, uint _tn, mat& _obs, urowvec _ind);	// Forward variable computation with partial observation

	private:
			std::vector<GMM_Model> durationGMMs;// Vector of GMM for adaptive duration state probabilities
			GMR 		*condDurationProb;					// GMR for computing conditional duration probabilities
			GMM_Model *condPd;									// Resulting conditional duration probability obtained from GMR
			uint 			nVARSDUR;									// Number of variables for duration probabilities
			uint 			nSTATESDUR;								// Number of states for duration probabilities. It is assumed that every GMM has the same number.

			void computeConditionedDurationProbs(uint _i, colvec u, urowvec _varIn, urowvec _varOut);	// Computes the conditional duration probabilities for a given external input
			void initializeRecursiveForwardVariable(uint T);																					// Initializes the forward variable
			void initializeCondDurProb(mat durGMMpriors, mat durGMMmu, mat durGMMsigma);							// Initialize duration state probabilities

};
} // PbdLib namespace

#endif
