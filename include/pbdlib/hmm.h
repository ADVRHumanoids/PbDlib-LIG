/**
Copyright (C) 2015, Martijn Zeestraten 


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
\brief HMM class
The class HMM implements learning for a classical hidden Markov model.

\todo Implement the whole class.

*/

#ifndef HMM_H
#define HMM_H

#include "pbdlib/datapoints.h"
#include "pbdlib/demonstration.h"
#include "pbdlib/gmm.h"
#include "pbdlib/mvn.h"
#include "armadillo"
#include <math.h>

using namespace arma;

namespace pbdlib{
class HMM: public GMM_Model
{
	protected:
		colvec alpha; // Forward variable

		// HMM Extension:
		mat TransitionMatrix;

	private:
		// EM Functions (overrides of GMM parts)
		void learnKMEANS(); 
		uint EM(double likelihood); 
		bool EM_isfinished(double _old, double l_new); 

	public:
		HMM(uint _nSTATES, uint _nVARS);
		HMM(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path,const std::string &transition_path);
		~HMM(){}

		mat& getTRANSITION(){return TransitionMatrix;}

		
		uint EM_Learn(); // Override function
		//double getProbability(const colvec& sample);  // override
		//double getLikelihood(const mat& SAMPLES); // override
		void setTRANSITION(const mat&);
};

} // End pbdlib namespace

#endif

