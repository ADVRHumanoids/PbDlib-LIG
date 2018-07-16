/**
 
 
Copyright (C) 2015, Martijn Zeestraten, Leonel Rozo, Ioannis Havoutis

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
\brief HSMM class
The class HSMM implements a Hidden semi-Markov model where the model states are represented by a Gaussian mixture model,
and the training phase is carried out by a HMM-based approach.

\author Martijn Zeestraten
\author Leonel Rozo
\author Ioannis Havoutis
*/

#ifndef HSMM_H
#define HSMM_H

#include "pbdlib/datapoints.h"
#include "pbdlib/hmm.h"
#include "pbdlib/demonstration.h"
#include "pbdlib/mvn.h"
#include "armadillo"
#include <math.h>

using namespace arma;

namespace pbdlib 
{

class HSMM:public HMM
{
	protected:
		// -> Variables for computing forward variable recursively (readable):
		colvec 	scalingFtr;	// Scaling factor for recursive computations
		mat 	 	recAlpha; 	// Recursive alpha variable

		// Variables for state duration
		mat   	Pd;     // Matrix with precomputed probability values for the duration;
		uint 		PdSize;	// Maximum maximum duration step

		// -> Functions for recursive forward variable
		void initializeRecursiveForwardVariable(uint T, bool initDur = true);


	private:
		std::vector<GaussianDistribution> DurationCOMPONENTS;

		// -> Variables for calculating Forward Variable efficiently:
		GaussianDistribution* Gtmp;
		mat 		Atmp1, Atmp2; // Help variables for the calculation of alpha variables
		mat    	ALPHA;     		// Variable to keep track of past alphas
		colvec 	bmx;       		// variable used to keep track of observation probability P(x|model)
		colvec 	btmp;      		// Unnormalized observation probabilities
		colvec 	S;         		// Probability of a state starting at t+1 given a partial observation

		// Initialization flags
		bool Initialized,tmpInit;	// Variables used to check if alpha calculation is initialized

		// For predictions:
		colvec 	alphatmp;  // Temp variable to hold in the loop alpha
		mat 		ALPHAtmp;  // to store ALPHA in the prediction loop
		mat    	AlphaPred; // Matrix to hold predictions
		colvec 	Stmp;      // to store S in the prediction loop

		//For demonstration integration
		mat hsmm_transition; //= zeros(init_nStates, init_nStates);
		mat hsmm_transition_ticks; //= zeros(init_nStates, init_nStates);
		rowvec hsmm_priors; //= zeros(init_nStates);
		rowvec hsmm_priors_ticks; // = zeros(init_nStates);

		std::vector <running_stat<double> > hsmm_duration_stats;


//		void initializeFwdCalculation();
		void updateBtmp(colvec& _Btmp,colvec& _obs);
		void updateBtmp(colvec& _Btmp,colvec& _obs,urowvec& _ind);
		void updateBmx(colvec& _bmx, mat& ALPHA, colvec& _btmp);
		void updateALPHA(mat& ALPHA, colvec& S, colvec& _bmx);
		void updateALPHA(mat& ALPHA, colvec& S);
		void updateS(colvec& _S, mat& ALPHA, colvec& _bmx);
		void updateS(colvec& _S, mat& ALPHA);
		void updateAlpha(colvec& _alpha, mat& ALPHA, colvec& btmp);
		void updateAlpha(colvec& _alpha, mat& ALPHA);
	
		// -> Forward variable Functions for efficient computation
		void initializeForwardVariable(); // initialization without observation
		void initializeForwardVariable(colvec& obs); // initialization with observation
		void initializeForwardVariable(colvec& obs, urowvec& ind); // initialization with observation
		void lstepForwardVariable(); // Forward step without observation
		void lstepForwardVariable(colvec& obs); // Forward step with observation
		void lstepForwardVariable(colvec& obs,urowvec& ind); // Forward step with partial observation
		
	public:
		// Constructors
		HSMM(uint _nSTATES, uint _nVARS);		
		HSMM(const std::string &priors_path,
			   const std::string &mu_path, 
			   const std::string &sigma_path,
			   const std::string &transition_path, 
			   const std::string &durMu_path, 
			   const std::string &durSigma_path);
		HSMM(const std::string &priors_path,
				 const std::string &mu_path,
				 const std::string &sigma_path,
				 const std::string &transition_path,
				 const std::string &durMu_path,
				 const std::string &durSigma_path,
				 uint T, bool initDur = true);
		// Destructor
		~HSMM(){}
		
		// -> Forward Variables for efficient computation:
		void stepForwardVariable(); 												// Forward step without observation
		void stepForwardVariable(colvec& obs); 							// Forward step with observation
		void stepForwardVariable(colvec& obs,urowvec& ind); // Forward step with partial observation
		colvec& getForwardVariable(){return alpha;} // Get current forward variable
		void resetForwardVariable(){initializeFwdCalculation();}
		mat& predictForwardVariable(uint N); // Predict forward variable without predicted observations
		void predictForwardVariable(mat& _AlphaPred); // Predict forward variable without predicted observations (implementation for real-time)

		void clear();
		void predictForwardVariableDeterministic(mat& _AlphaPred, int startState);
		void predictForwardVariableStochastic(mat& _AlphaPred);
		void predictForwardVariableStochasticStart(mat& _AlphaPred, int startState);

		void integrateDemonstration(Demonstration demo);
		void get_state_seq(uword state_seq[], mat pred);
		int getClosestState(const colvec P);

		void initializeFwdCalculation();
		// -> Forward variables for classic (readable) computation:
		void stepRecursiveForwardVariable(uint tn);						// Forward variable computation without observation
		void stepRecursiveForwardVariable(uint tn, mat& obs);	// Forward variable computation with full observation
		void stepRecursiveForwardVariable(uint tn, mat& obs, urowvec ind);	// Forward variable computation with partial observation
		mat& getRecursiveForwardVariable(){return recAlpha;};	// Returning matrix containing forward variable values
		void resetRecursiveForwardVariable();									// Set to zero the matrix containing the forward variable values

		// State duration probabilities (set and get)
		void setDurationCOMPONENTS(const std::vector<GaussianDistribution>& components);
		std::vector<GaussianDistribution>& getDurationCOMPONENTS();
		GaussianDistribution& getDurationCOMPONENTS(uint);
		void setDurMU(int, const colvec&);
		colvec& getDurMU(uint);
		void setDurSIGMA(int, const mat&);
		mat& getDurSIGMA(uint);
		void setMaxDuration(uint maxPd);
		uint getMaxDuration(){return PdSize;};

		void saveInFiles(std::string path);
		void loadFromFiles(std::string path);
	};

} // PbdLib namespace

#endif

