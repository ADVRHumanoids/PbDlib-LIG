/**
Copyright (C) 2014, Davide De Tommaso, Milad Malekzadeh

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

/*! \file gmr.h
\brief GMM_model class
The class GMM_model allows to use a Gaussian Mixture Model and to learn the parameters from task demonstrations

\author Davide De Tommaso, Milad Malekzadeh
\bug No known bugs.
\todo Make variables/methods to be inherited by HMM-related classes protected members (e.g., nVARS, nSTATES, PRIORS, learnKMEANS(), etc).
\todo EM() may be a virtual function to force re-implementation in child classes such as HMM.
\todo getProbability may be implemented as a template function so that both "mat" and "colvec" can be passed/returned to/by the function.

*/

#ifndef GMM_H
#define GMM_H

#include "pbdlib/datapoints.h"
#include "pbdlib/demonstration.h"
#include "pbdlib/mvn.h"
#include "armadillo"
#include <math.h>

using namespace arma;

namespace pbdlib
{

class GMM_Model
{
	private:

	protected:
    uint nVARS, nSTATES;
    rowvec PRIORS;
    mat gamma;

    std::vector<GaussianDistribution> COMPONENTS;
    std::vector<Demonstration> DEMONSTRATIONS;
    std::vector<std::string> vars_names;

    void learnKMEANS(double regularization = 0.);
    bool EM_isfinished(double l_old, double l_new);
    uint EM(double likelihood,double regularization);

	public:
		GMM_Model(std::vector<Demonstration> &demos, uint _nSTATES);
		GMM_Model(uint _nSTATES, uint _nVARS);
		GMM_Model(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path,const std::string &vars_path);
		GMM_Model(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path);
		~GMM_Model(){}

    void onlineEMDP(int N,colvec P,double lambda, double minSigma);
    uint EM_learn(double regularization = 1E-5);
		std::vector<std::string>& getVARSNames();
		std::string getVARSNames(int);
		uint getIndexOfVARName(const std::string& varname);
		uint getNumVARS();
		uint getNumSTATES();
		rowvec& getPRIORS();
		double getPRIORS(int);
		std::vector<GaussianDistribution>& getCOMPONENTS();
		GaussianDistribution& getCOMPONENTS(int);
		colvec& getMU(int);
		mat& getSIGMA(int);
		mat& getLAMBDA(int);
		mat& getGamma();

		double getProbability(const colvec& sample, bool usePriors = true);
		colvec getProbability(const mat& samples, bool usePriors = true);
		double getLikelihood(const mat& SAMPLES);  
		bool addDemo(Demonstration &demo);

		void setDEMONSTRATIONS(std::vector<Demonstration> demons);
		void setPRIORS(const rowvec& priors);
		void setCOMPONENTS(const std::vector<GaussianDistribution>& components);
		void setMU(int, const colvec&);
		void setSIGMA(int, const mat&);
		void setLAMBDA(int,const mat&);
		void setVARSNames(const std::vector<std::string>& vars);
		void saveInFiles();
		void clear();
	};

} //end of pbdlib namespace

#endif


