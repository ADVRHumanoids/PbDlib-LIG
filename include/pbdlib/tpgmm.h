/*
Copyright (C) 2014, Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso, João Silvério, Sylvain Calinon


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

/*! \file pgmm.h
\brief pgmm class
This class is the implementation of the task parameterized Gaussian mixture models (TP-GMM) as described in:

@inproceedings{Calinon13IROS,
	author="Calinon, S. and Alizadeh, T. and Caldwell, D. G.",
	title="On improving the extrapolation capability of task-parameterized movement models",
	booktitle="Proc. {IEEE/RSJ} Intl Conf. on Intelligent Robots and Systems ({IROS})",
	year="2013", month="November", address="Tokyo, Japan",
	pages="610--616",
}

\author Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso, João Silvério, Sylvain Calinon
\bug No known bugs.
*/
#ifndef TPGMM_H
#define TPGMM_H

#include "armadillo"
#include "pbdlib/gmm.h"
#include "pbdlib/gmr.h"
#include "pbdlib/taskparameters.h"
#include "pbdlib/tpdemonstration.h" 

using namespace arma;

namespace pbdlib
{

#define REALMIN 2.2251e-200
#define REALMAX 1.7977e200
enum EMatrixType {INVERTIBLE, SYMMETRIC, ORTHONORMAL, OTHER};

class TPGMM
{
	private:
	/*!
	nVARS: Number of variables
	nSATAES: Number of Gaussian components of the model (This sould be set by the user)
	nFRAMES: Number of frames (This sould be set by the user)
	GMMS: A vector of GMM models containing for each task parameter, the center and covariance of the TPGMM
	PRIORS: A vector containing the mixing coefficients of the Gaussian components
	*/
	uint nVARS, nSTATES, nTaskParameters;
	std::vector<GMM_Model> GMMS; 
	mat Data;
	rowvec PRIORS;
	std::vector<std::string> vars_names;

	// Auxiliary variables for product of Gaussian:
	GMM_Model* prodGMM;
	colvec accMu, tmpMu;
	mat accLambda, tmpLambda;

	std::vector<mat> invA;
	std::vector<mat> invAT;

	void setAuxiliaryVarsProdGauss();
	public:
	
		/*!
		initTensorGMM_timeBased(): Initialization of the centers and the covariances of GMMS using equally spaced time splits
		*/
		void initTensorGMM_timeBased();

		/*!
		initTensorGMM_kmeans(): Initialization of the centers and the covariances of GMMS using kmeans clustering
		*/
		void initTensorGMM_kmeans();

		/*!
		EM_tensorGMM(): The implementation of the Expectation-Maximization algorithm for pgmm
		*/
		void EM_tensorGMM();

		/*!
		estimateTensorGMM(): Calls initTensorGMM_timeBased() & EM_tensorGMM()
		*/
		uint estimateTensorGMM(); //EM algorithm with the last formulation (Tensor)

		TPGMM(std::vector<TPDemonstration>& demos, uint nSTATES);
		TPGMM(uint nVARS, uint nSTATES, uint _nTPs);
		~TPGMM(){}

		/*!
		loadPGMMfromMATLAB(): From a set of given .txt files reads the PRIORS, MU and SIGMA of the GMMS
		Currently, the MU and SIGMA for every task parameter are saved in different files in the following format:
		'MuFileName_P' nbFrame '.txt'
		'SigmaFileName_P' nbFrame '.txt'
		Note that in the MuFileName*.txt there is a matrix containing centers like this: [MU_1 MU_2 MU_3 ... MU_j]; j=1:nSTATES.
		Similarly, in the SigmaFileName*.txt there is matrix containing the covariance matrices like this: [SIGMA_1 SIGMA_2 ... SIGMA_j] ; j=1:nSTATES
		*/
		void loadTPGMMfromMATLAB(std::string PriorsFileName, std::string VarsFileName, std::string MuFileName, std::string SigmaFileName);
	
		/*!
		getNumVARS(): Gives nVARS
		*/
		uint getNumVARS();
		/*!
		getNumSTATES(): Gives nSTATES
		*/
		uint getNumSTATES();
		/*!
		getNumFRAMES(): Gives nFRAMES
		*/
		uint getNumFRAMES();
		/*!
		getPRIORS(): Gives the PRIORS
		*/
		
		rowvec& getPRIORS();
		double getPRIORS(uint);

		/*!
		getGMMS(): Gives GMMS
		*/
		std::vector<GMM_Model>& getGMMS();
		GMM_Model& getGMMS(uint id);

		GMM_Model* getTransformedGMM(std::vector<TaskParameter> _TPs, EMatrixType AType = OTHER );// [Calinon, 2012] Eq. (3) - Prod. of linearly transformed Gaussians
		GMM_Model* getTransformedGMM(TaskParameters _TPs, EMatrixType AType = OTHER );// [Calinon, 2012] Eq. (3) - Prod. of linearly transformed Gaussians

		GaussianDistribution TPGMR(std::vector<TaskParameter> _TPs, Datapoints* DataIn, double diagRegularizationFactor = 1E-4); // Computes GMR within the frames and multiplies the resulting Gaussians (TODO: Implement options with AType)
		GaussianDistribution TPGMR(TaskParameters _TPs, Datapoints* DataIn, double diagRegularizationFactor = 1E-4);

		void setVARSNames(const std::vector<std::string>& vars);
		std::vector<std::string>& getVARSNames();
};

} // end of namespace pbdlib

#endif

