/**
 *
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

#include "pbdlib/hmm.h"

namespace pbdlib
{
HMM::HMM(uint _nSTATES, uint _nVARS):
		GMM_Model(_nSTATES,_nVARS)
{
}

HMM::HMM(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path, const std::string &transition_path)
	:GMM_Model(priors_path,mu_path,sigma_path)
{
//	mat priors, mu, sigma; // GMM components
	mat transition; // HMM addition

	// Load from File
//	priors.load(priors_path, raw_ascii);
//	mu.load(mu_path, raw_ascii);
//	sigma.load(sigma_path, raw_ascii);
	transition.load(transition_path,raw_ascii);

	// Determine dimensions:
	/*
	nVARS = mu.n_rows;
	nSTATES = priors.n_elem;


	for (uint i = 0;i<nSTATES;i++)
	{
		this->setMU(i,mu.col(i));
		this->setSIGMA(i,sigma.cols(i*nVARS,(i+1)*nVARS-1));
	}
	this->setPRIORS(priors);
	*/

	// Load transition components:
	this->setTRANSITION(transition);
}

void HMM::setTRANSITION(const mat& _Transition)
{
	TransitionMatrix = _Transition;
}

/*
double HMM::getProbability(const colvec& sample)
{
	double P = 0.0;
	// See Eq. (2.0.2) in doc/TechnicalReport.pdf
	for(uint k=0; k<nSTATES; k++)
		//P += alpha[k] * as_scalar( this->gmm->getCOMPONENTS(k).getPDFValue(sample) ); // See Eq. (2.0.3) for "getPDFValue" in doc/TechnicalReport.pdf
	return P;
}



void HMM::learnKMEANS()
{
	mat DemosTmp = DEMONSTRATIONS[0].getDatapoints().getData();

	for(int i=1; i<DEMONSTRATIONS.size(); i++)
		DemosTmp.insert_cols(DEMONSTRATIONS[0].getDatapoints().getNumPOINTS(), DEMONSTRATIONS[i].getDatapoints().getData());

	//Criterion to stop the EM iterative update
	double cumdist_threshold = 1e-10;
	uint maxIter = 100;

	//Initialization of the parameters
	uint nbVar = DemosTmp.n_rows;
	uint nbData = DemosTmp.n_cols;
	double cumdist_old = -std::numeric_limits<double>::max();
	uint nbStep = 0;

	//srand (time(NULL));

	urowvec idTmp = sort_index(randn<urowvec>(nbData));

	uvec allrows = linspace<uvec>(0, nbVar-1, nbVar);

	mat Mu = DemosTmp.submat(allrows, idTmp.cols(0,nSTATES-1));
	cube Sigma = zeros(nbVar,nbVar,nSTATES);


	//k-means iterations
	while(true){
		mat distTmp = zeros(nSTATES,nbData);

		//E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for (uint i=0; i<nSTATES; i++){
			//Compute distances
			distTmp.row(i) = sum((DemosTmp - repmat(Mu.col(i), 1, nbData)) % (DemosTmp - repmat(Mu.col(i), 1, nbData)));
		}

		vec vTmp = zeros<vec>(nbData,1);
		vec idList = zeros<vec>(nbData,1);
		uword index;
		vec tempVec;
		double tempId;
		for (uint i=0; i<nbData; i++){
			tempVec = distTmp.col(i);
			vTmp.row(i) = min(tempVec);
			tempId = tempVec.min(index);
			idList.row(i) = index;
		}

		double cumdist = sum(vTmp);

		//M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		uvec idTmp1;
		rowvec priors = zeros(1,nSTATES);
		mat sigTmp = zeros(nbVar,nbVar);
		for (uint i=0; i<nSTATES; i++){
			idTmp1 = find(idList == i);
			priors(i) = idTmp1.n_elem;
			Mu.col(i) = mean(DemosTmp.submat(allrows, idTmp1), 1);
			sigTmp = cov( trans(join_rows(DemosTmp.submat(allrows, idTmp1),DemosTmp.submat(allrows, idTmp1))) ,0 );
			Sigma.slice(i) = sigTmp;
		}
		priors = priors/nbData;

		//Stopping criterion %%%%%%%%%%%%%%%%%%%%
		if (fabs(cumdist-cumdist_old) < cumdist_threshold){
			cout << "%%%%%%%%%%% KMEANS %%%%%%%%%%%%%%%%%%" << endl;
			cout << endl << "Priors_kmeans :" << endl << priors << endl;
			cout << endl << "Mu_kmeans :" << endl << Mu << endl;
			cout << endl << "Sigma_kmeans :" << endl << Sigma << endl;
			cout << "%%%%%%% KMEANS FINISHED %%%%%%%%%%%%%" << endl;

			GaussianDistribution *GD;
			colvec mu;

			for(int i=0; i<nSTATES; i++){
				PRIORS(i) = priors(i);
				mu = Mu.col(i);
				COMPONENTS[i].setMU(mu);
				COMPONENTS[i].setSIGMA(Sigma.slice(i));
			}

			return;
		}
		cumdist_old = cumdist;
		nbStep = nbStep + 1;

		//cout << endl << "nbStep   " << nbStep << endl;
	}
}


double HMM::getLikelihood(const mat &SAMPLES)
{
	double L=0.0;
	colvec s(SAMPLES.n_rows);

	// See Eq. (2.0.4) in doc/TechnicalReport.pdf
	for(uint j=0; j<SAMPLES.n_cols; j++){
		s = SAMPLES.col(j);
		L += log( getProbability( s ) );
	}
	return L;
}
*/




} // End pbdlib namespace
