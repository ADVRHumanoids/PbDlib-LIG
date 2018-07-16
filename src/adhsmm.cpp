/**
Copyright (C) 2015,	Leonel Rozo

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

#include "pbdlib/adhsmm.h"

namespace pbdlib
{
ADHSMM::ADHSMM(const std::string &priors_path,
		const std::string &mu_path,
		const std::string &sigma_path,
		const std::string &transition_path,
		const std::string &durPriors_path,
		const std::string &durMu_path,
		const std::string &durSigma_path,
		uint T):HSMM(priors_path,mu_path,sigma_path,transition_path,durMu_path,durSigma_path, T, false)
{
	// Setting duration state probabilities
	mat tmpDurPriors, tmpDurMu, tmpDurSigma;
	// Loading duration probability models
	tmpDurPriors.load(durPriors_path,raw_ascii);
	tmpDurMu.load(durMu_path,raw_ascii);
	tmpDurSigma.load(durSigma_path, raw_ascii);
	this->nSTATESDUR = tmpDurPriors.n_rows;
	this->nVARSDUR = tmpDurMu.n_rows;

	// Setting duration GMMs
	initializeCondDurProb(tmpDurPriors, tmpDurMu, tmpDurSigma);
	this->condDurationProb = new GMR();
}

void ADHSMM::initializeRecursiveForwardVariable(uint T)
{

}

void ADHSMM::initializeCondDurProb(mat tmpDurPriors, mat durGMMmu, mat durGMMsigma)
{
	GMM_Model _tmpDurationProb(this->nSTATESDUR, this->nVARSDUR);
	std::vector<GaussianDistribution> _tmpComponents;
	uint _indMu, _ind1Sigma, _ind2Sigma;
	colvec _tmpMu;
	mat _tmpSigma;

	// Assuming that the means and covariances matrices are ordered like this:
	// durGMMmu = [ mu_1,1 mu_1,2 ... mu_1,nSTATESDUR 	mu_2,1 mu_2,2, ... mu_2,nSTATESDUR 	... mu_nSTATES,1 mu_nSTATE,2 ... mu_nSTATES,nSTATESDUR ]
	// durGMMsigma = [ Sigma_1,1 Sigma_1,2 ... Sigma_1,nSTATESDUR   Sigma_2,1 Sigma_2,2 ... Sigma_2,nSTATESDUR   ... Sigma_nSTATES,1 Sigma_nSTATES,2 ... Sigma_nSTATES,nSTATESDUR]
	for(uint i=0 ; i<this->nSTATES ; i++)
	{
		for(uint j=0 ; j<this->nSTATESDUR ; j++)
		{
			// Temporally storing \mu
			_indMu = j+i*this->nSTATESDUR;
			_tmpMu = durGMMmu.col(_indMu);

			// Temporally storing \Sigma
			_ind1Sigma = (j*this->nVARSDUR)+i*this->nVARSDUR*this->nSTATESDUR;
			_ind2Sigma = (j*this->nVARSDUR+this->nVARSDUR-1)+i*this->nVARSDUR*this->nSTATESDUR;
			_tmpSigma = durGMMsigma.cols(_ind1Sigma, _ind2Sigma);

			// Saving current Gaussian distribution
			_tmpComponents.push_back(GaussianDistribution(_tmpMu, _tmpSigma));
		}
		// Saving GMM encoding duration probability for state "i"
		_tmpDurationProb.setCOMPONENTS(_tmpComponents);
		_tmpDurationProb.setPRIORS(tmpDurPriors.col(i).t());

		// Saving in vector of GMMs
		this->durationGMMs.push_back(_tmpDurationProb);

		// Clearing temporary saving
		_tmpComponents.clear();
	}
}

void ADHSMM::computeConditionedDurationProbs(uint _i, colvec _u, urowvec _varIn, urowvec _varOut)
{
	// Set Pd and PdSize for current external input value
	this->condDurationProb->setGMMModel(&this->durationGMMs[_i]);
	// Evaluating GMR for the given external input _u
	this->condPd = this->condDurationProb->regression(_u, _varIn, _varOut);
	// Setting maximum duration based on the resulting \mu and \Sigma
	this->PdSize = round(this->condPd->getMU(0)(0) + 2*this->condPd->getSIGMA(0)(0,0));
}

void ADHSMM::setNumSTATESDUR(uint _nSTATESDUR)
{
	this->nSTATESDUR = _nSTATESDUR;
}

void ADHSMM::setNumVARSDUR(uint _nVARSDUR)
{
	this->nVARSDUR = _nVARSDUR;
}

void ADHSMM::setDurationGMMs(std::vector<GMM_Model> _durGMMs)
{
	if(_durGMMs.size() == this->nSTATES)
			this->durationGMMs = _durGMMs;
	else
		std::cout << "\n [ERROR]::ADHSMM_Model::setDurationGMMs if(_durGMMs.size() == this->nSTATES) ... else .";
}

std::vector<GMM_Model>& ADHSMM::getDurationGMMs()
{
	return this->durationGMMs;
}

GMM_Model& ADHSMM::getDurationGMMs(uint _idGMM)
{
	return this->durationGMMs[_idGMM];
}

void ADHSMM::stepRecursiveForwardVariable(colvec _u, urowvec _varIn, urowvec _varOut, uint _tn)
{
	colvec colTn(1), colD(1);
	colTn(0) = _tn+1;	// This needs to be done for evaluating conditional duration probability. There is no "zero duration"

	for(uint i = 0 ; i < this->getNumSTATES() ; i++)
	{
		// Setting PdSize and Pd based on external input (GMR for current GMM encoding duration probability)
		this->computeConditionedDurationProbs(i, _u, _varIn, _varOut);

		// \alpha_i,t = \Pi_i * Pd_i(t|u) + ...
		if(_tn < PdSize)
			recAlpha(i,_tn) = this->getPRIORS()(i) * this->condPd->getCOMPONENTS(0).getPDFValue(colTn)(0);

		// \alpha_i,t = \sum_d \alpha_j,t-d * a_j,i * Pd_i(d|u)
		for(uint d = 0 ; d < std::min(_tn,PdSize) ; d++)
		{
			colD(0) = d+1;	// This needs to be done for evaluating conditional duration probability. There is no "zero duration"
			recAlpha(i,_tn) = recAlpha(i,_tn) + as_scalar(recAlpha.col(_tn-d-1).t() * this->getTRANSITION().col(i)) *
					this->condPd->getCOMPONENTS(0).getPDFValue(colD)(0);
		}
	}
}

void ADHSMM::stepRecursiveForwardVariable(colvec _u, urowvec _varIn, urowvec _varOut, uint _tn, mat& _obs)
{
	double tmpObsProb;
	colvec colTn(1), colD(1);
	colTn(0) = _tn+1;	// This needs to be done for evaluating conditional duration probability. There is no "zero duration"

	for(uint i = 0 ; i < this->getNumSTATES() ; i++)
	{
		// Setting PdSize and Pd (precomputed duration prob.)
		this->computeConditionedDurationProbs(i, _u, _varIn, _varOut);

		// \alpha_i,t = \Pi_i * Pd_i(t|u) * \prod_s N(\ksi_s) + ...
		if(_tn < PdSize)
		{
		 	tmpObsProb = arma::prod( scalingFtr.subvec(0, _tn) % this->getCOMPONENTS(i).getPDFValue(_obs) );
			recAlpha(i,_tn) = this->getPRIORS()(i) * this->condPd->getCOMPONENTS(0).getPDFValue(colTn)(0) * tmpObsProb;
		}

		// \alpha_i,t = \sum_d \alpha_j,t-d * a_j,i * Pd_i(d|u) * \prod_s N(\ksi_s)
		for(uint d = 0 ; d < std::min(_tn,PdSize) ; d++)
		{
			tmpObsProb = arma::prod( scalingFtr.subvec(_tn-d, _tn) % this->getCOMPONENTS(i).getPDFValue(_obs.cols(_tn-d, _tn)) );
			colD(0) = d+1;	// This needs to be done for evaluating conditional duration probability. There is no "zero duration"
			recAlpha(i,_tn) = recAlpha(i,_tn) + as_scalar(recAlpha.col(_tn-d-1).t() * this->getTRANSITION().col(i)) *
					this->condPd->getCOMPONENTS(0).getPDFValue(colD)(0) * tmpObsProb;
		}
	}

	// Computing scaling factor
	if(_tn < scalingFtr.n_rows-1)
		scalingFtr(_tn+1) = 1/sum(recAlpha.col(_tn));
}

void ADHSMM::stepRecursiveForwardVariable(colvec _u, urowvec _varIn, urowvec _varOut, uint _tn, mat& _obs, urowvec _ind)
{
	double tmpObsProb;
	colvec colTn(1), colD(1);
	colTn(0) = _tn+1;	// This needs to be done for evaluating conditional duration probability. There is no "zero duration"

	for(uint i = 0 ; i < this->getNumSTATES() ; i++)
	{
		// Setting PdSize and Pd (precomputed duration prob.)
		this->computeConditionedDurationProbs(i, _u, _varIn, _varOut);

		// \alpha_i,t = \Pi_i * Pd_i(t|u) * \prod_s N(\ksi_s) + ...
		if(_tn < PdSize)
		{
		 	tmpObsProb = arma::prod( scalingFtr.subvec(0, _tn) % this->getCOMPONENTS(i).getPDFValue(_obs, _ind) );
			recAlpha(i,_tn) = this->getPRIORS()(i) * this->condPd->getCOMPONENTS(0).getPDFValue(colTn)(0) * tmpObsProb;
		}

		// \alpha_i,t = \sum_d \alpha_j,t-d * a_j,i * Pd_i(d|u) * \prod_s N(\ksi_s)
		for(uint d = 0 ; d < std::min(_tn,PdSize) ; d++)
		{
			tmpObsProb = arma::prod( scalingFtr.subvec(_tn-d, _tn) % this->getCOMPONENTS(i).getPDFValue(_obs.cols(_tn-d, _tn), _ind) );
			colD(0) = d+1;	// This needs to be done for evaluating conditional duration probability. There is no "zero duration"
			recAlpha(i,_tn) = recAlpha(i,_tn) + as_scalar(recAlpha.col(_tn-d-1).t() * this->getTRANSITION().col(i)) *
					this->condPd->getCOMPONENTS(0).getPDFValue(colD)(0) * tmpObsProb;
		}
	}

	// Computing scaling factor
	if(_tn < scalingFtr.n_rows-1)
		scalingFtr(_tn+1) = 1/sum(recAlpha.col(_tn));
}

} // End pbdlib namespace
