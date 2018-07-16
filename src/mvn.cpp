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

#include "pbdlib/mvn.h"

namespace pbdlib
{

GaussianDistribution::GaussianDistribution(colvec& _MU, mat& _SIGMA)
{
	nVARS = _MU.n_rows;
	MU = _MU;
	SIGMA = _SIGMA;
	LAMBDA= _SIGMA.i();
}

GaussianDistribution::GaussianDistribution(uint _nVARS)
{
	nVARS = _nVARS;
}

uint GaussianDistribution::getNumVARS()
{
	return nVARS;
}


colvec &GaussianDistribution::getMU()
{
	return MU;
}


mat& GaussianDistribution::getSIGMA()
{
	return SIGMA;
}

mat& GaussianDistribution::getLAMBDA()
{
	return LAMBDA;
}



void GaussianDistribution::setMU(const colvec &_MU)
{
	MU = _MU;
}


void GaussianDistribution::setSIGMA(const mat &_SIGMA)
{
	SIGMA = _SIGMA;
	LAMBDA = SIGMA.i();
}

void GaussianDistribution::setLAMBDA(const mat &_LAMBDA)
{
	LAMBDA = _LAMBDA;
	SIGMA = _LAMBDA.i();
}


void GaussianDistribution::setNumVARS(uint numvars)
{
	nVARS = numvars;
}


colvec GaussianDistribution::getPDFValue(const mat &SAMPLES)
{
	// Allocate space for solution and temp variable:
	colvec Probs(SAMPLES.n_cols);
//	mat D_Tmp;
//	D_Tmp.set_size(MU.n_rows,SAMPLES.n_cols);

	// Calculate probabilities:
	getPDFValue(Probs,SAMPLES);

	return Probs;
}

colvec GaussianDistribution::getPDFValue(const mat &SAMPLES, urowvec ind)
{
	// Allocate space for solution and temp variable:
	colvec Probs(SAMPLES.n_cols);
//	mat D_Tmp;
//	D_Tmp.set_size(MU.n_rows,SAMPLES.n_cols);

	// Calculate probabilities:
	getPDFValue(Probs,SAMPLES, ind);

	return Probs;
}

void GaussianDistribution::getPDFValue(colvec& _probs, mat _SAMPLES)
{
	// Implementation for real time execution: doesnt require memory allocation:
	// The user should supply:
	// - _probs  : colvec of length nbSamples
	// - _SAMPLES: Matrix with samples (nbVar x nbSamples, also used as tmp variable)

	// calculate difference (x-mu)
	_SAMPLES= trans(_SAMPLES) - repmat(trans(MU),_SAMPLES.n_cols,1);

	// calculate exponential (x-mu)^T*inv(sigma)*(x-mu):
	_probs = sum((_SAMPLES* LAMBDA) % _SAMPLES, 1);

	// calculate Exponential:
	_probs = sqrt(fabs(det(LAMBDA))/pow(2*PI,LAMBDA.n_cols))*exp(-0.5*arma::abs(_probs));

	// Implementation with determinant based on SIGMA instead of Lambda:
	//	Probs = exp(-0.5*arma::abs(Probs)) / sqrt(pow((2*PI),SIGMA.n_cols) * (fabs(det(SIGMA)) + THRESHOLD_MIN));

}

void GaussianDistribution::getPDFValue(colvec& _probs, mat _SAMPLES, urowvec ind)
{
	// Implementation for real time execution: doesnt require memory allocation:
	// The user should supply:
	// - _probs  : colvec of length nbSamples
	// - _SAMPLES: Matrix with samples (nbVar x nbSamples, also used as tmp variable)

	urowvec colInd = linspace<urowvec>(0,_SAMPLES.n_cols-1,_SAMPLES.n_cols);

	// calculate difference (x-mu)
	_SAMPLES= trans(_SAMPLES.submat(ind,colInd)) - repmat(trans(MU.elem(ind)),_SAMPLES.n_cols,1);

	// calculate exponential (x-mu)^T*inv(sigma)*(x-mu):
	mat tmpLAMBDA = SIGMA.submat(ind,ind).i();
	_probs = sum((_SAMPLES* tmpLAMBDA) % _SAMPLES, 1);

	// calculate Exponential:
	_probs = sqrt(fabs(det(tmpLAMBDA))/pow(2*PI,ind.n_elem))*exp(-0.5*arma::abs(_probs));

	// Implementation with determinant based on SIGMA instead of Lambda:
	//	Probs = exp(-0.5*arma::abs(Probs)) / sqrt(pow((2*PI),SIGMA.n_cols) * (fabs(det(SIGMA)) + THRESHOLD_MIN));

}

mat GaussianDistribution::stochasticSampling(uint nbS)
{
	mat noisySamples;

	noisySamples = repmat(MU, 1, nbS) + sqrtm(SIGMA) * randn(nVARS, nbS);

	return noisySamples;
}

mat GaussianDistribution::sqrtm(const mat SIGMA)
{
	mat sqrtmSigma=zeros(nVARS, nVARS);
	vec eigval;
	mat eigvec;
	eig_sym(eigval,eigvec,SIGMA);
	for(int i=0; i<eigval.n_elem;i++)
	{
		sqrtmSigma += eigvec.col(i)*trans(eigvec.col(i)) * sqrt(eigval(i));
	}
	return sqrtmSigma;
}

void GaussianDistribution::setParamsFromData(const mat SAMPLES, const rowvec REWARD, const uint nbImportanceSampling)
{
	mat pTmp, eTmp, rTmpDiag;
	rowvec rTmp, rSrt;

	// Keep the nbImportanceSampling points with highest rewards
	urowvec idSrt = sort_index(REWARD, 1);
	rSrt = sort(REWARD, 1);
	//cout << endl << "rSrt:\n" << rSrt << endl;
	//cout << endl << "indices:\n" << idSrt << endl;
	uint nbP = idSrt.n_elem;

	if ((idSrt.n_elem > nbImportanceSampling) && (!nbImportanceSampling == 0))
		nbP = nbImportanceSampling;
	//cout << endl << "nbImportanceSampling: \n" << nbP << endl;

	pTmp = SAMPLES.cols(idSrt.cols(0,nbP-1));
	//cout << endl << "pTmp: \n" << pTmp << endl;

	rTmp = rSrt.cols(0,nbP-1);
	//cout << endl << "rTmp: \n" << rTmp << endl;

	//Compute error term
	eTmp = pTmp - repmat(MU, 1, pTmp.n_cols);
	//cout << endl << "eTmp: \n" << eTmp << endl;

	//Udpate the current best SAMPLES
	MU = MU + eTmp * trans(rTmp) /sum(rTmp);
	//cout << endl << "Mu: \n" << MU << endl;

	//Update the exploration noise (covariance matrix)
	rTmpDiag = zeros(rTmp.n_cols, rTmp.n_cols);
	rTmpDiag.diag() = rTmp;
	SIGMA = (eTmp * rTmpDiag * trans(eTmp)) / sum(rTmp);

	// Also calculate LAMBDA
	LAMBDA = SIGMA.i();
	//cout << endl << "SIGMA:\n" << SIGMA << endl;
}

void GaussianDistribution::setParamsFromData( const mat SAMPLES, const rowvec REWARD)
{
	// Reward weighted estimation without important sampling
	uint nbImportanceSampling = SAMPLES.n_cols;
	setParamsFromData(SAMPLES, REWARD, nbImportanceSampling);
}

void GaussianDistribution::setParamsFromData( const mat SAMPLES)
{
	// Calculation of MU and SIGMA (without weight and important sampling)
	MU = mean(SAMPLES, 1);
	SIGMA = cov(trans(SAMPLES), 1);
	LAMBDA = SIGMA.i();
	// OR
	//rowvec REWARD = ones(1, SAMPLES.n_cols);
	//uint nbImportanceSampling = SAMPLES.n_cols;
	//setParamsFromData(SAMPLES, REWARD, nbImportanceSampling);
}

} //end of pbdlib namespace

