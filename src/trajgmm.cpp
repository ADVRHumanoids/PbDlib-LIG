/**
Copyright (C) 2015, Leonel Rozo, Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon



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

#include "pbdlib/trajgmm.h"

namespace pbdlib
{
TrajGMM::TrajGMM(std::vector<Demonstration> &demos, uint _nSTATES,
		uint _nVARPOS, uint _nDERIV, double _dt):GMM_Model(demos, _nSTATES)
{
	// Initializing Trajectory GMM variables
	this->nVARPOS = _nVARPOS;
	this->nDERIV = _nDERIV;
	this->dt = _dt;
}

TrajGMM::TrajGMM(uint _nSTATES, uint _nVARPOS, uint _nDERIV, double _dt)
	:GMM_Model(_nSTATES, _nVARPOS*_nDERIV)
{
	// Initializing Trajectory GMM variables
	this->nVARPOS = _nVARPOS;
	this->nDERIV = _nDERIV;
	this->dt = _dt;
}

TrajGMM::TrajGMM(const std::string &priors_path,
		const std::string &mu_path, const std::string &sigma_path,
		const std::string &vars_path, const std::string &modelInfo_path)
		:GMM_Model(priors_path, mu_path, sigma_path, vars_path)
{
	// Setting Trajectory GMM variables
	mat modelInfo;
	modelInfo.load(modelInfo_path);
	this->nVARPOS = modelInfo(0);
	this->nDERIV = modelInfo(1);
	this->dt = modelInfo(2);
}

TrajGMM::TrajGMM(const std::string &priors_path,
		const std::string &mu_path, const std::string &sigma_path, const std::string &modelInfo_path)
		:GMM_Model(priors_path, mu_path, sigma_path)
{
	// Setting Trajectory GMM variables
	mat modelInfo;
	modelInfo.load(modelInfo_path);
	this->nVARPOS = modelInfo(0);
	this->nDERIV = modelInfo(1);
	this->dt = modelInfo(2);
}

void TrajGMM::init_TrajGMM_timeBased(double _regTerm)
{
	// Auxiliar variables
	mat 		_tmpData;
	rowvec 	_tmpPriors;
	uvec 		_idTmp;

	// Storing 1st demonstration in matrix format
	mat _Data = DEMONSTRATIONS[0].getDatapoints().getData();
	cout << "Size(data) = " << _Data.n_rows << "," << _Data.n_cols << endl;
	// Time dimension
	rowvec _timeDim = linspace<rowvec>(1, _Data.n_cols, _Data.n_cols);
	cout << "Size(_timeDim) = " << _timeDim.n_rows << "," << _timeDim.n_cols << endl;
	// Adding time dimension to first demo
	_Data.insert_rows(0,_timeDim);
	cout << "Size(data) = " << _Data.n_rows << "," << _Data.n_cols << endl;

	// Storing remaining demos
	for(int i=1; i<DEMONSTRATIONS.size(); i++)
	{
		_tmpData = DEMONSTRATIONS[i].getDatapoints().getData();	// Temporary storing
		_tmpData.insert_rows(0,_timeDim);		// Adding time dimension
		_Data = join_rows(_Data,_tmpData);	// Adding demo
		cout << "Size(data) = " << _Data.n_rows << "," << _Data.n_cols << endl;
	}

	// Linear Space Gaussian in time:
	mat _TimingSep = linspace(min(_Data.row(0)), max(_Data.row(0)),this->nSTATES+1).t();

	// Setting Mu, Sigma and Priors
	_tmpPriors = zeros(1,nSTATES);
	for (uint i=0; i<nSTATES; i++)
	{
		// Find All indices of data that
		_idTmp = find(_Data.row(0) >= _TimingSep(i));
		_tmpData = _Data.cols(_idTmp);
		_idTmp = find(_tmpData.row(0) <= _TimingSep(i+1));

		// Calculate Properties of GMM;
		this->setMU(i,mean(_tmpData.cols(_idTmp),1));
		this->setSIGMA(i, cov(_tmpData.cols(_idTmp).t()) + eye(this->nVARS, this->nVARS)*_regTerm);
		_tmpPriors(i) = _idTmp.n_elem;
	}
	this->setPRIORS(_tmpPriors/sum(_tmpPriors));
}

uint TrajGMM::EM_learn(double _regTerm)
{
	// Setting new (temporary) values for nVARS
	this->nVARS += 1;
	init_TrajGMM_timeBased(_regTerm);

	// Setting back nVARS value
	this->nVARS -= 1;

	// Removing time dimension from model parameters
	for(uint i = 0 ; i < this->nSTATES ; i++)
	{
		this->setMU(i, this->getMU(i).rows(1, this->nVARS));
		this->setSIGMA(i, this->getSIGMA(i).submat(1,1,this->nVARS, this->nVARS));
	}

	// EM learning
	return EM(std::numeric_limits<double>::min(),_regTerm);
}

void TrajGMM::constructPHI(uint _nData, uint _nSamples)
{
	// Initializing and setting a one dimensional operator for one datapoint
	mat oneDimOperator = zeros<mat>(this->nDERIV,this->nDERIV);
	oneDimOperator(0, oneDimOperator.n_cols-1) = 1;

	for(uint i = 1 ; i < this->nDERIV ; i++)
		oneDimOperator.row(i) = (oneDimOperator.row(i-1) - circshift(arma::conv_to<mat>::from(oneDimOperator.row(i-1)), 0, -1)) / this->dt;

	// Initializing and setting a one dimensional operator for all the datapoints of a dataset (e.g., demonstrations)
	sp_mat bigOperator = zeros<sp_mat>(_nData*this->nDERIV, _nData);
	sp_mat PHI0 = zeros<sp_mat>(_nData*this->nDERIV, _nData);
	bigOperator.submat((this->nDERIV-1)*this->nDERIV, 0, this->nDERIV*this->nDERIV-1, this->nDERIV-1) = oneDimOperator;

	// Setting matrix PHI for one dataset
	for(int t = 0 ; t <= _nData-this->nDERIV ; t++)
		PHI0 = PHI0 + circshift(bigOperator, this->nDERIV*t, t);

	// Handling of borders
	for(int i = 1 ; i < this->nDERIV ; i++)
	{
		bigOperator.row(this->nDERIV*this->nDERIV-i).zeros();
		bigOperator.col(i-1).zeros();
		PHI0 = PHI0 + circshift(bigOperator, -i*this->nDERIV, -i);
	}

	// Application to multiple dimensions and multiple datasets (e.g., demonstrations)
	this->PHI1 = conv_to<sp_mat>::from( kron(conv_to<mat>::from(PHI0), arma::eye(this->nVARPOS, this->nVARPOS)) );
	this->PHI  = conv_to<sp_mat>::from( kron(arma::eye(_nSamples, _nSamples), conv_to<mat>::from(this->PHI1)) );
}

template<class armaType>
armaType TrajGMM::circshift(const armaType& _matrix, int _shX, int _shY)
{
	int iOutputInd, iInputInd, ii, jj;
	int ydim = _matrix.n_cols;
	int xdim = _matrix.n_rows;
	armaType outMat = _matrix;

	for (int j = 0 ; j < ydim ; j++)
	{
		// Setting index for columns shift (vertical shifting in the input matrix)
		jj = (j + _shY) % ydim;	// Here operator % is computing the modulo

		if (jj <0)
			jj = jj + ydim;

	  for(int i = 0; i < xdim ; i++)
	  {
	  	// Setting index for rows shift (horizontal shifting in the input matrix)
	  	ii = (i + _shX) % xdim;	// Here operator % is computing the modulo

	  	if (ii <0)
	  		ii = ii + xdim;

	  	// Shifting
	  	outMat[(jj * xdim) + ii] = _matrix[(j * xdim) + i];
	  }
	}

	return outMat;
}

colvec TrajGMM::getMuQ(vec _q)
{
	colvec muQ = zeros<colvec>(this->nVARS*_q.n_elem);

	for(uint i = 0 ; i < _q.n_elem ; i++)
		muQ.rows(i*this->nVARS, (i+1)*this->nVARS-1) = this->getMU(_q(i));

	return muQ;
}

mat TrajGMM::getSigmaQ(vec _q)
{
	mat sigmaQ = zeros<mat>(this->nVARS*_q.n_elem, this->nVARS*_q.n_elem);
	uint ind1, ind2;

	for(uint i = 0 ; i < _q.n_elem ; i++)
	{
		ind1 = i*this->nVARS;
		ind2 = (i+1)*this->nVARS-1;
		sigmaQ.submat(ind1, ind1, ind2, ind2) = this->getSIGMA(_q(i));
	}
	return sigmaQ;
}

mat TrajGMM::leastSquaresOptimalData(colvec _q)
{
	mat optimalData;

	// Get big vector Mu and big matrix Sigma
	colvec muQ = getMuQ(_q);
	mat sigmaQ = getSigmaQ(_q);

	// Least squares computation using a readable code (optimized using sp_mat type and spsolve() function for solving the WLS problem)
	if(this->PHI1.is_empty())
		this->constructPHI(_q.n_elem,1);

	sp_mat PHI1invSigmaQ = this->PHI1.t() * conv_to<sp_mat>::from(sigmaQ.i());	// PHI1^\trsp * Sigma_Q
	sp_mat Rq = PHI1invSigmaQ * this->PHI1;	// (PHI1^\trsp * Sigma_Q)^-1 * PHI1
	mat rq = PHI1invSigmaQ * muQ;	// (PHI1^\trsp * Sigma_Q)^-1 * Mu_Q
	mat zeta = spsolve(Rq,rq,"lapack"); // Solving least squares problem
	optimalData = reshape(zeta, this->nVARPOS,_q.n_elem);

	// Least squares computation using Cholesky and QR decompositions
	// -> THE CODE BELOW WORKS BUT IT IS NOT OPTIMIZED FOR HANDLING SPARSE MATRICES AND THEREFORE IT IS MUCH SLOWER THAN THE PIECE
	// 		OF CODE ABOVE.
	/*
	mat T = chol(sigmaQ).t();
	mat Tphi = T.i() * this->PHI1;
	mat Tmu = T.i() * muQ;
	mat Q, R;
	qr_econ(Q, R, Tphi);
	mat z = Q.t() * Tmu;
	mat zeta1 = R.i() * z;
	optimalData = reshape(zeta1, this->nVARPOS, _q.n_elem);
	*/

	return optimalData;
}

std::vector<GaussianDistribution> TrajGMM::leastSquaresOptimalProb(colvec _q)
{
	std::vector<GaussianDistribution> outputProb;

	// Computing optimal path
	mat optimalData;

	// Get big vector Mu and big matrix Sigma
	colvec muQ = getMuQ(_q);
	mat sigmaQ = getSigmaQ(_q);

	// Least squares computation using a readable code (optimized using sp_mat type and spsolve() function for solving the WLS problem)
	if(this->PHI1.is_empty())
		this->constructPHI(_q.n_elem,1);

	sp_mat PHI1invSigmaQ = this->PHI1.t() * conv_to<sp_mat>::from(sigmaQ.i());
	sp_mat Rq = PHI1invSigmaQ * this->PHI1;
	mat rq = PHI1invSigmaQ * muQ;
	mat zeta = spsolve(Rq,rq,"lapack");
	optimalData = reshape(zeta, this->nVARPOS,_q.n_elem);

	// Covariance Matrix computation of ordinary least squares estimate
	mat mse = (muQ.t()*sigmaQ.i()*muQ - rq.t()*inv(conv_to<mat>::from(Rq))*rq) / ((this->nVARS-this->nVARPOS)*_q.n_elem);
	mat S = inv(conv_to<mat>::from(Rq)) * conv_to<double>::from(mse);

	// Saving optimal data and corresponding covariance matrix in a vector of Gaussian distributions
	uint 		id1, id2;
	colvec 	_tmpMu;
	mat 		_tmpSigma;
	for(uint t = 0 ; t < _q.n_elem ; t++)
	{
		_tmpMu = optimalData.col(t);
		id1 = t*this->nVARPOS;
		id2 = (t+1)*this->nVARPOS-1;
		_tmpSigma = S.submat(id1,id1,id2,id2) * _q.n_elem;

		// Saving current Gaussian distribution
		outputProb.push_back(GaussianDistribution(_tmpMu, _tmpSigma));
	}

	return outputProb;
}
}	//end of pbdlib namespace




