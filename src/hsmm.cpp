/**
Copyright (C) 2015,	Martijn Zeestraten, Leonel Rozo, Ioannis Havoutis

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

#include "pbdlib/hsmm.h"
#include <algorithm>    // std::max
namespace pbdlib
{

HSMM::HSMM(uint _nSTATES, uint _nVARS):
				HMM(_nSTATES, _nVARS)
{
}

HSMM::HSMM(const std::string &priors_path, 
			const std::string &mu_path, 
			const std::string &sigma_path, 
			const std::string &transition_path, 
			const std::string &durMu_path, 
			const std::string &durSigma_path):
	HMM(priors_path,mu_path,sigma_path,transition_path)
{
	mat durMu, durSigma;
	
	// Load state duration components
	durMu.load(durMu_path, raw_ascii);
	durSigma.load(durSigma_path, raw_ascii);

	// Set state duration components:
	mat _SIGMA = zeros(1,1);
	colvec _MU = zeros(1,1);
	std::vector<GaussianDistribution> components;

	for(uint i=0; i<this->getNumSTATES(); i++){
		_MU(0,0) = durMu(0,i);
		_SIGMA(0,0) = durSigma(0,i);
		// for debugging
		components.push_back(GaussianDistribution(_MU, _SIGMA));
	}
	setDurationCOMPONENTS(components);

	// Forward variable calculation initialization:
	initializeFwdCalculation();
}

HSMM::HSMM(const std::string &priors_path,
			const std::string &mu_path,
			const std::string &sigma_path,
			const std::string &transition_path,
			const std::string &durMu_path,
			const std::string &durSigma_path,
			uint T, bool initDur):
	HMM(priors_path,mu_path,sigma_path,transition_path)
{
	mat durMu, durSigma;
	mat _SIGMA = zeros(1,1);
	colvec _MU = zeros(1,1);
	std::vector<GaussianDistribution> components;

	if(initDur)
	{
		// Load state duration components
		durMu.load(durMu_path, raw_ascii);
		durSigma.load(durSigma_path, raw_ascii);

		// Set state duration components:
		for(uint i=0; i<this->getNumSTATES(); i++){
			_MU(0,0) = durMu(0,i);
			_SIGMA(0,0) = durSigma(0,i);
			// for debugging
			components.push_back(GaussianDistribution(_MU, _SIGMA));
		}
	}
	else
	{
		// Set state duration components to zero:
		for(uint i=0; i<this->getNumSTATES(); i++)
			components.push_back(GaussianDistribution(_MU, _SIGMA));
	}
	setDurationCOMPONENTS(components);

	// Forward variable calculation initialization:
	initializeRecursiveForwardVariable(T,initDur);
}

void HSMM::setDurationCOMPONENTS(const std::vector<GaussianDistribution>& components)
{
	if(components.size() == this->getNumSTATES())
		DurationCOMPONENTS = components;
	else
		std::cout << "\n [ERROR]::HMM::setCOMPONENTS if(components.size() == nSTATES) ... else .";
}

void HSMM::setDurMU(int id, const colvec& _MU)
{
	DurationCOMPONENTS[id].setMU(_MU);
}

void HSMM::setDurSIGMA(int id, const mat& _SIGMA)
{
	DurationCOMPONENTS[id].setSIGMA(_SIGMA);
}

void HSMM::setMaxDuration(uint maxPd)
{
	this->PdSize = maxPd;

	// Initializing duration probabilities
	Pd.set_size(this->getNumSTATES(),PdSize);

	// Duration input vector
	rowvec dur = linspace<rowvec>(1,PdSize, PdSize);

	// Pre-Calculate Pd (duration probabilities)
	for (uint i = 0;i<this->getNumSTATES();i++)
		Pd.row(i) = getDurationCOMPONENTS(i).getPDFValue(dur).t(); // getPDFValue accepts a rowvec of values, but returns a colvec....
}

std::vector<GaussianDistribution>& HSMM::getDurationCOMPONENTS()
{
	return DurationCOMPONENTS;
}

GaussianDistribution& HSMM::getDurationCOMPONENTS(uint id)
{
	return DurationCOMPONENTS[id];
}

colvec& HSMM::getDurMU(uint id)
{
	return DurationCOMPONENTS[id].getMU();
}

mat& HSMM::getDurSIGMA(uint id)
{
	return DurationCOMPONENTS[id].getSIGMA();
}

void HSMM::initializeFwdCalculation()
{
	// Calculate PD size:
//	cout << "Pd size : ";
	PdSize=0;
	for (uint i = 0;i<this->getNumSTATES();i++)
	{
		if (PdSize <accu(this->getDurMU(i)))
			PdSize = accu(this->getDurMU(i));
	}
	PdSize *=1.2; // Enlarge the Size of Pd for 'accuracy'
	//PdSize = 50;
//	cout << PdSize << endl;

	// This function is used to set the forward variable matrices to the appropriate size:
	ALPHA.zeros(this->getNumSTATES(), PdSize);
	Atmp1.zeros(this->getNumSTATES(),PdSize-1);
	Atmp2.zeros(this->getNumSTATES(),PdSize-1);
	alpha.zeros(this->getNumSTATES());
	bmx.zeros(this->getNumSTATES());
	btmp.zeros(this->getNumSTATES());
	S.zeros(this->getNumSTATES());
	Pd.zeros(this->getNumSTATES(),PdSize);
	

	// Duration input vector
	mat dur = linspace(1,PdSize, PdSize);
	// Pre-Calculate Pd
	for (uint i = 0;i<this->getNumSTATES();i++)
	{
		// Note that we need to transpose twice....
		// getPDFValue accepts a rowvec of values, but returns a colvec....
		Pd.row(i) = getDurationCOMPONENTS(i).getPDFValue(dur.t()).t();

	}

	// For forward variable:
	Gtmp = new GaussianDistribution(this->getNumVARS());

	Initialized = false;
}


// -------------------------------------------------------------------
/* 				Forward Variable Step calculations 					*/
// -------------------------------------------------------------------

void HSMM::stepForwardVariable()
{
	if (Initialized==false)
	{
		initializeForwardVariable();
		Initialized = true;
	}
	else
	{
		lstepForwardVariable();
	}
}

void HSMM::stepForwardVariable(colvec& obs)
{
	if (Initialized==false)
	{
		initializeForwardVariable(obs);
		Initialized = true;
	}
	else
	{
		lstepForwardVariable(obs);
	}
}

void HSMM::stepForwardVariable(colvec& obs, urowvec& ind)
{
	if (Initialized==false)
	{
		initializeForwardVariable(obs, ind);
		Initialized = true;
	}
	else
	{
		lstepForwardVariable(obs, ind);
	}
}
// -------------------------------------------------------------------
// ----------  INITIALIZATION of Forward Variable
void  HSMM::initializeForwardVariable()
{
	// ALPHA Variable
	ALPHA = Pd;
	ALPHA.each_col() %= this->getPRIORS().t(); // % is the element wise product

	// Update S
	updateS(S,ALPHA);

	// Update Alpha
	updateAlpha(alpha,ALPHA);
}

void HSMM::initializeForwardVariable(colvec& _obs)
{
	// Calculate the initial forward step
	updateBtmp(btmp,_obs);

	// ALPHA Variable
	ALPHA = Pd;
	ALPHA.each_col() %= this->getPRIORS().t(); // % is the element wise product

	// Update bmx:
	updateBmx(bmx, ALPHA, btmp);

	// Update S
	updateS(S,ALPHA,bmx);

	// Update Alpha
	updateAlpha(alpha,ALPHA,btmp);
}


void HSMM::initializeForwardVariable(colvec& _obs, urowvec& _ind)
{
	// Calculate the initial forward step
	updateBtmp(btmp,_obs, _ind);

	
	// ALPHA Variable
	ALPHA = Pd;
	ALPHA.each_col() %= this->getPRIORS().t(); // % is the element wise product

	// Update bmx:
	updateBmx(bmx, ALPHA, btmp);

	// Update S
	updateS(S,ALPHA,bmx);

	// Update Alpha
	updateAlpha(alpha,ALPHA,btmp);
}


// -------------------------------------------------------------------
// ---------- Step functions of Forward Variable
void HSMM::lstepForwardVariable()
{
	// No observation, assume all btmp is 1 (i.e. observation is equally likeli for all
	// states	
	
	// ALPHA Variable
//	cout << "ALPHA: " << endl;
	updateALPHA(ALPHA,S);
//	cout << ALPHA.t() << endl;

	// Update S
//	cout << "S   : ";
	updateS(S,ALPHA);
//	cout << S << endl;

	// Update Alpha
//	cout << "alpha : " ;
	updateAlpha(alpha,ALPHA);
//	cout << alpha.t() << endl;
}

void HSMM::lstepForwardVariable(colvec& _obs)
{
	// update Btmp:
	updateBtmp(btmp,_obs);
	//cout << "Btmp: " << endl << btmp << endl;
	
	// ALPHA Variable (using current bmx)
	updateALPHA(ALPHA,S,bmx);

	// Update bmx:
	updateBmx(bmx, ALPHA, btmp);
	//cout << "Bmx: " << endl << bmx << endl;

	// Update S
	updateS(S,ALPHA,bmx);

	// Update Alpha
	updateAlpha(alpha,ALPHA,btmp);
	//cout << "Alpha: " << endl << alpha << endl;
}

void HSMM::lstepForwardVariable(colvec& _obs, urowvec& _ind)
{
	// update Btmp:
	updateBtmp(btmp,_obs,_ind);
	
	// ALPHA Variable (using current bmx)
	updateALPHA(ALPHA,S,bmx);

	// Update bmx:
	updateBmx(bmx, ALPHA, btmp);

	// Update S
	updateS(S,ALPHA,bmx);

	// Update Alpha
	updateAlpha(alpha,ALPHA,btmp);
}



// -------------------------------------------------------------------
// ---------- Prediction Functions for Forward Variable
mat& HSMM::predictForwardVariable(uint _N)
{
	//cout << "Allocating Memory: ";
	AlphaPred.set_size(this->getNumSTATES(), _N);

	// Make Predictions:
	predictForwardVariable(AlphaPred);

	// Return Alpha Variable:
	return AlphaPred;
}

// Implementation for real-time state prediction (no checks on sizes done)
// We assume that _AlphaPred = nbStates x nbPred of size
void HSMM::predictForwardVariable(mat& _AlphaPred)
{
	
	tmpInit = Initialized;
	// Make copy of the current state of the system:
	ALPHAtmp = ALPHA;
	Stmp = S;

	
	for (uint i = 0;i<_AlphaPred.n_cols;i++)
	{
		// Alpha variable
		if (tmpInit==false)
		{
			// Initialize: 
			ALPHAtmp = Pd;
			ALPHAtmp.each_col() %= this->getPRIORS().t(); // % is the element wise product
			tmpInit = true;
		}
		else
		{
			updateALPHA(ALPHAtmp,Stmp);
		}	

		// Update S
		updateS(Stmp,ALPHAtmp);

		// Update Alpha
		updateAlpha(alphatmp,ALPHAtmp);

		// Save alpha
		_AlphaPred.col(i) = alphatmp;
	//	cout << "Done" << endl;
	}
}

void HSMM::clear(){
	hsmm_transition.zeros();
	hsmm_transition_ticks.zeros();
	hsmm_priors.clear();
	hsmm_priors_ticks.clear();
	GMM_Model::clear();
}

void HSMM::integrateDemonstration(Demonstration demo) {
	//Compute sequence of states
	int nPoints = demo.getDatapoints().getNumPOINTS();
	mat h(this->getNumSTATES(), nPoints);
	for (int i=0; i<(int)this->getNumSTATES(); i++){
		h.row(i) = trans(this->getCOMPONENTS(i).getPDFValue(demo.getDatapoints().getData()));
	}
	uword imax[demo.getDatapoints().getNumPOINTS()];
	get_state_seq(imax, h);

	int nStates = this->getNumSTATES();
	int nVars = this->getNumVARS();

	hsmm_transition.resize(nStates, nStates);
	hsmm_transition_ticks.resize(nStates, nStates);
	hsmm_priors.resize(nStates);
	hsmm_priors_ticks.resize(nStates);

	//update hsmm priors for stochastic sampling
	hsmm_priors_ticks(imax[0]) += 1;
	hsmm_priors =  hsmm_priors_ticks / accu(hsmm_priors_ticks);

	//update duration statistics
	hsmm_duration_stats.resize(nStates);

	uint s0 = imax[0], currStateDuration = 1;
	for (int t = 0; t < nPoints; t++) {
		if ( (imax[t] != s0) || (t+1 == nPoints) ) {
			if ( (t+1 == nPoints) ) { currStateDuration += 1; }
			hsmm_duration_stats[s0](currStateDuration);
			currStateDuration = 0.0;
			hsmm_transition(s0, imax[t]) += 1.0;
			hsmm_transition_ticks(s0, imax[t]) += 1.0;
			s0 = imax[t];
		}
		currStateDuration += 1;
	}

	// normalize the transition matrix
	for (int i = 0; i < hsmm_transition.n_rows ; i++){
		double row_sum = accu(hsmm_transition.row(i));
		if (row_sum > 1){
			double row_sum_ticks = accu(hsmm_transition_ticks.row(i));
			hsmm_transition.row(i) = hsmm_transition_ticks.row(i) / row_sum_ticks;
		}
	}

	this->setTRANSITION(hsmm_transition);

	// Set hsmmd durations:
	mat _SIGMA = zeros(1,1);
	colvec _MU = zeros(1,1);
	std::vector<GaussianDistribution> hsmm_components;

	for (uint i=0; i<nStates; i++) {
		_MU(0,0) = hsmm_duration_stats[i].mean();//state_duration_means(i);
		_SIGMA(0,0) = std::max( hsmm_duration_stats[i].var(), 1.0);
		hsmm_components.push_back(GaussianDistribution(_MU, _SIGMA));
	}

	this->setDurationCOMPONENTS(hsmm_components);

}

void HSMM::get_state_seq(uword state_seq[], mat pred){
	for (int t=0; t<pred.n_cols; t++){
		vec vTmp = pred.col(t);
		vTmp.max(state_seq[t]);
	}
}

int HSMM::getClosestState(const colvec P){
	/// distance to closest cluster
	double d = arma::norm(this->getMU(0)-P,2); //Initialized distance of P from first cluster
	int minIndex = 0;	//Index corresponding to current cluster
	//Find cluster corresponding to minimum distance
	for (int k=1; k<this->getNumSTATES(); k++) {
		double Dist = arma::norm(this->getMU(k) - P, 2);
		if (Dist < d){
			d = Dist;
			minIndex = k; //Index corresponding to minimum distance
		}
	}
	return minIndex;
}

void HSMM::predictForwardVariableDeterministic(mat& _AlphaPred, int startState) {
	rowvec tempPriors = this->getPRIORS();
	rowvec fixed_priors(this->getNumSTATES(), fill::zeros);
	fixed_priors(startState) = 1.0;
	this->setPRIORS(fixed_priors);
	this->initializeFwdCalculation();
	predictForwardVariable(_AlphaPred);
	this->setPRIORS(tempPriors);
}

void HSMM::predictForwardVariableStochasticStart(mat& _AlphaPred, int startState) {
	rowvec tempPriors = this->getPRIORS();
	this->setPRIORS(hsmm_priors);
	this->initializeFwdCalculation();
	int nbSt = 0;
	int currTime = 0;
	rowvec iList(1, fill::zeros);
	int nbData = _AlphaPred.n_cols;
	int nStates = this->getNumSTATES();
	mat h = zeros(nStates, nbData);
	rowvec h1, h2;

	int nbD = this->Pd.n_cols;

	uword Iminmax;
	while (currTime < nbData) {
		if (nbSt==0){
			iList(0) = startState;
			h1 = ones<rowvec>(nbData);
		}else{
			h1 = join_rows( join_rows(
					zeros<rowvec>(currTime), cumsum(this->Pd.row(iList(iList.n_elem-2))) ),
					ones<rowvec>( std::max( (nbData-currTime-nbD),0) ));
			currTime = currTime + round( this->getDurMU( (uint)iList(iList.n_elem-2))(0,0) );
		}
		h2 = join_rows( join_rows(
				ones<rowvec>(currTime),
				ones<rowvec>( this->Pd.row(iList(iList.n_elem-1)).n_elem )
				- cumsum(this->Pd.row(iList(iList.n_elem-1))) ),
				zeros<rowvec>( std::max( (nbData-currTime-nbD),0) ));

		h.row(iList(iList.n_elem-1)) = h.row(iList(iList.n_elem-1))
										+ min( join_cols(h1.cols(0,nbData-1),h2.cols(0,nbData-1)) );
		vec tmp= this->getTRANSITION().row(iList(iList.n_elem-1)).t() % randu(nStates);
		tmp.max(Iminmax);
		iList.resize(iList.n_elem+1);
		iList(iList.n_elem-1) = Iminmax;

		nbSt = nbSt+1;
	}
	h = h / repmat(sum(h,0) ,nStates,1);
	_AlphaPred = h;
	this->setPRIORS(tempPriors);
}

void HSMM::predictForwardVariableStochastic(mat& _AlphaPred) {
	rowvec tempPriors = this->getPRIORS();
	this->setPRIORS(hsmm_priors);
	this->initializeFwdCalculation();
	int nbSt = 0;
	int currTime = 0;
	rowvec iList(1, fill::zeros);
	int nbData = _AlphaPred.n_cols;
	int nStates = this->getNumSTATES();
	mat h = zeros(nStates, nbData);
	rowvec h1, h2;

	int nbD = this->Pd.n_cols;

	uword Iminmax;
	while (currTime < nbData) {
		if (nbSt==0){
			vec tmp = this->getPRIORS().t() % randu(nStates);
			tmp.max(Iminmax);
			iList(0) = Iminmax;
			h1 = ones<rowvec>(nbData);
		}else{
			h1 = join_rows( join_rows(
					zeros<rowvec>(currTime), cumsum(this->Pd.row(iList(iList.n_elem-2))) ),
					ones<rowvec>( std::max( (nbData-currTime-nbD),0) ));
			currTime = currTime + round( this->getDurMU( (uint)iList(iList.n_elem-2))(0,0) );
		}
		h2 = join_rows( join_rows(
				ones<rowvec>(currTime),
				ones<rowvec>( this->Pd.row(iList(iList.n_elem-1)).n_elem )
				- cumsum(this->Pd.row(iList(iList.n_elem-1))) ),
				zeros<rowvec>( std::max( (nbData-currTime-nbD),0) ));

		h.row(iList(iList.n_elem-1)) = h.row(iList(iList.n_elem-1))
										+ min( join_cols(h1.cols(0,nbData-1),h2.cols(0,nbData-1)) );
		vec tmp= this->getTRANSITION().row(iList(iList.n_elem-1)).t() % randu(nStates);
		tmp.max(Iminmax);
		iList.resize(iList.n_elem+1);
		iList(iList.n_elem-1) = Iminmax;

		nbSt = nbSt+1;
	}
	h = h / repmat(sum(h,0) ,nStates,1);
	_AlphaPred = h;
	this->setPRIORS(tempPriors);
}

/*		FUNCTIONS FACILITATING FORWARD VARIABLE CALCULATION
 *
 *
 */
void HSMM::updateBtmp(colvec& _btmp, colvec& _obs)
{
	// Calculate the initial forward step
	for (uint i =0;i<this->getNumSTATES();i++)
	{
		_btmp(i) = this->getCOMPONENTS(i).getPDFValue(_obs)(0)+1e-12;
	}

	_btmp = _btmp/accu(_btmp);

}
void HSMM::updateBtmp(colvec& _btmp, colvec& _obs, urowvec& _ind)
{
	Gtmp->setNumVARS(_ind.n_elem);
	// Calculate the initial forward step
	for (uint i =0;i<this->getNumSTATES();i++)
	{
		// Create temporary Gaussian that uses only the specified indices:
		Gtmp->setMU(this->getMU(i)(_ind));
		Gtmp->setSIGMA(this->getSIGMA(i)(_ind,_ind));
		
		// Evaluate the gaussian probability:
		_btmp(i) = Gtmp->getPDFValue(_obs.rows(_ind))(0)+1e-12;
	}
	_btmp = _btmp/accu(_btmp);
}


// Equation (12): ALPHA MATRIX without observation:
void HSMM::updateALPHA(mat& _ALPHA, colvec& _S)
{
	// Help variables for vector-wise multiplications:
	Atmp1 = Pd.cols(0,PdSize-2);
	Atmp1.each_col() %= _S;
	
	Atmp2 = _ALPHA.cols(1,PdSize-1);

	_ALPHA.cols(0,PdSize-2) = Atmp2  + Atmp1;
	_ALPHA.col(PdSize-1) = _S % Pd.col(PdSize-1);
}
// Equation (12): ALPHA matrix update
void HSMM::updateALPHA(mat& _ALPHA, colvec& _S, colvec& _bmx)
{
	// Help variables for vector-wise multiplications:
	Atmp1 = Pd.cols(0,PdSize-2);
	Atmp1.each_col() %= _S;

	Atmp2 = _ALPHA.cols(1,PdSize-1);
	Atmp2.each_col() %= _bmx;

	_ALPHA.cols(0,PdSize-2) = Atmp2 + Atmp1;
	_ALPHA.col(PdSize-1) = _S % Pd.col(PdSize-1);
}

void HSMM::updateBmx(colvec& _bmx, mat& _ALPHA, colvec& _Btmp)
{
	// Equation (2) & (3): bmx update
	_bmx = _Btmp/accu(_Btmp.t()*sum(_ALPHA,1));
}	

// Equation (6): Update S
void HSMM::updateS(colvec& _S, mat& _ALPHA)
{
	// Equations (5) & (6) calculate S:
	_S = this->getTRANSITION().t()*_ALPHA.col(0);
}

// Equation (6): Update S
void HSMM::updateS(colvec& _S, mat& _ALPHA, colvec& _bmx)
{
	// Equations (5) & (6) calculate S:
	_S = this->getTRANSITION().t()*(_bmx%_ALPHA.col(0));
}

void HSMM::updateAlpha(colvec& _alpha,  mat& _ALPHA)
{
	// Calculate forward varaible
	_alpha = sum(_ALPHA,1);
	// Normalize
//	_alpha = _alpha/accu(_alpha);
}

void HSMM::updateAlpha(colvec& _alpha, mat& _ALPHA, colvec& _Btmp)
{
	// Calculate forward varaible
	_alpha = _Btmp%sum(_ALPHA,1);
	// Normalize
//	_alpha = _alpha/accu(_alpha);
}

void HSMM::saveInFiles(std::string path) {
	int nSTATES = this->getNumSTATES();
	int nVARS = this->getNumVARS();
	mat priors(1, nSTATES);
	mat mu(nVARS, nSTATES);
	mat sigma(nVARS, nVARS*nSTATES);
	mat transition(nSTATES, nSTATES);
	mat durMu(1, nSTATES);
	mat durSigma(1, nSTATES);

	for(uint i=0; i<nSTATES; i++){
		priors(0,i) = getPRIORS()(i);
			durMu(0,i) = getDurMU(i)(0);
		durSigma(0,i) = getDurSIGMA(i)(0);
		for(uint j=0; j<nVARS; j++){
			mu(j,i) = getMU(i)(j);
			for(uint k=0; k<nVARS; k++)
				sigma(j,k + i*(nVARS)) = getSIGMA(i)(j,k);
		}
	}
	transition = this->getTRANSITION();

	cout << "Saving in path: "<< path << endl;
	priors.save(path + "priors.txt", raw_ascii);
	mu.save(path + "mu.txt", raw_ascii);
	sigma.save(path + "sigma.txt", raw_ascii);
	transition.save(path + "transition.txt", raw_ascii);
	durMu.save(path + "durMu.txt", raw_ascii);
	durSigma.save(path + "durSigma.txt", raw_ascii);

	hsmm_transition.save(path + "hsmm_transition.txt", raw_ascii);
	hsmm_transition_ticks.save(path + "hsmm_transition_ticks.txt", raw_ascii);
	hsmm_priors.save(path + "hsmm_priors.txt", raw_ascii);
	hsmm_priors_ticks.save(path + "hsmm_priors_ticks.txt", raw_ascii);
}

void HSMM::loadFromFiles(std::string path) {
	this->clear();
	std::string priors_path(path + "priors.txt");
	std::string mu_path(path + "mu.txt");
	std::string sigma_path(path + "sigma.txt");
	std::string transition_path(path + "transition.txt");
	std::string durMu_path(path + "durMu.txt");
	std::string durSigma_path(path + "durSigma.txt");
	*this = HSMM(priors_path,
			mu_path,
			sigma_path,
			transition_path,
			durMu_path,
			durSigma_path);
	hsmm_transition.load(path + "hsmm_transition.txt", raw_ascii);
	hsmm_transition_ticks.load(path + "hsmm_transition_ticks.txt", raw_ascii);
	hsmm_priors.load(path + "hsmm_priors.txt", raw_ascii);
	hsmm_priors_ticks.load(path + "hsmm_priors_ticks.txt", raw_ascii);
}

// ---------------------------------------------------------------------
/* Forward Variable Step calculations - Recursive (readable) approach	*/
// ---------------------------------------------------------------------
void HSMM::initializeRecursiveForwardVariable(uint T, bool initDur)
{
	// Initializing scaling factor
	scalingFtr.set_size(T);
	scalingFtr.zeros();
	scalingFtr(0) = 1;
	// Initializing recursive forward variable
	recAlpha.set_size(this->getNumSTATES(),T);
	recAlpha.zeros();
	// For forward variable:
	Gtmp = new GaussianDistribution(this->getNumVARS());

	// Initializing d
	if(initDur)
	{
		cout << "Pd size : ";
		PdSize=0;
		for (uint i = 0;i<this->getNumSTATES();i++)
		{
			if (PdSize <accu(this->getDurMU(i)))
				PdSize = accu(this->getDurMU(i));
		}
		PdSize *=1.2; // Enlarge the Size of Pd for 'accuracy'
		cout << PdSize << endl;

		// Initializing duration probabilities
		Pd.set_size(this->getNumSTATES(),PdSize);

		// Duration input vector
		rowvec dur = linspace<rowvec>(1,PdSize, PdSize);

		// Pre-Calculate Pd (duration probabilities)
		for (uint i = 0;i<this->getNumSTATES();i++)
			Pd.row(i) = getDurationCOMPONENTS(i).getPDFValue(dur).t(); // getPDFValue accepts a rowvec of values, but returns a colvec....
	}
	else
	{
		cout << "[WARNING] PdSize and Pd variables were not initialized. Received flag was set to false." << endl;
	}
}

void HSMM::stepRecursiveForwardVariable(uint tn)
{
	for(uint i = 0 ; i < this->getNumSTATES() ; i++)
	{
		// \alpha_i,t = \Pi_i * Pd_i(t) + ...
		if(tn < PdSize)
			recAlpha(i,tn) = this->getPRIORS()(i) * Pd(i,tn);

		// \alpha_i,t = \sum_d \alpha_j,t-d * a_j,i * Pd_i(d)
		for(uint d = 0 ; d < std::min(tn,PdSize) ; d++)
			recAlpha(i,tn) = recAlpha(i,tn) + as_scalar(recAlpha.col(tn-d-1).t() * this->getTRANSITION().col(i)) * Pd(i,d);
	}
}

void HSMM::stepRecursiveForwardVariable(uint tn, mat& obs)
{
	double tmpObsProb;

	for(uint i = 0 ; i < this->getNumSTATES() ; i++)
	{
		// \alpha_i,t = \Pi_i * Pd_i(t) * \prod_s N(\ksi_s) + ...
		if(tn < PdSize)
		{
		 	tmpObsProb = arma::prod( scalingFtr.subvec(0, tn) % this->getCOMPONENTS(i).getPDFValue(obs) );
			recAlpha(i,tn) = this->getPRIORS()(i) * Pd(i,tn) * tmpObsProb;
		}

		// \alpha_i,t = \sum_d \alpha_j,t-d * a_j,i * Pd_i(d) * \prod_s N(\ksi_s)
		for(uint d = 0 ; d < std::min(tn,PdSize) ; d++)
		{
			tmpObsProb = arma::prod( scalingFtr.subvec(tn-d, tn) % this->getCOMPONENTS(i).getPDFValue(obs.cols(tn-d,tn)) );
			recAlpha(i,tn) = recAlpha(i,tn) + as_scalar(recAlpha.col(tn-d-1).t() * this->getTRANSITION().col(i)) * Pd(i,d) * tmpObsProb;
		}
	}

	// Computing scaling factor
	if(tn < scalingFtr.n_rows-1)
		scalingFtr(tn+1) = 1/sum(recAlpha.col(tn));
}

void HSMM::stepRecursiveForwardVariable(uint tn, mat& obs, urowvec ind)
{
	double tmpObsProb;

	for(uint i = 0 ; i < this->getNumSTATES() ; i++)
		{
		// \alpha_i,t = \Pi_i * Pd_i(t) * \prod_s N(\ksi_s) + ...
			if(tn < PdSize)
			{
				tmpObsProb = arma::prod( scalingFtr.subvec(0, tn) % this->getCOMPONENTS(i).getPDFValue(obs,ind) );
			 	recAlpha(i,tn) = this->getPRIORS()(i) * Pd(i,tn) * tmpObsProb;
			}

			// \alpha_i,t = \sum_d \alpha_j,t-d * a_j,i * Pd_i(d) * \prod_s N(\ksi_s)
			for(uint d = 0 ; d < std::min(tn,PdSize) ; d++)
			{
				tmpObsProb = arma::prod( scalingFtr.subvec(tn-d, tn) % this->getCOMPONENTS(i).getPDFValue(obs.cols(tn-d,tn),ind) );
				recAlpha(i,tn) = recAlpha(i,tn) + as_scalar(recAlpha.col(tn-d-1).t() * this->getTRANSITION().col(i)) * Pd(i,d) * tmpObsProb;
			}
		}

		if(tn < scalingFtr.n_rows-1)
			scalingFtr(tn+1) = 1/sum(recAlpha.col(tn));
}

void HSMM::resetRecursiveForwardVariable()
{
	recAlpha.zeros();
}
} // End pbdlib namespace
