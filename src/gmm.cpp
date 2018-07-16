/**
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

#include "pbdlib/gmm.h"

namespace pbdlib
{

GMM_Model::GMM_Model(std::vector<Demonstration> &demos, uint _nSTATES)
{
	DEMONSTRATIONS = demos;
	nVARS = demos[0].getDatapoints().getNumVARS();
	nSTATES = _nSTATES;
	PRIORS = rowvec(_nSTATES);
	setVARSNames(demos[0].getDatapoints().getVarNames());

	// Initialize Components
	COMPONENTS.resize(nSTATES,GaussianDistribution(nVARS));
}

GMM_Model::GMM_Model(uint _nSTATES, uint _nVARS)
{
	nVARS = _nVARS;
	nSTATES = _nSTATES;
	PRIORS = rowvec(_nSTATES);
	
	// Initialize Components
	COMPONENTS.resize(nSTATES,GaussianDistribution(nVARS));
}

GMM_Model::GMM_Model(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path,const std::string &vars_path)
{
	mat priors, mu, sigma;

	priors.load(priors_path, raw_ascii);
	mu.load(mu_path, raw_ascii);
	sigma.load(sigma_path, raw_ascii);

	nVARS = mu.n_rows;
	nSTATES = priors.n_elem;

	// Load Var names from file.
	std::ifstream varsfile(vars_path.c_str()); // Open file stream
	std::vector<std::string> vars; 			   // Vector to hold variable names
	std::string varsnames;                     // tmp string to read line from file
	std::string buf; 						   // buffer string
	// We assume that the variable names are on the first line, separated by spaces
  	getline(varsfile,varsnames);
	std::stringstream ss(varsnames);
	// Read all strings:
	while (ss>>buf)
	{
		vars.push_back(buf);
	}

	setVARSNames(vars);

	// Load Sigma and MU
	mat _SIGMA = zeros(nVARS, nVARS);
	colvec _MU = zeros(nVARS,1);
	rowvec _PRIORS(1, nSTATES);
	std::vector<GaussianDistribution> components;

	for(uint i=0; i<nSTATES; i++){
		_PRIORS(0,i) = priors(0,i);
		for(uint j=0; j<nVARS; j++){
			_MU(j,0) = mu(j,i);
			for(uint k=0; k<nVARS; k++)
				_SIGMA(j,k) = sigma(j,k + i*(nVARS));
		}
		components.push_back(GaussianDistribution(_MU, _SIGMA));
	}

	setPRIORS(_PRIORS);
	setCOMPONENTS(components);
}

GMM_Model::GMM_Model(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path)
{
	// Constructor which does not require to specify variable names
	mat priors, mu, sigma;
	
	priors.load(priors_path, raw_ascii);
	mu.load(mu_path, raw_ascii);
	sigma.load(sigma_path, raw_ascii);

	nVARS = mu.n_rows;
	nSTATES = priors.n_elem;

	mat _SIGMA = zeros(nVARS, nVARS);
	colvec _MU = zeros(nVARS,1);
	rowvec _PRIORS(1, nSTATES);
	std::vector<GaussianDistribution> components;

	for(uint i=0; i<nSTATES; i++){
		_PRIORS(0,i) = priors(0,i);
		for(uint j=0; j<nVARS; j++){
			_MU(j,0) = mu(j,i);
			for(uint k=0; k<nVARS; k++)
				_SIGMA(j,k) = sigma(j,k + i*(nVARS));
		}
		components.push_back(GaussianDistribution(_MU, _SIGMA));
	}

	setPRIORS(_PRIORS);
	setCOMPONENTS(components);
}

void GMM_Model::saveInFiles()
{
	mat priors(1, nSTATES);
	mat mu(nVARS, nSTATES);
	mat sigma(nVARS, nVARS*nSTATES);
	std::ofstream varsfile ("data/gmm/GMM_vars.txt");

	for(uint i=0; i<nSTATES; i++){
		priors(0,i) = getPRIORS()(i);
		for(uint j=0; j<nVARS; j++){
			mu(j,i) = getMU(i)(j);
			for(uint k=0; k<nVARS; k++)
				sigma(j,k + i*(nSTATES+1)) = getSIGMA(i)(j,k);
		}
	}

	for(uint i=0; i<nVARS; i++){
		varsfile << vars_names[i];
		if(i<nVARS-1)
			varsfile << " ";
	}

	priors.save("data/gmm/GMM_priors.txt", raw_ascii);
	mu.save("data/gmm/GMM_mu.txt", raw_ascii);
	sigma.save("data/gmm/GMM_sigma.txt", raw_ascii);

	varsfile.close();
}

void GMM_Model::setDEMONSTRATIONS(std::vector<Demonstration> demons)
{
	DEMONSTRATIONS = demons;
}

void GMM_Model::setPRIORS(const rowvec& priors)
{
	if(priors.n_elem == nSTATES)
		PRIORS = priors;
	else

		std::cout << "\n [ERROR]::GMM_Model::setPRIORS if(priors.n_elem == nSTATES) ... else .";
}


void GMM_Model::setCOMPONENTS(const std::vector<GaussianDistribution>& components)
{
	if(components.size() == nSTATES)
		COMPONENTS = components;
	else
		std::cout << "\n [ERROR]::GMM_Model::setCOMPONENTS if(components.size() == nSTATES) ... else .";
}

void GMM_Model::setMU(int id, const colvec& _MU)
{
	COMPONENTS[id].setMU(_MU);
}

void GMM_Model::setSIGMA(int id, const mat& _SIGMA)
{
	COMPONENTS[id].setSIGMA(_SIGMA);
}

void GMM_Model::setLAMBDA(int id, const mat& _LAMBDA)
{
	COMPONENTS[id].setLAMBDA(_LAMBDA);
}

uint GMM_Model::getIndexOfVARName(const std::string& varname)
{
	std::vector<std::string>::iterator it;
	it = std::find (vars_names.begin(), vars_names.end(), varname);
	return std::distance(vars_names.begin(), it);
}

uint GMM_Model::getNumVARS()
{
	return nVARS;
}

uint GMM_Model::getNumSTATES()
{
	return nSTATES;
}

rowvec& GMM_Model::getPRIORS()
{
	return PRIORS;
}

double GMM_Model::getPRIORS(int id)
{
	return PRIORS[id];
}

std::vector<GaussianDistribution>& GMM_Model::getCOMPONENTS()
{
	return COMPONENTS;
}

GaussianDistribution& GMM_Model::getCOMPONENTS(int id)
{
	return COMPONENTS[id];
}

colvec& GMM_Model::getMU(int id)
{
	return COMPONENTS[id].getMU();
}

mat& GMM_Model::getSIGMA(int id)
{
	return COMPONENTS[id].getSIGMA();
}
mat& GMM_Model::getLAMBDA(int id)
{
	return COMPONENTS[id].getLAMBDA();
}

double GMM_Model::getProbability(const colvec& sample, bool usePriors)
{
	double P = 0.0;
	// See Eq. (2.0.2) in doc/TechnicalReport.pdf
	for(uint k=0; k<nSTATES; k++)
	{
		if(usePriors)
			P += PRIORS[k] * as_scalar( COMPONENTS[k].getPDFValue(sample) ); // See Eq. (2.0.3) for "getPDFValue" in doc/TechnicalReport.pdf
		else
			P += as_scalar( COMPONENTS[k].getPDFValue(sample) ); // See Eq. (2.0.3) for "getPDFValue" in doc/TechnicalReport.pdf
	}
	return P;
}

colvec GMM_Model::getProbability(const mat& samples, bool usePriors)
{
	colvec probs(samples.n_cols);
	colvec s(samples.n_rows);

	// See Eq. (2.0.4) in doc/TechnicalReport.pdf
	for(uint j=0 ; j<samples.n_cols ; j++){
		s = samples.col(j);
		probs(j) = log( getProbability(s, usePriors) );
	}
	return probs;
}

void GMM_Model::onlineEMDP(int N,colvec P,double lambda, double minSigma)
{	
	/*
	 Performs online GMM clustering by using an updated DP-Means algorithms
	Ref:
		Kulis, B. and Jordan, M. I. (2012). Revisiting k-means: New algorithms via bayesian nonparametrics. In Proc. Intl Conf. on Machine Learning (ICML)

	Input:
	N: Number of points processed
	P: current point being added to GMM
	lambda: splitting distance
	nimSigma: minimum covariance for regularization
	 */
	double d = arma::norm(this->getMU(0)-P,2); //Initialized distance of P from first cluster
	int minIndex = 0;	//Index corresponding to current cluster
	//Find cluster corresponding to minimum distance
	for(int k=1;k<nSTATES;k++){
		double Dist = arma::norm(this->getMU(k) - P,2);
		if (Dist < d){
				d = Dist;
				minIndex = k; //Index corresponding to minimum distance
		}
	}
	//Allocate new cluster if distance of each component higher than lambda
	if (lambda < d){
		minIndex = nSTATES;
		mat SigmaTmp = minSigma*eye(nVARS,nVARS); //Sigma of new component is the minimum cov
		GaussianDistribution newComponent(P,SigmaTmp); //Mean of new component is P
		COMPONENTS.push_back(newComponent); //add new component
		rowvec priorsTmp = zeros(1,nSTATES+1);
		priorsTmp.cols(0,nSTATES-1) = this->getPRIORS();
		priorsTmp(nSTATES) = 1./N; //prior for new component inversely proportional to #P
		priorsTmp = priorsTmp/arma::norm(priorsTmp,1); //evaluate new priors
		this->nSTATES = nSTATES +1; //update number of states
		this->setPRIORS(priorsTmp);
	}
	else{	
		/*
		*Update components belonging to P by using MAP estimate
		Ref:
		Gauvain, J.-L. and Lee, C.-H. (1994). Maximum a pos- teriori estimation for multivariate gaussian mixture observations of markov chians. IEE Transactions on Speech and Audio Processing, 2(2).  */
		double PriorsTmp = 1./N + PRIORS[minIndex];
		colvec MuTmp = 	1./PriorsTmp*(PRIORS[minIndex]*this->getMU(minIndex)+P/N);
		mat SigmaTmp = PRIORS[minIndex]/PriorsTmp*(this->getSIGMA(minIndex)+(this->getMU(minIndex)-MuTmp)*(this->getMU(minIndex)-MuTmp).t()) + 1./(N*PriorsTmp)*(minSigma*eye(nVARS,nVARS)+(P-MuTmp)*(P-MuTmp).t());
		this->setMU(minIndex,MuTmp);
		this->setSIGMA(minIndex,SigmaTmp);
		rowvec priors = this->getPRIORS();
		priors[minIndex] = PriorsTmp;
		priors = priors/arma::norm(priors,1);
		this->setPRIORS(priors);
	}
}

static uvec randperm( int n )
{
	//return sort_index(randn<urowvec>(nbData));
	return sort_index(randu<colvec>(n)); // was randn, should work equally well?
}

void GMM_Model::learnKMEANS(double regularization)
{

	mat Data = DEMONSTRATIONS[0].getDatapoints().getData();

	for(int i=1; i<DEMONSTRATIONS.size(); i++)
		Data.insert_cols(DEMONSTRATIONS[0].getDatapoints().getNumPOINTS(), DEMONSTRATIONS[i].getDatapoints().getData());

	//Criterion to stop the EM iterative update
	double cumdist_threshold = 1e-10;
	uint maxIter = 100;

	//Initialization of the parameters
	uint nbVar = Data.n_rows;
	uint nbData = Data.n_cols;

	double cumdist_old = -std::numeric_limits<double>::max();
	uint nbStep = 0;

	//srand (time(NULL));
	// random permutation 
	cout << "rand perm" << endl;
	uvec idTmp = randperm(nbData);

	cout << "linspace" << endl;
	uvec allrows = linspace<uvec>(0, nbVar-1, nbVar);

	// nSTATES means
	cout << "submat" << endl;
	mat Mu = Data.submat(allrows, idTmp.subvec(0,nSTATES-1));

	//k-means iterations
	while(true)
	{
		mat distTmp = zeros(nbData,nSTATES);
		
		//E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for (uint i=0; i<nSTATES; i++){
			//Compute distances
			distTmp.col(i) = trans( sum( pow(Data - repmat(Mu.col(i), 1, nbData), 2.0) ) );
		}

		vec vTmp = zeros<vec>(nbData,1);
		uvec idList = zeros<uvec>(nbData,1);

		for( uint i = 0; i < nbData; i++ )
		{
			// there is a chance that here two elements will be as close and will skip one?s
			vTmp[i] = ((rowvec)distTmp.row(i)).min(idList[i]);
		}

		double cumdist = sum(vTmp);

		//M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for (uint i=0; i<nSTATES; i++)
		{
			idTmp = find(idList == i);
			Mu.col(i) = mean(Data.submat(allrows, idTmp), 1);			
		}
		
		//Stopping criterion %%%%%%%%%%%%%%%%%%%%
		if (fabs(cumdist-cumdist_old) < cumdist_threshold){
			cout << "%%%%%%%%%%% KMEANS %%%%%%%%%%%%%%%%%%" << endl;
			cout << nbStep << "Steps" << endl;
			cout << endl << "Mu_kmeans :" << endl << Mu << endl;
			

			GaussianDistribution *GD;
			colvec mu;

			for(int i=0; i<nSTATES; i++){
				idTmp = find(idList == i);
				COMPONENTS[i].setMU( Mu.col(i) );

				PRIORS(i) = idTmp.n_elem; 

				COMPONENTS[i].setSIGMA( 
					cov( trans(join_rows(Data.submat(allrows, idTmp),Data.submat(allrows, idTmp))) ,0 )
					+ eye(nbVar, nbVar)*regularization
				);
			}

			PRIORS = PRIORS / nbData;
			cout << "PRIORS:" << endl << PRIORS << endl;
			
			cout << "%%%%%%% KMEANS FINISHED %%%%%%%%%%%%%" << endl;
			return;
		}
		cumdist_old = cumdist;
		nbStep = nbStep + 1;

	}
}

uint GMM_Model::EM_learn(double regularization)
{
    learnKMEANS(regularization);
    return EM(std::numeric_limits<double>::min(),regularization);
}

double GMM_Model::getLikelihood(const mat &SAMPLES)
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

bool GMM_Model::EM_isfinished(double l_old, double l_new)
{
	//cout << "\n error=" <<  fabs((l_new/l_old) - 1.0);
	if ( fabs((l_new/l_old) - 1.0)  <= THRESHOLD_MIN)
		return true;
	return false;
}


uint GMM_Model::EM(double likelihood,double regularization)
{
	double likelihood_new;

	uint k,i;
	rowvec E;
	colvec mu_tmp;
	mat Pxi, Pix, Pix_tmp, DataTmp1, sigma_tmp;

	mat demos = DEMONSTRATIONS[0].getDatapoints().getData();

	for(int i=1; i<DEMONSTRATIONS.size(); i++)
		demos.insert_cols(DEMONSTRATIONS[0].getDatapoints().getNumPOINTS(), DEMONSTRATIONS[i].getDatapoints().getData());

	//E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// See Eq. (2.0.5) in doc/TechnicalReport.pdf
	rowvec prior = PRIORS;

	Pxi = zeros(demos.n_cols, nSTATES);
	for(k=0; k<nSTATES; k++){
		Pxi.col(k) = COMPONENTS[k].getPDFValue( demos );
	}
	// compute posterior probability p(i|x)
	Pix_tmp = repmat(prior, demos.n_cols, 1) % Pxi;
	Pix = Pix_tmp / repmat(sum(Pix_tmp, 1), 1, nSTATES);
	this->gamma = Pix;
	// compute cumulated posterior probability
	E = sum(Pix);

	//M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// See Eq. (2.0.6),(2.0.7) and (2.0.8) in doc/TechnicalReport.pdf
	for (i=0; i<nSTATES; i++){
		// update the priors
		PRIORS(i) = E(i) / demos.n_cols;
		// update the centers
		mu_tmp = demos * Pix.col(i) / E(i);
		COMPONENTS[i].setMU(mu_tmp);
		// update the covariance matrices
		DataTmp1 = demos - repmat( COMPONENTS[i].getMU(), 1, demos.n_cols);
		sigma_tmp = (repmat(trans(Pix.col(i)), nVARS, 1) % DataTmp1 * trans(DataTmp1)) / E(i);
        sigma_tmp = sigma_tmp  + regularization * eye(nVARS,nVARS);
		COMPONENTS[i].setSIGMA( sigma_tmp );
	}


	likelihood_new = getLikelihood( demos ); // See Eq. (2.0.4) in doc/TechnicalReport.pdf

	if( EM_isfinished(likelihood, likelihood_new) )
		return 1;

    return  EM(likelihood_new,regularization)+1;
}

void GMM_Model::setVARSNames(const std::vector<std::string>& vars)
{
	if(vars.size() == nVARS)
		vars_names = vars;
	//else
	//	std::cout << "\n [ERROR]::GMM_Model::setVARSNames if(vars.size() == nVARS) ... else .";
}

std::vector<std::string>& GMM_Model::getVARSNames()
{
	return vars_names;
}

std::string GMM_Model::getVARSNames(int id)
{
	return vars_names[id];
}


bool GMM_Model::addDemo(Demonstration& demo)
{
	if(demo.getDatapoints().getNumVARS() == nVARS)
		DEMONSTRATIONS.push_back(demo);
	else
		return false;

	return true;
}

void GMM_Model::clear(){
	COMPONENTS.clear();
	nSTATES=1;
	setPRIORS(ones<vec>(1));
}

mat& GMM_Model::getGamma()
{
	return this->gamma;
}
} //end of pbdlib namespace

