/**
Copyright (C) 2014, Davide De Tommaso, Sylvain Calinon

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

#include "pbdlib/gmr.h"

namespace pbdlib
{

GMR::GMR(GMM_Model* model)
{
	gmm = NULL;
	setGMMModel(model);
}
GMM_Model* GMR::regression(Datapoints* data_in)
{
	// Change note 4/03/2015, Martijn Zeestraten
	// This regression function is the initial version of the 'regression', it uses variable names
	// from the data_in object to select the input and output variables. This function was quite cumbersome
	// and now an overload is implemented which more closely resembles the Matlab version
	//
	// The actual regression is performed by regression(data_in, in, out) while this function only
	// determines the input and output variables.
	uint i,j,k;
	
	// Check if variable names exist:
	if(data_in->getVarNames().size() == 0){
		std::cout << "\n [ERROR]::GMR::regression() if(data_in->getVarNames().size() == 0) ... else .";
		//return 0;
	}

	// Extract the input and output variable indices
	urowvec in,out;		// row vectors to hold the input and output indices
	bool found = false;	// 

	in.zeros(1, data_in->getNumVARS());	// Size of In is the number of vars in the input
	out.zeros(1, gmm->getNumVARS() - data_in->getNumVARS());// Size of Out is the diff between input and GMM size

	// Determine the input variables
	for(j=0; j<data_in->getNumVARS(); j++)
		// Check which index in GMM has varname
		for(i=0; i<gmm->getNumVARS(); i++)
			if ( data_in->getVarName(j).compare( gmm->getVARSNames(i) ) == 0)
				in(j) = i; // found, store variable index
	// Remainder of the variables is output variable
	for(i=0, k=0; i<gmm->getNumVARS(); i++){
		for(j=0; j<data_in->getNumVARS(); j++)
			if ( data_in->getVarName(j).compare( gmm->getVARSNames(i) ) == 0)
				found = true;
			if(!found)
				out(k++) = i;
		found = false;
	}


	// Extract the input data and put it into a mat shape
	mat x = data_in->getData();
	// Perform actual regression and return results
	return regression(x,in,out);
}

GMM_Model* GMR::regression( mat data_in, urowvec in, urowvec out)
{
	// This regression performs regression of the model and returns
	// a new GMM_Model object.
	
	// Create the GMM_model
	GMM_Model* gmmOut = new GMM_Model(data_in.n_cols, out.n_elem);
	// Perform regression
	regression(gmmOut,data_in,in,out);
	// Return the model
	return gmmOut;

	

}


void  GMR::regression(GMM_Model* gmmOut, mat data_in, urowvec in, urowvec out)
{
	// Base regression function. Use this one to be as fast as possible. 
	// One need to supply a gmmOut object of appropriate size (nbStates = data_in.n_cols, nbVar = size(out))
	
	
	Pxi = zeros(data_in.n_cols,gmm->getNumSTATES());
	//Compute the influence of each GMM component, given input x
	// See Eq. (3.0.5) in doc/TechnicalReport.pdf
	for(i=0; i<gmm->getNumSTATES(); i++){
		Mu = gmm->getMU(i)(in);
		Sigma = gmm->getSIGMA(i)(in,in);
		//Mu.print("Mu(i) = ");			// DEBUGGING
		//Sigma.print("Sigma(i) = ");	// DEBUGGING
		Gtmp->setMU(Mu);
		Gtmp->setSIGMA(Sigma);
		Pxi.col(i) = gmm->getPRIORS(i) * (Gtmp->getPDFValue(data_in));
	}
	//Priors.print("Priors = ");	// DEBUGGING
	//data_in->getData().print("t = ");	// DEBUGGING
	//Pxi.print("Pxi = ");	// DEBUGGING
	beta = Pxi / repmat(sum(Pxi,1),1,gmm->getNumSTATES()); // See Eq. (3.0.5) in doc/TechnicalReport.pdf

	//beta.print("beta = ");	// DEBUGGING
	for(t=0; t<data_in.n_cols; t++){
		MuOut.zeros(out.n_elem);
		SigmaOut.zeros(out.n_elem,out.n_elem);
		for(i=0; i<gmm->getNumSTATES(); i++){
			Mu = gmm->getMU(i);
			Sigma = gmm->getSIGMA(i);
			// Pre compute inverse, use inv_sympd to improve performance
			InvSigmaInIn = (Sigma(in,in).i());

			MuOutTmp = Mu(out) + Sigma(out,in) * InvSigmaInIn * (data_in.col(t)-Mu(in)); // See Eq. (3.0.3) in doc/TechnicalReport.pdf
			SigmaOutTmp = Sigma(out,out) - Sigma(out,in) * InvSigmaInIn * Sigma(in,out); // See Eq. (3.0.4) in doc/TechnicalReport.pdf
			MuOut = MuOut + beta(t,i) * MuOutTmp; // See Eq. (3.0.2) in doc/TechnicalReport.pdf
			SigmaOut = SigmaOut + beta(t,i) * SigmaOutTmp; // See Eq. (3.0.2) in doc/TechnicalReport.pdf
		}
//		COMPONENTS.push_back( GaussianDistribution(MuOut, SigmaOut) );
		gmmOut->getCOMPONENTS(t).setMU(MuOut);
		gmmOut->getCOMPONENTS(t).setSIGMA(SigmaOut);
			
	}

//	GMM_Model* gmmOut;
//	gmmOut->setCOMPONENTS(COMPONENTS);
//	return gmmOut;

}

void GMR::setGMMModel(GMM_Model* gmmmodel)
{
	gmm = gmmmodel;
	
	// Initial settings for regression 
	Gtmp= new GaussianDistribution(gmm->getNumSTATES());
}

} //end of pbdlib namespace

