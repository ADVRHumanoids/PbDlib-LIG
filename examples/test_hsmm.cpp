/**
Copyright (C) 2015, Martijn Zeestraten,Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon

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

/*! \file test_hsmm.cpp
\brief Learning GMM model
Learning a GMM model from a demonstration saved in the file data_txyz.txt

\author Martijn Zeestraten, Sylvain Calinon
\bug No known bugs.
*/

#include "pbdlib/hsmm.h"
#include "armadillo"

using namespace arma;
using namespace pbdlib;
using namespace std;

int main()
{

	// Create Object
	HSMM* myHSMM;
	
	// Load data from file
	string muPath = "data/hsmm/HSMM_test_mu.txt";
	string sigmaPath = "data/hsmm/HSMM_test_sigma.txt";
	string priorsPath = "data/hsmm/HSMM_test_priors.txt";
	string transPath = "data/hsmm/HSMM_test_trans.txt";
	string durMu = "data/hsmm/HSMM_test_durMu.txt";
	string durSigma = "data/hsmm/HSMM_test_durSigma.txt";
	
	myHSMM = new HSMM(priorsPath,muPath,sigmaPath,transPath,durMu,durSigma);
	

	colvec X = myHSMM->getMU(0);
	// Display Data
	
	cout << "--- HSMM components: " << endl;
	cout << "nStates " << myHSMM->getNumSTATES() <<endl;
	for (uint i=0;i<myHSMM->getNumSTATES();i++)
	{
		cout << "Priors" << endl;
		cout << myHSMM->getPRIORS(i) << endl;
		cout << "Mu" << endl;
		cout << myHSMM->getMU(i) << endl;
		cout << "Sigma" << endl;
		cout << myHSMM->getSIGMA(i) << endl;
		cout << "DurMu" << endl;
		cout << myHSMM->getDurMU(i) << endl;
		cout << "DurSigma" << endl;
		cout << myHSMM->getDurSIGMA(i) << endl;
	}

	// Test alpha variable (Let the model step based on observations)
 	cout << "--> Forward variable steps: " << endl;
	urowvec in;
	in << 0 << 1;
	cout << "Observations (constant): " << X.t() << endl;

	cout << "Alpha values:" << endl;
	for (uint i = 1;i<100;i++)
	{
		myHSMM->stepForwardVariable();

		if (i%10 ==0)
		{
		cout << "Alpha " << i << " : [" <<
		   	myHSMM->getForwardVariable()(0) << ", " <<
		   	myHSMM->getForwardVariable()(1) << ", " <<
		   	myHSMM->getForwardVariable()(2) << "]" << endl; 
		}
	}


	// Test Alpha Predictions
	// (Let the model generate predictions (i.e. state of the model is not adjusted)
	mat pred;

	myHSMM->resetForwardVariable(); // Reset the forward variable
	
	// Offline method (function allocates matrix itself)
	pred = myHSMM->predictForwardVariable(200);
	
	// Online method (we pre-allocate matrix ourselves)
	pred.set_size(3,100);
	myHSMM->predictForwardVariable(pred);

	cout << "Alpha Predictions" << endl << pred.t() << endl;

}

