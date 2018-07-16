/**
Copyright (C) 2015, Leonel Rozo, Martijn Zeestraten, Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon

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

/*! \file test_hsmm2.cpp
\brief Learning HSMM model
Testing the computation of the forward variable of an HSMM learned from demonstrations. Forward variable is computed using the classical (readable) formulation.

\author Leonel Rozo
\bug No known bugs.
*/

#include "pbdlib/hsmm.h"
#include "armadillo"

using namespace arma;
using namespace pbdlib;
using namespace std;

int main()
{

	// Create HSMM object
	HSMM* myHSMM;

	// Load data from files
	string muPath = "data/hsmm/HSMM_test2_mu.txt";
	string sigmaPath = "data/hsmm/HSMM_test2_sigma.txt";
	string priorsPath = "data/hsmm/HSMM_test2_priors.txt";
	string transPath = "data/hsmm/HSMM_test2_trans.txt";
	string durMu = "data/hsmm/HSMM_test2_durMu.txt";
	string durSigma = "data/hsmm/HSMM_test2_durSigma.txt";
	string testData = "data/hsmm/HSMM_test2_testData.txt";

	// Loading test data
	mat test;
	test.load(testData, raw_ascii);

	// Creating HSMM
	myHSMM = new HSMM(priorsPath,muPath,sigmaPath,transPath,durMu,durSigma,test.n_cols);
	myHSMM->setMaxDuration(80);

	// Display HSMM parameters
	cout << "--- HSMM components: " << endl;
	cout << "nStates " << myHSMM->getNumSTATES() <<endl;
	for (uint i=0 ; i<myHSMM->getNumSTATES() ; i++)
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


	// --> Testing duration-based forward variable
	cout << "Press [ENTER] to test duration-based forward variable...";
	cin.get();

	uint T = test.n_cols;
	// Loop for T samples
	for(uint t=0 ; t<T ; t++)
		myHSMM->stepRecursiveForwardVariable(t);

	// Getting alpha variable and normalizing it
	mat alpha = myHSMM->getRecursiveForwardVariable();	// Getting forward variable
	alpha = alpha / repmat(sum(alpha,0),myHSMM->getNumSTATES(),1);	// Normalizing forward variable

	// Printing results
	for(uint t=0 ; t<T ; t++)
	{
		if (t%10 == 0)
		{
			cout << "alpha(:," << t << ") = " << endl;
			alpha.col(t).print();
		}
	}
	// ---> Saving alpha for comparison purposes (in MATLAB)
	//alpha.save("Alpha1.mat", raw_ascii);

	// --> Testing forward variable computation with observation+duration
	cout << "Press [ENTER] to test forward variable using full observations";
	cin.get();

	myHSMM->resetRecursiveForwardVariable();
	mat tmpObs;
	// Loop for the observations saved in matrix "test"
	for(uint t=0 ; t<T ; t++)
	{
		tmpObs = test.cols(0,t);
		myHSMM->stepRecursiveForwardVariable(t,tmpObs);
	}
	alpha = myHSMM->getRecursiveForwardVariable();	// Getting forward variable
	alpha = alpha / repmat(sum(alpha,0),myHSMM->getNumSTATES(),1);	// Normalizing forward variable

	for(uint t=0 ; t<T ; t++)
	{
		if (t%10 == 0)
		{
			cout << "alpha(:," << t << ") = " << endl;
			alpha.col(t).print();
		}
	}
	// ---> Saving alpha for comparison purposes (in MATLAB)
	//alpha.save("Alpha2.mat", raw_ascii);

	// --> Testing forward variable computation with partial observation + duration
	cout << "Press [ENTER] to test forward variable using partial observations";
	cin.get();
	myHSMM->resetRecursiveForwardVariable();
	urowvec index;
	index << 0 ;
	// Loop for the observations saved in matrix "test"
	for(uint t=0 ; t<T ; t++)
	{
		tmpObs = test.cols(0,t);
		myHSMM->stepRecursiveForwardVariable(t,tmpObs,index);
	}
	alpha = myHSMM->getRecursiveForwardVariable();	// Getting forward variable
	alpha = alpha / repmat(sum(alpha,0),myHSMM->getNumSTATES(),1);	// Normalizing forward variable

	for(uint t=0 ; t<T ; t++)
	{
		if (t%10 == 0)
		{
			cout << "alpha(:," << t << ") = " << endl;
			alpha.col(t).print();
		}
	}
}

