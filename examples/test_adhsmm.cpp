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

/*! \file test_adhsmm.cpp
\brief Learning ADHSMM model
Example showing the computation of the forward variable of an ADHSMM for a given external input that modifies the duration probabilities of the model.

\author Leonel Rozo
\bug No known bugs.
*/

#include "pbdlib/adhsmm.h"
#include "armadillo"

using namespace arma;
using namespace pbdlib;
using namespace std;

int main()
{

	// ---> Create Object
	ADHSMM* myADHSMM;

	// ---> Load data from file
	string muPath 	  =	"data/adhsmm/ADHSMM_test_mu.txt";
	string sigmaPath  = "data/adhsmm/ADHSMM_test_sigma.txt";
	string priorsPath = "data/adhsmm/ADHSMM_test_priors.txt";
	string transPath  = "data/adhsmm/ADHSMM_test_trans.txt";
	string durPriors  = "data/adhsmm/ADHSMM_test_durPriors.txt";
	string durMu 	  = "data/adhsmm/ADHSMM_test_durMu.txt";
	string durSigma   = "data/adhsmm/ADHSMM_test_durSigma.txt";
	string testData   = "data/adhsmm/ADHSMM_test_testData.txt";

	// ---> Loading test data
	mat test;
	test.load(testData, raw_ascii);
	uint T = test.n_cols;

	// ---> Setting artificial external input
	//mat  u = zeros<mat>(1,T);
	mat u1 = zeros<mat>(1,45);
	mat u2 = ones<mat>(1,45);
	mat u3 = zeros<mat>(1,110);
	mat u = arma::join_horiz(u1,u2);
	u = arma::join_horiz(u,u3);
	urowvec varIn, varOut;
	varIn  << 0; // The input variables for regression of the duration probability. (external input)
	varOut << 1; // The output variables for regression of the duration probability. (duration)

	myADHSMM = new ADHSMM(priorsPath,muPath,sigmaPath,transPath,durPriors,durMu,durSigma,T);


	// ---> Display Data
	cout << "--- ADHSMM components: " << endl;
	cout << "nStates " << myADHSMM->getNumSTATES() <<endl;
	for (uint i=0;i<myADHSMM->getNumSTATES();i++)
	{
		cout << "Priors" << endl;
		cout << myADHSMM->getPRIORS(i) << endl;
		cout << "Mu" << endl;
		cout << myADHSMM->getMU(i) << endl;
		cout << "Sigma" << endl;
		cout << myADHSMM->getSIGMA(i) << endl;

		// Duration GMM
		cout << "Duration prob. [" << i << "]" << endl;
		myADHSMM->getDurationGMMs(i).getPRIORS().print("durPriors");
		for(uint j=0 ; j<myADHSMM->getNumSTATESDUR() ; j++)
		{
			myADHSMM->getDurationGMMs(i).getMU(j).print("durMu:");
			myADHSMM->getDurationGMMs(i).getSIGMA(j).print("durSigma:");
		}
	}
	cout << "Transition matrix = " << endl;
	myADHSMM->getTRANSITION().print();


	// ---> Testing duration-based forward variable
	cout << "Press [ENTER] to test duration-based forward variable...";
	cin.get();

	// Loop for T samples
	for(uint t=0 ; t<T ; t++)
		myADHSMM->stepRecursiveForwardVariable(u.col(t), varIn, varOut, t);

	// Getting alpha variable and normalizing it
	mat alpha = myADHSMM->getRecursiveForwardVariable();	// Getting forward variable
	alpha = alpha / repmat(sum(alpha,0),myADHSMM->getNumSTATES(),1);	// Normalizing forward variable

	// Printing results
	for(uint t=0 ; t<T ; t++)
	{
		if (t%4 == 0)
		{
			cout << "alpha(:," << t << ") = " << endl;
			alpha.col(t).print();
		}
		//cin.get();
	}
	// ---> Saving alpha for comparison purposes (in MATLAB)
	//alpha.save("Alpha1.mat", raw_ascii);


	// ---> Testing forward variable computation with observation+duration
	cout << "Press [ENTER] to test forward variable using full observations";
	cin.get();

	myADHSMM->resetRecursiveForwardVariable();
	mat tmpObs;
	// Loop for the observations saved in matrix "test"
	tmpObs = test.col(0);
	for(uint t=0 ; t<T ; t++)
	{
		myADHSMM->stepRecursiveForwardVariable(u.col(t), varIn, varOut, t,tmpObs);

		// Simulating the effects of the external input on the system.
		// If the external input is "1", then the systems entirely stops
		if(!u.col(t)(0))
			tmpObs = arma::join_horiz(tmpObs,test.col(min(t+1,T-1)));
		else
			tmpObs = arma::join_horiz(tmpObs,tmpObs.col(tmpObs.n_cols-1));
	}
	alpha = myADHSMM->getRecursiveForwardVariable();	// Getting forward variable
	alpha = alpha / repmat(sum(alpha,0),myADHSMM->getNumSTATES(),1);	// Normalizing forward variable

	for(uint t=0 ; t<T ; t++)
	{
		if (t%4 == 0)
		{
			cout << "alpha(:," << t << ") = " << endl;
			alpha.col(t).print();
		}
		//cin.get();
	}

	// ---> Saving alpha for comparison purposes (in MATLAB)
	//alpha.save("Alpha2.mat", raw_ascii);
}
