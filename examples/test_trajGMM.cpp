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

/*! \file test_trajGMM.cpp
\brief Learning trajectory GMM model
Learning a trajectory GMM model from demonstrations saved in the file data_xy.txt

\author Leonel Rozo
\bug No known bugs.
*/

#include "pbdlib/trajgmm.h"
#include <sstream>

#define nbStates 5
#define nbVarPos 2
#define nbDeriv 3
#define dt 0.01
#define nbData 100
#define nbSamples 4

using namespace pbdlib;

int main(int argc, char **argv)
{
	/**********************************************************************************************************/
	/*													Test TrajGMM learning from demos saved in txt files														*/
	/**********************************************************************************************************/
	// Demonstrations variables
	std::vector<Demonstration> demos;
	Demonstration demo =  Demonstration(nbVarPos,nbData*nbSamples);
	std::vector<std::string> vars;
	vars.push_back("x1");	vars.push_back("x2");
	vars.push_back("v1");	vars.push_back("v2");
	vars.push_back("a1");	vars.push_back("a2");

	// Trajectory model initialization
	TrajGMM *myTrajGMM1;
	myTrajGMM1 = new TrajGMM(nbStates, nbVarPos, nbDeriv, dt);
	myTrajGMM1->setVARSNames(vars);
	myTrajGMM1->constructPHI(nbData,nbSamples);

	// Files variables
	mat posData;
	std::string posDataPath;
	std::string commonPath =	"data/trajgmm/TrajGMM_posData_Sshape";
	std::stringstream convert; // stringstream used for the conversion

	// Loading and storing demos
	for(uint i = 1 ; i <= nbSamples ; i++)
	{
		// Loading i-th demo
		convert.str("");
		convert << i;
		posDataPath = commonPath + convert.str() + ".txt";
		posData.load(posDataPath);

		// Reshaping demonstration data for computing static and dynamic features
		posData.reshape(nbVarPos*nbData,1);
		mat trainingData = conv_to<mat>::from(myTrajGMM1->getPHI1()) * (posData * 1E2); //1E2 is used to avoid numerical computation problem
		trainingData.reshape(nbVarPos*myTrajGMM1->getNumDERIV(), nbData);

		// Storing transformed demo
		demo.getDatapoints().setData(trainingData);
		demo.getDatapoints().setVarNames(vars);
		demos.push_back(demo);
	}
	// Setting training data
	myTrajGMM1->setDEMONSTRATIONS(demos);

	// Learning and printing results
	cout << "\n Number of EM iterations: " << myTrajGMM1->EM_learn(1E-4);
	for (int i=0;i<nbStates;i++){
		cout << "\n Mu_" << i << " = \n" << myTrajGMM1->getMU(i);
		cout << "\n Sigma_" << i << " = \n" << myTrajGMM1->getSIGMA(i);
	}
	cout << "\nPress enter to test TrajGMM using a model learned in MATLAB..." << endl; getchar();


	/**********************************************************************************************************/
	/*													Test TrajGMM reproduction using a pre-trained model														*/
	/**********************************************************************************************************/
	// ---> Load data from file (TrajGMM trained with Sshape AMARSI data)
	std::string muPath 		=	"data/trajgmm/TrajGMM_test_mu.txt";
	std::string sigmaPath = "data/trajgmm/TrajGMM_test_sigma.txt";
	std::string priorsPath = "data/trajgmm/TrajGMM_test_priors.txt";
	std::string varNamesPath = "data/trajgmm/TrajGMM_test_varnames.txt";
	std::string modelInfoPath = "data/trajgmm/TrajGMM_test_info.txt";
	std::string testGammaPath = "data/trajgmm/TrajGMM_test_Gamma.txt";

	TrajGMM *myTrajGMM2;
	myTrajGMM2 = new TrajGMM(priorsPath,muPath,sigmaPath,varNamesPath,modelInfoPath);

	// Loading a precomputed state sequence (from MATLAB)
	mat q;
	q.load(testGammaPath);
	q = q - 1;	// Changing state format to C++

	// Computing optimal path from TrajGMM given the state sequence
	mat data1 = myTrajGMM2->leastSquaresOptimalData(arma::conv_to<colvec>::from(q));
	// Computing optimal path and related covariance stored in a vector of Gaussian distributions
	std::vector<GaussianDistribution> optimalDist;
	optimalDist = myTrajGMM2->leastSquaresOptimalProb(arma::conv_to<colvec>::from(q));

	// Printing and saving data for possible comparisons
	mat tmpMu = zeros<mat>(myTrajGMM2->getNumVARSPOS(), q.n_elem);
	mat tmpSigma = zeros<mat>(myTrajGMM2->getNumVARSPOS(), myTrajGMM2->getNumVARSPOS()*q.n_elem);
	for(uint t = 0 ; t < q.n_elem ; t++)
	{
		data1.col(t).print("Mu_1 = ");

		tmpMu.col(t) = optimalDist.at(t).getMU();
		tmpSigma.cols(t*myTrajGMM2->getNumVARSPOS(),(t+1)*myTrajGMM2->getNumVARSPOS()-1) = optimalDist.at(t).getSIGMA();

		tmpMu.col(t).print("Mu_2 = ");
		tmpSigma.cols(t*myTrajGMM2->getNumVARSPOS(),(t+1)*myTrajGMM2->getNumVARSPOS()-1).print("Sigma = ");
		//std::cin.get();
	}
	// ---> Saving data for comparison purposes (in MATLAB)
	//data1.save("data1TrajGMM.mat", raw_ascii);
	//tmpMu.save("data2muTrajGMM.mat", raw_ascii);
	//tmpSigma.save("data2sigmaTrajGMM.mat", raw_ascii);

	return 0;
}
