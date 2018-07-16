/*
* Copyright (c) 2014
* - Davide De Tommaso @ dtmdvd[at]gmail[dot]com
* - Milad Malekzadeh @ milad[dot]malekzadeh[at]gmail[dot]com
* - Leonel Rozo @ leonel[dot]rozo[at]iit[dot]it
* - Sylvain Calinon @ sylvain[dot]calinon[at]gmail[dot]com
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY <copyright holder> ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*
*
*	----------------- NOTE ---------------------------
*	This example is a modified version of the 'old'  example pgmm_test. This
*	new example uses the tpgmm class instead of pgmm. Changes that need to be 
*	made when swithcing from pgmm to tpgmm are:
*	- Delaration PGMM_Model -> TPGMM
*	- Product of Gaussian: prodGMM -> transformed
*
*
*
*
*/

#include "pbdlib/tpgmm.h"
#include "pbdlib/gmm.h" 
#include "pbdlib/gmr.h"
#include "pbdlib/datapoints.h"
#include "pbdlib/tpdemonstration.h" 
#include "pbdlib/taskparameters.h"
#include <sstream>

using namespace pbdlib;

int main(int argc, char **argv)
{
	// Model variables
	uint nFrames = 2;
	uint nStates = 3;
	uint nVars = 3; //t,x,y
	uint nDemos = 4; //Number of Demonstrations
	uint nData = 200; //Number of datapoints in a demonstration
	//
	// time is input (t) and the attractor position (2D vector [x y]) is the output.
	std::vector<std::string> varNames;
	varNames.push_back("t");	
	varNames.push_back("x");
	varNames.push_back("y");

	//..................................................................................
	//......... loading the demonstrations and task parameters from the txt files ......

	std::vector<TPDemonstration> demos; // Vector to hold demonstrations
	TPDemonstration tmpDemo; // Temp variable to load individual demonstrations

	char Datafilename[256];
	char TPfilename[256];
	cout << "Loading the demonstrations and the task parameters ..." << endl;
	for (uint m=0; m<nDemos; m++){   //Loading demos in the loop
		// Form filenames:
		sprintf(Datafilename, "../data/pgmm/Data0%d.txt",m+1);
		sprintf(TPfilename, "../data/pgmm/Param0%d.txt",m+1);

		// Load files in demo:
		// note we need to transpose the data (last argument)
		// because the data in the file a row for each data point
		// and a column for each variable. The agreement is
		// the other way around
		tmpDemo.loadFromFiles(Datafilename,TPfilename, true);

		// push back demonstration:
		demos.push_back(tmpDemo);
	}
	
	cout << "Demonstrations Loaded succesfully." << endl;
	cout<<"Press any key to continue..."<<endl; getchar();

	//..................................................................................
	//......... Learning the PGMM model parameters using the EM algorithm ......
	cout<<"Learning model parameters using EM..." << endl;

	TPGMM *tpgmm;
	tpgmm = new TPGMM(demos, nStates);

	cout << "\n Number of EM iterations: " << tpgmm->estimateTensorGMM() << endl;

	// saving the learned model parameters
	mat MuTmp = zeros(nVars, nStates*nFrames);
	mat SigmaTmp = zeros(nVars, nVars*nStates*nFrames);
	mat Priors = zeros(1,nStates);
	std::vector<GMM_Model> GMMSresult;
	GMMSresult = tpgmm->getGMMS();
	for(uint m=0; m<nFrames; m++){
		for (uint i=0;i<nStates; i++){
			MuTmp.col(m*nStates+i) = GMMSresult[m].getMU(i);
			SigmaTmp.cols(m*nVars*nStates+i*nVars, m*nVars*nStates+i*nVars+nVars-1) = GMMSresult[m].getSIGMA(i);
		}
	}
	sprintf(Datafilename, "../../data/pgmm/Mu0%d.txt",2);
	MuTmp.save(Datafilename, raw_ascii);
	sprintf(Datafilename, "../../data/pgmm/Sigma0%d.txt",2);
	SigmaTmp.save(Datafilename, raw_ascii);
	sprintf(Datafilename, "../../data/pgmm/Priors0%d.txt",2);
	Priors = tpgmm->getPRIORS();
	Priors.save(Datafilename, raw_ascii);
	cout<<"the learned model parameters are written into text files succesfully."<<endl; cout<<"Press any key to continue..."<<endl; getchar();

	
	// Reproduction variables (GMR)
	GMR	*gmr = new GMR();
	Datapoints	*repro = new Datapoints(1,1);  // Only one datapoint passed at each time step
	mat	reproData, auxYgmr;
	std::vector<std::string> reproVarNames;
	reproVarNames.push_back("t");              // time is the query point
	repro->setVarNames(reproVarNames);

	// loading new parameters for the reproduction from the text file
	// In this case we have considered that the parameters are the same all along the reproduction. They can be varying.
	TaskParameters ReproTPs(nVars,nFrames);
	sprintf(TPfilename, "../../data/pgmm/ParamRepro0%d.txt",1);
	ReproTPs.loadFromFile(TPfilename);

	GMM_Model* gmm;
	for(uint tn=1; tn<=nData; tn++){
		// Computing the resulting GMM given the set of parameters {A,b}
		gmm = tpgmm->getTransformedGMM(ReproTPs);
		gmm->setVARSNames(varNames);
		
		// Printing model components
		cout << "Resulting GMM given the set of parameters 'A' and 'b'" << endl;
		cout << "nbStates : " << gmm->getNumSTATES() << endl;
		cout << "nbVar    : " << gmm->getNumVARS() << endl;
		cout << "VarNames " ;
		for (uint i=0;i<gmm->getNumVARS();i++)
			cout << varNames[i] ;
		cout << endl;
	
		cout << "Priors   : "  << gmm->getPRIORS() << endl;
		for(uint i=0; i<nStates; i++){
			cout << "State #" << i << ":"<< endl;
			gmm->getMU(i).print("Mu = ");
			gmm->getSIGMA(i).print("Sigma = ");
		}

		// Computing GMR
		gmr->setGMMModel(gmm);
		mat repDataPt;
		repDataPt << (tn * 0.01);	// time step as input
		repDataPt.print("tn = ");
		repro->setData(repDataPt);
		GMM_Model* gmmOut = gmr->regression(repro);
		gmmOut->getMU(0).print("MuOut = ");

		cout << "Please press [ENTER] to continue or [CTRL+C] to finish" << endl;
		char key;
		std::cin.ignore(1);
	}
	return 0;
}




