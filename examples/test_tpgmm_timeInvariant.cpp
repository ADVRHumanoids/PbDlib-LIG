/*
 * Copyright (c) 2014
 * - João Silvério @ joao[dot]silverio[at]iit[dot]it
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
 */

#include "pbdlib/tpgmm.h"
#include "pbdlib/gmm.h" 
#include "pbdlib/gmr.h"
#include "pbdlib/datapoints.h"
#include "pbdlib/taskparameters.h"
#include <sstream>
#include <math.h>
#include <iostream>

#ifdef __APPLE__
	#include <sys/time.h>
#endif

using namespace std;
using namespace pbdlib;

int main(int argc, char **argv)
{
	// Filenames that contain model parameters
	std::string fParamsTPGMM("data/pgmm/timeInvModel/pgmm_timeInv.txt");
	std::string fVarsTPGMM  ("data/pgmm/timeInvModel/vars_timeInv.txt");
	std::string fPriorsTPGMM("data/pgmm/timeInvModel/priors_timeInv.txt");
	std::string fMuTPGMM    ("data/pgmm/timeInvModel/Zmu_timeInv");		// Initialize without file extension -> "_P#" will be added depending on nParams
	std::string fSigmaTPGMM ("data/pgmm/timeInvModel/Zsigma_timeInv");	// Initialize without file extension -> "_P#" will be added depending on nParams

	// Loading files and setting pgmm variables
	vec ModelParams;// The model parameters are loaded from the PGMM model file
	ModelParams.load(fParamsTPGMM);
	uint nVars = ModelParams[0];
	uint nFrames = ModelParams[1];
	uint nStates = ModelParams[2];

	// According to GMM class, the name of the variables should also come in a separate .txt file.
	// I'll use a toy file for now, during object initialization, and keep it like below for the PGMM object construction.
	std::vector<std::string> varNames;
	varNames.push_back("x_in");
	varNames.push_back("y_in");
	varNames.push_back("z_in");
	varNames.push_back("x_out");
	varNames.push_back("y_out");
	varNames.push_back("z_out");

	TPGMM* tpgmm;
	tpgmm = new TPGMM(nVars, nStates, nFrames);
	tpgmm->loadTPGMMfromMATLAB(fPriorsTPGMM, fVarsTPGMM, fMuTPGMM, fSigmaTPGMM);

	// Reproduction variables (GMR)
	GMR *gmr = new GMR();
	Datapoints *repro = new Datapoints(3,1);				// Only one datapoint passed at each time step

	std::vector<std::string> reproVarNames;
	reproVarNames.push_back("x_in");		// time is the query point
	reproVarNames.push_back("y_in");
	reproVarNames.push_back("z_in");
	repro->setVarNames(reproVarNames);

	
	// Defining "artificial" parameters (New Implementation TPGMM)
	TaskParameters TPs(nVars, nFrames);
	TaskParameter auxTP;

	// First Frame:
	auxTP.A = eye<mat>(nVars,nVars);
	auxTP.b = zeros(nVars,1);
	TPs.setTaskParameters(0,auxTP);

	// Second and third frame
	for (uint i=0;i<2;i++)
	{
		auxTP.b.at(3) = tpgmm->getGMMS(i).getMU(0)[0];		// the object positions are set in the .txt file
		auxTP.b.at(4) = tpgmm->getGMMS(i).getMU(0)[1];
		auxTP.b.at(5) = tpgmm->getGMMS(i).getMU(0)[2];
		TPs.setTaskParameters(i,auxTP);
	}

	/* OLD IMPLEMENTATION PGMM class, is left for comparison:
	std::vector<mat> A_tmp;
	mat auxA;
	std::vector<colvec> b_tmp;
	colvec auxb;

	auxA = eye<mat>(nVars, nVars);
	A_tmp.push_back(auxA);
	A_tmp.push_back(auxA);	// P=2 frames

	auxb = zeros(nVars, 1);
	tpgmm->getGMMS(0).getMU(0).print();

	auxb.at(3) = tpgmm->getGMMS(0).getMU(0)[0];		// the object positions are set in the .txt file
	auxb.at(4) = tpgmm->getGMMS(0).getMU(0)[1];
	auxb.at(5) = tpgmm->getGMMS(0).getMU(0)[2];
	b_tmp.push_back(auxb);
	auxb.at(3) = tpgmm->getGMMS(1).getMU(0)[0];
	auxb.at(4) = tpgmm->getGMMS(1).getMU(0)[1];
	auxb.at(5) = tpgmm->getGMMS(1).getMU(0)[2];
	b_tmp.push_back(auxb);
	*/

	/* for tracking the time */
	struct timeval tim;
	double t0, t;

	ofstream myfile;
	myfile.open ("output_log.txt");

	mat repDataPt;
	repDataPt = zeros(3, 1);

	for(int tn = -400 ; tn <= 400 ; tn++){

		// Record time before computing GMM for new {A,b} and GMR
		gettimeofday(&tim, NULL);
		t0=tim.tv_sec*1000+(tim.tv_usec/1000.0); // tv_usec comes in u-sec -> divide by 1M for sec.; 1k for m-sec

		// Computing the resulting gmm given the set of parameters {A,b}
		GMM_Model* gmm;
		gmm = tpgmm->getTransformedGMM(TPs, ORTHONORMAL); // New Implementation of TPGMM
		//gmm = tpgmm->getTransformedProdGMM(A_tmp,b_tmp, ORTHONORMAL); // Old implementatin of pgmm
		gmm->setVARSNames(varNames);
		// Shaving off the "aux" frame - not needed when using a single frame
		// A_tmp.erase(A_tmp.begin()+1);
		// b_tmp.erase(b_tmp.begin()+1);

		
		// Printing model components
		cout << "Resulting GMM given the set of parameters 'A' and 'b'" << endl;
		for(uint i = 0 ; i < nStates ; i++){
			cout << "State #" << i << ":"<< endl;
			gmm->getCOMPONENTS(i).getMU().print("Mu = ");
			gmm->getCOMPONENTS(i).getSIGMA().print("Sigma = ");
		}

		// Computing GMR
		gmr->setGMMModel(gmm);				//<-problem
		repDataPt[1] = tn*0.001;
//		repDataPt.at(0) = 0.005;
//		repDataPt.at(1) = 0.4950;
//		repDataPt.at(2) = 0.2;
		repDataPt.print("x = ");
		repro->setData(repDataPt);
		GMM_Model* gmmOut = gmr->regression(repro);
		
		gmmOut->getMU(0).print("MuOut = ");

		// Calculating elapsed time from the beginning of GMM recomputation to GMR
		gettimeofday(&tim, NULL);
		t=tim.tv_sec*1000+(tim.tv_usec/1000.0);
		cout << "GMR computational time: " << t-t0 << " ms" << endl;

		cout << "Please press [ENTER] to continue or [CTRL+C] to finish" << endl;
		char key;
		std::cin.ignore(1);
		
	}
	myfile.close();
	return 0;
}
