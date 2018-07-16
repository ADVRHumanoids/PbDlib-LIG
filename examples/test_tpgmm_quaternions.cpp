/**
Copyright (C) 2015, João Silvério

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

/*! \file test_tpgmm_quaternions.cpp
\brief An example of TP-GMM using quaternions.

This source code is given for free! In exchange, I would be grateful if you cite
the following reference in any academic publication that uses this code or part of it:

@inproceedings{Silverio15IROS,
	author="Silv\'erio, J. and Rozo, L. and Calinon, S. and Caldwell, D. G.",
	title="Learning bimanual end-effector poses from demonstrations using task-parameterized dynamical systems",
	booktitle='Proc. {IEEE/RSJ} Intl Conf. on Intelligent Robots and Systems ({IROS})",
	year="2015",
	month="Sept.-Oct.",
	address="Hamburg, Germany",
	pages=""
}

\author João Silvério
\bug No known bugs.
*/

#include "pbdlib/tpgmm.h"
#include "pbdlib/gmm.h" 
#include "pbdlib/gmr.h"
#include "pbdlib/datapoints.h"
#include "pbdlib/taskparameters.h"
#include "pbdlib/quaternion.h"
#include <sstream>
#include <math.h>
#include <iostream>

using namespace std;
using namespace pbdlib;

int main(int argc, char **argv)
{
	// Filenames that contain model parameters
	std::string fParamsTPGMM("data/pgmm/quaternionModel/tpgmm_quaternions.txt");
	std::string fVarsTPGMM  ("data/pgmm/quaternionModel/vars_quaternions.txt");
	std::string fPriorsTPGMM("data/pgmm/quaternionModel/priors_quaternions.txt");
	std::string fMuTPGMM    ("data/pgmm/quaternionModel/Zmu_quaternions"); // Initialize without file extension -> "_P#" will be added depending on nParams
	std::string fSigmaTPGMM ("data/pgmm/quaternionModel/Zsigma_quaternions");	// Initialize without file extension -> "_P#" will be added depending on nParams

	// Loading files and setting tpgmm variables
	vec ModelParams; // The model parameters are loaded from the TPGMM model file
	ModelParams.load(fParamsTPGMM);
	uint nVars = ModelParams[0];
	uint nFrames = ModelParams[1];
	uint nStates = ModelParams[2];

	// According to GMM class, the name of the variables should also come in a separate .txt file.
	// I'll use a toy file for now, during object initialization, and keep it like below for the PGMM object construction.
	std::vector<std::string> varNames;
	varNames.push_back("t"); // time driven regression
	varNames.push_back("q0_out"); // quaternion variables - the order depends on the convention with which the quaternion is initialized
	varNames.push_back("q1_out"); // -> SCALAR_VEC means that q=[u vec], aka [w x y z] (the default is VEC_SCALAR)
	varNames.push_back("q2_out");
	varNames.push_back("q3_out");

	TPGMM* tpgmm;
	tpgmm = new TPGMM(nVars, nStates, nFrames);
	tpgmm->loadTPGMMfromMATLAB(fPriorsTPGMM, fVarsTPGMM, fMuTPGMM, fSigmaTPGMM);

	// Reproduction variables (GMR)
	GMR *gmr = new GMR();

	std::vector<std::string> inputVarNames;
	inputVarNames.push_back("t");			 // time-driven movement
	Datapoints *repro = new Datapoints(1,1); // Only one datapoint passed at each time step (the example uses moving frames, so regression can't be done in batch)
	repro->setVarNames(inputVarNames);

	mat repDataPt;				  // stores the value of the regression variable (used by 'repro' object)
	repDataPt = zeros(1, 1);	  // Keeping it as matrix to remain generic since time-driven movements => size(repDataPt)=1
	Quaternion q_out(SCALAR_VEC); // a quaternion object to store the output of GMR
								  //    -> quatOrder = SCALAR_VEC means that q=[u vec], aka [w x y z] (the default is VEC_SCALAR)

	// Declare Task Parameter objects
	TaskParameters TPs(nVars, nFrames);
	TaskParameter auxTP;

	// The frame will encode the orientation of the end-effector with respect to an object
	// -> the object's orientation will be initialized as the identity but it might change during the task.
	Quaternion q_frame(SCALAR_VEC);	// default constructor -> q_frame=[1 0 0 0]

	// Use these angle-axis variables to define new orientations of the frame during reproductions
	double angle = M_PI/5;
	vec3 axis;
	axis << 0 << 0 << 1 << endr;	// [pi/2 , [0 0 1]] -> rotation of 90deg around z-axis

	// First Frame	(this example considers one single frame)
	auxTP.A = eye<mat>(nVars,nVars);
	auxTP.A.submat(1,1,4,4) = q_frame.matrix();
	auxTP.b = zeros(nVars,1);
	TPs.setTaskParameters(0,auxTP);

	for(int tn = 0 ; tn <= 80 ; tn++){	// the model encodes an 8s movement

		// Adjust angle, axis according to new frame orientations
		q_frame = Quaternion(angle,axis,SCALAR_VEC);

		//Update task-parameters
		auxTP.A.submat(1,1,4,4) = q_frame.matrix();
		TPs.setTaskParameters(0,auxTP);

		// Computing the resulting gmm given the set of parameters {A,b}
		GMM_Model* gmm;
		gmm = tpgmm->getTransformedGMM(TPs, ORTHONORMAL); // New Implementation of TPGMM
		gmm->setVARSNames(varNames);

		// Printing model components
		cout << "Resulting GMM given the set of parameters 'A' and 'b'" << endl;
		for(uint i = 0 ; i < nStates ; i++){
			cout << "State #" << i << ":"<< endl;
			gmm->getCOMPONENTS(i).getMU().print("Mu = ");
			gmm->getCOMPONENTS(i).getSIGMA().print("Sigma = ");
		}

		// Computing GMR
		gmr->setGMMModel(gmm);
		repDataPt[0] = tn*0.1;
		repDataPt.print("t = ");
		repro->setData(repDataPt);
		GMM_Model* gmmOut = gmr->regression(repro);
		q_out = Quaternion(gmmOut->getMU(0),SCALAR_VEC);
		q_out.normalize();

		cout << "\nFrame orientation: q_frame = " ;
		q_frame.print();

		cout << "\nNew end-effector orientation: q_new = " ;
		q_out.print();

		cout << "\nEnd-effector orientation with respect to frame: q_rel = " ;
		(q_out*q_frame.conjugate()).print();

		cout << "\nPress [ENTER] to continue or [CTRL+C] to finish" << endl;
		std::cin.ignore(1);

	}
	return 0;
}
