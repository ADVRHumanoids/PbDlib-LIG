
/**
Copyright (C) 2015, Martijn Zeestraten

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

/*! \file test_mpc_
\brief Testing mpc 

\author Martijn Zeestraten 
\bug No known bugs.
*/

#include <iostream>
#include <sstream>
#include "pbdlib/mpc.h"
#include "armadillo"
#include <vector>


using namespace pbdlib;
using namespace arma;
using namespace std;

int main(int argc, char **argv)
{
	// Settings
	int Np      = 20;     // Prediction horizon
	int Nc      = 10;     // Control horizon
	int nbVarPos= 1;      // Use 1D system
	int nbDeriv = 2;      // Number of derivatives (position and velocity)
	float dt    = 0.1;    // time step
	float alpha = 0.05;   // control cost factor
	int nbData  = 250;    // Number of data points
	colvec Xinit;         // Initial state of the system 
	Xinit << 0 << 0;
	colvec Target;        // Final desired final value for the system
	Target << 10 << 0; 


	//-----------  Create System Matrices
	// Assistive Matrices for System construction:	
	
	cout << "Creating system matrices:"<< endl;
	mat A,B,C,Ad,Bd;

	// A matrix:
	A << 0 << 1 << endr << 0 << 0;

	// B Matrix:
	B << 0 << endr <<1<< endr;

	// C Matrix:
	C << 1 << 0; 

	// Discretize system matrices (Euler discretization):
	Ad = A*dt+eye<mat>(2,2);
	Bd = B*dt;

	cout <<"Ad: "<< endl << Ad << endl;
	cout <<"Bd: "<< endl << Bd << endl;
	
	// ----------- Create Target:
	// We create a simple target, to let the system go to zero
	colvec MuQ= repmat(Target,Np,1); // Np predictions of dimension 2 (position and velocity)
	cout << "Size MuQ: " << MuQ.n_rows<< "x"  << MuQ.n_cols << endl;

	// Create tracking cost matrix:
	mat QTmp = zeros(2,2);
	QTmp(0,0) = 1;
	mat Q = kron(QTmp,eye<mat>(Np,Np));
	cout << "Size Q: "  << Q.n_rows << "x" << Q.n_cols << endl;

	mat R = eye(nbVarPos*Nc,nbVarPos*Nc)*(alpha*alpha);
	cout << "Size R: "  << R.n_rows << "x" << R.n_cols << endl;

	// ----- Perform MPC:
	cout << "Performing MPC" << endl;
	
	// Create MPC object:
	MPC myMPC(Ad,Bd,C,Np,Nc);

	//Output data
	colvec Xtmp = Xinit;
	colvec U;
	for (int t=0;t<nbData;t++){
		// Compute control command:
		U = myMPC.computeControlCommand(Xtmp,MuQ,Q,R);
		// Perform iteration:
		Xtmp = Ad*Xtmp+Bd*U;

		// Display Result:
		if (t%10==0 || t==0 || t==nbData)
		{
			std::cout << "X("<< t << "): [" << Xtmp(0) << ", "<< Xtmp(1) <<
		   		"]\t U(" << t << "): " << U(0) << endl;
		}
	}

	return 0;
}
