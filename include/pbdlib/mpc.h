/**
Copyright (C) 2015, Martijn Zeestraten

This file is part of PbDLib (Programming-by-demonstration C++ Library).

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

/*! \file mpc.h
\brief lwr class
	This class is the implementation of a minimal intervention controller with 
	linear quadratic regulation, as described in:

\author Martijn Zeestraten
\bug No known bugs.
*/


#ifndef MPC_H_
#define MPC_H_

#include "armadillo"

using namespace arma;
namespace pbdlib {
class MPC
{
	public:

	// Constructor
	// A: [nbVar x nbVar] System matrix of discretized system
	// B: [nbVar x nbContr] Input Matrix of the discretized system
	// C: [nbVar x nbVarOutput] Output Matrix 
	// Np: Number of state predictions made during MPC
	// Nc: Number of control actions predicted during MPC
	MPC(mat& _Ad, mat& _Bd, mat& _C, uint _Np, uint _Nc);

	// Function that computes the control command
	// X:   [nbVar x 1]            Input vector size nbVar
	// MuQ: [nbVar*Np x 1]         Vector of state reference on the prediction horizon
	// Q:   [nbVar*Np x nbVar*Np]  Matrix representing the tracking costs
	// R:   [nbContr*Nc x nbContr*Nc] Matrix representing the control cost
	colvec&  computeControlCommand(colvec& X, colvec& Muq, mat& Q, mat& R);


	// Function to set the system matrices:
	// A: [nbVar x nbVar] System matrix of discretized system
	// B: [nbVar x nbContr] Input Matrix of the discretized system
	// C: [nbVar x nbVarOutput] Output Matrix 
	void setSystem(mat& _A, mat& _B, mat& _C);

	// Function to set the prediction horizon:
	void setPredictionHorizon(uint _Np);

	// Function to set the control horizon:
	void setControlHorizon(uint _Nc);

	// Getters:
 	uint getNp(){return Np;}
	uint getNc(){return Nc;}

	private:
	mat F;   // dynamic matrix
	mat Phi; // Input matrix
	mat CT;  // C M
	mat I;   // Output matrix (to select first control command)

	mat PhiQPhiR;// help variables for the least square solver
	vec PhiQMuqFx;
	uint nbVar; // Number of tracking variables
	uint nbContr; // Number of control variables
	colvec U;
	
	mat A, B, C; // System matrices

	void constructDynamicMatrices(mat& A, mat& B, mat& C,uint _Np, uint _Nc);
	void prepareForLeastSquares();

	protected:
	uint Np, Nc; // Prediction and control Horizon

};
} // end of pbdlib namespace

#endif
