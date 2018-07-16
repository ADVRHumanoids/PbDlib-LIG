// May 2015, Martijn Zeestraten


#ifndef PTRAJMPC_H_
#define PTRAJMPC_H_


#include "armadillo"
#include "mpc.h"
#include "pbdlib/mvn.h"
#include "pbdlib/hsmm.h"

using namespace arma;
using namespace std;

namespace pbdlib
{
class TrajMPC 
{
	private: // Variables:
		// variables:
		uint nbVarPos; // Number of position variables
		uint nbDeriv;  // Number of derivatives used
		uint nbVar;    // Total number of variables
		urowvec in;    // in values

		MPC* myMPC;
		HSMM* hsmm;

		// System Matrices:
		mat A, B, C; // System Matrices
		mat Ad, Bd;  // Discrete system matrices

		// For state sequence:
		mat AlphaPred; // Matrix to hold alpha predictions
		mat maxAlpha;  // Matrix to hold maximum predicted values:
		mat tmpMax;
		mat tmpAlphaCol;
		uvec q;        // State sequence

		mat Q;         // Tracking cost Matrix
		mat Lambda;    // Precision Matrices 
		colvec Muq;    // Reference 
		mat R;         // Control cost matrix

		colvec U;
	public:
		// Constructor 1:
		// hsmm     : Hiddem semi-markov Model
		// nbVarPos : Number of position variables
		// nbDeriv  : Number of derivative considered (1:Position, 2: Position+Velocity 3: ...)
		// Np       : Prediction horizon
		// Nc       : Control horizon 
		// Alpha    : control cost alpha
		TrajMPC(HSMM* _hsmm, double dt, uint nbVarPos, uint nbDeriv, uint _Np, uint _Nc,double _alpha); 

		// Constructor 2:
		// hsmm  : Hiddem semi-markov Model
		// A     : [nbVar x nbVar] System matrix of disretized system
		// B     : [nbVar x nbVarIn] Input matrix of discretized system 
		// C     : [nbVarOut x nbVar] Output matrix
		TrajMPC(HSMM* _hsmm, mat& Ad, mat& Bd, mat& C, uint _Np, uint _Nc, double _alpha);

		// Function that computes the current control command
		colvec& computeControlCommand(colvec& X);

		void reset(){hsmm->resetForwardVariable();} // reset HSMM
		mat getSystemDynamics(){return Ad;}
		mat getInputDynamics(){return Bd;}
		uvec getq(){return q;}
		colvec getCurrentAttractor();
		uint getNumVARS(){return hsmm->getNumVARS();}
		uint getNumVARSPos(){return nbVarPos;}


	private:
		// functions:
		void updateReference(uvec& q);
		void updateCostMatrix(double alpha);
		void createSystemMatrices(double dt, uint nbVarPos, uint nbDeriv);
		void createPrecisionMatrices(HSMM* _hsmm);
		void allocateReference(uint nbVar, uint nbVarPos, uint Np, uint Nc);
		void assignHSMM(HSMM* hsmm, uint nbVarPos, uint nbDeriv);


};
} // End namespace
#endif

