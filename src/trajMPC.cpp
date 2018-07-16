// Martijn Zeestraten, June 2015
//
#include "pbdlib/trajMPC.h"

namespace pbdlib
{

TrajMPC::TrajMPC(HSMM* _hsmm, double _dt, uint _nbVarPos, uint _nbDeriv, uint _Np, uint _Nc, double _alpha) 
{

	nbVarPos = _nbVarPos;           // Assume fully actuated system
	nbDeriv  = _nbDeriv;   
	nbVar    = _nbVarPos*_nbDeriv;  // Total number of variables 

	// Check HSMM size and copy pointer:
	assignHSMM(_hsmm,_nbVarPos,_nbDeriv);

	// Invert all covariance matrices:
	createPrecisionMatrices( _hsmm);

	// Create system matrices:
	createSystemMatrices(_dt, _nbVarPos, _nbDeriv);

	// Create MPC object
	myMPC = new MPC(Ad,Bd,C,_Np,_Nc);

	// Initialize R and Q matrices for appropriate size
	allocateReference(nbVar, nbVarPos, _Np, _Nc);

	// Update Cost Matrix
	updateCostMatrix(_alpha);

	// Create vector that selects the variables used for observation probability.
	// We assume that we only want to condition on the position information
	in.set_size(nbVarPos);
	for (uint i=0;i<nbVarPos;i++)
		in(i) = i;

	// Allocate memory for state sequence generation
	maxAlpha.set_size(1,_Np);
	tmpMax.set_size(_hsmm->getNumVARS(),1);
	tmpAlphaCol.set_size(_hsmm->getNumVARS(),1);

	AlphaPred.set_size(_hsmm->getNumSTATES(),_Np);
}

TrajMPC::TrajMPC(HSMM* _hsmm,mat& _Ad, mat& _Bd, mat& _C, uint _Np, uint _Nc, double _alpha)
{
	// Get variable sizes:
	nbVarPos = _Bd.n_cols;      // Assume fully actuated system
	nbVar    = _Ad.n_rows;      // Total number of variables 
	nbDeriv  = nbVar/nbVarPos;  // 

	// Check hsmm vars
	assignHSMM(_hsmm,nbVarPos, nbDeriv);

	// Invert all covariance matrices:
	createPrecisionMatrices( _hsmm);

	// Create MPC object
	myMPC = new MPC(_Ad,_Bd,_C,_Np,_Nc);

	// Initialize R and Q matrices for appropriate size
	allocateReference(nbVar, nbVarPos, _Np, _Nc);
	
	// Update Cost Matrix
	updateCostMatrix(_alpha);

	// Create vector that selects the variables used for observation probability.
	// We assume that we only want to condition on the position information
	in.set_size(nbVarPos);
	for (uint i=0;i<nbVarPos;i++)
		in(i) = i;

	// Allocate memory for state sequence generation
	maxAlpha.set_size(1,_Np);
	tmpMax.set_size(_hsmm->getNumVARS(),1);
	tmpAlphaCol.set_size(_hsmm->getNumVARS(),1);

	AlphaPred.set_size(_hsmm->getNumSTATES(),_Np);
}

colvec TrajMPC::getCurrentAttractor()
{

	// Make prediction of first most likeli state:
	mat AlphaPred = hsmm->predictForwardVariable(1);

	// Find corresponding corresponding indices (i.e. argmax())
	maxAlpha = max(AlphaPred,0);
	tmpAlphaCol = AlphaPred;
	tmpMax   = ones(hsmm->getNumSTATES(),1)*maxAlpha;
	uint qAt = accu(find(tmpAlphaCol==tmpMax, 1,"first"));

	// Return attractor of most likeli state:
	return  hsmm->getMU(qAt);

}

void TrajMPC::createPrecisionMatrices(HSMM* _hsmm)
{
	// Get Relevant variables:
	uint _nbStates = _hsmm->getNumSTATES();
	uint _nbVar = _hsmm->getNumVARS();

	// Resize Lambda
	Lambda.set_size(_nbVar,_nbStates*_nbVar);

	// Copy Data
	for (uint i = 0;i<_nbStates;i++)
		Lambda.cols(i*_nbVar,_nbVar*(i+1)-1) = _hsmm->getSIGMA(i).i();

}

void TrajMPC::assignHSMM(HSMM* _hsmm, uint _nbVarPos, uint _nbDeriv)
{
	if (_hsmm->getNumVARS() != _nbVarPos*_nbDeriv)
	{
		cout << "The number of varialbles in the HSMM ("<< _hsmm->getNumVARS()
			<<	") does not correspond with the specified number of position variables and derivatives (" 
			<< _nbVarPos << "*" << _nbDeriv << " = " << _nbVarPos*_nbDeriv<<".";
		throw std::invalid_argument("Invalid Parameters");
	}
	else
	{
		hsmm = _hsmm;
	}

}

void TrajMPC::allocateReference(uint _nbVar, uint _nbVarPos, uint _Np, uint _Nc)
{
	R = zeros(_nbVarPos*_Nc,_nbVarPos*_Nc);
	Q = zeros(_nbVar*_Np, _nbVar*_Np);
	Muq = zeros(_nbVar*_Np);
	// Allocate space for q
	q = uvec(_Np);
}

void TrajMPC::createSystemMatrices(double _dt, uint _nbVarPos, uint _nbDeriv)
{
	// Create System Matrices
	// Assistive Matrices for System construction:	
	uint _nbVar = _nbVarPos*_nbDeriv;
	mat IvarPos, Idisc; 
	IvarPos = eye(_nbVarPos,_nbVarPos);
	Idisc = eye(_nbVar,_nbVar);
	
	// A matrix:
	A = zeros(_nbDeriv,_nbDeriv);
	A.submat(0,1,_nbDeriv-2,_nbDeriv-1) = eye(_nbDeriv-1,_nbDeriv-1);
	A = kron(A,IvarPos);

	// B Matrix:
	B = zeros(_nbDeriv,1);
	B(_nbDeriv-1,0)=1;
	B = kron(B,IvarPos);

	// C Matrix:
	C = zeros(1,_nbDeriv);
	C(0,0) = 1;
	C = kron(C,IvarPos);	

	// Discretize system matrices (Euler discretization):
	Ad = A*_dt+Idisc;
	Bd = B*_dt;


	cout <<"Ad: "<< endl << Ad << endl;
	cout <<"Bd: "<< endl << Bd << endl;

}

colvec&  TrajMPC::computeControlCommand(colvec& X)
{
	// Update hidden markov model
	hsmm->stepForwardVariable(X, in);

	// Create State Sequence:
	// Alpha Predictions
	hsmm->predictForwardVariable(AlphaPred);

	// Calculate the largest alpha values
	maxAlpha = max(AlphaPred,0);

	// Find corresponding corresponding indices (i.e. argmax())
	for (uint i = 0;i<myMPC->getNp();i++)
	{
		tmpAlphaCol = AlphaPred.col(i);
		tmpMax = ones(hsmm->getNumSTATES(),1)*maxAlpha(i);
		q(i)= accu(find(tmpAlphaCol==tmpMax, 1,"first"));
	}

	// Update update Reference Q and Muq:
	updateReference(q);

	// Compute control command:
	return myMPC->computeControlCommand(X,this->Muq,this->Q,this->R);
	
}

void TrajMPC::updateReference(uvec& _q)
{
	for (uint i = 0; i<_q.n_rows;i++)
	{
		// Construct Precision matrix based on state sequence
		Q.submat(i*nbVar,
			   	i*nbVar, 
				(i+1)*nbVar-1,
				(i+1)*nbVar-1) =
			Lambda.cols(_q(i)*nbVar,(_q(i)+1)*nbVar-1);

		// Construct Mean vector
		Muq.rows(i*nbVar,(i+1)*nbVar-1) =hsmm->getMU(_q(i));
	}

}

void TrajMPC::updateCostMatrix(double _alpha)
{
	// Construct Cost Matrix 
	R = eye(nbVarPos*myMPC->getNc(),nbVarPos*myMPC->getNc())*(_alpha*_alpha);
}

} // end namespace
