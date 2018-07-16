// May 2015, Martijn Zeestraten
//
//

#include "pbdlib/mpc.h"

/* MPC Constructor:
 * Ad : Discretized system matrix A
 * Bd : Discretized input matrix B
 * C  ; Ouput matrix
 * Np : Prediction horizon
 * Nc : Control Horizon
*/
namespace pbdlib{
MPC::MPC(mat& _Ad, mat& _Bd, mat& _C, uint _Np, uint _Nc)
{

	A = _Ad;
	B = _Bd;
	C = _C;

	Np = _Np;
	Nc = _Nc;

	// Construct dynamic matrices required to perform MPC
	constructDynamicMatrices(_Ad,_Bd,_C,_Np,_Nc);

	// Allocate appropriate size for Least squares:
	prepareForLeastSquares();
}

void MPC::constructDynamicMatrices(mat& _A, mat& _B, mat& _C, uint _Np, uint _Nc)
{
	// Determine Sizes:
	int An = _A.n_rows;
	int Am = _A.n_cols;
	int Bm = _B.n_cols;
  
	
	// Construct F:
	//     [ CA       ]
	//     [ CA^2     ]
	//F =  [ ..       ]
	//     [ CA^{Np-1}]
	//     [ CA^Np    ] 
	//
	F.reshape(An*_Np,Am);
	mat C1 = mat(An*_Np,Bm);

	F.rows(0,An-1) = _A;
	C1.rows(0,An-1) = _B;
	for (uint i = 1;i<_Np;i++)
	{
		F.rows(i*An,(i+1)*An-1) = 
			F.rows((i-1)*An,i*An-1)*_A;
		C1.rows(i*An,(i+1)*An-1) = 
			F.rows((i-1)*An,i*An-1)*_B;
	}
//	cout << "F created " << endl;
//	cout << F << endl;

	// Construct Phi
	//        [ CB     0     0 
	//        [ CAB    CB    
	//        [ CA^2B  CAB
	//
	// Phi =  [ ..      ..  ..  ..
	//        [ CA^(Np-1)
	//        [ 
	Phi.resize(An*_Np,Bm*_Nc);
	for (uint i=0;i<_Nc;i++)
	{

		Phi.submat(i*An, i*Bm, _Np*An-1, (i+1)*Bm-1) = 
			C1.rows(0,(_Np-i)*An-1);
	}
//	cout << "Phi Created" << endl;
//	cout << Phi << endl;

	// Construct I Matrix (to select first control command
	I.resize(Bm, Bm*Nc);
	I.cols(0,Bm-1) = eye(Bm,Bm);

	CT = kron(eye(_Np,_Np),_C);
}
void MPC::prepareForLeastSquares()
{
	// Function is used to set all matrices to the appropriate size
	
	// Get sizes
	nbVar = C.n_rows;   // Number of tracking variables
	nbContr = B.n_cols; // Number of control variables

	PhiQPhiR.resize(nbContr*Nc,nbContr*Nc);
	PhiQMuqFx.resize(nbContr*Nc);
}


colvec& MPC::computeControlCommand(colvec& X,colvec& ref, mat& Q, mat& R)
{
	// Solves a damped least-squares problem
	
	// Define the two parts:
	PhiQPhiR= Phi.t()*Q*Phi+R;
	PhiQMuqFx= Phi.t()*Q*(ref-F*X);

	// Solve the least squares and return the first control command:
	U =I*solve(PhiQPhiR,PhiQMuqFx);
	return U;
}

void MPC::setSystem(mat& _A, mat& _B, mat& _C)
{
	// Copy matrices to local variables
	A = _A;
	B = _B;
	C = _C; 

	// Recompute dynamic matrices
	constructDynamicMatrices(A,B,C,Np,Nc);

	// Allocate appropriate size for Least squares:
	prepareForLeastSquares();
}

void MPC::setPredictionHorizon(uint _Np)
{
	// Copy variable
	Np = _Np;

	// Recompute dynamic matrices
	constructDynamicMatrices(A,B,C,Np,Nc);

	// Allocate appropriate size for Least squares:
	prepareForLeastSquares();
}

void MPC::setControlHorizon(uint _Nc)
{

	// Copy variable
	Nc = _Nc;

	// Recompute dynamic matrices
	constructDynamicMatrices(A,B,C,Np,Nc);
	
	// Allocate appropriate size for Least squares:
	prepareForLeastSquares();
}

}
