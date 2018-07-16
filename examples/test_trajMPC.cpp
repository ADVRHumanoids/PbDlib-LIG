#include <time.h>
#include "armadillo"
// Pbdlib dependencies:
#include "pbdlib/hsmm.h"
#include "pbdlib/taskparameters.h"
#include "pbdlib/gmm.h"
#include "pbdlib/trajMPC.h"


using namespace arma;
using namespace pbdlib;
using namespace std;


int main(int argc, char* argv[])
{
		
cout << "Setting Variables: ";

// Model Settings:
uint nVarsPos = 2; // Number of position variables 
uint nDeriv   = 2; // Number of derivatives considered

// TrajMPC Settings:
uint Np       = 70;  // Number of predictions
uint Nc       = 30;  // Number of control predictions
uint N        = 600; // Number of iterations of simulation
double dt     = 0.01;  // time step considered
double alpha  = 0.001; // control cost variable

// Initial state of the system:
colvec Xinit;
Xinit << -27.0 << 42.0<< 0.0 << 0.0;

string prefix = "HSMM_test";

// Source Files for the HSMM Model:
string muPath        = "data/hsmm/" + prefix + "_mu.txt";
string sigmaPath     = "data/hsmm/" + prefix + "_sigma.txt";
string priorsPath    = "data/hsmm/" + prefix + "_priors.txt";
string transPath     = "data/hsmm/" + prefix + "_trans.txt";
string durMu         = "data/hsmm/" + prefix + "_durMu.txt";
string durSigma      = "data/hsmm/" + prefix + "_durSigma.txt";

// Create pbdlib objects:
cout << "Create pbdlib objects ";
HSMM*    myHSMM;        // HSMM to handle state sequences
TrajMPC* myTrajMPC;   
cout << "done " << endl;

// Load HSMM:
cout << "load HSMM ";
myHSMM = new HSMM(priorsPath,muPath,sigmaPath,transPath,durMu,durSigma);
cout << "done " << endl;
cout << "HSMM nbStates: " << myHSMM->getNumSTATES() << endl;
cout << "HSMM nbVars  : " << myHSMM->getNumVARS() << endl;

// Load TrajMPC:
cout << "load TrajMPC ";
myTrajMPC = new TrajMPC(myHSMM,dt,nVarsPos,nDeriv,Np,Nc,alpha);
cout << "done " << endl;
cout << "TrajMPC nbVars  : " << myTrajMPC->getNumVARS() << endl;

// **********************************************************************

// Make sure that the TrajMPC model is correclty initialized by resetting it before starting
myTrajMPC->reset();

// Simulate reproduction:
colvec U;
mat Ad = myTrajMPC->getSystemDynamics();
mat Bd = myTrajMPC->getInputDynamics();

//colvec curPos = myTrajMPC->getCurrentAttractor();
colvec curPos = Xinit;
cout << "Current Attractor: " << endl <<curPos<< endl;
mat TestData= zeros(myHSMM->getNumVARS(),N);
mat Alphas = zeros(myHSMM->getNumSTATES(),N);
mat ControlData = zeros(nVarsPos,N);
umat qPredict;
qPredict.set_size(Np,N);

Alphas.col(0) = myHSMM->getForwardVariable();
TestData.col(0) = curPos;

int start = clock();
cout << "Simulation of the reproduction: " << endl;
for (uint it =1;it<N;it++)
{
	// Make step:
	U = myTrajMPC->computeControlCommand(curPos);
	
	// Simulate movement of attractor:
	curPos = Ad*curPos + Bd*U;
	
	// Display current Pos
//	cout << "it: " << it << ": " << curPos.t() << endl;

	// Save Data
	TestData.col(it)    = curPos;
	Alphas.col(it)      = myHSMM->getForwardVariable();
	ControlData.col(it) = U.subvec(0,nVarsPos-1);
	qPredict.col(it)    = myTrajMPC->getq();

}
int end = clock();
std::cout << "it took " << end - start << " ticks, or " << ((float)end - start)/CLOCKS_PER_SEC << " seconds." << std::endl;
std::cout << "freq:   " << float(N)/(((float)end - start)/CLOCKS_PER_SEC ) << endl;


cout << "Writing algorithm output to files... " ;
TestData.save("trajMPC_states.txt",raw_ascii);
Alphas.save("trajMPC_alpha.txt", raw_ascii);
qPredict.save("trajMPC_q.txt",raw_ascii);
ControlData.save("trajMPC_u.txt",raw_ascii);
cout << "Done" << endl;

}

