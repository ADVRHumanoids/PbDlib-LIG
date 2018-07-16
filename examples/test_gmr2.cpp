/**
Copyright (C) 2015, Martijn Zeestraten, Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon

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

/*! \file repro_gmr.cpp
\brief Learning GMM model and testing reproductions with GMR
Learning a GMM model from a demonstration saved in the file data_txyz.txt and after


In this example of GMR we set the input and output variables of the regression
by index and not by name. This operation is faster, and allows for more flexibility when
performing online regression. In addition, we show that we can perform regression on only part
of the GMM (i.e we use x and y and input and z as output. t is not used)

\author Martijn Zeestraten, Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon

\bug No known bugs.
*/

#include "pbdlib/gmm.h"
#include "pbdlib/gmr.h"
#include <sstream>

#define nbStates 3
#define nbVar 4
#define nbData 200

using namespace pbdlib;

int main(int argc, char **argv)
{
	// Get Demonstration data
	std::vector<Demonstration> demos;
	Demonstration demo =  Demonstration(nbVar,nbData);
	std::vector<std::string> vars,vars2;
	vars.push_back("t");
	vars.push_back("x");
	vars.push_back("y");
	vars.push_back("z");

	demo.getDatapoints().loadFromFile("data/data_txyz.txt");
	demo.getDatapoints().setVarNames(vars);
	
	cout << endl << "Data is: " << demo.getDatapoints().getData() << endl;
	cout << endl << "nbVar is: " << demo.getDatapoints().getNumVARS() << endl;
	cout << endl << "nbDatapoints is: " << demo.getDatapoints().getNumPOINTS() << endl;
	cout << "\nPress enter to continue..." << endl; getchar();

	demos.push_back(demo);
	// Create GMM
	GMM_Model* gmm;
	gmm = new GMM_Model(demos,nbStates);
	gmm->setVARSNames(vars);
	cout << "\n Number of EM iterations: " << gmm->EM_learn();
	
	for (int i=0;i<nbStates;i++){
		cout << "\n Mu_" << i << " = \n" << gmm->getMU(i);
		cout << "\n Sigma_" << i << " = \n" << gmm->getSIGMA(i);
	}
	cout << "\nPress enter to continue..." << endl; getchar();

	// Perform GMR
	// In this example of GMR we set the input and output variables of the regression
	// by index and not by name. This operation is faster, and allows for more flexibility when
	// performing online regression. In addition, we show that we can perform regression on only part
	// of the GMM (i.e we use x and y and input and z as output. t is not used)
	
	// Input data: We select x, and y as input data
	cout<< "Get Data" << endl;
	mat dataIn = demo.getDatapoints().getData().submat(1,0,2,nbData-1);
	cout << "Data In" << endl << dataIn << endl;
	GMR* reg = new GMR(gmm);

	// Perform Regression
	urowvec in, out;
	in << 1 << 2; // The input variables for regression (x, y)
	out << 3;     // The output variables for regression (z) 
	GMM_Model* gmmOut = reg->regression(dataIn,in,out);

	cout << "\n\n Reproduction with GMR, output: z, input : x,y :\n";

	for (int t=0;t<nbData;t+=40){
		cout << "\n MU_" << t << " = \n"  << gmmOut->getMU(t);
		cout << "\n SIGMA_" << t << " = \n"  << gmmOut->getSIGMA(t);
	}

	return 0;
}

