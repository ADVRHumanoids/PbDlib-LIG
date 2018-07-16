/**
Copyright (C) 2014, Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon

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

/*! \file learn_gmm.cpp
\brief Learning GMM model
Learning a GMM model from a demonstration saved in the file data_txyz.txt

\author Davide De Tommaso, Milad Malekzadeh, Sylvain Calinon
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
	std::vector<Demonstration> demos;
	Demonstration demo =  Demonstration(nbVar,nbData);
	std::vector<std::string> vars;
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

	GMM_Model *gmm;
	gmm = new GMM_Model(demos,nbStates);
	gmm->setVARSNames(vars);

	cout << "\n Number of EM iterations: " << gmm->EM_learn();
	for (int i=0;i<nbStates;i++){
		cout << "\n Mu_" << i << " = \n" << gmm->getMU(i);
		cout << "\n Sigma_" << i << " = \n" << gmm->getSIGMA(i);
	}
	//Save GMM in files GMM_priors.txt, GMM_mu.txt, GMM_sigma.txt and GMM_vars.txt
	gmm->saveInFiles();
	cout << "\nPress enter to continue..." << endl; getchar();

	//Loading the GMM model from files...  
	GMM_Model *gmm2 = new GMM_Model("data/gmm/GMM_priors.txt",
									"data/gmm/GMM_mu.txt",
									"data/gmm/GMM_sigma.txt",
									"data/gmm/GMM_vars.txt");

	for (int i=0;i<nbStates;i++)
		gmm2->setMU(i, gmm2->getMU(i)*10.0);
	for (int i=0;i<nbStates;i++)
		cout << "\n Mu_" << i << " = \n" << gmm2->getMU(i);
	return 0;
}

