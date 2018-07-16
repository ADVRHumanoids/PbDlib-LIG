
/**
Copyright (C) 2014, Danilo Bruno

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

/*! \file test_onlineDP
\brief Testing online Dirichlet Process clustering

\author Danilo Bruno
\bug No known bugs.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include "pbdlib/gmm.h"
#include "armadillo"
#include <vector>
#include <string>


using namespace pbdlib;
using namespace arma;

int main(int argc, char **argv)
{

    //load data

    std::ifstream dataFile("data/data.csv");

    std::string line;

    std::vector<colvec> Data;
    colvec PTmp = zeros(2,1);


    while(std::getline(dataFile,line)){
        std::stringstream lineStream(line);
        std::string charTmp;
        std::getline(lineStream,charTmp,',');
        PTmp(0) = std::atof(charTmp.c_str());
        std::getline(lineStream,charTmp,',');
        PTmp(1) = std::atof(charTmp.c_str());
        Data.push_back(PTmp);
    }

    GMM_Model gmm(1,2);
    double minSigma = 1E-5;
    double lambda = 0.05;

    mat minSIGMA = minSigma*eye(2,2);
    std::vector<GaussianDistribution> comps;
    GaussianDistribution componentTmp(Data[0],minSIGMA);

    comps.push_back(componentTmp);
    gmm.setCOMPONENTS(comps);

    rowvec priorsTmp = zeros(1,1);
    priorsTmp(0) = 1;
    gmm.setPRIORS(priorsTmp);
    int N = 1;
    for (int i=1;i<Data.size();i++){
        N++;
        colvec P = Data.at(i);
        gmm.onlineEMDP(N,P,lambda,minSigma);
    }

    for (int i=0;i<gmm.getNumSTATES();i++)
    {
        std::cout << gmm.getMU(i) << std::endl;
        std::cout << gmm.getSIGMA(i) << std::endl;
        std::cout << gmm.getPRIORS() << std::endl;
    }


	return 0;
}
