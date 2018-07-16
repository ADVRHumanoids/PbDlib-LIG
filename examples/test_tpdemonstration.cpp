/**
Copyright (C) 2015,Martijn Zeestraten, Tohid Alizadeh, Davide De Tommaso

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

/*! \file test_demonstration.cpp
\brief testing Demonstration and Parameters classes
Testing basic features of the Demonstration and Parameters classes

\author Martijn Zeestraten 
\bug No known bugs.
*/

#include "pbdlib/tpdemonstration.h"
#include "pbdlib/taskparameters.h"
#include <sstream>

using namespace pbdlib;
using namespace arma;

int main()
{
	uint nDemos = 4; // load from 4 files

	// loading the demonstrations and task parameters from the txt files:
	std::vector<TPDemonstration> demos; // Vector to hold demonstrations
	TPDemonstration tmpDemo; // Temp variable to load individual demonstrations

	char Datafilename[256];
	char TPfilename[256];
	cout << "Loading the demonstrations and the task parameters ..." << endl;
	for (uint m=0; m<nDemos; m++){   //Loading demos in the loop
		// Form filenames:
		sprintf(Datafilename, "../../data/pgmm/Data0%d.txt",m+1);
		sprintf(TPfilename, "../../data/pgmm/Param0%d.txt",m+1);

		// Load files in demo:
		// note we need to transpose the data (last argument)
		// because the data in the file a row for each data point
		// and a column for each variable. The agreement is
		// the other way around
		tmpDemo.loadFromFiles(Datafilename,TPfilename, true);

		// push back demonstration:
		demos.push_back(tmpDemo);
	}
	
	cout << "Demonstrations Loaded succesfully." << endl;

	// Show Data:
	for (uint m =0;m<nDemos ;m++)
	{
		cout << "----- DEMONSTRATION : " << m << endl;
		cout << "            NbVar: " << demos[m].getNumVARS() << endl;
		cout << "         NbPoints: " << demos[m].getNumPOINTS() << endl;
		cout << "NbTask Parameters: " << demos[m].getNumTASKPARAMETERS() << endl;
		cout << "----- Global Data.cols(0,10): " << endl;
		cout << demos[m].getDataInGlobalSpace().getData().cols(0,10).t() << "..." << endl;
		cout << "[Enter] To continue to Task Parameter Info" << endl;
		getchar();
		for (uint i=0;i<demos[m].getNumTASKPARAMETERS();i++)
		{
			cout << " Info Task Parameter " << i << endl;	
			// For datapoint 0 show task parameter i:
			cout << "b^T: " << 
				demos[m].getDataPointTPs(0).getTaskParameters(i).b.t() 
				<< endl;
			cout << "A: " << endl <<
				demos[m].getDataPointTPs(0).getTaskParameters(i).A.t() 
				<< endl;
			cout << "Data.cols(0,10) in Task Parameter " << i << " space" << endl;
			cout << demos[m].getDataInTaskParameters(i).cols(0,10).t() << "...." << endl;
			cout << "[Enter] To continue to next Task Parameter" << endl;
			getchar();
		}
	}
			
	// Replacing Task Parameters:
	cout << "Replacing Task variables for all Data...[Enter] " << endl;
	getchar();

	// Create new Idendity Parameters:
	TaskParameters newTPs(demos[0].getNumVARS(), 1);
	newTPs.getTaskParameters(0).A = eye<mat>(newTPs.getNumVARS(),newTPs.getNumVARS());
	newTPs.getTaskParameters(0).b = zeros(newTPs.getNumVARS(),1);

	// Apply parameters:
	
	for (uint m=0;m<nDemos;m++)
	{
		// Set the same parameters for all datapoints:
		demos[m].setDataPointTPs(newTPs);
		
		// Display new info:
		for (uint i=0;i<demos[m].getNumTASKPARAMETERS();i++)
		{
			cout << " Info Task Parameter " << i << endl;	
			// For datapoint 0 show task parameter i:
			cout << "b^T: " << 
				demos[m].getDataPointTPs(0).getTaskParameters(i).b.t() 
				<< endl;
			cout << "A: " << endl <<
				demos[m].getDataPointTPs(0).getTaskParameters(i).A.t() 
				<< endl;
			cout << "Data.cols(0,10) in Task Parameter " << i << " space" << endl;
			cout << demos[m].getDataInTaskParameters(i).cols(0,10).t() << "...." << endl;
			cout << "[Enter] To continue to next Task Parameter" << endl;
			getchar();
		}

	}
	return 0;

}
