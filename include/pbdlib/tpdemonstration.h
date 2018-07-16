/**
Copyright (C) 2014, Davide De Tommaso

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

/*! \file demonstration.h
\brief Demonstration class
The class Demonstration model a single demonstration in PbD framework

\author Davide De Tommaso
\bug No known bugs.
*/


#ifndef TPDEMONSTRATION_H
#define TPDEMONSTRATION_H

#include "armadillo"
#include <fstream>

#include "pbdlib/datapoints.h"
#include "pbdlib/taskparameters.h"

using namespace arma;

namespace pbdlib
{

class TPDemonstration
{
    private:

		// We store the 'global' demonstration data
        Datapoints Data;

		// We store unique task parameters for each datapoint:
		std::vector<TaskParameters> PointTPs;
		bool uniqueTPs; // flag that indicates if all TPs have unique TPs (used for saving)

    public:
        TPDemonstration(){}
		TPDemonstration(mat& _data, std::vector<TaskParameters> VTPs);
		TPDemonstration(mat& _data, TaskParameters TPs);
		TPDemonstration(std::string data_path, std::string tp_path, bool DataTrans= false);
		



		// Function used to set the data of the demonstration in the global reference
		void setGlobalData(mat& Data);

		// --------------- 	PARAMETER SET FUNCTIONS ---------------
		// Internally each data point has a unique value Task parameters object:
		// Each datapoint thus has nTaskParameters task parameters A and b.
		// This is done since during demonstration the task parameters can change 
		// from timestep to timestep. By the following function we can
		// set the unique task parameters corresponding to the datapoints:
		void setDataPointTPs(std::vector<TaskParameters> _TPs);
	
		// If you just want to set the same set of task parameters for each
		// data point you can use the following function. Here you need to pass the
		// Parameters
		// Eeach data point will then have the same task parameters. 
		void setDataPointTPs(TaskParameters _TPs);
		void setDataPointTPs(uint i, TaskParameters _Tps);
		void setDataPointTPs(std::vector<TaskParameter> _TPs);
		void setDataPointTPs(std::vector<mat> _A, std::vector<colvec> _b);
		
		// Functions To Load data from Files:
		void loadFromFiles(std::string data_path, std::string tp_path, bool Transpose =false);
		void loadTaskParametersFromFile(std::string tp_path);
		void loadDataFromFile(std::string data_path, bool DataTrans=false);
		
        void saveInFiles(std::string prefix,std::vector<std::string> par_names,
			   	std::vector<std::string> TPnames);
		// ----------------------- DATA GET FUNCTIONS ---------------------------
		// The data is stored in the global space. These functions are used to linear transform the data
		// according to the specified task parameters of each frame
		mat getDataInTaskParameters(uint i);
		mat getDataInTaskParameters();
		std::vector<TaskParameters>& getDataPointTPs(){return PointTPs;}
		TaskParameters& getDataPointTPs(uint pIndex){return PointTPs[pIndex];}
		uint getNumVARS(){return this->Data.getNumVARS();}
		uint getNumPOINTS(){return this->Data.getNumPOINTS();}
		uint getNumTASKPARAMETERS();
		Datapoints getDataInGlobalSpace(){return this->Data;}
		


};

} //end of pbdlib namespace

#endif
