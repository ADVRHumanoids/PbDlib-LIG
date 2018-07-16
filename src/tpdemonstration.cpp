/**
 * 
Copyright (C) 2015, Martijn Zeestraten 

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

#include "pbdlib/tpdemonstration.h"

namespace pbdlib
{

// -----------------   CONSTRUCTORS -----------------------------------
TPDemonstration::TPDemonstration(std::string data_path, std::string tp_path, bool DataTrans)
{
	// Load Data:
	loadDataFromFile(data_path, DataTrans);

	// Load taskparameters:
	loadTaskParametersFromFile(tp_path);

}
TPDemonstration::TPDemonstration(mat& _data, std::vector<TaskParameters> _VTPs)
{
	// Set Data:
	this->Data.setData(_data);
	// Set Task Parameters:
	this->setDataPointTPs(_VTPs);
}
TPDemonstration::TPDemonstration(mat& _data, TaskParameters _TPs)
{
	// Set Data:
	this->Data.setData(_data);
	// Set Task Parmaters:
	this->setDataPointTPs(_TPs);
}







// ------------------ SET TASK PARAMETERS --------------------------------
void TPDemonstration::setDataPointTPs(std::vector<TaskParameters> _PointTPs)
{
	// input
	// _PointTPs: A vector containing [nPOINTS] objects TaskParameters each containing _nTaskParameters.
	//
	// Task-parameters are thus specified for each individual data point. This is usefull
	// when the task-parameters were changing during the demonstration phase.
	
	// Check Size:
	if (_PointTPs.size()!=this->getNumPOINTS())
	{
		throw std::invalid_argument("[ERROR]TPDemonstration::setDataPointTPs( ... ): \nDimensions of _PointTPs are not consistent with class initialization.");
	}
	
	// Different task-parameters for each Data point are given, assign to PointTPs:
	PointTPs = _PointTPs;

	// Set flag for unique Task Parameters
	uniqueTPs=true;

}

void TPDemonstration::setDataPointTPs(TaskParameters _TPs)
{
	// Input:
	// _TPs: A object containing _nTaskParameters number of task-parameters
	//
	// Different frames for different task parameters are not given. Therefore
	// we assume that the taskparameters are the same for each datapoint
	
	// Assign the same task-parameters for each datapoint:

	PointTPs.clear();

	for (uint i = 0;i<this->getNumPOINTS();i++)
		PointTPs.push_back(_TPs);

	// reset flag for unique TPs:
	uniqueTPs = false;
}

void TPDemonstration::setDataPointTPs(uint pIndex, TaskParameters _TPs)
{
	// This function writes a set of task parameters to a specific data point
	//
	// Input:
	// pIndex : Point index of the datapoint we are writing to
	// _TPs: An object containing _nTaskParameters number of task-parameters
	

	// Check index
	if (pIndex>=this->getNumPOINTS())
		throw std::invalid_argument("TPDemonstration::setDataPointTPs(...): \nIndex out of bounds.");

	// Check Parameter dimensions:
	if (_TPs.getNumTASKPARAMETERS() !=this->getNumTASKPARAMETERS()
			|| _TPs.getNumVARS()!= this->getNumVARS())
		throw std::invalid_argument("TPDemonstration::setDataPointTPs(...):\nTask parameter dimensions are incompatible with the number of dimensions of the parameter object.");

	// Write Parameters:
	PointTPs[pIndex] = _TPs;

	// Set flag for unique task parameters:
	uniqueTPs = true;
}

void TPDemonstration::setDataPointTPs(std::vector<TaskParameter> _VTP)
{
	// Input _VTP : A vector of the object TaskParameter 
	//                 Each element in the vector represents a parameter set A, b
	//                 
	// No task parameters are specified for each datapoint individually. Therefore
	// we assume that the taskparameters are the same for each datapoint
	
	// Copy parameters into parameters object:
	uint nVARS = this->getNumVARS();

	TaskParameters _TPs(nVARS, nVARS);
	for (uint i=0;i<_VTP.size();i++)
	{
		// Check Dimensions:
		if (_VTP[i].A.n_elem !=nVARS*nVARS 
				|| _VTP[i].b.n_elem != nVARS)
			throw std::invalid_argument("[ERROR]TPDemonstration::setDataPointTPs(...) \n : Task parameter dimensions are incompatible with the number of variables in the data.");
		else
			_TPs.setTaskParameters(i,_VTP[i]);
	}

	// Set the same for all data points:
	setDataPointTPs(_TPs);

	// reset flag for unique task parameters:
	uniqueTPs = false;
	
}

void TPDemonstration::setDataPointTPs(std::vector<mat> _A, std::vector<colvec> _b)
{
	// Function puts the Transformation Matrices containd in _A with 
	// corresponding translation matrices contained in _b into a TaskParameters
	// object. The taskparameters will hold for all data points
	// in the data of this class.
	//
	// Inputs:
	// _A : Vector [nTPs] of matrices [nVARS x nVARS];
	// _b : Vector [nTps] of colvecs [nVars x 1];


	// Check dimensions of _A and _b
	if (_A.size() != _b.size())
		throw std::invalid_argument("[ERROR]TPDemonstration::setDataPointTPs(...): \nSize of Vectors _A and _b are not consistent.");

	uint nVARS = this->getNumVARS();
	uint nTPs  = _A.size();

	// Copy parameters into TaskParametere structure:
	TaskParameters tmpTPs(this->getNumVARS(),nTPs);
	for (uint m=0;m<nTPs;m++)
	{
		// Check dimensions
		if (_A[m].n_elem !=nVARS*nVARS || _b[m].n_elem !=nVARS)
			throw std::invalid_argument("[ERROR]TPDemonstration::setDataPointTPs(...): \nTask parameter dimensions are incompatible with the number of variables in the data.");
		else
		{
			// Write A and b to temp parameters
			tmpTPs.getTaskParameters(m).A = _A[m];
			tmpTPs.getTaskParameters(m).b = _b[m];
		}
	}

	// Set the same for all data points:
	setDataPointTPs(tmpTPs);

	// reset flag for unique task parameters:
	uniqueTPs = false;
}





// ---------------------------- SET GLOBAL DATA -----------------------------------
void TPDemonstration::setGlobalData(mat& _Data)
{
	// Check Size:
	if (_Data.n_cols != this->getNumPOINTS()
			&& _Data.n_cols != this->getNumVARS())
	{
			throw std::invalid_argument("[ERROR] TPDemonstration::setGlobalData(...): \nDimensions of _Data are not consistent with class initialization.");
	}
	this->Data.setData(_Data);
}





// ---------------------- LOAD DATA FROM FILES -----------------------------------
void TPDemonstration::loadFromFiles(std::string data_path, std::string tp_path, bool DataTrans)
{
	// Function loads the demonstration data and corresponding
	// task parameters from the files specified.
	
	// First load data
	this->Data.loadFromFile(data_path, DataTrans);

	// Set dimensions:
	PointTPs.reserve(this->getNumPOINTS());

	// Load Parameters:
	loadTaskParametersFromFile(tp_path);
}

void TPDemonstration::loadDataFromFile(std::string data_path, bool DataTrans)
{

	mat dataRaw;
	if( !dataRaw.load(data_path, raw_ascii) )
	{
		throw std::invalid_argument( "[ERROR]::TaskParameters::readParamsFromTxtFile(...) if( paramTmp.load(path, raw_ascii) )   ... else .");
	}
	if (DataTrans)
		dataRaw = dataRaw.t();

	// Set global data:
	setGlobalData(dataRaw);
}







// --------------------- LOAD TASK PARAMETERS -----------------------------------
void TPDemonstration::loadTaskParametersFromFile(std::string tp_path)
{
	// Function loads Taskparameters from file.
	// The file is assumed to have the following structure:
	//
	//             <---- Data Points ---->
	//      ^    [b11^T, b12^T, ..., b1p^T]
	//      |    [A11  , A12  , ..., A1p  ]
	//      |    [b21^T, b22^T, ..., b2p^T]
	//     TPs   [A21  , A22  , ..., A2p  ]
	//      |    [...  , ...  , ..., ...  ]
	//      |    [bm1^T, bm2^T, ..., bmp^T]
	//      v    [Am1  , Am2  , ..., Amp  ]
	//
	// Where p represents the number of datapoints nPOINTS and 
	// m the number of taskparameters nTaskParameters. When the file only
	// contains one column of task parmeters. It is assumed that all datapoints
	// where recorded using the same task parameters.
	//
	// Input:
	// Path to data file

	// Load data from file:
	mat TPraw;
	uint nVARS = this->getNumVARS();

	if( !TPraw.load(tp_path, raw_ascii) )
	{
		throw std::invalid_argument("[ERROR]TPDemonstration::loadTaskParametersFromFile(...): \n Failure loading parameter file.");
	}
	
	
	// check if the task parameters have the correct dimensions:
	// Number of cols consistent|| Number of rows consistent
	if ((TPraw.n_cols % nVARS)!=0 || (TPraw.n_rows % (nVARS+1)) !=0)
	{
		throw std::invalid_argument("[Error]TPDemonstration::loadTaskParametersFromFile(...):\n The size of the parameter file is inconsistent with the data dimensions." );
	}
	uint nTPs = TPraw.n_rows/(nVARS+1); // Number of TPs per column


	// Check if the number TPs is equal to 1 (meaning the same set of TPs 
	// for each datapoint), or equal to NumPOINTS or equal to NumPOINTS
	// (meaning unique TPs for each datapoint)
	uint _nPoints = TPraw.n_cols/nVARS; // Number of columns with TPs
	if (_nPoints!=this->getNumPOINTS()&& _nPoints!=1) 
		throw std::invalid_argument("[Error]TPDemonstration::loadTaskParametersFromFile(...): \n The size of the data in the file is inconsistent with the specified	parameter dimensions." );

	// Create a vector of Taskparameters object
	std::vector<TaskParameters> tmpVTPs;

	cout << "loading " << endl;
	cout << "TPraw" << endl << TPraw << endl;
	cout << "nVars: " << nVARS << endl;
	cout << "nTPs : " << nTPs << endl;
	cout << "_nPoints: " << _nPoints << endl;

	// Create a tmp TaskParameters Object
	uint TPind;
	TaskParameters tmpTPs(nVARS,nTPs);
	// For each column of Task Parameters
	for(uint p=0;p<_nPoints;p++)
	{
		// For each Task Parameter in the column
		for(uint m=0; m<nTPs; m++)
		{
			TPind= m*(nVARS+1);
			tmpTPs.getTaskParameters(m).A =
				TPraw.submat(TPind+1 ,p*nVARS, TPind+nVARS,(p+1)*nVARS-1);
			tmpTPs.getTaskParameters(m).b =
				trans(TPraw.submat(TPind, p*nVARS,TPind,(p+1)*nVARS-1));
		}
		// Add to vector of TPs:
		cout << "Push back ";
		tmpVTPs.push_back(tmpTPs); 
	}
		cout << "done" << endl;

	// Set task parameters to object:
	if(_nPoints==1)
	{
		// Set the same task parameters for each data point
		cout << "Set Data Points1 " << endl;
		cout << tmpVTPs.size() << endl;
		setDataPointTPs(tmpVTPs[0]);
	}
	else
	{
		// Set unique task parameters for each data point:
		cout << "Set Data Points2 " ;
		setDataPointTPs(tmpVTPs);
	}
	cout << "Done " << endl;
}






void TPDemonstration::saveInFiles(std::string prefix,
		std::vector<std::string> varnames, 
		std::vector<std::string> TPnames)
{
	// Function Saves data and task parameters into files:
	// - prefix_data.txt
	// - prefix_TaskParameters.txt
	// - prefix_VarNames.txt
	// - prefix_TPNamex.txt
	//
	// The file of the task parameters will have the following
	// structure: 
	//
	//             <---- Data Points ---->
	//      ^    [b11^T, b12^T, ..., b1p^T]
	//      |    [A11  , A12  , ..., A1p  ]
	//      |    [b21^T, b22^T, ..., b2p^T]
	//     TPs   [A21  , A22  , ..., A2p  ]
	//      |    [...  , ...  , ..., ...  ]
	//      |    [bm1^T, bm2^T, ..., bmp^T]
	//      v    [Am1  , Am2  , ..., Amp  ]
	//
	// Where p represents the number of datapoints nPOINTS and 
	// m the number of taskparameters nTaskParameters. 
	//
	// When the file only contains unique task parameters then
	// the task parameter file will have the following structure:
	//
	//      ^    [b11^T]
	//      |    [A11  ]
	//      |    [b21^T]
	//     TPs   [A21  ]
	//      |    [...  ]
	//      |    [bm1^T]
	//      v    [Am1  ]

	// Check size of varnames:
	if (varnames.size()!=this->getNumVARS())
		throw std::invalid_argument("[ERROR]TPDemonstration::saveInFiles(...): \n The number of provided filenames does not correspond to the number of variables.");

	// Check size of TaskParameter names
	if (TPnames.size()!=this->getNumTASKPARAMETERS())
		throw std::invalid_argument("[ERROR]TPDemonstration::saveInFiles(...): \n The number of provided Task Parameter names does not correspond to the number Task Parameters.");

	
	std::string data_fname = prefix + "_data.txt";
	std::string TPs_fname  = prefix + "_TaskParameters.txt";
	std::string varNames_fname = prefix + "_VarNames.txt";
	std::string TPnames_fname = prefix + "_TPNames.txt";


	// Write Data
	this->Data.setVarNames(varnames);
	this->Data.saveInFile(data_fname);
	
	// Save variable names:
	std::ofstream fVarNames(varNames_fname.c_str());  
	for(unsigned int i=0; i<varnames.size(); i++)    
		fVarNames<< varnames[i] << " ";
	fVarNames.close();

	// Put Task Parameters in Matrix:
	if (uniqueTPs)
	{
		
		TaskParameters tmpTPs(this->getNumVARS(),this->getNumTASKPARAMETERS());
		TaskParameter tmpTP;
		mat tmpMat1;
		mat tmpMat2;
		for (uint p=0;p<this->getNumPOINTS();p++)
		{
			tmpTPs = this->getDataPointTPs(p);
			for (uint m=0;m<this->getNumTASKPARAMETERS();m++)
			{
				// Concatenate TP in sets per data point:
				tmpTP = tmpTPs.getTaskParameters(m);
				if (m==0)
					tmpMat1 = join_cols(tmpTP.b.t(),tmpTP.A);
				else
					tmpMat1 = join_cols(tmpMat1,join_cols(tmpTP.b.t(),tmpTP.A));
			}
			
			// Concatenate TPs of data points in columns:
			if (p==0)
				tmpMat2 = tmpMat1;
			else
				tmpMat2 = join_rows(tmpMat2,tmpMat1);
		}
		// Save to File:
		tmpMat2.save(TPs_fname.c_str(),raw_ascii);

		// Save Task Parameter names:
		std::ofstream fTPnames(TPnames_fname.c_str());  
		for(unsigned int i=0; i<TPnames.size(); i++)    
			fTPnames<< TPnames[i] << endl;
		fTPnames.close();
	}
	else
	{
		PointTPs[0].saveToFile(prefix, TPnames);
	}


	

}



// ----------------------- GET INT TASKPARAMETERS FUNCTIONS ---------------------------

// The data is stored in the global space. These functions are used to linear transform the data
// according to the specified task parameters of each frame
mat TPDemonstration::getDataInTaskParameters(uint _Pindex)
{
	// For each Data point
	uint nPOINTS = this->getNumPOINTS();
	mat dataOut = zeros(this->getNumVARS(),nPOINTS);
	TaskParameter tmpParam;
	for (uint i =0;i<nPOINTS;i++)
	{
		// Get parameters _Pindex corresponding to datapoint i:
		tmpParam = PointTPs[i].getTaskParameters(_Pindex);

		// Project data from global to local frame:  g = Al+b => l =A\(g-b) 
		dataOut.col(i) = solve(tmpParam.A,this->Data.getData(i)-tmpParam.b);
	}
	return dataOut;
}

mat TPDemonstration::getDataInTaskParameters()
{
	// Create Matrix of appropriate size:
	uint nVARS   = this->getNumVARS();
	uint nPOINTS = this->getNumPOINTS();
	mat dataOut  = zeros(nVARS*this->getNumTASKPARAMETERS(),nPOINTS);

	for (uint m=0;m<this->getNumTASKPARAMETERS();m++)
	{
		// Stack data of each frame:
		dataOut.rows(m*nVARS, (m+1)*nVARS-1) = getDataInTaskParameters(m);
	}
	
	return dataOut;


}


uint TPDemonstration::getNumTASKPARAMETERS()
{
	if (PointTPs.size()==0)
		return 0;
	else
		return this->PointTPs[0].getNumTASKPARAMETERS();
}

} // end of pbdlib namespace
