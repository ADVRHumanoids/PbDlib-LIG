/**
Copyright (C) 2014, Davide De Tommaso, Milad Malekzadeh

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


#include "pbdlib/taskparameters.h"

namespace pbdlib
{

TaskParameters::TaskParameters(uint _nVARS, uint nTaskParameters)
{
	nVARS = _nVARS;
	MyTPs.resize(nTaskParameters);
	// Initialize task parameters to zero:
	for (uint i = 0;i<nTaskParameters;i++)
	{
		MyTPs[i].A = zeros(_nVARS,_nVARS);
		MyTPs[i].b = zeros(_nVARS,1);
	}
}

std::vector<TaskParameter>& TaskParameters::getTaskParameters()
{
	return MyTPs;
}

TaskParameter& TaskParameters::getTaskParameters(uint _index)
{
	return MyTPs[_index];
}

void TaskParameters::setTaskParameters(std::vector<TaskParameter> _taskParameters)
{
	MyTPs= _taskParameters;
}

void TaskParameters::setTaskParameters(uint _index, TaskParameter _param)
{
	if (MyTPs.size()<=_index)
		throw std::invalid_argument("TaskParameters::setTaskParameters(uint _index, TaskParameter _param)\n _index cannot be larger or equal than number of frames in TPGMM.");
	else
		MyTPs[_index] = _param;
}

void TaskParameters::loadFromFile(std::string path)
{
	mat rawTPs;
	// Load Data
	if( !rawTPs.load(path, raw_ascii) )
	{
		throw std::invalid_argument( "\n [ERROR]::TaskParameters::loadFromFile(std::string path) if( !rawTPs.load(path, raw_ascii) )   ... else .");
		return;
	}
	
	// Check dimensions
	if (rawTPs.n_rows % (nVARS+1)!=0|| rawTPs.n_cols!=nVARS)
	{
		throw std::invalid_argument("[Error]::TaskParameters::loadFromFile(std::string path) the size of the data in the file is not consistent with the specified	parameter dimensions" );
	}

	
	// Allocate memory:
	uint _nTPs = rawTPs.n_rows/(nVARS+1);
	MyTPs.resize(_nTPs);

	// Write task parameters:
	uint m;
	uint index;

	for(m=0; m<MyTPs.size(); m++)
	{
		index = m*(nVARS+1);
		MyTPs[m].b = trans(rawTPs.row(index));
		MyTPs[m].A = rawTPs.rows(index+1 , index+nVARS);
	}

}

void TaskParameters::saveToFile(std::string prefix, std::vector<std::string> TPnames)
{

	// Check size of TaskParameter name
	if (TPnames.size()!=this->getNumTASKPARAMETERS())
		throw std::invalid_argument("[ERROR]TaskParameters::saveToFile(...): \n The number of provided Task Parameter names does not correspond to the number Task Parameters.");

	std::string data_fname = prefix + "_TaskParameters.txt";
	std::string TPnames_fname = prefix + "_TPNames.txt";

	mat tmpMat;	
	TaskParameter tmpTP;
	for (uint m=0;m<this->getNumTASKPARAMETERS();m++)
	{
		tmpTP = this->getTaskParameters(m);
		if (m==0)
			tmpMat = join_cols(tmpTP.b.t(),tmpTP.A);
		else
			tmpMat = join_cols(tmpMat,join_cols(tmpTP.b.t(),tmpTP.A));
	}
	
	// Save Frames:
	tmpMat.save(data_fname.c_str(),raw_ascii);

	// Save Task Parameter names:
	std::ofstream fTPnames(TPnames_fname.c_str());  
	for(unsigned int i=0; i<TPnames.size(); i++)    
    	fTPnames<< TPnames[i] << " ";
	fTPnames.close();


}

} //end of pbdlib namespace

