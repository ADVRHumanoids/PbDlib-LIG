/**
Copyright (C) 2014, Tohid Alizadeh, Milad Malekzadeh, Davide De Tommaso

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

/*! \file taskparameters.h
\brief TaskParameters class


A Task parameter represents a linear transformations Ax+b. A Task parameter can be used in two ways:
1) To transform demonstration data from a 'global' representation 'g'to the local representation 'l': l = A\(g-b)
2) To transform local representation of data(e.g. GMMs defined in local data) into global representations : g = Al+b
and vector b. 

The class TaskParameters represents a set of Task Parameters. By defining multiple task parameters global data can be
viewd from different local perspectives.

b: Position vector
A: Projection matrix
TaskParameter: struct of matrix A and vector b defining a linear transformation
nPARAMS: Number of task parameters - according to the task this should be determined before
nVARS: Number of variables
taskparams: The set of task parameters


\author Martijn Zeestraten, Tohid Alizadeh, Milad Malekzadeh, Davide De Tommaso
\bug No known bugs.
*/


#ifndef TASKPARAMETERS_H
#define TASKPARAMETERS_H

#include <fstream>
#include "armadillo"
#include "pbdlib/datapoints.h"


using namespace arma;

namespace pbdlib
{


struct TaskParameter 
{
    vec b;
    mat A;
};

class TaskParameters
{
    private:
        std::vector<TaskParameter> MyTPs; // Vector of task parameters
        uint nVARS; // Number of variables

    public:
        ~TaskParameters(){}
        TaskParameters(uint nVARS,uint nTaskParameters);

        std::vector<TaskParameter>& getTaskParameters();
		TaskParameter& getTaskParameters(uint index);
		uint getNumTASKPARAMETERS(){return MyTPs.size();}
		uint getNumVARS(){return nVARS;}

		void setTaskParameters(std::vector<TaskParameter> _TPs);
		void setTaskParameters(uint index, TaskParameter _TP);
        void loadFromFile(std::string path);
		void saveToFile(std::string filename, std::vector<std::string> TPnames);

};

} // end of pbdlib namespace

#endif // PARAMETERS_H
