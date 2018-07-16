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

#include "pbdlib/demonstration.h"

namespace pbdlib
{

Demonstration::Demonstration(uint _nVARS, uint _nPOINTS)
{
	data = Datapoints(_nVARS,_nPOINTS);
}

Datapoints& Demonstration::getDatapoints()
{
	return data;
}

void Demonstration::saveInFile(std::string path)
{
	data.getData().save(path, raw_ascii);
}

} // end of pbdlib namespace
