/**
Copyright (C) 2015, João Silvério

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

/*! \file test_quaternions.cpp
\brief testing Quaternion class

Testing basic features of the Quaternion class.

\author João Silvério
\bug No known bugs.
*/

#include "pbdlib/quaternion.h"
#include <sstream>

using namespace pbdlib;
using namespace arma;

int main(int argc, char **argv)
{
	SortConvention sort = SCALAR_VEC;

	cout << endl;
	cout << " ** Testing pbdlib::Quaternion class **" << endl;
	cout << endl;

	// Generate an identity quaternion
	cout << "-- Creating a quaternion using default constructor... " << endl;
	Quaternion q1(sort);
	cout << "q1 = ";
	q1.print();
	cout << "norm = " << q1.norm() << endl << endl;

	// Initialize a quaternion from a vec3 and a double
	cout << "-- Constructing quaternion from a vec3 and a double... " << endl;
	vec3 vecPart;
	double scalarPart;
	vecPart[0] = 0.5332;
	vecPart[1] = 0.5928;
	vecPart[2] = 0.0831;
	scalarPart = 0.5978;
	Quaternion q2(vecPart,scalarPart,sort);
	q2.normalize();
	cout << "q2 = ";
	q2.print();
	cout << "norm = " << q2.norm() << endl << endl;

	// Initialize a quaternion from a vec4
	cout << "-- Constructing quaternion from a vec4... " << endl ;
	vec4 quat;

//	quat[0] = scalarPart;
//	quat.subvec(1,3) = vecPart;
	quat.subvec(0,2) = vecPart;
	quat[3] = scalarPart;

	Quaternion q3(quat, sort);	// the order in quat is reverted on purpose, so that it is clear that 'sort' works
	q3.normalize();
	cout << "q3 = ";
	q3.print();
	cout << "norm = " << q3.norm() << endl << endl;

	// Initialize a quaternion from 4 elements
	cout << "-- Constructing quaternion from 4 double... " << endl;
	Quaternion q4(quat[0],quat[1],quat[2],quat[3], sort);
	q4.normalize();
	cout << "q4 = ";
	q4.print();
	cout << "norm = " << q4.norm() << endl << endl;

	// Initialize a quaternion from an angle and an axis
	cout << "-- Constructing quaternion from axis-angle pair " << endl;
	double angle = M_PI/2.0;
	vec3 axis = zeros(3);
	axis[0] = 1.0;
	cout << "Angle = " << angle << "rad, " << " Axis = "; axis.t().print();
	Quaternion qAxisAngle(angle, axis, sort);
//	qAxisAngle.normalize();
	cout << "\nqAxisAngle = ";
	qAxisAngle.print();
	cout << "norm = " << qAxisAngle.norm() << endl << endl;

	// Initialize a quaternion from an vec3
//	cout << "-- Constructing quaternion from rotation vector (norm gives the angle, direction gives the axis) " << endl;
//	vec3 rotVec;
//	rotVec << -0.8626 << -1.8820 << -2.0997;
//	cout << "Rotation vector = " << rotVec.t().print();
//	Quaternion qRotVec(rotVec, sort);
//	cout << "\nqRotVec = ";
//	qRotVec.print();
//	cout << "norm = " << qRotVec.norm() << endl << endl;

	// Generate matrix representation of quaternion
	cout << "-- Generating matrix representation of quaternion... " << endl;
	mat44 quatMatrix = q4.matrix();
	cout << "Q4 = " << endl;
	quatMatrix.print();
	cout << endl;

	// Test quaternion product computation
	cout << "-- Computing quaternion product q3*q2... " << endl;
	Quaternion q5 = q3*q2;
	cout << "q5 = ";
	q5.print();
	cout << "norm = " << q5.norm() << endl << endl;

	// Compute quaternion product using matrix reverse
	cout << "-- Computing quaternion product using matrix reverse... " << endl;
	vec4 quatProdOut;
    quatProdOut = q2.matrixReverse()*q3.getAllCoeffs();
    cout << "q5 = ";
    quatProdOut.t().print();
    cout << endl;

	// Getting the 4-elements of a quaternion into a vec4/colvec4
	cout << "-- Retrieving the 4 elements of the quaternion as a vec4... " << endl;
	vec4 coeffs;
	coeffs = q5.getAllCoeffs();
	cout << "q5 = ";
	coeffs.t().print();
	cout << endl;

	// Computing the logarithm of a quaternion
	cout << "-- Log of quaternion q5... " << endl;
	cout << "log(q5) = ";
	qLog(q5).t().print();
	cout << endl;

	// Computing the exponential of a
	cout << "-- Exponential of log(q5) (should be again q5)... " << endl;
	cout << "exp(log(q5)) = ";
	qExp(qLog(q5),q5.getSort()).getAllCoeffs().t().print();
	cout << endl;

	// Angular displacement between q2 and q3
	cout << "-- Angular displacement between q2 and q3 (= angular velocity that rotates q2 into q3 in the unit time)... " << endl;
	cout << "omega = ";
	angDiff(q3,q2).t().print();
	cout << endl;

	// Rotate a 3D point
	cout << "-- Rotation of a 3D point by q5... " << endl;
	vec3 testPoint;
	testPoint << 0.8147 << 0.9058 << 0.1270;
	testPoint.print("p = ");
	cout << endl;
	vec3 rotatedPoint = q5*testPoint;
	rotatedPoint.print("Rotated p = ");
	cout << endl;

	// Rotation matrix
	cout << "-- Rotation matrix from q5... " << endl;
	mat33 rotMat = q5.getRot();
	rotMat.print("Rotation Matrix = ");
	cout << endl;

	return 0;
}
