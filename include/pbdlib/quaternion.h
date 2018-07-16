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

/*! \file quaternion.h
\brief Defines a class for basic quaternion handling using armadillo functions and types.

It aims to facilitate quaternion modeling in task-parameterized models.

The code is adapted from http://www.cs.stanford.edu/~acoates/quaternion.h to:
	-> use armadillo data types
	-> allow the usage of quaternions as [u vec] or [vec u]
	-> compute the quaternion matrix that permits changing multiplication order

Functions will be implemented gradually on a need-to basis.

\author João Silvério
\bug No known bugs.
 */


#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>
#include <math.h>
#include <sstream>
#include "armadillo"

using namespace arma;

enum SortConvention { SCALAR_VEC , VEC_SCALAR };

namespace pbdlib{

class Quaternion {

	colvec4 coeffs;
	SortConvention sort;

public:

	// Constructors
	Quaternion(SortConvention sortInit = VEC_SCALAR); 	// Default constructor - initializes identity quaternion
	Quaternion(const vec3& vec, double scalar, SortConvention sortInit = VEC_SCALAR); 	// Initialize a quaternion from vec3 and scalar
	Quaternion(const vec4& vec, SortConvention sortInit = VEC_SCALAR); // Initialize a quaternion from a vec4
	Quaternion(const double* array, SortConvention sortInit = VEC_SCALAR); // Initialize a quaternion from an array of double
	Quaternion(double q1, double q2, double q3, double q4,  SortConvention sortInit = VEC_SCALAR);	// Initialize a quaternion from 4 scalars of type double
	Quaternion(double angle, vec3 axis,  SortConvention sortInit = VEC_SCALAR); 	// Initialize a quaternion from an angle and an axis
//	Quaternion(const vec3& angleAxis,  SortConvention sortInit = VEC_SCALAR); 	// Initialize a quaternion from a vec3 whose direction is the rotation axis and norm is the angle of rotation

	// Coefficients
	double x() const;	// Return the quaternion elements as x,y,z,w
	double y() const;	// -> using .x(),.y(),.z(),.w() may read more neatly though (it also avoids lots of 'if' clauses in member functions)
	double z() const;
	double w() const;
	vec3 vector() const; // Return the vector part of a quaternion as a vec3 type
	double scalar() const { return w();} // Return the scalar part of a quaternion
	colvec4 getAllCoeffs() const { return coeffs;}	// Return the quaternion elements as a colvec4
	SortConvention getSort() const {return sort;} // Return the sorting type

	// Properties
	Quaternion conjugate(void) const {return Quaternion(-vector(), scalar(), sort);} // Returns quaternion conjugate
	double norm() const { return sqrt(dot(coeffs,coeffs));} // Computes quaternion norm
	Quaternion inverse(void) const {return conjugate() / norm();} // Computes quaternion inverse (same thing as conjugate if norm==1, just slower)

	// Operations
	void normalize(){coeffs /= norm();} // Normalize the quaternion
	Quaternion product(const Quaternion& rhs) const; // Computes quaternion product
	Quaternion operator*(const Quaternion& rhs) const {return product(rhs);} // Quaternion product operator
	vec3 operator*(const vec3& rhs) const; // Quaternion rotation operator
	Quaternion operator-() const {return Quaternion(-(this->coeffs),sort);} 	// Unary negation
	Quaternion operator/(double s) const; 	// Quaternion scalar division operator.

	// Quaternion matrices
	mat44 matrix() const; // Quaternion matrix
	mat44 matrixReverse() const; // Quaternion matrix for switching quaternion order in product

	// SO(3) conversions
	mat33 getRot() const;	// Conversion to Rotation matrix

	// Stdout
	void print(){ cout << " " << w() << " <" << x() << ", " << y() << ", " << z() << ">" << endl;} // Prints the quaternion
};

// SO(3) conversions
double inner(const Quaternion& q1, const Quaternion& q2); // inner product between two quaternions
vec3 qLog(const Quaternion& q);	// Logarithm of a quaternion
Quaternion qExp(const vec3& aa, SortConvention sortType = VEC_SCALAR); // Exponential of an angular displacement (= a quaternion) [maybe the name of the function may be improved?]
vec3 angDiff(const Quaternion& qh, const Quaternion& q); // Rotation between two quaternions as a 3D vector (angular displacement)

// SLERP
Quaternion SLERP(const double& t, Quaternion qFrom, const Quaternion& qTo);

};

#endif /* QUATERNION_H */
