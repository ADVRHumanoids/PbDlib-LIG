/**
 *
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

#include "pbdlib/quaternion.h"

namespace pbdlib
{

// -- Constructors --
// Default constructor - initializes identity quaternion
Quaternion::Quaternion(SortConvention sortInit)
{
	sort = sortInit;
	if(sort==VEC_SCALAR)
	{
		coeffs.subvec(0,2) = zeros(3);
		coeffs(3)   = 1;
	}
	else
	{
		coeffs(0)   = 1;
		coeffs.subvec(1,3) = zeros(3);
	}
}

// Initialize a quaternion from vec3 and scalar
Quaternion::Quaternion(const vec3& vec, double scalar, SortConvention sortInit)
{
	sort = sortInit;
	if(sort==VEC_SCALAR)
	{
		coeffs.subvec(0,2) = vec;
		coeffs(3) = scalar;
	}
	else
	{
		coeffs(0) = scalar;
		coeffs.subvec(1,3) = vec;
	}
}

// Initialize a quaternion from a vec4
Quaternion::Quaternion(const vec4& vec, SortConvention sortInit)
{
	sort = sortInit;
	coeffs = vec;
}

// Initialize a quaternion from an array of double
Quaternion::Quaternion(const double* array, SortConvention sortInit)
{
	sort = sortInit;
	for (uint i = 0; i < 4; i++)
	{
		coeffs[i] = array[i];
	}
}

// Initialize a quaternion from 4 scalars of type double
Quaternion::Quaternion(double q1, double q2, double q3, double q4,  SortConvention sortInit)
{
	sort = sortInit;
	coeffs[0] = q1;
	coeffs[1] = q2;
	coeffs[2] = q3;
	coeffs[3] = q4;
}

// Initialize a quaternion from an angle and an axis
Quaternion::Quaternion(double angle, vec3 axis,  SortConvention sortInit)
{
	sort = sortInit;

	double scalar = cos(angle/2.0);
	vec3 vector = axis*sin(angle/2.0);

	if(sort==VEC_SCALAR)
	{
		coeffs.subvec(0,2) = vector;
		coeffs(3) = scalar;
	}
	else
	{
		coeffs(0) = scalar;
		coeffs.subvec(1,3) = vector;
	}
}

// Initialize a quaternion from a vec3 whose direction is the rotation axis and norm is the angle of rotation
// -> this is not working because apparently it conflicts with the one that uses vec4 as type of the first argument
//Quaternion::Quaternion(const vec3& angleAxis,  SortConvention sortInit)
//{
//	sort = sortInit;
//
//	double angle = arma::norm(angleAxis);
//	double scalar = cos(angle/2.0);
//	vec3 vector = arma::normalise(angleAxis)*sin(angle/2.0);
//
//	if(sort==VEC_SCALAR)
//	{
//		coeffs.subvec(0,2) = vector;
//		coeffs(3) = scalar;
//	}
//	else
//	{
//		coeffs(0) = scalar;
//		coeffs.subvec(1,3) = vector;
//	}
//}


// -- Coefficients --
// Return the quaternion elements as x,y,z,w
// -> using .x(),.y(),.z(),.w() may read more neatly though (it also avoids lots of 'if' clauses down the road)
double Quaternion::x() const
{
	if(sort==VEC_SCALAR) return coeffs[0];
	else return coeffs[1];
}
double Quaternion::y() const
{
	if(sort==VEC_SCALAR) return coeffs[1];
	else return coeffs[2];
}
double Quaternion::z() const
{
	if(sort==VEC_SCALAR) return coeffs[2];
	else return coeffs[3];
}
double Quaternion::w() const
{
	if(sort==VEC_SCALAR) return coeffs[3];
	else return coeffs[0];
}

// Return the vector part of a quaternion as a vec3 type
vec3 Quaternion::vector() const
{
	vec3 vPart;
	vPart[0] = x();
	vPart[1] = y();
	vPart[2] = z();
	return vPart;
}


//	-- Operations --
// Computes quaternion product
Quaternion Quaternion::product(const Quaternion& rhs) const
{
	if(sort!=rhs.sort)
	{
		cout << "WARNING: Product between two quaternions using different conventions on the order of their elements." <<endl;
	}
	return Quaternion(this->matrix()*rhs.coeffs,sort); // Note that the quaternion matrix was used here, since its easier to implement (also faster I think)
}

// Rotation of a 3D point using quaternion product
vec3 Quaternion::operator*(const vec3& rhs) const
{
	Quaternion qOut = (*this)*Quaternion(0,rhs[0],rhs[1],rhs[2],this->getSort())*this->conjugate(); // rotation operation: q*Quaternion([0 px py pz])*q.inv
	return qOut.vector();
}

//	/**
//	 * @brief Quaternion scalar product operator.
//	 * @param s A scalar by which to multiply all components
//	 * of this quaternion.
//	 * @return The quaternion (*this) * s.
//	 */
//	Quaternion operator*(double s) const {
//	}
//		return Quaternion(complex()*s, real()*s);
//
//	/**
//	 * @brief Produces the sum of this quaternion and rhs.
//	 */
//	Quaternion operator+(const Quaternion& rhs) const {
//		return Quaternion(x()+rhs.x(), y()+rhs.y(), z()+rhs.z(), w()+rhs.w());
//	}
//
//	/**
//	 * @brief Produces the difference of this quaternion and rhs.
//	 */
//	Quaternion operator-(const Quaternion& rhs) const {
//		return Quaternion(x()-rhs.x(), y()-rhs.y(), z()-rhs.z(), w()-rhs.w());
//	}

// Quaternion scalar division operator.
Quaternion Quaternion::operator/(double s) const
{
	if (s == 0) cout << "Dividing quaternion by 0." << std::endl;
	return Quaternion(vector()/s, scalar()/s, sort);
}


// -- Quaternion matrices --

// Returns a matrix representation of (*this) quaternion for implementing
// the quaternion product through matrix-vector product, using the same order
// of the quaternions:
//
// 		q_new = (*this) * q
// 			<=>
//		q_new = this->matrix * q
mat44 Quaternion::matrix() const
{
	mat44 m;
	if(sort==VEC_SCALAR)
	{
		//			double m[16]= {
		//					 w(),-z(), y(), x(),
		//				 	 z(), w(),-x(), y(),
		//					-y(), x(), w(), z(),
		//					-x(),-y(),-z(), w()
		//			};
		m <<  w() << -z() <<  y() << x() << endr
	      <<  z() <<  w() << -x() << y() << endr
		  << -y() <<  x() <<  w() << z() << endr
		  << -x() << -y() << -z() << w() << endr;
	}
	else
	{
		//			double m[16]= {
		//					 w(), x(), y(),	z(),
		//					-x(), w(),-z(),	y(),
		//					-y(), z(), w(),-x(),
		//					-z(),-y(), x(), w()
		//			};
		m <<  w() << -x() << -y() << -z() << endr
		  <<  x() <<  w() << -z() <<  y() << endr
		  <<  y() <<  z() <<  w() << -x() << endr
		  <<  z() << -y() <<  x() <<  w() << endr;
	}
	return m;

}

// Returns a matrix representation of (*this) quaternion for implementing
// the quaternion product through matrix-vector product, using the *reversed* order
// of the quaternions:
//
// 		q_new = q * (*this)
// 			<=>
//		q_new = this->matrix * q
//
mat44 Quaternion::matrixReverse() const
{
	mat44 m;
	if(sort==VEC_SCALAR)
	{
		m <<  w() <<  z() << -y() <<  x() << endr
		  << -z() <<  w() <<  x() <<  y() << endr
		  <<  y() << -x() <<  w() <<  z() << endr
		  << -x() << -y() << -z() <<  w() << endr;
	}
	else
	{
		m <<  w() << -x() << -y() << -z() << endr
		  <<  x() <<  w() <<  z() << -y() << endr
		  <<  y() << -z() <<  w() <<  x() << endr
		  <<  z() <<  y() << -x() <<  w() << endr;
	}
	return m;
}

mat33 Quaternion::getRot() const
{
	mat33 rotMatrix;
	rotMatrix << 	1-2*(y()*y() + z()*z())	   << 		2*(x()*y()-z()*w()) 	<< 		2*(x()*z()+y()*w()) 	 << endr
			  <<	 2*(x()*y()+z()*w())	   << 1-2*(pow(x(),2) + pow(z(),2)) << 		2*(y()*z()-x()*w()) 	 << endr
			  <<	 2*(x()*z()-y()*w())	   << 		2*(y()*z()+x()*w()) 	<< 1-2*(pow(x(),2) + pow(y(),2)) << endr;
	return rotMatrix;
}

// -- SO(3) conversions --
// inner product between two quaternions
double inner(const Quaternion& q1, const Quaternion& q2)
{
	return dot(q1.getAllCoeffs(),q2.getAllCoeffs());
};

// Logarithm of a quaternion
vec3 qLog(const Quaternion& q)
{
	return acos(q.scalar())*q.vector()/arma::norm(q.vector());
};

// Exponential of an angular displacement (= a quaternion) [maybe the name of the function may be improved?]
Quaternion qExp(const vec3& aa, SortConvention sortType)
{
	return Quaternion(2*arma::norm(aa),arma::normalise(aa),sortType);
};

// Rotation between two quaternions as a 3D vector (angular displacement)
vec3 angDiff(const Quaternion& qh, const Quaternion& q)
{
	// consider the closest to qh (q or -q)
	if(inner(qh,q)<0)
		return 2*qLog(qh*Quaternion(-q.getAllCoeffs(),q.getSort()).inverse());
	else
		return 2*qLog(qh*q.inverse());
}

// -- SLERP --
Quaternion SLERP(const double& t, Quaternion qFrom, const Quaternion& qTo)
{
	if(inner(qTo,qFrom)<0)
		qFrom = Quaternion(-qFrom.getAllCoeffs());

	vec3 axis = t*qLog(qTo*qFrom.inverse());
	return pbdlib::qExp(axis,qFrom.getSort())*qFrom;	//add a reference to this expression?
}
//	/**
//	 * @brief Returns a matrix representation of this
//	 * quaternion for right multiplication.
//	 *
//	 * Specifically this is the matrix such that:
//	 *
//	 * q.vector().transpose() * this->matrix() = (q *
//	 * (*this)).vector().transpose() for any quaternion q.
//	 *
//	 * Note that this is @e NOT the rotation matrix that may be
//	 * represented by a unit quaternion.
//	 */
//	TMatrix4 rightMatrix() const {
//		double m[16] = {
//				+w(), -z(),  y(), -x(),
//				+z(),  w(), -x(), -y(),
//				-y(),  x(),  w(), -z(),
//				+x(),  y(),  z(),  w()
//		};
//		return TMatrix4(m);
//	}
//
//	/**
//	 * @brief Returns this quaternion as a 4-vector.
//	 *
//	 * This is simply the vector [x y z w]<sup>T</sup>
//	 */
//	TVector4 vector() const { return TVector4(mData); }
//

//
//	/**
//	 * @brief Computes the rotation matrix represented by a unit
//	 * quaternion.
//	 *
//	 * @note This does not check that this quaternion is normalized.
//	 * It formulaically returns the matrix, which will not be a
//	 * rotation if the quaternion is non-unit.
//	 */
//	TMatrix3 rotationMatrix() const {
//		double m[9] = {
//				1-2*y()*y()-2*z()*z(), 2*x()*y() - 2*z()*w(), 2*x()*z() + 2*y()*w(),
//				2*x()*y() + 2*z()*w(), 1-2*x()*x()-2*z()*z(), 2*y()*z() - 2*x()*w(),
//				2*x()*z() - 2*y()*w(), 2*y()*z() + 2*x()*w(), 1-2*x()*x()-2*y()*y()
//		};
//		return TMatrix3(m);
//	}
//
//	/**
//	 * @brief Returns the scaled-axis representation of this
//	 * quaternion rotation.
//	 */
//	TVector3 scaledAxis(void) const {
//		double w[3];
//		HeliMath::scaled_axis_from_quaternion(w, mData);
//		return TVector3(w);
//	}
//
//	/**
//	 * @brief Sets quaternion to be same as rotation by scaled axis w.
//	 */
//	void scaledAxis(const TVector3& w) {
//		double theta = w.norm();
//		if (theta > 0.0001) {
//			double s = sin(theta / 2.0);
//			TVector3 W(w / theta * s);
//			mData[0] = W[0];
//			mData[1] = W[1];
//			mData[2] = W[2];
//			mData[3] = cos(theta / 2.0);
//		} else {
//			mData[0]=mData[1]=mData[2]=0;
//			mData[3]=1.0;
//		}
//	}
//
//	/**
//	 * @brief Returns a vector rotated by this quaternion.
//	 *
//	 * Functionally equivalent to:  (rotationMatrix() * v)
//	 * or (q * Quaternion(0, v) * q.inverse()).
//	 *
//	 * @warning conjugate() is used instead of inverse() for better
//	 * performance, when this quaternion must be normalized.
//	 */
//	TVector3 rotatedVector(const TVector3& v) const {
//		return (((*this) * Quaternion(v, 0)) * conjugate()).complex();
//	}
//
//
//
//	/**
//	 * @brief Computes the quaternion that is equivalent to a given
//	 * euler angle rotation.
//	 * @param euler A 3-vector in order:  roll-pitch-yaw.
//	 */
//	void euler(const TVector3& euler) {
//		double c1 = cos(euler[2] * 0.5);
//		double c2 = cos(euler[1] * 0.5);
//		double c3 = cos(euler[0] * 0.5);
//		double s1 = sin(euler[2] * 0.5);
//		double s2 = sin(euler[1] * 0.5);
//		double s3 = sin(euler[0] * 0.5);
//
//		mData[0] = c1*c2*s3 - s1*s2*c3;
//		mData[1] = c1*s2*c3 + s1*c2*s3;
//		mData[2] = s1*c2*c3 - c1*s2*s3;
//		mData[3] = c1*c2*c3 + s1*s2*s3;
//	}
//
//	/** @brief Returns an equivalent euler angle representation of
//	 * this quaternion.
//	 * @return Euler angles in roll-pitch-yaw order.
//	 */
//	TVector3 euler(void) const {
//		TVector3 euler;
//		const static double PI_OVER_2 = M_PI * 0.5;
//		const static double EPSILON = 1e-10;
//		double sqw, sqx, sqy, sqz;
//
//		// quick conversion to Euler angles to give tilt to user
//		sqw = mData[3]*mData[3];
//		sqx = mData[0]*mData[0];
//		sqy = mData[1]*mData[1];
//		sqz = mData[2]*mData[2];
//
//		euler[1] = asin(2.0 * (mData[3]*mData[1] - mData[0]*mData[2]));
//		if (PI_OVER_2 - fabs(euler[1]) > EPSILON) {
//			euler[2] = atan2(2.0 * (mData[0]*mData[1] + mData[3]*mData[2]),
//					sqx - sqy - sqz + sqw);
//			euler[0] = atan2(2.0 * (mData[3]*mData[0] + mData[1]*mData[2]),
//					sqw - sqx - sqy + sqz);
//		} else {
//			// compute heading from local 'down' vector
//			euler[2] = atan2(2*mData[1]*mData[2] - 2*mData[0]*mData[3],
//					2*mData[0]*mData[2] + 2*mData[1]*mData[3]);
//			euler[0] = 0.0;
//
//			// If facing down, reverse yaw
//			if (euler[1] < 0)
//				euler[2] = M_PI - euler[2];
//		}
//		return euler;
//	}
//
//	/**
//	 * @brief Computes a special representation that decouples the Z
//	 * rotation.
//	 *
//	 * The decoupled representation is two rotations, Qxy and Qz,
//	 * so that Q = Qxy * Qz.
//	 */
//	void decoupleZ(Quaternion* Qxy, Quaternion* Qz) const {
//		TVector3 ztt(0,0,1);
//		TVector3 zbt = this->rotatedVector(ztt);
//		TVector3 axis_xy = ztt.cross(zbt);
//		double axis_norm = axis_xy.norm();
//
//		double axis_theta = acos(HeliMath::saturate(zbt[2], -1,+1));
//		if (axis_norm > 0.00001) {
//			axis_xy = axis_xy * (axis_theta/axis_norm); // limit is *1
//		}
//
//		Qxy->scaledAxis(axis_xy);
//		*Qz = (Qxy->conjugate() * (*this));
//	}
//
//	/**
//	 * @brief Returns the quaternion slerped between this and q1 by fraction 0 <= t <= 1.
//	 */
//	Quaternion slerp(const Quaternion& q1, double t) {
//		return slerp(*this, q1, t);
//	}
//
//	/// Returns quaternion that is slerped by fraction 't' between q0 and q1.
//	static Quaternion slerp(const Quaternion& q0, const Quaternion& q1, double t) {
//
//		double omega = acos(HeliMath::saturate(q0.mData[0]*q1.mData[0] +
//				q0.mData[1]*q1.mData[1] +
//				q0.mData[2]*q1.mData[2] +
//				q0.mData[3]*q1.mData[3], -1,1));
//		if (fabs(omega) < 1e-10) {
//			omega = 1e-10;
//		}
//		double som = sin(omega);
//		double st0 = sin((1-t) * omega) / som;
//		double st1 = sin(t * omega) / som;
//
//		return Quaternion(q0.mData[0]*st0 + q1.mData[0]*st1,
//				q0.mData[1]*st0 + q1.mData[1]*st1,
//				q0.mData[2]*st0 + q1.mData[2]*st1,
//				q0.mData[3]*st0 + q1.mData[3]*st1);
//	}
//
//	/**
//	 * @brief Returns pointer to the internal array.
//	 *
//	 * Array is in order x,y,z,w.
//	 */
//	double* row(uint32_t i) { return mData + i; }
//	// Const version of the above.
//	const double* row(uint32_t i) const { return mData + i; }

};


