/**
Copyright (C) 2014, Danilo Bruno, Sylvain Calinon

This file is part of PbDLib (Programming-by-demonstration C++ Library).

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

#include "pbdlib/lqr.h"
#include <iostream>

namespace pbdlib
{

LQR::LQR(mat _A, mat _B, double _dt)
{
	this->A = _A;
	this->B = _B;
	this->dt = _dt;

	this->nbVar = A.n_rows;
	this->nbCtrl = B.n_cols;

	this->prob_set = 0;

}

bool LQR::evaluate_gains_finiteHorizon(mat SFinal,colvec dFinal){

	if (this->prob_set == 0) {return 1;};

	mat DTarget = LQR::diff(this->Target);

	this->S[this->nbData-1] = SFinal;
	this->d.col(this->nbData-1) = dFinal;

	mat invR = this->R.i();

	for (int t=this->nbData-2;t>=0;t--){
		this->S[t] =  this->S[t+1] + this->dt * (this->A.t()*this->S[t+1] + this->S[t+1]*this->A - this->S[t+1] * this->B * invR * this->B.t() * this->S[t+1] + this->Q[t+1]); //See Eq. (5.1.11) in doc/TechnicalReport.pdf
		this->d.col(t) =  this->d.col(t+1) - this->dt * (this->S[t+1]*this->B*invR * this->B.t() - this->A.t()) * this->d.col(t+1) + this->dt * this->S[t+1] * DTarget.col(t+1)  - this->dt * (this->S[t+1] * this->A * this->Target.col(t+1)); //See Eq. (5.1.29) in doc/TechnicalReport.pdf
	}
	for (int t=0;t<nbData;t++){
		this->L[t] = invR*this->B.t()*this->S[t]; //See Eq. (5.1.30) in doc/TechnicalReport.pdf
		this->M.col(t) = invR*this->B.t()*this->d.col(t); //See Eq. (5.1.30) in doc/TechnicalReport.pdf
	}
	return 0;
};

bool LQR::evaluate_gains_infiniteHorizon(){

	if (this->prob_set == 0) {return 1;};
	mat DTarget = LQR::diff(this->Target);
	mat invR = this->R.i();
	for (int t=0;t < this->nbData;t++){
		this->S[t] = this->solveAlgebraicRiccati(this->A,this->B,this->Q[t],this->R); //See Sec. (5.2) in doc/TechnicalReport.pdf
		this->d.col(t) = (this->S[t]*this->B*invR*this->B.t() - this->A.t()).i() * (this->S[t]*DTarget.col(t) - this->S[t]*this->A*this->Target.col(t)); //See Eq. (5.1.29) in doc/TechnicalReport.pdf
		this->L[t] = invR*this->B.t()*this->S[t]; //See Eq. (5.1.30) in doc/TechnicalReport.pdf
		//std::cout << this->L[t] << std::endl;
		this->M.col(t) = invR*this->B.t()*this->d.col(t); //See Eq. (5.1.30) in doc/TechnicalReport.pdf
	}
	return 0;
};

mat LQR::diff(mat V){

	mat DV(V.n_rows,V.n_cols);

	if (V.n_cols==1) { DV = zeros<mat>(V.n_rows,V.n_cols);}
	else
	{
		DV.col(0) = (V.col(1) - V.col(0))/this->dt;
		for (int i=1;i<V.n_cols-1;i++){
			DV.col(i) = (V.col(i+1)-V.col(i-i))/(2*this->dt);
		}
		DV.col(V.n_cols-1) =  (V.col(V.n_cols-1) - V.col(V.n_cols-2))/this->dt;
	}
	return DV;
}

mat LQR::solveAlgebraicRiccati(mat A, mat B, mat Q, mat R) //See Sec. (5.2) in doc/TechnicalReport.pdf
{
	int n = A.n_rows;

	mat G = B*R.i()*B.t();
	mat Z(2*n,2*n); //See Eq. (5.2.3) in doc/TechnicalReport.pdf
	Z(span(0,n-1),span(0,n-1)) = A;
	Z(span(n,2*n-1),span(0,n-1)) = -Q;
	Z(span(0,n-1),span(n,2*n-1)) = -G;
	Z(span(n,2*n-1),span(n,2*n-1)) = -A.t();

	// Using schur decomposition
	/*mat U(2*n,2*n);
	mat T(2*n,2*n);
	auxlib::schur_dec(U,T,Z);	//Missing ordered schur...
	//Ordered schur decomposition from lapack to add
	mat S = U(span(0,n-1),span(0,n-1)).t().i()*U(span(n,2*n-1),span(0,n-1)).t();
	*/

	//Using diagonalization (See Eq. (5.2.4) in doc/TechnicalReport.pdf)
	cx_mat U(2*n,n);
	cx_vec dd(2*n);
	cx_mat V(2*n,2*n);

	arma::eig_gen(dd,V,Z);
	int i = 0;
	for(int j=0;j<2*n;j++){
		if (real(dd(j)) < 0) {
			U.col(i) = V.col(j);
			i++;
		}
	}
	mat S1 = zeros(n,n);
	cx_mat Sc = U(span(0,n-1),span(0,n-1)).t().i()*U(span(n,2*n-1),span(0,n-1)).t(); //See Eq. (5.2.5) in doc/TechnicalReport.pdf
	mat S = real(Sc);
	return S;
};


bool LQR::setProblem(mat _R, std::vector<mat> _Q, mat _Target)
{
	this->R = _R;
	this->Q = _Q;
	this->Target = _Target;

	this->nbData = _Q.size();

	for (int i=0;i<this->nbData;i++){
		this->S.push_back(zeros(this->nbVar,this->nbVar));
		this->L.push_back(zeros(this->nbCtrl,this->nbVar));
	}
	this->d = zeros<mat>(this->nbVar,this->nbData);
	this->M = zeros<mat>(this->nbCtrl,this->nbData);

	this->prob_set = 1;
}


} //end of pbdlib namespace


