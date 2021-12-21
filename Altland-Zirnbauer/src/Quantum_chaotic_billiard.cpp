#include <iostream>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <ctime>
#include <chrono>
#include <random>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/QR>
#include "../include/Quantum_chaotic_billiard.h"
#include "../include/Auxiliary_Functions.h"
#define PI 3.14159265

using namespace std;

Quantum_chaotic_billiard::Quantum_chaotic_billiard(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2, int N1, int N2, int electron_hole_deg){
	Set_Setup(H, W, C1, C2, N1, N2, electron_hole_deg);
}

void Quantum_chaotic_billiard::Set_Setup(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2, int N1, int N2, int electron_hole_deg){
	
	_H = H;
	_W = W;
	_C1 = C1;
	_C2 = C2;
	_N1 = N1;
	_N2 = N2;
	_electron_hole_deg = electron_hole_deg;
}

void Quantum_chaotic_billiard::Calculate_Smatrix(){

	complex<double> number_2(2,0);
	complex<double> complex_identity(0,1);

	int N1 = _N1;
	int N2 = _N2;
	int n = N1+N2;

	MatrixXcd identityS = MatrixXcd::Identity(_W.cols(), _W.cols());

	MatrixXcd D(_H.rows(), _H.cols());

	D << (-_H + complex_identity*M_PI*_W*(_W.adjoint()));
	PartialPivLU<MatrixXcd> lu(D);
	MatrixXcd D_inv_W = lu.inverse()*_W;

	// Scatering Matrix //
	
	MatrixXcd S(_W.cols(), _W.cols());

	S << identityS - number_2*complex_identity*M_PI*(_W.adjoint())*D_inv_W;

//	MatrixXcd paulimatrix_y(2,2);
	
//	paulimatrix_y << 0, -complex_identity,
		         complex_identity, 0;

//	cout << "\nParticle-Hole symmetry constraint: \n" << (S-Kronecker_Product(MatrixXcd::Identity(n,n), paulimatrix_y)*S.conjugate()*Kronecker_Product(MatrixXcd::Identity(n,n), paulimatrix_y)).cwiseAbs() << endl;

//	cout << "\nTime-reversal symmetry constraint: (Only CI class) \n" << (S-S.transpose()).cwiseAbs() << endl;

	this -> _S = S;
}

void Quantum_chaotic_billiard::Calculate_G_and_P(){

	MatrixXcd ttdaga = _C1*_S*_C2*(_S.adjoint());

	MatrixXcd identityP = MatrixXcd::Identity(ttdaga.rows(), ttdaga.cols());

	_G = ttdaga.trace();
	_P = (ttdaga*(identityP - ttdaga)).trace();
}

complex<double> Quantum_chaotic_billiard::getG(){
	return this -> _G;
}

complex<double> Quantum_chaotic_billiard::getP(){
	return this -> _P;
}

void Quantum_chaotic_billiard::Calculate_Concurrence(){

	MatrixXcd t = _S.block(_electron_hole_deg * _N1, 0, _electron_hole_deg * _N2, _electron_hole_deg * _N1);

	MatrixXcd ttdaga = t*t.adjoint();

	VectorXcd eigenvalues_ttdaga = ttdaga.eigenvalues();

	double tau_1 = eigenvalues_ttdaga(0).real();
	double tau_2 = eigenvalues_ttdaga(_N1).real();

	_Concurrence = 2*(sqrt(tau_1*(1-tau_1)*tau_2*(1-tau_2))/(tau_1+tau_2-2*tau_1*tau_2));

	_Entanglement = -((1+sqrt(1-pow(_Concurrence,2)))/2)*log2((1+sqrt(1-pow(_Concurrence,2)))/2) - (1-(1+sqrt(1-pow(_Concurrence,2)))/2)*log2(1-(1+sqrt(1-pow(_Concurrence,2)))/2);
}

double Quantum_chaotic_billiard::getConcurrence(){
	return this -> _Concurrence;
}

double Quantum_chaotic_billiard::getEntanglement(){
	return this -> _Entanglement;
}

