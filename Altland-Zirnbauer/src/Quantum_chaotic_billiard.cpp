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
