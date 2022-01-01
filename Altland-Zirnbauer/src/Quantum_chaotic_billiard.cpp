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

MatrixXcd Create_Unitary_Random_Matrix();
complex<double> Calculate_Noise(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R);
complex<double> Calculate_Noise_Fixed_Base(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R);

void Quantum_chaotic_billiard::Calculate_Bell_Parameter(){

	MatrixXcd r, t, U_L(2,2), U_R(2,2), U_Lprime(2,2), U_Rprime(2,2), paulimatrix_z(2,2);
	complex<double> C_a_b, C_a_bprime, C_aprime_b, C_aprime_bprime;

	paulimatrix_z << 1, 0,
			 0, -1;

	t = _S.block(_electron_hole_deg * _N1, 0, _electron_hole_deg * _N2, _electron_hole_deg * _N1);
	r = _S.block(0, 0, _electron_hole_deg * _N1, _electron_hole_deg * _N1);

	U_L = Create_Unitary_Random_Matrix();
	U_R = Create_Unitary_Random_Matrix();
	U_Lprime = Create_Unitary_Random_Matrix();
	U_Rprime = Create_Unitary_Random_Matrix();

	C_a_b = Calculate_Noise(r, t, U_L, U_R);
	C_a_bprime = Calculate_Noise(r, t, U_L, U_Rprime);
	C_aprime_b = Calculate_Noise(r, t, U_Lprime, U_R);
	C_aprime_bprime = Calculate_Noise(r, t, U_Lprime, U_Rprime);

	_Bell_Parameter = abs(((C_a_b + C_aprime_b + C_a_bprime - C_aprime_bprime).real()));
	_Bell_Parameter_Dephase = 2*abs((paulimatrix_z*r*t.adjoint()*paulimatrix_z*t*r.adjoint()).trace())/((r.adjoint()*r*t.adjoint()*t).trace()).real();
}

void Quantum_chaotic_billiard::Calculate_Bell_Parameter_Fixed_Base(){

	complex<double> complex_identity(0,1);
	int chiral_deg = 2;
	
	MatrixXcd r, t, U_L(2,2), U_R(2,2), U_Lprime(2,2), U_Rprime(2,2), paulimatrix_x(2,2), paulimatrix_y(2,2), paulimatrix_z(2,2);
	complex<double> C_a_b, C_a_bprime, C_aprime_b, C_aprime_bprime;

	paulimatrix_x << 0, 1,
		      	 1, 0;
	
	paulimatrix_y << 0, -complex_identity,
		      	 complex_identity, 0;

	paulimatrix_z << 1, 0,
			 0, -1;

	t = _S.block(chiral_deg * _N1, 0, chiral_deg * _N2, chiral_deg * _N1);
	r = _S.block(0, 0, chiral_deg * _N1, chiral_deg * _N1);

	MatrixXd T_y_theta(3,3);
	
	double degree = 0;
	double theta = degree*PI/180;

	T_y_theta << cos(theta), 0, sin(theta),
	  	     0, 1, 0,
	   	     -sin(theta), 0, cos(theta);	     

	Vector3d Alice_vector(0, 0, 1);
	Vector3d Alice_vector_prime(1, 0, 0);
	Vector3d Bob_vector(-1/(sqrt(2)), 0, -1/(sqrt(2)));
	Vector3d Bob_vector_prime(-1/(sqrt(2)), 0, 1/(sqrt(2)));
	
	Alice_vector = T_y_theta*Alice_vector;
	Alice_vector_prime = T_y_theta*Alice_vector_prime;
	Bob_vector = T_y_theta*Bob_vector;
	Bob_vector_prime = T_y_theta*Bob_vector_prime;

	U_L = product_vector_matrix(Alice_vector, paulimatrix_x, paulimatrix_y, paulimatrix_z);
	U_Lprime = product_vector_matrix(Alice_vector_prime, paulimatrix_x, paulimatrix_y, paulimatrix_z);
	U_R = product_vector_matrix(Bob_vector, paulimatrix_x, paulimatrix_y, paulimatrix_z);
	U_Rprime = product_vector_matrix(Bob_vector_prime, paulimatrix_x, paulimatrix_y, paulimatrix_z);

	C_a_b = Calculate_Noise_Fixed_Base(r, t, U_L, U_R);
	C_a_bprime = Calculate_Noise_Fixed_Base(r, t, U_L, U_Rprime);
	C_aprime_b = Calculate_Noise_Fixed_Base(r, t, U_Lprime, U_R);
	C_aprime_bprime = Calculate_Noise_Fixed_Base(r, t, U_Lprime, U_Rprime);

	_Bell_Parameter = abs(((C_a_b + C_aprime_b + C_a_bprime - C_aprime_bprime).real()));
}

double Quantum_chaotic_billiard::getBell_Parameter(){
	return this -> _Bell_Parameter;
}

double Quantum_chaotic_billiard::getBell_Parameter_Dephase(){
	return this -> _Bell_Parameter_Dephase;
}


MatrixXcd Create_Unitary_Random_Matrix(){

	// Function to Create Unitary Random Matrix distributed with Haar Measure //

	complex<double> complex_identity(0,1);

	MatrixXcd Z(2,2), A(2,2), B(2,2), Q(2,2), R(2,2);
	MatrixXcd Diag_R, Delta;

	ColPivHouseholderQR<MatrixXcd> qr(Z.rows(), Z.cols());

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::normal_distribution<double> distribution(0.0,1.0);
	std::default_random_engine generator(seed);
	
	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < 3; i++){
			for (int j = 1; j < 3; j++){
			double aux = distribution(generator);
			B(i-1,j-1) = aux;
		}
	}

	Z = (1/sqrt(2))*(A+complex_identity*B);

	qr.compute(Z);

	Q = qr.householderQ().setLength(qr.nonzeroPivots());
	R = qr.matrixR().template triangularView<Upper>();

	Diag_R = R.diagonal().matrix().asDiagonal();

	Delta = Diag_R*Diag_R.cwiseAbs().inverse();

	Q = Q*Delta;
       
	return Q; // Unitary random matrix distributed with Haar measure //
}

complex<double> Calculate_Noise(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R){
	
	complex<double> C, Const_Norm;

	MatrixXcd paulimatrix_z(2,2), a_dot_sigma(2,2), b_dot_sigma(2,2);

	paulimatrix_z << 1, 0,
			 0, -1;

	a_dot_sigma = U_L.adjoint()*paulimatrix_z*U_L;
	b_dot_sigma = U_R.adjoint()*paulimatrix_z*U_R;

	Const_Norm = ((r*r.adjoint()).trace())*((t*t.adjoint()).trace())-(r*t.adjoint()*t*r.adjoint()).trace();

	C = ((((a_dot_sigma)*r*r.adjoint()).trace())*(((b_dot_sigma)*t*t.adjoint()).trace()) - (((a_dot_sigma)*r*t.adjoint()*(b_dot_sigma)*t*r.adjoint()).trace()))/Const_Norm;

	return C;
}

complex<double> Calculate_Noise_Fixed_Base(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R){
	
	complex<double> C, Const_Norm;

	MatrixXcd paulimatrix_z(2,2), a_dot_sigma(2,2), b_dot_sigma(2,2);

	paulimatrix_z << 1, 0,
			 0, -1;

	a_dot_sigma = U_L;
	b_dot_sigma = U_R;

	Const_Norm = ((r*r.adjoint()).trace())*((t*t.adjoint()).trace())-(r*t.adjoint()*t*r.adjoint()).trace();

	C = ((((a_dot_sigma)*r*r.adjoint()).trace())*(((b_dot_sigma)*t*t.adjoint()).trace()) - (((a_dot_sigma)*r*t.adjoint()*(b_dot_sigma)*t*r.adjoint()).trace()))/Const_Norm;

	return C;
}
