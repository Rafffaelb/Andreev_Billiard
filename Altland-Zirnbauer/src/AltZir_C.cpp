#include <iostream>
#include "../include/AltZir_C.h"
#include "../include/Auxiliary_Functions.h"
#include <cmath>
#include <complex>
#include <random>
#include <ctime>
#include <chrono>
#include <fstream>
#include <string>

using namespace std;

AltZir_C::AltZir_C(double lambda, int num_steps, int spin_deg, int electron_hole_deg){

	this -> _lambda = lambda;
	this -> _num_steps = num_steps;
	this -> _spin_deg = spin_deg;
	this -> _electron_hole_deg = electron_hole_deg;
}

AltZir_C::~AltZir_C() {}

void AltZir_C::Create_W(MatrixXcd* W_pointer, int ress, int N1, int N2, double lambda, double y){

	MatrixXcd W1(ress,N1);
	MatrixXcd W2(ress,N2);

	W1.setZero(); W2.setZero();

	for (int j=1; j < ress+1; j++ ){
		for (int k=1; k < N1+1; k++){
			if (j == k){
				std::complex<double> aux(y*sqrt(lambda/M_PI), 0);
				W1(j-1,k-1) = aux;
			}
		}
	}

	for (int j=1; j < ress+1; j++ ){
		for (int k=1; k < N2+1; k++){
			if (j==k){
				std::complex<double> aux(y*sqrt(lambda/M_PI), 0);
				W2(j+N1-1,k-1) = aux;
		
			}
		}
	}

	MatrixXcd W_aux(ress, (N1+N2));
	W_aux.setZero();
	W_aux << W1, W2;

	MatrixXcd W(2*ress, 2*(N1+N2));
	W.setZero();

	W << Kronecker_Product(W_aux, MatrixXcd::Identity(2,2));

	*W_pointer = W;
}

void AltZir_C::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);

	int n = N1 + N2;

	MatrixXcd C1_aux(n,n);
	MatrixXcd C2_aux(n,n);

	C1_aux.block(0, 0, N1, N1) << identity1; C1_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C1_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C1_aux.block(N1, N1, N2, N2) << MatrixXcd::Zero(N2, N2);

	C2_aux.block(0, 0, N1, N1) << MatrixXcd::Zero(N1, N1); C2_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C2_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C2_aux.block(N1, N1, N2, N2) << identity2;

	MatrixXcd C1(_electron_hole_deg * n, _electron_hole_deg * n);
	MatrixXcd C2(_electron_hole_deg * n, _electron_hole_deg * n);

	C1 << Kronecker_Product(C1_aux, MatrixXcd::Identity(_electron_hole_deg, _electron_hole_deg));
	C2 << Kronecker_Product(C2_aux, MatrixXcd::Identity(_electron_hole_deg, _electron_hole_deg));

	*C1_pointer << C1;
	*C2_pointer << C2;
}

void AltZir_C::Create_H(MatrixXcd* H_pointer, int ress, double _lambda){

	complex<double> complex_identity(0,1);

	MatrixXcd paulimatrix_x(2,2);
	MatrixXcd paulimatrix_y(2,2);
	MatrixXcd paulimatrix_z(2,2);

	paulimatrix_x << 0, 1,
		         1, 0;

	paulimatrix_y << 0, -complex_identity,
		         complex_identity, 0;

	paulimatrix_z << 1, 0,
		      	 0, -1;


	auto seed = std::chrono::system_clock::now().time_since_epoch().count();

	std::normal_distribution<double> distribution(0.0, 1.0);
	
	std::default_random_engine generator(seed);

	MatrixXcd A(ress, ress); MatrixXcd B(ress, ress);
	MatrixXcd C(ress, ress); MatrixXcd D(ress, ress);

	A.setZero(); B.setZero();
	C.setZero(); D.setZero();

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			B(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			C(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			D(i-1,j-1) = aux;
		}
	}

	MatrixXcd H_0_aux(ress, ress); MatrixXcd H_1_aux(ress, ress);
       	MatrixXcd H_2_aux(ress, ress); MatrixXcd H_3_aux(ress, ress);

	H_0_aux.setZero(); H_1_aux.setZero();
       	H_2_aux.setZero(); H_3_aux.setZero();

	for (int i = 1; i < ress + 1; i++){
		H_0_aux(i-1,i-1) = 0;
		H_1_aux(i-1,i-1) = (_lambda*(1/(2*sqrt(3*ress))))*B(i-1,i-1);
		H_2_aux(i-1,i-1) = (_lambda*(1/(2*sqrt(3*ress))))*C(i-1,i-1);
		H_3_aux(i-1,i-1) = (_lambda*(1/(2*sqrt(3*ress))))*D(i-1,i-1);
		for(int j = i + 1; j < ress + 1; j++){
			H_0_aux(i-1,j-1) = (_lambda*(1/(sqrt(4*ress))))*A(i-1,j-1);
			H_1_aux(i-1,j-1) = (_lambda*(1/(sqrt(4*ress))))*B(i-1,j-1);
			H_2_aux(i-1,j-1) = (_lambda*(1/(sqrt(4*ress))))*C(i-1,j-1);
			H_3_aux(i-1,j-1) = (_lambda*(1/(sqrt(4*ress))))*D(i-1,j-1);
		}
	}
	
	MatrixXcd H_0(ress, ress); MatrixXcd H_1(ress, ress);
       	MatrixXcd H_2(ress, ress); MatrixXcd H_3(ress, ress);

	H_0.setZero(); H_1.setZero();
	H_2.setZero(); H_3.setZero();

	H_0 << H_0_aux - H_0_aux.transpose();
	H_1 << H_1_aux + H_1_aux.transpose();
	H_2 << H_2_aux + H_2_aux.transpose();
	H_3 << H_3_aux + H_3_aux.transpose();


	MatrixXcd H(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * ress);
	H.setZero();

	H << complex_identity*Kronecker_Product(H_0, MatrixXcd::Identity(_electron_hole_deg, _electron_hole_deg)) + (Kronecker_Product(H_1, paulimatrix_x) + Kronecker_Product(H_2, paulimatrix_y) + Kronecker_Product(H_3, paulimatrix_z));
	*H_pointer = H;	
}

void AltZir_C::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps){
	std::ofstream output_G("Data_Analysis/Channel/Andreev_G_C_Channel.txt");
	std::ofstream output_P("Data_Analysis/Channel/Andreev_P_C_Channel.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 10; j++){
			if (j == 9){
				output_G << G(i,j).real() << std::endl;
				output_P << P(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
				output_P << P(i,j).real() << "\t";
			}
		}
	}
}
