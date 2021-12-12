#include <iostream>
#include "../include/AltZir_DIII.h"
#include "../include/Auxiliary_Functions.h"
#include <cmath>
#include <complex>
#include <random>
#include <ctime>
#include <chrono>
#include <fstream>
#include <string>

using namespace std;

AltZir_DIII::AltZir_DIII(double lambda, int num_steps, int spin_deg, int electron_hole_deg){

	this -> _lambda = lambda;
	this -> _num_steps = num_steps;
	this -> _spin_deg = spin_deg;
	this -> _electron_hole_deg = electron_hole_deg;
}

AltZir_DIII::~AltZir_DIII() {}

void AltZir_DIII::Create_W(MatrixXcd* W_pointer, int ress, int N1, int N2, double lambda, double y){

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

void AltZir_DIII::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);
	
	MatrixXcd C1_aux((N1+N2), (N1+N2));
	MatrixXcd C2_aux((N1+N2), (N1+N2));

	C1_aux.block(0, 0, N1, N1) << identity1; C1_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C1_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C1_aux.block(N1, N1, N2, N2) << MatrixXcd::Zero(N2, N2);

	C2_aux.block(0, 0, N1, N1) << MatrixXcd::Zero(N1, N1); C2_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C2_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C2_aux.block(N1, N1, N2, N2) << identity2;

	MatrixXcd C1_without_spin(_electron_hole_deg * (N1+N2), _electron_hole_deg * (N1+N2));
	MatrixXcd C2_without_spin(_electron_hole_deg * (N1+N2), _electron_hole_deg * (N1+N2));

	C1_without_spin << Kronecker_Product(C1_aux, MatrixXcd::Identity(_electron_hole_deg, _electron_hole_deg));
	C2_without_spin << Kronecker_Product(C2_aux, MatrixXcd::Identity(_electron_hole_deg, _electron_hole_deg));

	MatrixXcd C1(_electron_hole_deg * _spin_deg * (N1+N2), _electron_hole_deg * _spin_deg * (N1+N2));
	MatrixXcd C2(_electron_hole_deg * _spin_deg * (N1+N2), _electron_hole_deg * _spin_deg * (N1+N2));

	C1 << Kronecker_Product(C1_without_spin, MatrixXcd::Identity(_spin_deg, _spin_deg));
	C2 << Kronecker_Product(C2_without_spin, MatrixXcd::Identity(_spin_deg, _spin_deg));

	*C1_pointer << C1;
	*C2_pointer << C2;
}

void AltZir_DIII::Create_H(MatrixXcd* H_pointer, int ress, double _lambda){

	complex<double> complex_identity(0,1);

	MatrixXcd paulimatrix_x(2,2);
	MatrixXcd paulimatrix_z(2,2);

	paulimatrix_x << 0, 1,
		  	 1, 0;

	paulimatrix_z << 1, 0,
			 0, -1;

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::normal_distribution<double> distribution(0.0, 1.0);
	std::default_random_engine generator(seed);

	MatrixXcd A(ress, ress);
	A.setZero();

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	MatrixXcd H1_aux(ress, ress); MatrixXcd H3_aux(ress, ress);
	H1_aux.setZero(); H3_aux.setZero();

	for (int i = 1; i < ress + 1; i++){
		H1_aux(i-1,i-1) = 0;
		H3_aux(i-1,i-1) = 0;
		for (int j = i + 1; j < ress + 1; j++){
			H1_aux(i-1,j-1) = A(i-1,j-1)*(_lambda*(1/(sqrt(2*ress))));
			H3_aux(i-1,j-1) = A(j-1,i-1)*(_lambda*(1/(sqrt(2*ress))));
		}
	}

	MatrixXcd H1(ress, ress); MatrixXcd H3(ress, ress);
	H1.setZero(); H3.setZero();

	H1 << H1_aux - H1_aux.transpose();
	H3 << H3_aux - H3_aux.transpose();

	MatrixXcd H(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * ress);
	H.setZero();

	H << complex_identity*(Kronecker_Product(H1, paulimatrix_x) + Kronecker_Product(H3, paulimatrix_z));

	*H_pointer << H;
}

void AltZir_DIII::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps){

	std::ofstream output_G("Data_Analysis/Channel/Andreev_G_DIII_Channel.txt");
	std::ofstream output_P("Data_Analysis/Channel/Andreev_P_DIII_Channel.txt");

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
