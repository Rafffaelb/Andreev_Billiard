#include <iostream>
#include "../include/AltZir_D.h"
#include <cmath>
#include <complex>
#include <random>
#include <ctime>
#include <chrono>
#include <fstream>
#include <string>

using namespace std;

AltZir_D::AltZir_D(double lambda, int spin_deg, int electron_hole_deg){

	this -> _lambda = lambda;
	this -> _spin_deg = spin_deg;
	this -> _electron_hole_deg = electron_hole_deg;
}

AltZir_D::~AltZir_D() {}

void AltZir_D::Create_W(MatrixXcd* W_pointer, int ress, int N1, int N2, double lambda, double y){

	MatrixXcd W1(ress,N1); MatrixXcd W2(ress,N2);

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

	MatrixXcd W(ress, (N1+N2));
	W.setZero();
	W << W1, W2;

	*W_pointer = W;
}

void AltZir_D::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);

	int n = N1 + N2;

	MatrixXcd C1(n,n);
	MatrixXcd C2(n,n);

	C1.block(0, 0, N1, N1) << identity1; C1.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C1.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C1.block(N1, N1, N2, N2) << MatrixXcd::Zero(N2, N2);

	C2.block(0, 0, N1, N1) << MatrixXcd::Zero(N1, N1); C2.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C2.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C2.block(N1, N1, N2, N2) << identity2;

	*C1_pointer << C1;
	*C2_pointer << C2;
}

void AltZir_D::Create_H(MatrixXcd* H_pointer, int ress, double _lambda){

	complex<double> complex_identity(0,1);

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

	MatrixXcd H_aux(ress, ress);
	H_aux.setZero();

	for (int i = 1; i < ress + 1; i++){
		H_aux(i-1,i-1) = 0;
		for(int j = i + 1; j < ress + 1; j++){
			H_aux(i-1,j-1) = (_lambda*(1/(sqrt(ress))))*A(i-1,j-1);
		}
	}
	
	MatrixXcd Antisymmetric(ress, ress);
	Antisymmetric.setZero();

	Antisymmetric << H_aux - H_aux.transpose();

	MatrixXcd H(ress, ress);
	H.setZero();

	H << complex_identity*Antisymmetric;
	*H_pointer = H;	
}

void AltZir_D::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps){
	std::ofstream output_G("Data_Analysis/Channel/Andreev_G_D_Channel.txt");
	std::ofstream output_P("Data_Analysis/Channel/Andreev_P_D_Channel.txt");

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

void AltZir_D::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1){

	std::ofstream output_G("Data_Analysis/Gamma/Andreev_G_D_Gamma_N"+to_string(N1)+".txt");
	std::ofstream output_P("Data_Analysis/Gamma/Andreev_P_D_Gamma_N"+to_string(N1)+".txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 21; j++){
			if (j == 20){
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

void AltZir_D::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps, int N1) {}

void AltZir_D::Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps) {}

void AltZir_D::Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps) {}

void AltZir_D::Save_txt_files_Bell_Parameter_Fixed_Base(MatrixXd Bell_Parameter_Fixed_Base, int num_steps) {}

void AltZir_D::Save_txt_files_Energy(MatrixXcd G, int num_steps, int N1) {}

void AltZir_D::Save_txt_files_Energy_Gamma(MatrixXcd G, int num_steps, int N1, int gamma_idx) {}



