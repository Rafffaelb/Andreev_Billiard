#include <iostream>
#include <chrono>
#include "../../include/AltZir_AbstractClass_h/AltZir.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include "../../include/AltZir_AbstractClass_h/Run_Simulation_Conductance_Gamma.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

void AltZir::Run_Simulation_Conductance_Gamma(){

	auto start = chrono::system_clock::now();

	double Gamma, y, V;
       	int ress, N1, N2, n, _num_steps;

	ress = 100;
	_num_steps = 1000000;
	V = _lambda*_lambda/ress;

	for (int i = 1; i < 5; i++){

		if (i == 4){
			N1 = 10;
		}
		else{
			N1 = i;
		}

		N2 = N1;
		n = N1 + N2;

		// Create_ProjectionMatrices //

		MatrixXcd C1(_electron_hole_deg * _spin_deg * 2*N1, _electron_hole_deg * _spin_deg * 2*N1); MatrixXcd C2(_electron_hole_deg *_spin_deg * 2*N2, _electron_hole_deg * _spin_deg * 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
	
		MatrixXcd G(_num_steps, 21);
		MatrixXcd P(_num_steps, 21);

		G.setZero();
		P.setZero();

		for (int gamma_idx = 1; gamma_idx < 22; gamma_idx++){
	
			if (gamma_idx == 1){
				
				Gamma = 0.0001;	
			}
			else{

				Gamma = double(gamma_idx - 1) / double(20);
			}

			double y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));

			// Create W Matrices //

			MatrixXcd W(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * n);
			W.setZero();
			MatrixXcd *W_pointer = &W;

			Create_W(W_pointer, ress, N1, N2, _lambda, y);
		
			#pragma omp parallel for shared(W, C1, C2)
			for (int step = 1; step < _num_steps + 1; step++){
		
				// Generate Hamiltonian Matrix //

				MatrixXcd H(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * ress);
				H.setZero();
				MatrixXcd* H_pointer = &H;

				Create_H(H_pointer, ress, _lambda);

				// Create billiard setup //

				Quantum_chaotic_billiard billiard_setup(H, W, C1, C2, N1, N2, _electron_hole_deg);

				// Scattering Matrix //
		
				billiard_setup.Calculate_Smatrix(0);

				// Conductance (G) and Power Shot Noise (P) //
		
				billiard_setup.Calculate_G_and_P();

				G(step-1, gamma_idx-1) = billiard_setup.getG();
				P(step-1, gamma_idx-1) = billiard_setup.getP();
				
				if (step % _num_steps == 0){
					std::cout << "\nCurrent number of steps: " << step << "| Current index of Gamma: " << gamma_idx << "| Current number of open Channel N: " << N1 <<  std::endl;
				}
			}
	
			//Save G and P matrices as txt files //
			Save_txt_files_Gamma(G, P, _num_steps, N1);
		}
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}
