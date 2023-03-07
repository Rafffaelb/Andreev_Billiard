#include <iostream>
#include <chrono>
#include "../../include/AltZir_AbstractClass_h/AltZir.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include "../../include/AltZir_AbstractClass_h/Run_Simulation_Conductance_Energy.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

void AltZir::Run_Simulation_Conductance_Energy(){

	auto start = chrono::system_clock::now();

	double Gamma, Delta, y, V, Energy, small_gamma;
       	int N1, N2, n, _num_steps;

	Delta = 0.01;
	Gamma = 1;
	y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));
	int ress = 150;
	_num_steps = 5000; 

	for (int i = 1; i < 16; i++){

		MatrixXcd G(_num_steps, 101);

		G.setZero();

		N1 = i;

		N2 = N1;
		n = N1 + N2;

		// Create W Matrices //

		MatrixXcd W(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * n);
		W.setZero();
		MatrixXcd* W_pointer = &W;

		Create_W(W_pointer, ress, N1, N2, _lambda, y);

		// Create_ProjectionMatrices //
		
		MatrixXcd C1(_electron_hole_deg * _spin_deg * 2*N1, _electron_hole_deg * _spin_deg * 2*N1);
		MatrixXcd C2(_electron_hole_deg * _spin_deg * 2*N2, _electron_hole_deg * _spin_deg * 2*N2);

		C1.setZero();
		C2.setZero();

		MatrixXcd *C1_pointer = &C1;
		MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
		
		for (int step = 1; step < _num_steps + 1; step++){
		
			// Generate Hamiltonian Matrix //

			MatrixXcd H(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * ress);
			H.setZero();
			MatrixXcd* H_pointer = &H;

			Create_H(H_pointer, ress, _lambda);

			// Create billiard setup //

			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2, N1, N2, _electron_hole_deg);

			small_gamma = (N1*Gamma*Delta/(4*M_PI));

			#pragma omp parallel for shared(H, W, C1, C2) firstprivate(billiard_setup)
			for (int energy_idx = 1; energy_idx < 102; energy_idx++){

				Energy = small_gamma*((double)energy_idx-51);

				// Scattering Matrix //

				billiard_setup.Calculate_Smatrix(Energy);

				// Conductance (G) and Power Shot Noise (P) //
		
				billiard_setup.Calculate_G_and_P();

				G(step-1, energy_idx-1) = billiard_setup.getG();
			}

			if (step % _num_steps == 0){
				std::cout << "\nCurrent number of steps: " << step << "| Current number of open Channel N: " << N1 <<  std::endl;
			}
		}
		//Save G matrix as txt files //
		Save_txt_files_Energy(G, _num_steps, N1);

	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}
