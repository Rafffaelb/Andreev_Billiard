#include <iostream>
#include <chrono>
#include "../../include/AltZir_AbstractClass_h/AltZir.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include "../../include/AltZir_AbstractClass_h/Run_Simulation_Conductance_Energy_Gamma.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

void AltZir::Run_Simulation_Conductance_Energy_Gamma(){

        auto start = chrono::system_clock::now();

        double Gamma, Delta, y, Energy, small_gamma;
        int N1, N2, n, ress, _num_steps;

        _num_steps = 30000;

        Eigen::VectorXf N_vector(6); N_vector << 1, 2, 5, 10, 15, 30;
        Eigen::VectorXf ress_vector(6); ress_vector << 100, 100, 100, 100, 150, 300;
        Eigen::VectorXf Delta_vector(6); Delta_vector << 0.01, 0.01, 0.01, 0.01, 0.0066, 0.0033;

        for (int i = 0; i < 6; i++){

                MatrixXcd G(_num_steps, 101);
                G.setZero();

                N1 = N_vector[i];
                ress = ress_vector[i];

                Delta = Delta_vector[i];

                N2 = N1;
                n = N1 + N2;

                // Create_ProjectionMatrices //

                MatrixXcd C1(_electron_hole_deg * _spin_deg * 2*N1, _electron_hole_deg * _spin_deg * 2*N1);
	       	MatrixXcd C2(_electron_hole_deg * _spin_deg * 2*N2, _electron_hole_deg * _spin_deg * 2*N2);
                
		C1.setZero(); 
		C2.setZero();
                
		MatrixXcd *C1_pointer = &C1;
	       	MatrixXcd *C2_pointer = &C2;

                Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);

		for (int gamma_idx = 1; gamma_idx < 12; gamma_idx++){

                        if (gamma_idx == 1){

                                Gamma = 0.001;
                        }
                        else{

                                Gamma = double(gamma_idx - 1) / double(10);
                        }

                        y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));


                        // Create W Matrices //

                        MatrixXcd W(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * n);
                        W.setZero();
                        MatrixXcd *W_pointer = &W;

                        Create_W(W_pointer, ress, N1, N2, _lambda, y);
			
			for (int step = 1; step < _num_steps + 1; step++){
			
				// Generate Hamiltonian Matrix //
				
				MatrixXcd H(_electron_hole_deg * _spin_deg * ress, _electron_hole_deg * _spin_deg * ress);
				H.setZero();
				MatrixXcd* H_pointer = &H;
				
				Create_H(H_pointer, ress, _lambda);
				// Create billiard setup //
				
				Quantum_chaotic_billiard billiard_setup(H, W, C1, C2, N1, N2, _electron_hole_deg);
				
				if (N1 == 5 || N1 == 10){
					small_gamma = (gamma_idx < 4) ? (Delta/(4*M_PI)) : (N1*Gamma*Delta/ (4*M_PI));
				}
				else{
					if (N1 == 15 || N1 == 30){
						small_gamma = (gamma_idx < 2) ? (Delta/(4*M_PI)) : (N1*Gamma*Delta/ (4*M_PI));
					}
					else{
						small_gamma = (gamma_idx < 7) ? (Delta/(4*M_PI)) : (N1*Gamma*Delta/ (2*M_PI));
					}
				}

				#pragma omp parallel for shared(W, C1, C2, H) firstprivate(billiard_setup)
				for (int energy_idx = 1; energy_idx < 102; energy_idx++){
					
					Energy = small_gamma*((double)energy_idx - 51);
					
					// Scattering Matrix //
					billiard_setup.Calculate_Smatrix(Energy);
					
					// Conductance (G) and Power Shot Noise (P) //

					billiard_setup.Calculate_G_and_P();


					G(step-1, energy_idx-1) = billiard_setup.getG();
				}
				
				if (step % _num_steps == 0){
					std::cout << "\nCurrent number of steps: " << step << "| Current index of Gamma: " << gamma_idx << "| Current number of open Channel N: " << N1 << std::endl;
				}
			}
			//Save G matrix as txt files //
			
			Save_txt_files_Energy_Gamma(G, _num_steps, N1, gamma_idx);
		}
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end - start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}

