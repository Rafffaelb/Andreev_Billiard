#include <eigen3/Eigen/Dense>
#ifndef ALTZIR_H
#define ALTZIR_H

using namespace Eigen;

class AltZir{

	protected:

		int _num_steps;
		int _spin_deg;
		int _electron_hole_deg;
		double _lambda;

	public:

		AltZir();

		virtual ~AltZir();

		virtual void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double y) = 0;

		void Run_Simulation_Conductance_Channels();
		void Run_Simulation_Conductance_Gamma();
		void Run_Simulation_Concurrence_Gamma();		


		virtual void Create_H(MatrixXcd* H_pointer, int _ress, double _lambda) = 0;

		virtual void Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2) = 0;

		virtual void Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps) = 0;
		virtual void Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1) = 0;
		virtual void Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps, int N1) = 0;
};

#endif
